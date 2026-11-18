"""Integration coverage for figures, structural comparison, and trajectory tools."""

import csv
from importlib import import_module
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from pymol import cmd, publication, publication_recipe, structure_compare, trajectory_analysis


def fixture_trajectory():
    cmd.fab("ACDEFGHIK", "trajectory")
    initial = cmd.get_coordset("trajectory")
    rotation = np.asarray(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
    cmd.load_coordset(initial @ rotation + (3, 5, 8), "trajectory", state=2)
    changed = initial.copy()
    changed[0] += (4, 0, 0)
    cmd.load_coordset(changed, "trajectory", state=3)


class TestFigureUndo(unittest.TestCase):
    def test_only_the_last_style_is_undone(self):
        cmd.reinitialize()
        cmd.fab("ACDEF", "prot")
        fig = publication.PublicationFigure(cmd, track_undo=True)
        fig.apply_style(theme=publication.THEME_JOURNAL)
        cmd.set("sphere_scale", 0.8, "prot and name CA")
        cmd.set_bond("stick_radius", 0.6, "prot and name N", "prot and name CA")
        before = cmd.get_coords("prot").copy()
        fig.apply_style(theme=publication.THEME_ILLUSTRATION)
        cmd.translate([3, 0, 0], "prot", camera=0)
        cmd.fragment("benzene", "new_object")
        cmd.show("spheres", "new_object")
        self.assertTrue(fig.undo_style())
        self.assertEqual(cmd.get_color_tuple(cmd.get_setting_int("bg_rgb")), (1, 1, 1))
        np.testing.assert_allclose(cmd.get_coords("prot"), before + [3, 0, 0], atol=1e-6)
        self.assertEqual(cmd.count_atoms("new_object and rep spheres"), cmd.count_atoms("new_object"))
        atom_scales = []
        cmd.iterate("prot and name CA", "scales.append(s.sphere_scale)", space={"scales": atom_scales})
        np.testing.assert_allclose(atom_scales, 0.8)
        restored = [value for _, bonds in cmd.get_bond("stick_radius", "prot")
                    for _, _, value in bonds if value is not None]
        self.assertTrue(restored)
        np.testing.assert_allclose(restored, 0.6)
        self.assertFalse(fig.undo_style())
        self.assertFalse(any(name.startswith("_pubfig_undo_") for name in cmd.get_names("selections")))


class TestComparisonAndTrajectory(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()
        fixture_trajectory()

    def test_comparison_of_two_states_keeps_sources_and_renders_both_sides(self):
        before = cmd.get_coords("trajectory", state=2).copy()
        result = structure_compare.compare(cmd, "trajectory", "trajectory", 1, 2,
                                           "polymer.protein", "polymer.protein")
        self.assertEqual(result.matched, 9)
        self.assertLess(result.rmsd, 1e-5)
        np.testing.assert_allclose(cmd.get_coords(result.mobile), cmd.get_coords(result.reference), atol=1e-5)
        structure_compare.show_overlay(cmd, result)
        view = cmd.get_view()
        self.assertNotIn("trajectory", cmd.get_names("objects", enabled_only=1))
        with TemporaryDirectory() as folder:
            left, right = [Path(folder) / name for name in ("ref.png", "mob.png")]
            structure_compare.render_pair(cmd, result, left, right, width=240, height=240)
            self.assertTrue(left.read_bytes().startswith(b"\x89PNG"))
            self.assertTrue(right.read_bytes().startswith(b"\x89PNG"))
        self.assertEqual(cmd.get_view(), view)
        structure_compare.clear_comparison(cmd, result)
        self.assertEqual(cmd.get_names("objects", enabled_only=1), ["trajectory"])
        np.testing.assert_array_equal(cmd.get_coords("trajectory", state=2), before)

    def test_rmsd_distance_range_and_csv(self):
        rmsd = trajectory_analysis.analyze(cmd, "trajectory", "rmsd", 1, 3, 1, "all")
        self.assertLess(max(rmsd.values[:2]), 1e-5)
        self.assertGreater(rmsd.values[2], 0.1)
        distance = trajectory_analysis.analyze(cmd, "trajectory", "distance", 1, 3, 2,
                                               "index 1", "index 2")
        self.assertEqual(distance.frames, (1, 3))
        expected = [np.linalg.norm(cmd.get_coords("index 1", state=frame) - cmd.get_coords("index 2", state=frame))
                    for frame in distance.frames]
        np.testing.assert_allclose(distance.values, expected)
        with TemporaryDirectory() as folder:
            path = Path(folder) / "distance.csv"
            distance.write_csv(path)
            with path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([row["state"] for row in rows], ["1", "3"])
            np.testing.assert_allclose([float(row["distance_angstrom"]) for row in rows], expected)
        self.assertFalse(any(name.startswith("_analysis_atoms") for name in cmd.get_names("selections")))

    def test_vmd_scope_uses_strict_residue_bounds_in_comparison(self):
        scope = "resid < 10 and resid > 5"
        result = structure_compare.compare(cmd, "trajectory", "trajectory", 1, 2, scope, scope,
                                           syntax="vmd")
        self.addCleanup(structure_compare.clear_comparison, cmd, result)
        residues = set()
        cmd.iterate(result.reference, "residues.add(resi)", space={"residues": residues})
        self.assertEqual(residues, {"6", "7", "8", "9"})
        self.assertEqual(result.matched, 4)


try:
    from pymol.Qt import PYQT_NAME, QtCore, QtWidgets
    QtTest = import_module(f"{PYQT_NAME}.QtTest")
except ImportError:
    QtWidgets = None


@unittest.skipIf(QtWidgets is None, "Qt unavailable")
class TestWorkspaceQt(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        cmd.reinitialize()
        self.window = QtWidgets.QMainWindow()
        self.window.cmd = cmd
        self.addCleanup(self.window.close)

    def test_periodic_box_toggle_uses_visible_objects_and_current_frame(self):
        from pmg_qt.movie_speed import PeriodicBoxButton
        from pymol.constants import repres

        fixture_trajectory()
        cmd.set_symmetry("trajectory", 40, 41, 42, 90, 90, 90, state=1)
        cmd.set_symmetry("trajectory", 51, 52, 53, 80, 85, 95, state=2)
        cmd.pseudoatom("hidden_box")
        cmd.set_symmetry("hidden_box", 20, 21, 22, 90, 90, 90)
        cmd.disable("hidden_box")
        cmd.pseudoatom("no_box")
        messages = []
        button = PeriodicBoxButton(cmd, messages.append, self.window)
        box_rep = repres["cell"]

        button.click()
        self.assertEqual(button.text(), "关闭盒子")
        self.assertTrue(button.isChecked())
        self.assertIn(box_rep, cmd.get_vis()["trajectory"][2])
        self.assertNotIn(box_rep, cmd.get_vis()["hidden_box"][2])
        self.assertNotIn(box_rep, cmd.get_vis()["no_box"][2])
        cmd.frame(2)
        button.sync()
        self.assertEqual(cmd.get_symmetry("trajectory")[:3], [51, 52, 53])
        self.assertTrue(button.isChecked())

        button.click()
        self.assertEqual(button.text(), "显示盒子")
        self.assertNotIn(box_rep, cmd.get_vis()["trajectory"][2])
        cmd.show("cell", "trajectory")
        button.sync()
        self.assertTrue(button.isChecked())
        cmd.disable("trajectory")
        button.click()
        self.assertFalse(button.isChecked())
        self.assertEqual(messages, ["当前可见分子没有有效的周期盒数据。"])

    def test_recipe_reopens_two_chain_figure_with_one_ligand(self):
        from pmg_qt.publication_gui import PublicationDialog

        cmd.load(str(Path(__file__).resolve().parents[2] / "data" / "1oky.pdb.gz"), "complex")
        cmd.alter("complex and polymer.protein and resi 72-150", 'chain = "B"')
        cmd.create("second_ligand", "complex and resn STU")
        cmd.orient("complex")
        cmd.turn("y", 22)
        dialog = PublicationDialog(self.window, cmd)
        dialog.set_theme_ui(publication.THEME_ILLUSTRATION)
        dialog.color_controls.set_chain_color(("complex", "", "B"), (0.8, 0.1, 0.2))
        dialog.ligand_list.item(1).setCheckState(QtCore.Qt.CheckState.Unchecked)
        dialog.spin_w.setValue(720)
        dialog.chk_keep_appearance.setChecked(True)
        dialog.chk_keep_view.setChecked(True)
        cmd.bg_color("0xEAF5FF")
        dialog.on_preview()
        view = cmd.get_view()
        with TemporaryDirectory() as folder:
            path = Path(folder) / "figure.json"
            session = Path(folder) / "scene.pse"
            with patch("pmg_qt.publication_gui.getSaveFileNameWithExt", return_value=str(path)):
                dialog.on_save_recipe()
            dialog.fig.clear_undo()
            cmd.save(str(session))
            cmd.reinitialize()
            cmd.load(str(session))
            cmd.bg_color("red")
            cmd.turn("x", 45)
            reopened = PublicationDialog(self.window, cmd)
            data = publication_recipe.load_recipe(path, cmd, publication.detect_scene())
            self.assertEqual(data["unmatched"], [])
            reopened.apply_recipe(data, style=data["style"], restore_view=True)
            np.testing.assert_allclose(cmd.get_view(), view, atol=1e-5)
            self.assertEqual(cmd.get_setting_int("bg_rgb"), cmd.get_color_index("0xEAF5FF"))
            self.assertEqual(reopened.spin_w.value(), 720)
            self.assertEqual([lig.object for lig in reopened._current_ligands()], ["complex"])
            self.assertEqual(cmd.count_atoms("second_ligand and rep sticks"), 0)
            chain_color = cmd.get_color_index("0xCC1A33")
            self.assertEqual(cmd.count_atoms(f"chain B and polymer.protein and color {chain_color}"),
                             cmd.count_atoms("chain B and polymer.protein"))
            reopened.fig.clear_snapshot()
            cmd.alter("chain B", 'chain = "C"')
            data = publication_recipe.load_recipe(path, cmd, publication.detect_scene())
            self.assertIn("链 complex//B", data["unmatched"])

    def test_curve_click_jumps_to_the_sampled_frame(self):
        from pmg_qt.trajectory_analysis_gui import TrajectoryAnalysisDialog

        fixture_trajectory()
        dialog = TrajectoryAnalysisDialog(self.window)
        dialog.calculate()
        plot = dialog.plot
        plot.resize(500, 250)
        coords = plot.axes()
        point = QtCore.QPoint(round(coords.right() - 2), round(coords.center().y()))
        QtTest.QTest.mouseClick(plot, QtCore.Qt.MouseButton.LeftButton, pos=point)
        self.assertEqual(cmd.get_state(), 3)
        self.assertIn("第 3 帧", dialog.status.text())

    def test_vmd_atom_selection_drives_trajectory_distance(self):
        from pmg_qt.trajectory_analysis_gui import TrajectoryAnalysisDialog

        fixture_trajectory()
        dialog = TrajectoryAnalysisDialog(self.window)
        dialog.syntax.setCurrentIndex(dialog.syntax.findData("vmd"))
        dialog.mode.setCurrentIndex(dialog.mode.findData("distance"))
        dialog.selection_a.setText("resname GLY and resid < 10 and resid > 5 and name CA")
        dialog.selection_b.setText("index 0")
        dialog.calculate()
        trace = dialog.plot.trace
        self.assertIsNotNone(trace)
        self.assertEqual(trace.selection_syntax, "vmd")
        expected = [np.linalg.norm(cmd.get_coords("resi 6 and name CA", state=frame)
                                   - cmd.get_coords("index 1", state=frame)) for frame in trace.frames]
        np.testing.assert_allclose(trace.values, expected, rtol=1e-6)
