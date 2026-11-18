"""Color choices affect all selected ligands without recoloring heteroatoms."""

from pathlib import Path
import unittest
from unittest.mock import patch

from pymol import cmd, publication as pub


def atom_rgb(selection):
    colors = []
    cmd.iterate("first (%s)" % selection, "colors.append(color)", space={"colors": colors})
    return cmd.get_color_tuple(colors[0])


class ColorFixture:
    def setUp(self):
        cmd.reinitialize()
        cmd.load(str(Path(__file__).resolve().parents[2] / "data" / "1oky.pdb.gz"), "complex")
        cmd.create("second", "complex and resn STU")
        self.ligands = pub.detect_scene().ligands
        self.assertEqual(len(self.ligands), 2)

    def assertRGB(self, actual, expected):
        for value, target in zip(actual, expected):
            self.assertAlmostEqual(value, target, delta=0.004)

    def add_second_chain(self):
        # Different IDs in the same object, even with case-insensitive selectors.
        cmd.alter("polymer.protein and resi 72-150", 'chain = "a"')
        cmd.set("ignore_case_chain", 1)
        return ("polymer.protein and resi 72-150", "polymer.protein and resi 151-359")


class ColorBehavior(ColorFixture, unittest.TestCase):
    def test_default_chain_and_element_colors_are_independent(self):
        first, second = self.add_second_chain()
        fig = pub.PublicationFigure()
        complex_ligand = self.ligands[0]
        before = None
        for theme in (pub.THEME_JOURNAL, pub.THEME_ILLUSTRATION):
            fig.apply_style(theme=theme, ligands=[complex_ligand])
            colors = atom_rgb(first), atom_rgb(second)
            self.assertNotEqual(*colors)
            if before:
                self.assertNotEqual(colors[0], before[0])
                self.assertNotEqual(colors[1], before[1])
            before = colors
            ligand = f"({complex_ligand.selection})"
            self.assertGreater(cmd.count_atoms(f"{ligand} and rep sticks"), 0)
            self.assertEqual(cmd.count_atoms("second and rep sticks"), 0)
            self.assertRGB(atom_rgb(f"{ligand} and elem N"), cmd.get_color_tuple("nitrogen"))
            self.assertRGB(atom_rgb(f"{ligand} and elem O"), cmd.get_color_tuple("oxygen"))
            self.assertNotEqual(atom_rgb(f"{ligand} and elem C"), colors[0])
            self.assertNotEqual(atom_rgb(f"{ligand} and elem C"), colors[1])

    def test_all_themes_color_all_selected_ligands(self):
        fig = pub.PublicationFigure()
        for theme in pub.THEME_ORDER:
            fig.apply_style(theme=theme, ligands=self.ligands, keep_view=True,
                            protein_color=(0.2, 0.4, 0.6), ligand_color=(1, 0, 1))
            self.assertRGB(atom_rgb("polymer.protein and elem C"), (0.2, 0.4, 0.6))
            for ligand in self.ligands:
                self.assertRGB(atom_rgb("(%s) and elem C" % ligand.selection), (1, 0, 1))
                self.assertRGB(atom_rgb("(%s) and elem N" % ligand.selection), cmd.get_color_tuple("nitrogen"))

    def test_invalid_color_rejected_before_snapshot(self):
        fig = pub.PublicationFigure()
        before = cmd.get_setting_tuple("bg_rgb")
        for color in ((float("nan"), 0, 0), (-1, 0, 0), (0, 1), (0, 0, 2)):
            with self.assertRaises(ValueError):
                fig.apply_style(ligands=self.ligands, protein_color=color)
            self.assertFalse(fig.has_snapshot())
            self.assertEqual(cmd.get_setting_tuple("bg_rgb"), before)

    def test_color_does_not_change_disabled_object(self):
        cmd.fragment("benzene", "disabled")
        cmd.color("green", "disabled")
        cmd.disable("disabled")
        pub.PublicationFigure().apply_style(ligands=self.ligands, ligand_color=(1, 0, 1))
        self.assertRGB(atom_rgb("disabled and elem C"), cmd.get_color_tuple("green"))

    def test_deep_theme_reset_clears_atom_override(self):
        fig = pub.PublicationFigure()
        fig.apply_style(ligands=self.ligands, theme=pub.THEME_SURFACE)
        fig.apply_style(ligands=self.ligands, theme=pub.THEME_JOURNAL)
        values = []
        cmd.iterate("resn STU", "values.append(s.sphere_scale)", space={"values": values})
        self.assertTrue(all(abs(value - 1) < 0.001 for value in values))


class ColorControls(ColorFixture, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pymol.Qt import QtWidgets
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_ligand_color_is_shared_and_cancel_keeps_it(self):
        from pymol.Qt import QtWidgets, QtGui
        from pmg_qt.publication_color_controls import PublicationColorControls
        theme = ["journal"]
        controls = PublicationColorControls(cmd, lambda: theme[0])
        controls.set_color("ligand", (1, 0, 1), emit=False)
        theme[0] = "surface"
        self.assertEqual(controls.overrides()["ligand_color"], (1, 0, 1))
        self.assertEqual(controls.overrides()["ligand_mode"], "carbon")
        with patch.object(QtWidgets.QColorDialog, "getColor", return_value=QtGui.QColor()):
            controls.choose("ligand")
        self.assertEqual(controls.overrides()["ligand_color"], (1, 0, 1))
        controls.reset_current()
        self.assertIsNone(controls.overrides()["ligand_color"])
        self.assertEqual(controls.overrides()["ligand_mode"], "element")
        controls.deleteLater()

    def test_preview_and_export_keep_all_colors(self):
        from pymol.Qt import QtGui, QtWidgets
        from pmg_qt import publication_gui as pg
        first, second = self.add_second_chain()
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        dialog = pg.get_or_create_publication_dialog(parent)
        self.addCleanup(pg.clear_publication_dialog_ref)
        dialog._on_select_all()
        controls = dialog.color_controls
        self.assertEqual(controls.chain_combo.count(), 2)
        controls.chain_combo.setCurrentIndex(next(
            i for i in range(controls.chain_combo.count())
            if controls.chain_combo.itemData(i) == ("complex", "", "a")
        ))
        self.assertEqual(controls.chain_combo.currentData(), ("complex", "", "a"))
        with patch.object(QtWidgets.QColorDialog, "getColor", return_value=QtGui.QColor("#FF00FF")):
            controls.protein_button.click()
        dialog.on_preview()
        self.assertRGB(atom_rgb(first), (1, 0, 1))
        self.assertRGB(atom_rgb("second and elem N"), cmd.get_color_tuple("nitrogen"))
        other_chain_color = atom_rgb(second)

        self.assertTrue(pg.apply_smart_theme(parent, pub.THEME_ILLUSTRATION))
        self.assertEqual(controls.chain_combo.currentData(), ("complex", "", "a"))
        self.assertRGB(atom_rgb(first), (1, 0, 1))
        self.assertNotEqual(atom_rgb(second), other_chain_color)
        self.assertRGB(atom_rgb("second and elem N"), cmd.get_color_tuple("nitrogen"))

        seen = {}
        real_export = dialog.fig.export_png
        def export(*args, **kwargs):
            seen.update(kwargs)
            return real_export(*args, **kwargs)
        class InlineWorker:
            def __init__(self, function, finishslot):
                self.function, self.finishslot = function, finishslot
            def __call__(self):
                self.finishslot((self.function(), None))
            def isRunning(self):
                return False
        with patch.object(pg, "getSaveFileNameWithExt", return_value="ignored.png"), \
                patch.object(pg, "AsyncFunc", InlineWorker), \
                patch.object(dialog.fig, "export_png", side_effect=export), \
                patch.object(cmd, "png"):
            dialog.on_export()
        self.assertEqual(len(seen["ligands"]), 2)
        self.assertEqual(seen["theme"], pub.THEME_ILLUSTRATION)
        self.assertEqual(seen["protein_mode"], "chain")
        self.assertEqual(seen["ligand_mode"], "element")
        self.assertEqual(seen["chain_colors"], {("complex", "", "a"): (1, 0, 1)})
        self.assertRGB(atom_rgb(first), (1, 0, 1))
        parent.deleteLater()


if __name__ == "__main__":
    unittest.main()
