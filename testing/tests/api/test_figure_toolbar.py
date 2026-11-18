# ==================================================
# 功能说明：出图工具栏纯逻辑与可选 Qt 集成的行为测试。
# 使用方法：由仓库测试运行器或 python -m unittest 发现；禁止依赖 pytest。
# 依赖环境：已安装的 pymol；Qt 用例外跳过不可用环境。
# 生成时间：2026-11-14 13:10
# ==================================================

"""figure_toolbar 设置、表示与 keep_appearance 行为测试（标准 unittest）。"""

from __future__ import annotations

import unittest

from pymol import cmd
from pymol import publication as pub

from pmg_qt import figure_toolbar as ft


def _qt_available() -> bool:
    """检测是否可导入 pymol.Qt。"""
    try:
        from pymol.Qt import QtWidgets  # noqa: F401
        return True
    except Exception:
        return False


class TestFigureToolbarSettings(unittest.TestCase):
    """背景 / 光照 / 效果应写入真实 PyMOL 设置。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        ft.clear_appearance_customized(cmd)

    def test_background_disables_gradient_and_matches(self) -> None:
        """背景色块关闭渐变，且能匹配 white/grey/black。"""
        ft.set_background(cmd, "black")
        self.assertEqual(int(cmd.get_setting_int("bg_gradient")), 0)
        self.assertEqual(ft.match_background(cmd), "black")
        self.assertTrue(ft.is_appearance_customized(cmd))

        ft.set_background(cmd, "white")
        self.assertEqual(ft.match_background(cmd), "white")
        ft.set_background(cmd, "grey")
        self.assertEqual(ft.match_background(cmd), "grey")
        cmd.set("bg_gradient", 1)
        self.assertIsNone(ft.match_background(cmd))

    def test_lighting_presets_deterministic_and_independent(self) -> None:
        """光照预设为完整参数包；切换后与目标一致，且不改阴影/轮廓。"""
        cmd.set("ray_shadows", 1)
        cmd.set("ray_trace_mode", 1)
        ft.apply_light_preset(cmd, "平光")
        for name, val in ft.LIGHT_PRESETS["平光"].items():
            try:
                self.assertAlmostEqual(float(cmd.get_setting_float(name)), float(val), places=2)
            except Exception:
                if name == "ambient_occlusion_mode":
                    continue
                raise
        self.assertEqual(int(cmd.get_setting_int("ray_shadows")), 1)
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 1)

        ft.apply_light_preset(cmd, "立体")
        self.assertEqual(ft.match_light_preset(cmd), "立体")
        # 前一预设痕迹不应残留为“平光”
        self.assertNotAlmostEqual(
            float(cmd.get_setting_float("ambient")),
            float(ft.LIGHT_PRESETS["平光"]["ambient"]),
            places=2,
        )

    def test_effect_buttons_shadow_outline_fog(self) -> None:
        """效果开关对应 ray_shadows / ray_trace_mode / depth_cue。"""
        ft.set_effect_shadow(cmd, True)
        ft.set_effect_outline(cmd, True)
        ft.set_effect_fog(cmd, True)
        st = ft.get_effect_state(cmd)
        self.assertTrue(st["shadow"])
        self.assertTrue(st["outline"])
        self.assertTrue(st["fog"])

        ft.set_effect_shadow(cmd, False)
        ft.set_effect_outline(cmd, False)
        ft.set_effect_fog(cmd, False)
        st = ft.get_effect_state(cmd)
        self.assertFalse(st["shadow"])
        self.assertFalse(st["outline"])
        self.assertFalse(st["fog"])


class TestFigureToolbarRepresentations(unittest.TestCase):
    """结构表示作用域与对象保全。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        ft.clear_appearance_customized(cmd)

    def test_empty_scene_noop(self) -> None:
        """空场景调用表示切换不应抛错、不造对象。"""
        before = cmd.get_names("all")
        ft.apply_representation(cmd, "cartoon")
        self.assertEqual(cmd.get_names("all"), before)

    def test_representation_excludes_solvent_and_preserves_disabled(self) -> None:
        """仅改启用分子的结构表示；禁用对象与溶剂原子数不变。"""
        cmd.fab("ACDEFGHIK", "prot")
        cmd.fragment("benzene", "offmol")
        cmd.disable("offmol")
        cmd.show("lines", "offmol")
        # 溶剂分子对象
        cmd.read_pdbstr(
            "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O\nEND\n",
            "wat",
        )
        cmd.show("nb_spheres", "wat")
        n_wat = cmd.count_atoms("wat")
        n_off = cmd.count_atoms("offmol")
        off_lines = cmd.count_atoms("offmol and rep lines")

        ft.apply_representation(cmd, "cartoon")
        self.assertGreater(cmd.count_atoms("prot and rep cartoon"), 0)
        self.assertEqual(cmd.count_atoms("prot and not solvent and rep sticks"), 0)
        self.assertEqual(cmd.count_atoms("wat"), n_wat)
        self.assertEqual(cmd.count_atoms("offmol"), n_off)
        self.assertEqual(cmd.count_atoms("offmol and rep lines"), off_lines)
        self.assertNotIn("offmol", cmd.get_names("objects", enabled_only=1))

        ft.apply_representation(cmd, "spheres")
        self.assertGreater(cmd.count_atoms("prot and not solvent and rep spheres"), 0)
        self.assertEqual(cmd.count_atoms("prot and not solvent and rep cartoon"), 0)

    def test_cartoon_keeps_ligand_visible_and_hides_initial_lines(self):
        cmd.fab("AAAA", "prot")
        cmd.fragment("benzene", "lig")
        cmd.show("lines", "prot or lig")
        ft.apply_representation(cmd, "cartoon")
        self.assertGreater(cmd.count_atoms("prot and rep cartoon"), 0)
        self.assertGreater(cmd.count_atoms("lig and rep sticks"), 0)
        self.assertEqual(cmd.count_atoms("prot and rep lines"), 0)

    def test_labels_and_measurement_untouched(self) -> None:
        """结构切换不清除标签表示，也不删除测量对象。"""
        cmd.fab("AAAA", "prot")
        cmd.label("prot and name CA", '"A"')
        cmd.distance("d1", "prot and name CA and resi 1", "prot and name CA and resi 2")
        names_before = set(cmd.get_names("all"))
        n_label = cmd.count_atoms("prot and rep labels")

        ft.apply_representation(cmd, "sticks")
        self.assertEqual(set(cmd.get_names("all")), names_before)
        self.assertIn("d1", cmd.get_names("all"))
        self.assertEqual(cmd.count_atoms("prot and rep labels"), n_label)


class TestPublicationKeepAppearance(unittest.TestCase):
    """出图 keep_appearance 保留背景与光照，仍改表示。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        ft.clear_appearance_customized(cmd)

    def test_keep_appearance_preserves_bg_light_shadow(self) -> None:
        """keep_appearance=True 时不强制白底/默认光照/关阴影。"""
        cmd.fab("ACDEFGHIK", "prot")
        ft.set_background(cmd, "black")
        ft.apply_light_preset(cmd, "立体")
        ft.set_effect_shadow(cmd, True)
        amb = float(cmd.get_setting_float("ambient"))

        fig = pub.PublicationFigure(cmd)
        fig.apply_style(style=pub.STYLE_OVERVIEW, ligand=None, keep_view=True, keep_appearance=True)

        self.assertEqual(ft.match_background(cmd), "black")
        self.assertAlmostEqual(float(cmd.get_setting_float("ambient")), amb, places=2)
        self.assertEqual(int(cmd.get_setting_int("ray_shadows")), 1)
        # 表示样式仍应施加
        self.assertGreater(cmd.count_atoms("prot and polymer.protein and rep cartoon"), 0)

    def test_default_appearance_still_white(self) -> None:
        """未 keep_appearance 时仍为白底默认出版外观。"""
        cmd.fab("ACDEFGHIK", "prot")
        ft.set_background(cmd, "black")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(style=pub.STYLE_OVERVIEW, ligand=None, keep_view=True, keep_appearance=False)
        self.assertEqual(ft.match_background(cmd), "white")

    def test_default_resets_outline_and_lighting(self):
        cmd.fragment("benzene", "lig")
        ligand = pub.detect_scene(cmd).ligands[0]
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(ligand=ligand)
        default_reflect = cmd.get_setting_float("reflect")
        ft.apply_light_preset(cmd, "立体")
        ft.set_effect_outline(cmd, True)
        fig.apply_style(ligand=ligand)
        self.assertEqual(cmd.get_setting_int("ray_trace_mode"), 0)
        self.assertEqual(cmd.get_setting_int("ambient_occlusion_mode"), 0)
        self.assertAlmostEqual(cmd.get_setting_float("reflect"), default_reflect)

    def test_export_transparency_is_per_export(self):
        from unittest.mock import patch
        cmd.fragment("benzene", "lig")
        ft.set_background(cmd, "grey")
        fig = pub.PublicationFigure(cmd)
        ligand = pub.detect_scene(cmd).ligands[0]
        with patch.object(cmd, "png"):
            fig.export_png("ignored.png", ligand=ligand, transparent=True, keep_appearance=True)
            self.assertEqual(cmd.get_setting_int("ray_opaque_background"), 0)
            fig.export_png("ignored.png", ligand=ligand, transparent=False, keep_appearance=True)
            self.assertEqual(cmd.get_setting_int("ray_opaque_background"), 1)
        self.assertEqual(ft.match_background(cmd), "grey")

    def test_keep_appearance_with_theme_honors_toolbar(self):
        """keep_appearance=True 时主题不得覆盖工具栏背景/光照/轮廓/阴影。"""
        cmd.fab("ACDEFGHIK", "prot")
        ft.set_background(cmd, "black")
        ft.apply_light_preset(cmd, "立体")
        ft.set_effect_shadow(cmd, True)
        ft.set_effect_outline(cmd, True)
        amb = float(cmd.get_setting_float("ambient"))
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(
            theme=pub.THEME_ILLUSTRATION,
            ligand=None,
            keep_view=True,
            keep_appearance=True,
        )
        self.assertEqual(ft.match_background(cmd), "black")
        self.assertAlmostEqual(float(cmd.get_setting_float("ambient")), amb, places=2)
        self.assertEqual(int(cmd.get_setting_int("ray_shadows")), 1)
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 1)

    def test_theme_false_keep_appearance_sets_theme_look(self):
        """keep_appearance=False 时由所选主题决定外观。"""
        cmd.fab("ACDEFGHIK", "prot")
        ft.set_background(cmd, "black")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(
            theme=pub.THEME_ILLUSTRATION,
            ligand=None,
            keep_view=True,
            keep_appearance=False,
        )
        self.assertEqual(cmd.get_setting_int("bg_rgb"), cmd.get_color_index("0xF3EDE3"))
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 1)


@unittest.skipUnless(_qt_available(), "pymol.Qt 不可用")
class TestFigureToolbarQt(unittest.TestCase):
    """可选：工具栏控件与勾选同步的轻量 Qt 检查。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        ft.clear_appearance_customized(cmd)
        from pymol.Qt import QtWidgets
        app = QtWidgets.QApplication.instance()
        if app is None:
            self.app = QtWidgets.QApplication([])
        else:
            self.app = app

    def test_clicking_smart_styles_changes_the_viewport_settings(self) -> None:
        from pymol.Qt import QtWidgets
        from pmg_qt.publication_gui import clear_publication_dialog_ref

        self.addCleanup(clear_publication_dialog_ref)
        cmd.fab("ACDEFGHIK", "prot")
        win = QtWidgets.QMainWindow()
        win.cmd = cmd
        self.addCleanup(win.close)
        bar = ft.install_figure_toolbar(win, cmd)

        bar._theme_btns[pub.THEME_JOURNAL].click()
        self.assertEqual(ft.match_background(cmd), "white")
        journal_radius = cmd.get_setting_float("cartoon_loop_radius", "prot")

        bar._theme_btns[pub.THEME_ILLUSTRATION].click()
        self.assertEqual(cmd.get_setting_int("bg_rgb"), cmd.get_color_index("0xF3EDE3"))
        self.assertGreater(cmd.get_setting_float("cartoon_loop_radius", "prot"), journal_radius)

        bar._theme_btns[pub.THEME_JOURNAL].click()
        self.assertEqual(ft.match_background(cmd), "white")
        self.assertAlmostEqual(cmd.get_setting_float("cartoon_loop_radius", "prot"), journal_radius)

    def test_toolbar_buttons_drive_settings(self) -> None:
        """点击背景/效果按钮应改变设置，并在 refresh 后勾选正确。"""
        from pymol.Qt import QtWidgets

        win = QtWidgets.QMainWindow()
        # 最小反馈定时器，模拟主窗口挂载点
        win.feedback_timer = __import__("pymol.Qt", fromlist=["QtCore"]).QtCore.QTimer(win)
        win.feedback_timer.setSingleShot(True)
        win._quit_requested = False
        win._shutting_down = False
        win.menudict = {"Display": win.menuBar().addMenu("显示")}

        bar = ft.install_figure_toolbar(win, cmd)
        self.assertIs(win.figure_toolbar, bar)

        bar.btn_bg_black.click()
        self.assertEqual(ft.match_background(cmd), "black")
        bar.refresh_checked_states()
        self.assertTrue(bar.btn_bg_black.isChecked())

        bar.btn_shadow.setChecked(True)
        bar._on_shadow()
        self.assertEqual(int(cmd.get_setting_int("ray_shadows")), 1)

        # 外部改设置后刷新，勾选不得陈旧
        cmd.set("ray_shadows", 0)
        cmd.set("depth_cue", 0)
        bar.refresh_checked_states()
        self.assertFalse(bar.btn_shadow.isChecked())
        self.assertFalse(bar.btn_fog.isChecked())

        # 外部自定义状态不匹配预设时，独占组应清除旧勾选。
        cmd.bg_color("red")
        ft.apply_light_preset(cmd, "柔和")
        bar.refresh_checked_states()
        cmd.set("ambient", 0.91)
        bar.refresh_checked_states()
        self.assertFalse(any(b.isChecked() for b in (bar.btn_bg_white, bar.btn_bg_grey, bar.btn_bg_black)))
        self.assertFalse(any(b.isChecked() for b in bar._light_btns.values()))

        # 显示菜单含显隐动作
        acts = [a.text() for a in win.menudict["Display"].actions()]
        self.assertIn("出图工具栏", acts)

        win.close()


if __name__ == "__main__":
    unittest.main()
