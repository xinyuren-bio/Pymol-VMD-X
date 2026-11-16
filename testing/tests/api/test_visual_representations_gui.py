# ==================================================
# 功能说明：选择/表示 Qt 窗口与色板的 offscreen 动作、取消与关闭保护测试。
# 使用方法：由仓库测试运行器或 python -m unittest 发现；使用 QT_QPA_PLATFORM=offscreen。
# 依赖环境：pymol + PySide6（pymol.Qt）；无显示环境可跑。
# 生成时间：2026-11-14 16:50
# ==================================================

"""visual_representations_gui / color_palette 测试（标准 unittest）。"""

from __future__ import annotations

import os
import unittest
from unittest.mock import patch

from chempy import Atom
from chempy.models import Indexed

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _qt_available() -> bool:
    """检测 pymol.Qt 是否可用。"""
    try:
        from pymol.Qt import QtWidgets  # noqa: F401
        return True
    except Exception:
        return False


def _bead(name, resn, chain, resi, x, y, z):
    """构造珠。"""
    at = Atom()
    at.name = name
    at.resn = resn
    at.chain = chain
    at.resi = str(resi)
    at.symbol = "C"
    at.coord = [float(x), float(y), float(z)]
    return at


@unittest.skipUnless(_qt_available(), "需要 pymol.Qt")
class TestColorPalette(unittest.TestCase):
    """色板：取消不修改、确认才改、引用清理。"""

    @classmethod
    def setUpClass(cls) -> None:
        from pymol.Qt import QtWidgets

        app = QtWidgets.QApplication.instance()
        if app is None:
            cls.app = QtWidgets.QApplication([])
        else:
            cls.app = app

    def setUp(self) -> None:
        from pymol import cmd

        cmd.reinitialize()
        self.cmd = cmd
        model = Indexed()
        model.add_atom(_bead("CA", "ALA", "A", 1, 0, 0, 0))
        cmd.load_model(model, "palobj", zoom=0)

    def tearDown(self) -> None:
        from pmg_qt.color_palette import clear_color_palette_refs

        clear_color_palette_refs()

    def test_more_colors_cancel_no_change(self) -> None:
        """QColorDialog 取消不修改（模拟 invalid）。"""
        from pymol.Qt import QtGui, QtWidgets
        from pmg_qt import color_palette as cp

        parent = QtWidgets.QWidget()
        parent.cmd = self.cmd
        seen = []

        pop = cp.ColorPalettePopup(parent, current_rgb=(1.0, 0.0, 0.0))
        pop.colorSelected.connect(lambda r, g, b: seen.append((r, g, b)))

        # 猴子补丁：取消
        pop._on_more()
        pop._color_dialog.setCurrentColor(QtGui.QColor(0,255,0))
        pop._color_dialog.reject()
        self.assertEqual(seen, [])

    def test_swatch_emits_and_truecolor_api(self) -> None:
        """色块发出 RGB；truecolor 着色 API。"""
        from pymol.visual_representations import color_selection_truecolor, rgb_to_hex_color
        from pmg_qt.color_palette import ColorSwatchButton
        from pymol.Qt import QtWidgets

        parent = QtWidgets.QWidget()
        got = []
        btn = ColorSwatchButton((0.2, 0.4, 0.6), parent)
        btn.colorChosen.connect(lambda r, g, b: got.append((r, g, b)))
        btn.click()
        self.assertEqual(got, [(0.2, 0.4, 0.6)])
        hexcol = color_selection_truecolor(self.cmd, "palobj", (1.0, 0.0, 0.0))
        self.assertEqual(hexcol, rgb_to_hex_color((1.0, 0.0, 0.0)))

    def test_clear_refs(self) -> None:
        """清理全局引用。"""
        from pymol.Qt import QtWidgets
        from pmg_qt import color_palette as cp

        parent = QtWidgets.QWidget()
        parent.cmd = self.cmd
        cp.open_object_palette(parent, target="palobj")
        cp.clear_color_palette_refs()
        self.assertIsNone(cp._object_palette_ref)

    def test_property_buttons_color_only_the_target(self) -> None:
        from pymol.Qt import QtWidgets
        from pmg_qt.color_palette import open_object_palette

        cmd = self.cmd
        target = "sample-1"
        selection = "'sample-1'"
        data = (
            ("C1", "1", "ALA", "H", 10),
            ("C2", "01", "ALA", "S", 50),
            ("C3", "1", "GLY", "L", 90),
            ("C4", "01", "GLY", "H", 10),
            ("C5", "1", "ALA", "", 120),
        )
        for idx, (name, chain, resn, ss, b) in enumerate(data, 1):
            cmd.pseudoatom(target, name=name, chain=chain, resn=resn,
                           resi=str(idx), b=b, pos=[idx, 0, 0])
        cmd.alter(selection, "ss = structures[name]",
                  space={"structures": {name: ss for name, _, _, ss, _ in data}})
        cmd.select("palette_selection", f"{selection} and not name C5")

        def colors(selection):
            result = {}
            cmd.iterate(selection, "result[name] = color", space={"result": result})
            return result

        parent = QtWidgets.QWidget()
        parent.cmd = cmd
        unchanged = cmd.get_color_index("gray40")
        for scope in (target, "palette_selection"):
            for label in ("按链", "按残基名", "按二级结构", "按 B 因子"):
                with self.subTest(scope=scope, label=label):
                    cmd.color("gray40", "all")
                    palette = open_object_palette(parent, scope)
                    button = next(btn for btn in palette._popup.findChildren(QtWidgets.QPushButton)
                                  if btn.text() == label)
                    button.click()
                    self.assertFalse(palette._popup.isVisible())
                    self.assertEqual(colors("palobj"), {"CA": unchanged})
                    result = colors(selection)
                    if scope == "palette_selection":
                        self.assertEqual(result["C5"], unchanged)
                    else:
                        self.assertNotEqual(result["C5"], unchanged)
                    c1, c2, c3, c4 = (result[f"C{i}"] for i in range(1, 5))
                    if label == "按链":
                        self.assertEqual((c1, c2), (c3, c4))
                        self.assertNotEqual(c1, c2)
                    elif label == "按残基名":
                        self.assertEqual((c1, c3), (c2, c4))
                        self.assertNotEqual(c1, c3)
                    elif label == "按二级结构":
                        self.assertEqual(c1, c4)
                        self.assertEqual(len({c1, c2, c3}), 3)
                    else:
                        self.assertEqual(c1, c4)
                        self.assertEqual(len({c1, c2, c3}), 3)
                        low, high = (cmd.get_color_tuple(color) for color in (c1, c3))
                        self.assertGreater(low[2], low[0])
                        self.assertGreater(high[0], high[2])

    def test_residue_legend_tracks_colors_in_the_selected_scope(self) -> None:
        from pymol.Qt import QtWidgets
        from pmg_qt.color_palette import open_object_palette, rgb_to_hex

        cmd = self.cmd
        for i, resn in enumerate(("LIG", "ALA", "LIG", "GLY"), 1):
            cmd.pseudoatom("sample", name=f"C{i}", resn=resn, resi=i, pos=[i, 0, 0])
        cmd.select("legend_scope", "sample and not name C4")
        parent = QtWidgets.QWidget()
        parent.cmd = cmd

        palette = open_object_palette(parent, "legend_scope")
        buttons = {button.text(): button for button in palette._popup.findChildren(QtWidgets.QPushButton)}
        buttons["按残基名"].click()
        legend = palette._residue_legend
        self.assertTrue(legend.isVisible())
        self.assertFalse(palette._popup.isVisible())
        self.assertIn("legend_scope", legend.windowTitle())
        self.assertIn("整个选区", legend.scope_label.text())
        self.assertIn("隐藏或被遮挡", legend.scope_label.text())

        expected = {}
        cmd.iterate("legend_scope", "expected[resn] = color", space={"expected": expected})
        self.assertEqual(legend.table.topLevelItemCount(), 2)
        for i, (resn, count) in enumerate((('ALA', '1'), ('LIG', '2'))):
            row = legend.table.topLevelItem(i)
            self.assertEqual(row.text(0), resn)
            self.assertEqual(row.text(1), rgb_to_hex(cmd.get_color_tuple(expected[resn])))
            self.assertEqual(row.text(2), count)
            self.assertFalse(row.icon(1).isNull())

        legend.close()
        cmd.color("red", "sample and name C1")
        cmd.color("green", "sample and name C3")
        open_object_palette(parent, "legend_scope")
        buttons["残基颜色图例…"].click()
        lig = legend.table.topLevelItem(1)
        self.assertIn("多色", lig.text(1))
        self.assertEqual({lig.child(i).text(1) for i in range(lig.childCount())},
                         {"#FF0000", "#00FF00"})

        cmd.color("yellow", "sample and resn LIG")
        open_object_palette(parent, "palobj")._popup.close()
        legend.refresh()
        self.assertEqual(legend.table.topLevelItem(1).text(1), "#FFFF00")
        self.assertIn("legend_scope", legend.windowTitle())


@unittest.skipUnless(_qt_available(), "需要 pymol.Qt")
class TestVisualRepGui(unittest.TestCase):
    """表示窗口：无效输入、关闭保护、基本动作。"""

    @classmethod
    def setUpClass(cls) -> None:
        from pymol.Qt import QtWidgets

        app = QtWidgets.QApplication.instance()
        if app is None:
            cls.app = QtWidgets.QApplication([])
        else:
            cls.app = app

    def setUp(self) -> None:
        from pymol import cmd
        from pymol.Qt import QtWidgets

        cmd.reinitialize()
        self.cmd = cmd
        atoms = []
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A"]):
            atoms.append(_bead(nm, "DPPC", "L", 1, float(i) * 2, 0, 0))
        for i, nm in enumerate(["BB1", "BB2", "BB3", "SC1"]):
            atoms.append(_bead(nm, "U", "A", 1, 20 + float(i), 0, 0))
        model = Indexed()
        for at in atoms:
            model.add_atom(at)
        cmd.load_model(model, "guiobj", zoom=0)

        self.parent = QtWidgets.QWidget()
        self.parent.cmd = cmd
        self.parent._quit_requested = False
        self.parent._shutting_down = False

    def tearDown(self) -> None:
        from pmg_qt.visual_representations_gui import clear_visual_representations_ref
        from pmg_qt.color_palette import clear_color_palette_refs

        clear_visual_representations_ref()
        clear_color_palette_refs()

    def test_open_and_select_target(self) -> None:
        """打开并选中 target。"""
        from pmg_qt.visual_representations_gui import (
            open_visual_representations,
            clear_visual_representations_ref,
        )

        dlg = open_visual_representations(self.parent, target="guiobj")
        self.assertIsNotNone(dlg)
        self.assertEqual(dlg._current_source(), "guiobj")
        clear_visual_representations_ref()

    def test_invalid_selection_keeps_scene(self) -> None:
        """无效输入友好拒绝且不丢前次显示。"""
        from pmg_qt.visual_representations_gui import open_visual_representations
        from pymol import visual_representations as vr

        dlg = open_visual_representations(self.parent, target="guiobj")
        # 先成功应用一层
        dlg.ctrl.set_layers(
            [
                vr.RepresentationLayer(
                    layer_id=vr.new_layer_id(),
                    style=vr.STYLE_SPHERES,
                    color_mode=vr.COLOR_SOLID,
                    color_rgb=(1, 0, 0),
                    selection="resname DPPC",
                )
            ]
        )
        dlg.apply_layers()
        self.assertTrue(dlg.ctrl.owned)
        owned = list(dlg.ctrl.owned)

        dlg.ctrl.set_layers(
            [
                vr.RepresentationLayer(
                    layer_id=vr.new_layer_id(),
                    style=vr.STYLE_SPHERES,
                    color_mode=vr.COLOR_SOLID,
                    color_rgb=(0, 1, 0),
                    selection="all; cmd.color",
                )
            ]
        )
        dlg.apply_layers()
        self.assertEqual(dlg.ctrl.owned, owned)
        self.assertIn("失败", dlg.status.text())

    def test_close_keeps_display_and_blocks_engine(self) -> None:
        """关闭保留显示；关闭后 engine_ok 为假。"""
        from pmg_qt.visual_representations_gui import open_visual_representations
        from pymol import visual_representations as vr

        dlg = open_visual_representations(self.parent, target="guiobj")
        dlg.ctrl.set_layers(
            [
                vr.RepresentationLayer(
                    layer_id=vr.new_layer_id(),
                    style=vr.STYLE_SPHERES,
                    color_mode=vr.COLOR_BYRESN,
                    selection="nucleic",
                )
            ]
        )
        dlg.apply_layers()
        owned = list(dlg.ctrl.owned)
        dlg.close()
        self.assertTrue(dlg._closed)
        self.assertFalse(dlg._engine_ok())
        # 关闭不清理副本
        for n in owned:
            self.assertIn(n, self.cmd.get_names("all") or [])

    def test_preview_and_component_layers(self) -> None:
        """预览与按组分分层动作。"""
        from pmg_qt.visual_representations_gui import open_visual_representations

        dlg = open_visual_representations(self.parent, target="guiobj")
        dlg.ed_sele.setText("resname DPPC")
        dlg._preview()
        self.assertIn("4", dlg.lbl_preview.text().replace(",", ""))
        dlg._layers_by_component()
        self.assertGreaterEqual(len(dlg.ctrl.layers), 2)

    def test_quit_flag_blocks(self) -> None:
        """父窗口退出标志阻止引擎操作。"""
        from pmg_qt.visual_representations_gui import open_visual_representations

        dlg = open_visual_representations(self.parent, target="guiobj")
        self.parent._quit_requested = True
        self.assertFalse(dlg._engine_ok())


@unittest.skipUnless(_qt_available(), "需要 pymol.Qt")
class TestGuiHooksShape(unittest.TestCase):
    """主 GUI 钩子形状：信号与方法存在（不启动完整 execapp）。"""

    def test_signal_and_methods_declared(self) -> None:
        """PyMOLQtGUI 声明 visual/palette 信号与 open 方法。"""
        from pmg_qt.pymol_qt_gui import PyMOLQtGUI

        self.assertTrue(hasattr(PyMOLQtGUI, "visualsignal"))
        self.assertTrue(hasattr(PyMOLQtGUI, "palettesignal"))
        self.assertTrue(callable(getattr(PyMOLQtGUI, "open_visual_representations")))
        self.assertTrue(callable(getattr(PyMOLQtGUI, "open_color_palette")))


if __name__ == "__main__":
    unittest.main()
