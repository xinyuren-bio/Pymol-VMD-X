# ==================================================
# 功能说明：可复用中文非模态色板弹窗（色块 + 任意 RGB），供原生 C 与表示层着色。
# 使用方法：open_object_palette(parent, target='all')；退出清理 clear_color_palette_refs()。
# 依赖环境：PyMOL + PySide6（pymol.Qt）；着色走 truecolor，不覆盖全局颜色名。
# 生成时间：2026-11-14 16:50
# ==================================================

"""中文非模态色板：色块、分组着色与自选颜色。"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Callable, List, Optional, Sequence, Tuple

from pymol.Qt import QtCore, QtGui, QtWidgets

Qt = QtCore.Qt

# 模块级引用：parent 持有弹窗防 GC；清理时置空
_palette_refs: List[QtWidgets.QWidget] = []
_object_palette_ref: Optional["ObjectColorPalette"] = None

# 约 24 个常用色（0..1 RGB）
SWATCH_COLORS: List[Tuple[float, float, float]] = [
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 1.0, 0.0),
    (1.0, 0.0, 1.0),
    (0.0, 1.0, 1.0),
    (1.0, 0.5, 0.0),
    (0.5, 0.0, 1.0),
    (0.80, 0.40, 0.47),
    (0.27, 0.47, 0.67),
    (0.53, 0.34, 0.66),
    (0.72, 0.60, 0.22),
    (0.27, 0.67, 0.60),
    (0.90, 0.62, 0.00),
    (0.60, 0.60, 0.60),
    (0.20, 0.20, 0.20),
    (1.0, 1.0, 1.0),
    (0.85, 0.85, 0.85),
    (0.55, 0.27, 0.07),
    (0.13, 0.55, 0.13),
    (0.00, 0.50, 0.50),
    (0.70, 0.13, 0.13),
    (0.25, 0.41, 0.88),
    (1.0, 0.75, 0.80),
]

_PROPERTY_ACTIONS = (
    ("chain", "按链", "相同链 ID 使用相同颜色。"),
    ("resn", "按残基名", "相同 resname 使用相同颜色，着色后打开图例。"),
    ("ss", "按二级结构", "使用现有的二级结构标记：螺旋珊瑚色，折叠金色，其他灰色。"),
    ("b", "按 B 因子", "按当前范围的 B-factor 着色：低值蓝色，中间白色，高值红色。"),
)


def _parent_alive(parent: Optional[QtWidgets.QWidget]) -> bool:
    """判断父窗口是否仍可用。"""
    if parent is None:
        return False
    try:
        from pymol.Qt import sip
    except Exception:
        try:
            import shiboken6 as sip  # type: ignore
        except Exception:
            sip = None  # type: ignore
    if sip is not None:
        try:
            if sip.isdeleted(parent):
                return False
        except Exception:
            pass
    try:
        parent.objectName()
        return True
    except RuntimeError:
        return False


def _engine_ok(parent: Optional[QtWidgets.QWidget], cmd: Any) -> bool:
    """实例/引擎仍活跃时才允许着色。"""
    if not _parent_alive(parent):
        return False
    if getattr(parent, "_quit_requested", False) or getattr(parent, "_shutting_down", False):
        return False
    if cmd is None or getattr(cmd, "_COb", None) is None:
        return False
    try:
        from pmg_qt.publication_gui import publication_export_running

        if publication_export_running():
            return False
    except Exception:
        pass
    try:
        cmd.get_names("public_objects", enabled_only=0)
        return True
    except Exception:
        return False


def rgb_to_hex(rgb: Sequence[float]) -> str:
    """0..1 RGB → #RRGGBB。"""
    r, g, b = [max(0.0, min(1.0, float(x))) for x in rgb[:3]]
    return "#%02X%02X%02X" % (
        int(round(r * 255)),
        int(round(g * 255)),
        int(round(b * 255)),
    )


def make_swatch_icon(rgb: Sequence[float], size: int = 28) -> QtGui.QIcon:
    """生成色块图标。"""
    color = QtGui.QColor.fromRgbF(
        float(rgb[0]), float(rgb[1]), float(rgb[2])
    )
    pm = QtGui.QPixmap(size, size)
    pm.fill(color)
    painter = QtGui.QPainter(pm)
    painter.setPen(QtGui.QColor(80, 80, 80))
    painter.drawRect(0, 0, size - 1, size - 1)
    painter.end()
    return QtGui.QIcon(pm)


class ColorSwatchButton(QtWidgets.QToolButton):
    """色块按钮：tooltip 显示 RGB/HEX。"""

    colorChosen = QtCore.Signal(float, float, float)

    def __init__(
        self,
        rgb: Tuple[float, float, float],
        parent: Optional[QtWidgets.QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self._rgb = (float(rgb[0]), float(rgb[1]), float(rgb[2]))
        self.setIcon(make_swatch_icon(self._rgb))
        self.setIconSize(QtCore.QSize(22, 22))
        self.setAutoRaise(True)
        hexv = rgb_to_hex(self._rgb)
        self.setToolTip(
            "RGB (%.2f, %.2f, %.2f) / %s"
            % (self._rgb[0], self._rgb[1], self._rgb[2], hexv)
        )
        self.clicked.connect(self._emit)

    def _emit(self) -> None:
        """点击发出颜色。"""
        r, g, b = self._rgb
        self.colorChosen.emit(r, g, b)


class ColorPalettePopup(QtWidgets.QWidget):
    """
    非模态色板弹窗。

    确认更多颜色才修改；取消不修改。可选原色/按元素动作。
    """

    colorSelected = QtCore.Signal(float, float, float)
    actionOriginal = QtCore.Signal()
    actionByElement = QtCore.Signal()
    actionByProperty = QtCore.Signal(str)
    residueLegendRequested = QtCore.Signal()

    def __init__(
        self,
        parent: Optional[QtWidgets.QWidget] = None,
        current_rgb: Optional[Tuple[float, float, float]] = None,
        show_legacy_actions: bool = True,
        *,
        show_property_actions: bool = False,
    ) -> None:
        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.Window)
        self.setWindowTitle("颜色")
        self._current = current_rgb or (0.8, 0.2, 0.2)
        self._closed = False
        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        tip = QtWidgets.QLabel("选择颜色")
        tip.setStyleSheet("color:#444;font-size:11px;")
        root.addWidget(tip)

        cur_row = QtWidgets.QHBoxLayout()
        cur_row.addWidget(QtWidgets.QLabel("当前："))
        self.btn_current = ColorSwatchButton(self._current, self)
        self.btn_current.setEnabled(False)
        cur_row.addWidget(self.btn_current)
        self.lbl_current = QtWidgets.QLabel(rgb_to_hex(self._current))
        cur_row.addWidget(self.lbl_current)
        cur_row.addStretch(1)
        root.addLayout(cur_row)

        grid = QtWidgets.QGridLayout()
        grid.setSpacing(4)
        cols = 6
        for i, rgb in enumerate(SWATCH_COLORS):
            btn = ColorSwatchButton(rgb, self)
            btn.colorChosen.connect(self._on_swatch)
            grid.addWidget(btn, i // cols, i % cols)
        root.addLayout(grid)

        if show_legacy_actions:
            row = QtWidgets.QHBoxLayout()
            self.btn_orig = QtWidgets.QPushButton("原色")
            self.btn_elem = QtWidgets.QPushButton("按元素")
            self.btn_orig.setToolTip("恢复按元素/默认着色（legacy 可选）")
            self.btn_elem.setToolTip("按元素着色（cbag）")
            self.btn_orig.clicked.connect(self._on_original)
            self.btn_elem.clicked.connect(self._on_element)
            row.addWidget(self.btn_orig)
            row.addWidget(self.btn_elem)
            row.addStretch(1)
            root.addLayout(row)

        if show_property_actions:
            actions = QtWidgets.QGridLayout()
            for i, (field, label, tip) in enumerate(_PROPERTY_ACTIONS):
                button = QtWidgets.QPushButton(label)
                button.setToolTip(tip)
                button.clicked.connect(lambda checked=False, field=field: self._on_property(field))
                actions.addWidget(button, i // 2, i % 2)
            root.addLayout(actions)
            legend_button = QtWidgets.QPushButton("残基颜色图例…")
            legend_button.setToolTip("查看当前对象或选区中的残基名与原子颜色。")
            legend_button.clicked.connect(self._on_residue_legend)
            root.addWidget(legend_button)

        self.btn_more = QtWidgets.QPushButton("更多颜色…")
        self.btn_more.clicked.connect(self._on_more)
        root.addWidget(self.btn_more)

    def set_current_rgb(self, rgb: Sequence[float]) -> None:
        """更新当前色显示。"""
        self._current = (float(rgb[0]), float(rgb[1]), float(rgb[2]))
        self.btn_current.setIcon(make_swatch_icon(self._current))
        self.lbl_current.setText(rgb_to_hex(self._current))

    def _on_swatch(self, r: float, g: float, b: float) -> None:
        """色块确认。"""
        self.colorSelected.emit(r, g, b)
        self.close()

    def _on_original(self) -> None:
        """原色动作。"""
        self.actionOriginal.emit()
        self.close()

    def _on_element(self) -> None:
        """按元素动作。"""
        self.actionByElement.emit()
        self.close()

    def _on_property(self, field: str) -> None:
        self.close()
        self.actionByProperty.emit(field)

    def _on_residue_legend(self) -> None:
        self.close()
        self.residueLegendRequested.emit()

    def _on_more(self) -> None:
        """QColorDialog：确认才修改，取消不改。"""
        cur = QtGui.QColor.fromRgbF(*self._current)
        parent = self.parentWidget()
        dialog_parent = parent.window() if parent is not None else None
        # Release the popup's mouse grab before entering a modal color dialog.
        self.close()
        dialog = QtWidgets.QColorDialog(cur, dialog_parent)
        dialog.setWindowTitle("选择任意颜色")
        dialog.setOption(QtWidgets.QColorDialog.ColorDialogOption.DontUseNativeDialog, True)
        dialog.setWindowModality(Qt.WindowModality.WindowModal)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        def accept_color(color):
            if color.isValid():
                r, g, b, _a = color.getRgbF()
                self.colorSelected.emit(float(r), float(g), float(b))
        dialog.colorSelected.connect(accept_color)
        self._color_dialog = dialog
        _palette_refs.append(dialog)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def closeEvent(self, event: Any) -> None:
        """关闭标记。"""
        self._closed = True
        super().closeEvent(event)


class ResidueColorLegend(QtWidgets.QDialog):
    """按残基名汇总指定范围内的原子颜色。"""

    def __init__(self, parent: QtWidgets.QWidget) -> None:
        super().__init__(parent)
        self.resize(390, 380)
        layout = QtWidgets.QVBoxLayout(self)
        self.scope_label = QtWidgets.QLabel()
        self.scope_label.setWordWrap(True)
        layout.addWidget(self.scope_label)
        self.table = QtWidgets.QTreeWidget()
        self.table.setHeaderLabels(["残基名", "原子颜色", "总原子／珠数"])
        self.table.setUniformRowHeights(True)
        self.table.header().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.setColumnWidth(0, 104)
        self.table.setColumnWidth(2, 110)
        layout.addWidget(self.table)

        self.status = QtWidgets.QLabel()
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        buttons = QtWidgets.QHBoxLayout()
        refresh = QtWidgets.QPushButton("刷新")
        close = QtWidgets.QPushButton("关闭")
        refresh.clicked.connect(self.refresh)
        close.clicked.connect(self.close)
        buttons.addWidget(refresh)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout.addLayout(buttons)

    def show_for_target(self, cmd: Any, target: str) -> None:
        self.cmd = cmd
        self.target = target
        self.setWindowTitle(f"残基颜色图例 · {target}")
        self.refresh()
        self.show()
        self.raise_()
        self.activateWindow()

    def refresh(self) -> None:
        if not _engine_ok(self.parentWidget(), self.cmd):
            return
        from pymol.visual_representations import color_target_selection

        self.table.clear()
        counts = defaultdict(Counter)
        try:
            selection = color_target_selection(self.cmd, self.target)
            if self.target == "all":
                scope = "所有对象"
            elif self.target in self.cmd.get_names("selections"):
                scope = "整个选区"
            else:
                scope = "整个对象"
            self.scope_label.setText(f"统计范围：{scope}，包含隐藏或被遮挡的原子／珠。")
            self.cmd.iterate(selection, "counts[resn][color] += 1", space={"counts": counts})
        except Exception as exc:
            self.status.setText(f"读取颜色失败：{exc}")
            self.status.show()
            return

        color_labels = {}
        icons = {}
        for resn, atom_counts in sorted(counts.items(), key=lambda item: item[0].casefold()):
            colors = Counter()
            for index, count in atom_counts.items():
                if index not in color_labels:
                    rgb = self.cmd.get_color_tuple(index) if index >= 0 else None
                    label = rgb_to_hex(rgb) if rgb is not None else "动态颜色"
                    color_labels[index] = label
                    if rgb is not None and label not in icons:
                        icons[label] = make_swatch_icon(rgb, size=16)
                colors[color_labels[index]] += count

            row = QtWidgets.QTreeWidgetItem([resn or "（未命名）", "", f"{sum(colors.values()):,}"])
            self.table.addTopLevelItem(row)
            if len(colors) == 1:
                label = next(iter(colors))
                row.setText(1, label)
                if label in icons:
                    row.setIcon(1, icons[label])
            else:
                row.setText(1, f"多色（{len(colors)} 种）")
                for label, count in colors.most_common():
                    child = QtWidgets.QTreeWidgetItem(row, ["", label, f"{count:,}"])
                    if label in icons:
                        child.setIcon(1, icons[label])

        self.status.setText("当前范围内没有原子。" if not counts else "")
        self.status.setVisible(not counts)


class ObjectColorPalette(QtWidgets.QWidget):
    """
    面向原生 C 按钮的对象/选择着色入口（非模态）。

    实际弹出 ColorPalettePopup；parent 持有引用。
    """

    def __init__(self, parent: QtWidgets.QWidget, cmd: Any, target: str = "all") -> None:
        super().__init__(parent)
        self.setVisible(False)
        self.cmd = cmd
        self.target = target or "all"
        self._popup: Optional[ColorPalettePopup] = None
        self._residue_legend: Optional[ResidueColorLegend] = None
        self._closed = False

    def set_target(self, target: str) -> None:
        """更新着色目标。"""
        self.target = target or "all"

    def show_palette(self) -> None:
        """弹出色板。"""
        if self._closed or not _engine_ok(self.parent(), self.cmd):
            return
        if self._popup is None:
            self._popup = ColorPalettePopup(self, show_property_actions=True)
            # Original colors are a layer mode.
            self._popup.btn_orig.hide()
            self._popup.colorSelected.connect(self._apply_rgb)
            self._popup.actionOriginal.connect(self._apply_original)
            self._popup.actionByElement.connect(self._apply_element)
            self._popup.actionByProperty.connect(self._apply_property)
            self._popup.residueLegendRequested.connect(self._show_residue_legend)
            _palette_refs.append(self._popup)
        try:
            from pymol.visual_representations import color_target_selection
            colors = []
            self.cmd.iterate("first (%s)" % color_target_selection(self.cmd, self.target),
                             "colors.append(color)", space={"colors": colors})
            rgb = self.cmd.get_color_tuple(colors[0]) if colors else (0.7, 0.7, 0.7)
            self._popup.set_current_rgb(rgb)
        except Exception:
            self._popup.set_current_rgb((0.7, 0.7, 0.7))
        pos = QtGui.QCursor.pos()
        screen = QtGui.QGuiApplication.screenAt(pos) or self.parentWidget().screen()
        if screen is not None:
            bounds = screen.availableGeometry().adjusted(4, 4, -4, -4)
            self._popup.adjustSize()
            size = self._popup.frameSize()
            pos.setX(max(bounds.left(), min(pos.x(), bounds.right() - size.width() + 1)))
            pos.setY(max(bounds.top(), min(pos.y(), bounds.bottom() - size.height() + 1)))
        self._popup.move(pos)
        self._popup.show()
        self._popup.raise_()

    def _apply_rgb(self, r: float, g: float, b: float) -> None:
        """truecolor 着色。"""
        if not _engine_ok(self.parent(), self.cmd):
            return
        try:
            from pymol.visual_representations import color_selection_truecolor

            color_selection_truecolor(self.cmd, self.target, (r, g, b))
        except Exception as exc:
            print("着色失败：%s" % exc)

    def _show_residue_legend(self) -> None:
        if not _engine_ok(self.parent(), self.cmd):
            return
        if self._residue_legend is None:
            self._residue_legend = ResidueColorLegend(self.parentWidget())
            _palette_refs.append(self._residue_legend)
        self._residue_legend.show_for_target(self.cmd, self.target)

    def _apply_original(self) -> None:
        """原色：按元素回退。"""
        self._apply_element()

    def _apply_element(self) -> None:
        """按元素着色。"""
        if not _engine_ok(self.parent(), self.cmd):
            return
        try:
            self.cmd.util.cbag(self.target)
        except Exception:
            try:
                self.cmd.color("atomic", self.target)
            except Exception as exc:
                print("按元素着色失败：%s" % exc)

    def _apply_property(self, field: str) -> None:
        if not _engine_ok(self.parent(), self.cmd):
            return
        from pymol import util
        from pymol.visual_representations import color_target_selection

        try:
            selection = color_target_selection(self.cmd, self.target)
            if field == "ss":
                util.cbss(selection, helix_color="0xE76F51", sheet_color="0xE9C46A",
                          loop_color="0x9099A4", _self=self.cmd)
            else:
                palette = "blue_white_red" if field == "b" else "rainbow"
                self.cmd.spectrum(field, palette, selection)
        except Exception as exc:
            print("着色失败：%s" % exc)
            return
        if field == "resn":
            self._show_residue_legend()


def clear_color_palette_refs() -> None:
    """清除模块级色板引用（应用退出清理）。"""
    global _object_palette_ref, _palette_refs
    _object_palette_ref = None
    for w in list(_palette_refs):
        try:
            w.close()
        except Exception:
            pass
    _palette_refs = []


def open_color_palette_popup(
    parent: QtWidgets.QWidget,
    current_rgb: Optional[Tuple[float, float, float]] = None,
    on_color: Optional[Callable[[float, float, float], None]] = None,
    on_original: Optional[Callable[[], None]] = None,
    on_element: Optional[Callable[[], None]] = None,
    show_legacy_actions: bool = True,
) -> ColorPalettePopup:
    """
    打开可复用色板弹窗（表示层编辑颜色等）。

    parent 必须持有返回值以防 GC。
    """
    pop = ColorPalettePopup(
        parent, current_rgb=current_rgb, show_legacy_actions=show_legacy_actions
    )
    if on_color is not None:
        pop.colorSelected.connect(lambda r, g, b: on_color(r, g, b))
    if on_original is not None:
        pop.actionOriginal.connect(on_original)
    if on_element is not None:
        pop.actionByElement.connect(on_element)
    _palette_refs.append(pop)
    pos = QtGui.QCursor.pos()
    pop.move(pos)
    pop.show()
    pop.raise_()
    return pop


def open_object_palette(
    parent: QtWidgets.QWidget, target: str = "all"
) -> ObjectColorPalette:
    """打开（或激活）原生 C 对象色板；固定入口供主 GUI / 测试调用。"""
    global _object_palette_ref
    cmd = getattr(parent, "cmd", None)
    if cmd is None:
        from pymol import cmd as cmd  # noqa: N813
    if _object_palette_ref is None or _object_palette_ref.parent() is not parent:
        _object_palette_ref = ObjectColorPalette(parent, cmd, target=target)
    else:
        _object_palette_ref._closed = False
        _object_palette_ref.set_target(target)
        _object_palette_ref.cmd = cmd
    _object_palette_ref.show_palette()
    return _object_palette_ref
