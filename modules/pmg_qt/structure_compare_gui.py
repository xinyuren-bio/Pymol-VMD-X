"""Desktop entry point for structural overlays and paired figures."""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from pymol import structure_compare
from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.Qt.utils import PopupOnException, getSaveFileNameWithExt


class StructureCompareDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle("结构对照")
        self.cmd = parent.cmd
        self.comparison = None
        self.labels = ()
        root = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        self.reference = self._inputs(form, "参考结构")
        self.mobile = self._inputs(form, "比较结构")
        self.syntax = QtWidgets.QComboBox()
        self.syntax.addItem("PyMOL", "pymol")
        self.syntax.addItem("VMD", "vmd")
        form.addRow("选择语法", self.syntax)
        root.addLayout(form)
        help_label = QtWidgets.QLabel("范围留空时全选蛋白。用匹配的蛋白 Cα 对齐。")
        help_label.setWordWrap(True)
        root.addWidget(help_label)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        actions = QtWidgets.QHBoxLayout()
        for text, callback in (
            ("刷新对象", self.refresh), ("对齐并叠合", self.overlay),
            ("同尺度并排图…", self.preview_pair), ("清除对照", self.clear),
        ):
            button = QtWidgets.QPushButton(text)
            button.clicked.connect(callback)
            actions.addWidget(button)
        root.addLayout(actions)
        self.syntax.currentIndexChanged.connect(self._select_syntax)
        self._select_syntax()
        self.refresh()

    def _inputs(self, form, label):
        row = QtWidgets.QHBoxLayout()
        objects = QtWidgets.QComboBox()
        states = QtWidgets.QSpinBox()
        states.setRange(1, 1)
        states.setPrefix("状态 ")
        scope = QtWidgets.QLineEdit()
        scope.setAccessibleName(f"{label}的选择范围")
        row.addWidget(objects)
        row.addWidget(states)
        row.addWidget(scope, 1)
        form.addRow(label, row)
        objects.currentTextChanged.connect(
            lambda obj: states.setMaximum(max(1, self.cmd.count_states(obj)) if obj else 1))
        return objects, states, scope

    def _select_syntax(self, *_):
        example = ("例如 protein and resid 5 to 10" if self.syntax.currentData() == "vmd"
                   else "例如 polymer.protein and resi 5-10")
        for _, _, scope in (self.reference, self.mobile):
            scope.setPlaceholderText(example)

    def refresh(self):
        exclude = set(self.comparison.objects if self.comparison else ())
        names = [obj for obj in self.cmd.get_names("public_objects")
                 if obj not in exclude and self.cmd.get_type(obj) == "object:molecule"]
        for i, (objects, states, _) in enumerate((self.reference, self.mobile)):
            old, state = objects.currentText(), states.value()
            objects.clear()
            objects.addItems(names)
            if old in names:
                objects.setCurrentText(old)
                states.setValue(state)
            elif names:
                objects.setCurrentIndex(min(i, len(names) - 1))

    @PopupOnException.decorator
    def overlay(self):
        ref, ref_state, ref_scope = self.reference
        mob, mob_state, mob_scope = self.mobile
        if not ref.currentText() or not mob.currentText():
            self.status.setText("请选择两个结构或一个对象的两个坐标状态。")
            return
        result = structure_compare.compare(
            self.cmd, ref.currentText(), mob.currentText(), ref_state.value(), mob_state.value(),
            ref_scope.text(), mob_scope.text(),
            syntax=self.syntax.currentData(),
        )
        if self.comparison:
            sources = (ref.currentText(), mob.currentText())
            was_visible = set(self.comparison.visible_sources) | set(result.visible_sources)
            visible = tuple(dict.fromkeys(source for source in sources if source in was_visible))
            result = replace(result, visible_sources=visible)
            self.clear()
        self.comparison = result
        structure_compare.show_overlay(self.cmd, self.comparison)
        self.labels = (f"{ref.currentText()} · 状态 {ref_state.value()}",
                       f"{mob.currentText()} · 状态 {mob_state.value()}")
        self.status.setText(f"匹配 {self.comparison.matched} 个 Cα 残基，"
                            f"对齐后 RMSD {self.comparison.rmsd:.3f} Å。参考为蓝色，比较为橙色。")

    @PopupOnException.decorator
    def preview_pair(self):
        if not self.comparison:
            self.status.setText("请先点击「对齐并叠合」。")
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            with TemporaryDirectory(prefix="pymolx-compare-") as folder:
                left, right = (Path(folder) / name for name in ("reference.png", "comparison.png"))
                structure_compare.render_pair(self.cmd, self.comparison, left, right)
                panels = QtGui.QImage(str(left)), QtGui.QImage(str(right))
                if any(panel.isNull() for panel in panels):
                    raise RuntimeError("无法读取并排图。")
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

        width, height = panels[0].width(), panels[0].height()
        gap, heading = 20, 40
        image = QtGui.QImage(2 * width + gap, height + heading, QtGui.QImage.Format.Format_RGB32)
        image.fill(QtGui.QColor("white"))
        painter = QtGui.QPainter(image)
        try:
            for i, (title, panel) in enumerate(zip(self.labels, panels)):
                x = i * (width + gap)
                painter.drawText(x + 16, 26, title)
                painter.drawImage(x, heading, panel)
        finally:
            painter.end()

        dialog = QtWidgets.QDialog(self)
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.setWindowTitle("同尺度并排图")
        layout = QtWidgets.QVBoxLayout(dialog)
        view = QtWidgets.QLabel()
        view.setPixmap(QtGui.QPixmap.fromImage(image).scaledToWidth(
            800, QtCore.Qt.TransformationMode.SmoothTransformation))
        layout.addWidget(view)
        save = QtWidgets.QPushButton("导出并排 PNG…")
        layout.addWidget(save)

        @PopupOnException.decorator
        def export(_checked=False):
            path = getSaveFileNameWithExt(dialog, "导出并排 PNG", "", filter="PNG 图像 (*.png)")
            if path and not image.save(path, "PNG"):
                raise OSError(f"无法写入 {path}")

        save.clicked.connect(export)
        dialog.show()

    def clear(self):
        if self.comparison:
            structure_compare.clear_comparison(self.cmd, self.comparison)
            self.comparison = None
            self.status.clear()
