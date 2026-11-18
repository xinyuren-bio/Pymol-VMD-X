"""Trajectory measurements and click-to-frame plots."""

from bisect import bisect_left

from pymol import trajectory_analysis
from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.Qt.utils import PopupOnException, getSaveFileNameWithExt


class TracePlot(QtWidgets.QWidget):
    frame_picked = QtCore.Signal(int)

    def __init__(self):
        super().__init__()
        self.setMinimumSize(400, 220)
        self.setAccessibleName("轨迹曲线")
        self.trace = None
        self.selected = -1

    def set_trace(self, trace):
        self.trace = trace
        self.selected = -1
        self.update()

    def axes(self):
        return QtCore.QRectF(66, 26, self.width() - 88, self.height() - 65)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), self.palette().base())
        painter.setPen(self.palette().text().color())
        if self.trace is None:
            painter.drawText(self.rect(), QtCore.Qt.AlignmentFlag.AlignCenter, "计算后显示轨迹曲线")
            return
        plot = self.axes()
        frames, values = self.trace.frames, self.trace.values
        low, high = min(values), max(values)
        if high - low < 1e-8:
            low, high = max(0, low - 0.1), high + 0.1
        span = frames[-1] - frames[0]
        painter.drawLine(plot.topLeft(), plot.bottomLeft())
        painter.drawLine(plot.bottomLeft(), plot.bottomRight())
        painter.drawText(QtCore.QPointF(10, 18), "Å")
        painter.drawText(QtCore.QPointF(8, plot.top() + 5), f"{high:.3g}")
        painter.drawText(QtCore.QPointF(8, plot.bottom()), f"{low:.3g}")
        painter.drawText(QtCore.QPointF(plot.left(), plot.bottom() + 19), str(frames[0]))
        painter.drawText(QtCore.QRectF(plot.right() - 55, plot.bottom() + 2, 55, 22),
                         QtCore.Qt.AlignmentFlag.AlignRight, str(frames[-1]))
        painter.drawText(QtCore.QRectF(plot.left(), plot.bottom() + 2, plot.width(), 25),
                         QtCore.Qt.AlignmentFlag.AlignHCenter, "帧")

        def point(i):
            x = plot.left() + plot.width() * ((frames[i] - frames[0]) / span if span else 0.5)
            y = plot.bottom() - plot.height() * (values[i] - low) / (high - low)
            return QtCore.QPointF(x, y)

        painter.setPen(QtGui.QPen(QtGui.QColor("#3C82C4"), 1.8))
        painter.setBrush(QtGui.QColor("#3C82C4"))
        if len(frames) > 1:
            painter.drawPolyline(QtGui.QPolygonF([point(i) for i in range(len(frames))]))
        else:
            painter.drawEllipse(point(0), 3, 3)
        if self.selected >= 0:
            painter.setPen(QtGui.QPen(QtGui.QColor("#E3874B"), 2))
            painter.setBrush(QtGui.QColor("#E3874B"))
            painter.drawEllipse(point(self.selected), 4, 4)

    def mousePressEvent(self, event):
        plot = self.axes()
        if (self.trace is None or event.button() != QtCore.Qt.MouseButton.LeftButton
                or not plot.contains(event.position())):
            return
        frames = self.trace.frames
        fraction = (event.position().x() - plot.left()) / plot.width()
        frame = frames[0] + fraction * (frames[-1] - frames[0])
        i = min(len(frames) - 1, bisect_left(frames, frame))
        if i > 0 and frame - frames[i - 1] < frames[i] - frame:
            i -= 1
        self.selected = i
        self.update()
        self.frame_picked.emit(i)


class TrajectoryAnalysisDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle("轨迹分析")
        self.cmd = parent.cmd
        root = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        self.object = QtWidgets.QComboBox()
        form.addRow("轨迹对象", self.object)
        self.mode = QtWidgets.QComboBox()
        self.mode.addItem("对齐后 RMSD", "rmsd")
        self.mode.addItem("两原子距离", "distance")
        form.addRow("指标", self.mode)
        self.start, self.end, self.step, self.reference = (QtWidgets.QSpinBox() for _ in range(4))
        row = QtWidgets.QHBoxLayout()
        for label, widget in (("从", self.start), ("到", self.end), ("步长", self.step)):
            widget.setMinimum(1)
            row.addWidget(QtWidgets.QLabel(label))
            row.addWidget(widget)
        form.addRow("帧范围", row)
        self.reference.setMinimum(1)
        form.addRow("参考帧", self.reference)
        self.syntax = QtWidgets.QComboBox()
        self.syntax.addItem("PyMOL", "pymol")
        self.syntax.addItem("VMD", "vmd")
        form.addRow("选择语法", self.syntax)
        self.selection_a = QtWidgets.QLineEdit("polymer.protein and name CA")
        self.selection_b = QtWidgets.QLineEdit()
        self.selection_b.setPlaceholderText("例如 index 20")
        form.addRow("原子选择 A", self.selection_a)
        form.addRow("原子选择 B", self.selection_b)
        root.addLayout(form)
        self.selection_tip = QtWidgets.QLabel()
        self.selection_tip.setWordWrap(True)
        root.addWidget(self.selection_tip)
        self.plot = TracePlot()
        root.addWidget(self.plot, 1)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        actions = QtWidgets.QHBoxLayout()
        for text, callback in (("刷新对象", self.refresh), ("计算", self.calculate)):
            button = QtWidgets.QPushButton(text)
            button.clicked.connect(callback)
            actions.addWidget(button)
        actions.addStretch()
        self.csv = QtWidgets.QPushButton("导出 CSV…")
        self.csv.setEnabled(False)
        self.csv.clicked.connect(self.export)
        actions.addWidget(self.csv)
        root.addLayout(actions)
        self.object.currentTextChanged.connect(self._select_object)
        self.mode.currentIndexChanged.connect(self._select_mode)
        self.syntax.currentIndexChanged.connect(self._select_syntax)
        for widget in (self.start, self.end, self.step, self.reference):
            widget.valueChanged.connect(self._clear_result)
        for widget in (self.selection_a, self.selection_b):
            widget.textChanged.connect(self._clear_result)
        self.plot.frame_picked.connect(self.goto_frame)
        self._select_syntax()
        self._select_mode()
        self.refresh()

    def refresh(self):
        old = self.object.currentText()
        objects = [name for name in self.cmd.get_names("public_objects")
                   if self.cmd.get_type(name) == "object:molecule" and self.cmd.count_states(name) > 1]
        with QtCore.QSignalBlocker(self.object):
            self.object.clear()
            self.object.addItems(objects)
            if old in objects:
                self.object.setCurrentText(old)
        self._select_object(self.object.currentText())

    def _select_object(self, obj):
        count = self.cmd.count_states(obj) if obj else 1
        for widget in (self.start, self.end, self.step, self.reference):
            widget.setMaximum(count)
        self.start.setValue(1)
        self.end.setValue(count)
        self._clear_result()
        if not obj:
            self.status.setText("请先载入有多个坐标状态的轨迹，再刷新。")

    def _select_mode(self, *_):
        rmsd = self.mode.currentData() == "rmsd"
        self.reference.setEnabled(rmsd)
        self.selection_b.setEnabled(not rmsd)
        self._clear_result()

    def _select_syntax(self, *_):
        vmd = self.syntax.currentData() == "vmd"
        if self.selection_a.text().strip() in ("protein and name CA", "polymer.protein and name CA"):
            self.selection_a.setText("protein and name CA" if vmd else "polymer.protein and name CA")
        if vmd:
            self.selection_tip.setText("例如 resname GLY 或 resid < 10 and resid > 5。VMD 的 index 从 0 起算。")
        else:
            self.selection_tip.setText("帧号和 PyMOL 的 index 均从 1 起算。")
        self._clear_result()

    def _clear_result(self, *_):
        self.plot.set_trace(None)
        self.status.clear()
        self.csv.setEnabled(False)

    @PopupOnException.decorator
    def calculate(self):
        if not self.object.currentText():
            self.status.setText("请选择轨迹对象。")
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            trace = trajectory_analysis.analyze(
                self.cmd, self.object.currentText(), self.mode.currentData(),
                self.start.value(), self.end.value(), self.step.value(),
                self.selection_a.text(), self.selection_b.text(), self.reference.value(),
                syntax=self.syntax.currentData(),
            )
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.plot.set_trace(trace)
        self.csv.setEnabled(True)
        self.status.setText(f"{len(trace.frames)} 个采样帧，范围 {min(trace.values):.3f}–"
                            f"{max(trace.values):.3f} Å。点击曲线跳到最近的采样帧。")

    def goto_frame(self, index):
        trace = self.plot.trace
        if trace is None:
            return
        frame = trace.frames[index]
        if trace.object not in self.cmd.get_names("objects") or self.cmd.count_states(trace.object) < frame:
            self._clear_result()
            self.status.setText("轨迹已变化，请刷新后重新计算。")
            return
        self.cmd.mstop()
        self.cmd.frame(frame, trigger=0)
        if self.cmd.get_state() != frame:
            self.cmd.set("state", frame)
        if self.cmd.get_setting_int("state", trace.object) != frame:
            self.cmd.unset("state", trace.object)
        if self.cmd.get_setting_boolean("all_states", trace.object):
            self.cmd.set("all_states", 0, trace.object)
        self.cmd.enable(trace.object)
        self.status.setText(f"第 {frame} 帧 · {trace.values[index]:.4f} Å")

    @PopupOnException.decorator
    def export(self):
        trace = self.plot.trace
        if trace is None:
            return
        path = getSaveFileNameWithExt(self, "导出轨迹 CSV", "", filter="CSV 文件 (*.csv)")
        if path:
            trace.write_csv(path)
            self.status.setText(f"已导出 {len(trace.frames)} 帧到 CSV。")
