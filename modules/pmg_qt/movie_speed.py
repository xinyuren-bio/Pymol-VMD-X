"""Controls for the engine's movie playback rate."""

from bisect import bisect_left, bisect_right
import math

from pymol.Qt import QtGui, QtWidgets


class MovieSpeedControl(QtWidgets.QWidget):
    rates = (0.3, 1, 2, 5, 10, 15, 24, 30, 45, 60)

    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.setObjectName('movie_speed_control')
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(QtWidgets.QLabel('播放速度'))
        layout.addStretch()

        self.slower = QtWidgets.QToolButton()
        self.slower.setText('−')
        self.slower.setToolTip('降低播放速度')
        self.slower.setAccessibleName('降低播放速度')
        self.slower.clicked.connect(lambda: self.step(-1))
        layout.addWidget(self.slower)

        self.current = QtWidgets.QToolButton()
        self.current.setObjectName('movie_speed_current')
        self.current.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.current.setAccessibleName('播放速度')
        self.current.setToolTip('选择目标播放帧率，实际帧率取决于渲染速度。')
        layout.addWidget(self.current)
        menu = QtWidgets.QMenu(self.current)
        self.current.setMenu(menu)
        group = QtGui.QActionGroup(menu)
        self.actions = {}
        for fps in self.rates + (0, -1):
            label = f'{fps:g} 帧/秒' if fps > 0 else ('按帧延迟' if fps == 0 else '不限速')
            action = menu.addAction(label)
            action.setCheckable(True)
            group.addAction(action)
            action.triggered.connect(lambda checked=False, rate=fps: self.set_rate(rate))
            self.actions[fps] = action
        menu.aboutToShow.connect(self.sync)

        self.faster = QtWidgets.QToolButton()
        self.faster.setText('+')
        self.faster.setToolTip('提高播放速度')
        self.faster.setAccessibleName('提高播放速度')
        self.faster.clicked.connect(lambda: self.step(1))
        layout.addWidget(self.faster)
        self.sync()

    def _rate(self):
        configured = self.cmd.get_setting_float('movie_fps')
        if configured == 0:
            delay = self.cmd.get_setting_float('movie_delay')
            effective = 1000 / delay if delay > 0 else math.inf
        else:
            effective = configured if configured > 0 else math.inf
        return configured, effective

    def sync(self, *_):
        configured, effective = self._rate()
        label = f'{effective:.4g} 帧/秒' if math.isfinite(effective) else '不限速'
        self.current.setText(label)
        self.current.setAccessibleDescription(label)
        for rate, action in self.actions.items():
            action.setChecked(math.isclose(configured, rate, rel_tol=1e-6)
                              or (rate == -1 and configured < 0))
        effective = round(effective, 6)
        self.slower.setEnabled(effective > self.rates[0])
        self.faster.setEnabled(effective < self.rates[-1])

    def set_rate(self, fps):
        self.cmd.set('movie_fps', fps)
        self.sync()

    def step(self, direction):
        _, effective = self._rate()
        effective = round(effective, 6)
        index = (bisect_left(self.rates, effective) - 1 if direction < 0
                 else bisect_right(self.rates, effective))
        if 0 <= index < len(self.rates):
            self.set_rate(self.rates[index])
