"""Independent protein and ligand color controls for publication figures."""

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.publication_colors import (
    ELEMENT_CARBON_COLOR, LIGAND_CARBON_COLORS, chain_color, normalize_color,
)


class PublicationColorControls(QtWidgets.QWidget):
    changed = QtCore.Signal()

    def __init__(self, cmd, current_theme, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.current_theme = current_theme
        self._colors = {}
        self._chain_colors = {}
        self._chains = []
        layout = QtWidgets.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.protein_mode = QtWidgets.QComboBox()
        self.protein_mode.addItem("按链", "chain")
        self.protein_mode.addItem("单色", "single")
        self.chain_combo = QtWidgets.QComboBox()
        self.chain_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.chain_combo.setMinimumContentsLength(8)
        self.chain_combo.setToolTip("选择要单独调整颜色的蛋白链。")
        self.protein_button = QtWidgets.QPushButton("选色…")
        layout.addWidget(QtWidgets.QLabel("蛋白"), 0, 0)
        layout.addWidget(self.protein_mode, 0, 1)
        layout.addWidget(self.chain_combo, 0, 2)
        layout.addWidget(self.protein_button, 0, 3)

        self.ligand_mode = QtWidgets.QComboBox()
        self.ligand_mode.addItem("按元素", "element")
        self.ligand_mode.addItem("碳色可调", "carbon")
        self.ligand_button = QtWidgets.QPushButton("配体碳色…")
        self.ligand_button.setToolTip("在“碳色可调”模式下调整碳原子颜色，其他原子保留元素色。")
        self.reset_button = QtWidgets.QPushButton("恢复默认配色")
        layout.addWidget(QtWidgets.QLabel("配体"), 1, 0)
        layout.addWidget(self.ligand_mode, 1, 1)
        layout.addWidget(self.ligand_button, 1, 2)
        layout.addWidget(self.reset_button, 1, 3)
        layout.setColumnStretch(2, 1)

        self.protein_mode.currentIndexChanged.connect(self._mode_changed)
        self.ligand_mode.currentIndexChanged.connect(self._mode_changed)
        self.chain_combo.currentIndexChanged.connect(self.refresh)
        self.protein_button.clicked.connect(lambda: self.choose("protein"))
        self.ligand_button.clicked.connect(lambda: self.choose("ligand"))
        self.reset_button.clicked.connect(self.reset_current)
        self.refresh()

    def set_chains(self, chains):
        selected = self.chain_combo.currentData()
        self._chains = list(chains)
        with QtCore.QSignalBlocker(self.chain_combo):
            self.chain_combo.clear()
            for chain in self._chains:
                self.chain_combo.addItem(chain.label, chain.key)
            index = next((i for i, chain in enumerate(self._chains) if chain.key == selected), 0)
            self.chain_combo.setCurrentIndex(index)
        self.refresh()

    def overrides(self):
        return {"protein_color": self._colors.get("protein"), "ligand_color": self._colors.get("ligand"),
                "protein_mode": self.protein_mode.currentData(),
                "ligand_mode": self.ligand_mode.currentData(),
                "chain_colors": dict(self._chain_colors)}

    def load_overrides(self, colors):
        self._colors = {
            role: colors[f"{role}_color"] for role in ("protein", "ligand")
            if colors[f"{role}_color"] is not None
        }
        self._chain_colors = dict(colors["chain_colors"])
        for combo, mode in ((self.protein_mode, colors["protein_mode"]),
                            (self.ligand_mode, colors["ligand_mode"])):
            with QtCore.QSignalBlocker(combo):
                combo.setCurrentIndex(combo.findData(mode))
        self.refresh()

    def set_color(self, role, color, emit=True):
        if role not in ("protein", "ligand"):
            raise ValueError("未知配色目标。")
        rgb = normalize_color(color, self.cmd)
        if rgb is None:
            self._colors.pop(role, None)
        else:
            self._colors[role] = rgb
            combo = self.protein_mode if role == "protein" else self.ligand_mode
            with QtCore.QSignalBlocker(combo):
                combo.setCurrentIndex(1)
        self.refresh()
        if emit:
            self.changed.emit()

    def set_chain_color(self, key, value):
        rgb = normalize_color(value, self.cmd)
        if rgb is None:
            self._chain_colors.pop(key, None)
        else:
            self._chain_colors[key] = rgb
        with QtCore.QSignalBlocker(self.protein_mode):
            self.protein_mode.setCurrentIndex(0)
        self._mode_changed()

    def _mode_changed(self):
        self.refresh()
        self.changed.emit()

    def _display_color(self, role):
        theme = self.current_theme()
        if role == "protein" and self.protein_mode.currentData() == "chain":
            key = self.chain_combo.currentData()
            return self._chain_colors.get(key, chain_color(theme, max(0, self.chain_combo.currentIndex())))
        if role == "ligand" and self.ligand_mode.currentData() == "element":
            return normalize_color(ELEMENT_CARBON_COLOR, self.cmd)
        value = self._colors.get(role)
        if value is None:
            value = chain_color(theme, 0) if role == "protein" else LIGAND_CARBON_COLORS[theme]
        return normalize_color(value, self.cmd)

    def choose(self, role):
        if getattr(self.cmd, "_COb", None) is None:
            return
        if role == "ligand" and self.ligand_mode.currentData() == "element":
            return
        by_chain = role == "protein" and self.protein_mode.currentData() == "chain"
        if by_chain and self.chain_combo.currentIndex() < 0:
            return
        title = self.chain_combo.currentText() if by_chain else (
            "选择蛋白颜色" if role == "protein" else "选择配体碳颜色")
        color = QtWidgets.QColorDialog.getColor(QtGui.QColor.fromRgbF(*self._display_color(role)), self, title)
        if color.isValid():
            if by_chain:
                self.set_chain_color(self.chain_combo.currentData(), color.getRgbF()[:3])
            else:
                self.set_color(role, color.getRgbF()[:3])

    def reset_current(self):
        self._colors.clear()
        self._chain_colors.clear()
        for combo in (self.protein_mode, self.ligand_mode):
            with QtCore.QSignalBlocker(combo):
                combo.setCurrentIndex(0)
        self._mode_changed()

    def refresh(self):
        if getattr(self.cmd, "_COb", None) is None:
            return
        by_chain = self.protein_mode.currentData() == "chain"
        self.chain_combo.setVisible(by_chain)
        self.protein_button.setEnabled(bool(self._chains) or not by_chain)
        self.ligand_button.setEnabled(self.ligand_mode.currentData() == "carbon")
        for role, button in (("protein", self.protein_button), ("ligand", self.ligand_button)):
            swatch = QtGui.QPixmap(32, 32)
            swatch.setDevicePixelRatio(2)
            swatch.fill(QtGui.QColor.fromRgbF(*self._display_color(role)))
            button.setIcon(QtGui.QIcon(swatch))
            button.setIconSize(QtCore.QSize(16, 16))
        self.reset_button.setEnabled(bool(
            self._colors or self._chain_colors
            or self.protein_mode.currentIndex() or self.ligand_mode.currentIndex()))
