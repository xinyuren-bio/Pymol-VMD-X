"""Per-theme protein/ligand color choices for the figure dialog."""

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.publication_colors import THEME_COLORS, normalize_color


class PublicationColorControls(QtWidgets.QWidget):
    changed = QtCore.Signal()

    def __init__(self, cmd, current_theme, parent=None):
        super().__init__(parent)
        self.cmd = cmd
        self.current_theme = current_theme
        self._profiles = {}
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.protein_button = QtWidgets.QPushButton("蛋白颜色")
        self.ligand_button = QtWidgets.QPushButton("配体碳颜色")
        self.reset_button = QtWidgets.QPushButton("恢复风格配色")
        self.protein_button.clicked.connect(lambda: self.choose("protein"))
        self.ligand_button.clicked.connect(lambda: self.choose("ligand"))
        self.reset_button.clicked.connect(self.reset_current)
        for button in (self.protein_button, self.ligand_button, self.reset_button):
            row.addWidget(button)
        self.ligand_button.setToolTip("改变所选配体的碳原子颜色，保留氧、氮等元素配色。")
        self.refresh()

    def overrides(self, theme=None):
        profile = self._profiles.get(theme or self.current_theme(), {})
        return {"protein_color": profile.get("protein"),
                "ligand_color": profile.get("ligand")}

    def set_color(self, role, color, theme=None, emit=True):
        if role not in ("protein", "ligand"):
            raise ValueError("未知配色目标。")
        theme = theme or self.current_theme()
        rgb = normalize_color(color, self.cmd)
        profile = self._profiles.setdefault(theme, {})
        if rgb is None:
            profile.pop(role, None)
        else:
            profile[role] = rgb
        self.refresh()
        if emit:
            self.changed.emit()

    def _display_color(self, role):
        theme = self.current_theme()
        value = self._profiles.get(theme, {}).get(role)
        if value is None:
            value = THEME_COLORS[theme][role]
        return normalize_color(value, self.cmd)

    def choose(self, role):
        if getattr(self.cmd, "_COb", None) is None:
            return
        rgb = self._display_color(role)
        title = "选择蛋白颜色" if role == "protein" else "选择配体碳颜色"
        color = QtWidgets.QColorDialog.getColor(
            QtGui.QColor.fromRgbF(*rgb), self, title)
        if color.isValid():
            self.set_color(role, color.getRgbF()[:3])

    def reset_current(self):
        self._profiles.pop(self.current_theme(), None)
        self.refresh()
        self.changed.emit()

    def refresh(self):
        if getattr(self.cmd, "_COb", None) is None:
            return
        for role, button in (("protein", self.protein_button),
                             ("ligand", self.ligand_button)):
            color = QtGui.QColor.fromRgbF(*self._display_color(role))
            swatch = QtGui.QPixmap(32, 32)
            swatch.setDevicePixelRatio(2)
            swatch.fill(color)
            button.setIcon(QtGui.QIcon(swatch))
            button.setIconSize(QtCore.QSize(16, 16))
        self.reset_button.setEnabled(bool(self._profiles.get(self.current_theme())))
