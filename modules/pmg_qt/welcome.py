"""Welcome page for the empty Qt viewport."""

from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.branding import APP_NAME


Qt = QtCore.Qt


class _MolecularArtwork(QtWidgets.QWidget):
    def __init__(self, path, parent):
        super().__init__(parent)
        self.image = QtGui.QPixmap(path)

    def paintEvent(self, event):
        if self.image.isNull():
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform)
        side = min(self.width(), self.height() - 26)
        x = (self.width() - side) / 2
        y = (self.height() - side - 18) / 2
        painter.drawPixmap(QtCore.QRectF(x, y, side, side), self.image,
                           QtCore.QRectF(self.image.rect()))
        painter.setPen(QtGui.QColor('#a3bdc9'))
        caption = QtCore.QRectF(0, y + side, self.width(), 20)
        painter.drawText(caption, Qt.AlignmentFlag.AlignCenter, 'PDK1 · PDB 1OKY')


class WelcomePanel(QtWidgets.QWidget):
    def __init__(self, gui):
        super().__init__(gui.pymolwidget)
        self.gui = gui
        self.setObjectName('vmdx_welcome')
        self.setAcceptDrops(True)
        self._dismissed = False
        self._viewport_insets = (0, 0)
        self.parentWidget().installEventFilter(self)
        self.setStyleSheet('''
            #vmdx_welcome QLabel { color: #eaf4f4; background: transparent; }
            #vmdx_welcome QLabel#welcome_kicker { color: #8ce1cc; font-weight: 600; }
            #vmdx_welcome QLabel#welcome_hint { color: #a3bdc9; }
            #vmdx_welcome QPushButton {
                background: #203c49; color: #eaf4f4; border: 1px solid #40616c;
                border-radius: 7px; padding: 7px 12px;
            }
            #vmdx_welcome QPushButton:hover, #vmdx_welcome QPushButton:focus {
                background: #2a4b58; border-color: #86e3cc;
            }
            #vmdx_welcome QPushButton#welcome_open {
                background: #8ce1cc; color: #10282d; border-color: #8ce1cc; font-weight: 600;
            }
            #vmdx_welcome QPushButton#welcome_open:hover, #vmdx_welcome QPushButton#welcome_open:focus {
                background: #b4f1e0; border-color: #eaf4f4;
            }
            #vmdx_welcome QPushButton[flat="true"] {
                color: #b0cbd4; background: transparent; border: 1px solid transparent;
                padding: 4px 3px;
            }
            #vmdx_welcome QPushButton[flat="true"]:hover,
            #vmdx_welcome QPushButton[flat="true"]:focus { color: white; border-color: #40616c; }
        ''')

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(30, 15, 24, 15)
        layout.setSpacing(18)
        content = QtWidgets.QVBoxLayout()
        content.setSpacing(10)
        layout.addLayout(content, 5)
        content.addStretch()

        kicker = QtWidgets.QLabel('分子可视化工作台')
        kicker.setObjectName('welcome_kicker')
        content.addWidget(kicker)
        self.title = QtWidgets.QLabel(APP_NAME)
        self.title.setWordWrap(True)
        content.addWidget(self.title)
        caption = QtWidgets.QLabel('PyMOL 内核 · VMD 式选择 · ChimeraX 式工具栏')
        caption.setWordWrap(True)
        content.addWidget(caption)
        content.addSpacing(6)

        actions = QtWidgets.QHBoxLayout()
        open_button = QtWidgets.QPushButton('打开文件…')
        open_button.setObjectName('welcome_open')
        open_button.clicked.connect(gui.file_open)
        actions.addWidget(open_button)
        if gui.open_recent_menu is not None:
            recent_button = QtWidgets.QPushButton('最近文件')
            recent_button.setMenu(gui.open_recent_menu)
            actions.addWidget(recent_button)
        actions.addStretch()
        content.addLayout(actions)

        hint = QtWidgets.QLabel('也可拖入结构或会话文件')
        hint.setObjectName('welcome_hint')
        hint.setWordWrap(True)
        content.addWidget(hint)
        content.addStretch()

        footer = QtWidgets.QHBoxLayout()
        about = QtWidgets.QPushButton('关于与许可', flat=True)
        about.clicked.connect(gui.show_about)
        footer.addWidget(about)
        dismiss = QtWidgets.QPushButton('进入空白视口', flat=True)
        dismiss.clicked.connect(self.dismiss)
        footer.addWidget(dismiss)
        footer.addStretch()
        content.addLayout(footer)

        art_path = gui.cmd.exp_path('$PYMOL_DATA/pymol/vmdx_welcome_structure.png')
        self.artwork = _MolecularArtwork(art_path, self)
        layout.addWidget(self.artwork, 4)
        self.hide()

    def refresh(self):
        """Leave the first loaded scene and the native text view unobstructed."""
        if self._dismissed:
            return
        cmd = self.gui.cmd
        if cmd.get_names('all') or cmd.get_setting_boolean('text'):
            self.dismiss()
            return
        native_scale = cmd.get_setting_int('display_scale_factor')
        scale = max(1, native_scale) / self.parentWidget().devicePixelRatioF()
        panel = cmd.get_setting_int('internal_gui_width') if cmd.get_setting_boolean('internal_gui') else 0
        lines = cmd.get_setting_int('internal_feedback')
        self._viewport_insets = (round(panel * scale), round((18 + 12 * (lines - 1)) * scale) if lines else 0)
        self._fit_viewport()

    def dismiss(self):
        self._dismissed = True
        self.hide()

    def _fit_viewport(self):
        if self._dismissed:
            return
        right, bottom = self._viewport_insets
        w = self.parentWidget().width() - right
        h = self.parentWidget().height() - bottom
        if w < 260 or h < 230:
            self.hide()
            return
        self.setGeometry(0, 0, w, h)
        self.show()
        self.raise_()

    def eventFilter(self, watched, event):
        if watched is self.parentWidget() and event.type() == QtCore.QEvent.Type.Resize:
            self._fit_viewport()
        return False

    def resizeEvent(self, event):
        wide = self.width() >= 560 and self.height() >= 260
        self.artwork.setVisible(wide and not self.artwork.image.isNull())
        font = self.title.font()
        font.setPixelSize(43 if wide else 33)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        self.title.setFont(font)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        background = QtGui.QLinearGradient(0, 0, self.width(), self.height())
        background.setColorAt(0, QtGui.QColor('#0b1927'))
        background.setColorAt(1, QtGui.QColor('#153646'))
        painter.fillRect(self.rect(), background)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            self.parentWidget().dropEvent(event)
            self.refresh()
