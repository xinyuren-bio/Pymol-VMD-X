"""
Contains main class for PyMOL QT GUI
"""


from collections import defaultdict
import os
import re
import sys

import pymol
import pymol._gui
from pymol import colorprinting, save_shortcut
from pymol.branding import APP_NAME, COMMAND_PROMPT, local_build_info
from pymol.zh_ui import zh

from pymol.Qt import QtGui, QtCore, QtWidgets
from pymol.Qt.utils import (getSaveFileNameWithExt, UpdateLock, WidgetMenu,
        MainThreadCaller,
        PopupOnException,
        connectFontContextMenu, getMonospaceFont)

from .pymol_gl_widget import PyMOLGLWidget
from . import keymapping

from pmg_qt import properties_dialog, file_dialogs

Qt = QtCore.Qt
QFileDialog = QtWidgets.QFileDialog
getOpenFileNames = QFileDialog.getOpenFileNames


class PyMOLQtGUI(QtWidgets.QMainWindow, pymol._gui.PyMOLDesktopGUI):
    '''
    PyMOL QMainWindow GUI
    '''

    from pmg_qt.file_dialogs import (
            load_dialog,
            load_mae_dialog,
            file_fetch_pdb,
            file_save_png,
            file_save_mpeg,
            file_save_map,
            file_save_aln,
            file_save
    )

    _ext_window_visible = True
    _initialdir = ''

    # 任意线程可请求退出；Queued 投递到 GUI 线程，避免 closeEvent 内直接 PExit
    quitsignal = QtCore.Signal(int)
    # for thread-safe viewport command
    viewportsignal = QtCore.Signal(int, int)
    # 原生 V/C：Queued 投递，避免在 C++/OpenGL/菜单栈内直接建窗
    visualsignal = QtCore.Signal(str)
    palettesignal = QtCore.Signal(str)

    def keyPressEvent(self, ev):
        args = keymapping.keyPressEventToPyMOLButtonArgs(ev)

        if args is not None:
            self.pymolwidget.pymol.button(*args)

    def closeEvent(self, event):
        # 关闭窗口：接受销毁仅在正式关机阶段；平时忽略并排队优雅退出
        if getattr(self, '_shutting_down', False):
            event.accept()
            return
        event.ignore()
        self.request_quit(0)

    def pymolviewport(self, w, h):
        cw, ch = self.cmd.get_viewport()
        pw = self.pymolwidget
        scale = pw.fb_scale

        # maintain aspect ratio
        if h < 1:
            if w < 1:
                pw.pymol.reshape(int(scale * pw.width()),
                                 int(scale * pw.height()), True)
                return
            h = (w * ch) / cw
        if w < 1:
            w = (h * cw) / ch

        win_size = self.size()
        delta = QtCore.QSize(w - cw, h - ch) / scale

        # window resize
        self.resize(delta + win_size)

    def get_view(self):
        self.cmd.get_view(2, quiet=0)
        QtWidgets.QApplication.clipboard().setText(self.cmd.get_view(3))
        print(" get_view: matrix copied to clipboard.")

    def __init__(self):  # noqa
        QtWidgets.QMainWindow.__init__(self)
        self.setDockOptions(QtWidgets.QMainWindow.DockOption.AllowTabbedDocks |
                            QtWidgets.QMainWindow.DockOption.AllowNestedDocks)

        # 退出状态：幂等请求 + 关机阶段标志
        self._quit_requested = False
        self._shutting_down = False
        self._shutdown_done = False
        self._exit_code = 0
        self.quitsignal.connect(
            self._handle_quit_request,
            Qt.ConnectionType.QueuedConnection)

        # resize Window before it is shown
        options = pymol.invocation.options
        self.resize(
            options.win_x + (220 if options.internal_gui else 0),
            options.win_y + (246 if options.external_gui else 18))

        # for thread-safe viewport command
        self.viewportsignal.connect(self.pymolviewport)
        self.visualsignal.connect(
            self.open_visual_representations,
            Qt.ConnectionType.QueuedConnection)
        self.palettesignal.connect(
            self.open_color_palette,
            Qt.ConnectionType.QueuedConnection)

        # reusable dialogs
        self.dialog_png = None
        self.advanced_settings_dialog = None
        self.props_panel = None
        self.builder = None
        self.shortcut_menu_filter_dialog = None
        self.scene_panel_dialog = None

        # setting index -> callable
        self.setting_callbacks = defaultdict(list)

        # "session_file" setting in window title
        self.setting_callbacks[440].append(
            lambda v: self.setWindowTitle(f'{APP_NAME} ({os.path.basename(v)})' if v else APP_NAME)
        )

        # "External" Command Line and Loggin Widget
        self._setup_history()
        self.lineedit = CommandLineEdit()
        self.lineedit.setObjectName("command_line")
        self.browser = QtWidgets.QPlainTextEdit()
        self.browser.setObjectName("feedback_browser")
        self.browser.setReadOnly(True)

        # convenience: clicking into feedback browser gives focus to command
        # line. Drawback: Copying with CTRL+C doesn't work in feedback
        # browser -> clear focus proxy while text selected
        self.browser.setFocusProxy(self.lineedit)

        @self.browser.copyAvailable.connect
        def _(yes):
            self.browser.setFocusProxy(None if yes else self.lineedit)
            self.browser.setFocus()

        # Font
        self.browser.setFont(getMonospaceFont())
        connectFontContextMenu(self.browser)

        lineeditlayout = QtWidgets.QHBoxLayout()
        command_label = QtWidgets.QLabel(COMMAND_PROMPT)
        command_label.setObjectName("command_label")
        lineeditlayout.addWidget(command_label)
        lineeditlayout.addWidget(self.lineedit)
        # 命令栏右侧紧凑 V：打开选择/表示；不改变命令行解释与历史
        self.btn_visual_rep = QtWidgets.QToolButton()
        self.btn_visual_rep.setText("V")
        self.btn_visual_rep.setToolTip("选择/显示（V）")
        self.btn_visual_rep.setAutoRaise(True)
        self.btn_visual_rep.setFixedWidth(28)
        self.btn_visual_rep.clicked.connect(
            lambda: self.open_visual_representations("all"))
        lineeditlayout.addWidget(self.btn_visual_rep)
        self.lineedit.setToolTip(f'''Command Input Area

Get the list of commands by hitting <TAB>

Get the list of arguments for one command with a question mark:
{COMMAND_PROMPT} color ?

Read the online help for a command with "help":
{COMMAND_PROMPT} help color

Get autocompletion for many arguments by hitting <TAB>
{COMMAND_PROMPT} color ye<TAB>    (will autocomplete "yellow")
''')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.browser)
        layout.addLayout(lineeditlayout)

        quickbuttonslayout = QtWidgets.QVBoxLayout()
        quickbuttonslayout.setSpacing(2)

        extguilayout = QtWidgets.QBoxLayout(QtWidgets.QBoxLayout.Direction.LeftToRight)
        extguilayout.setContentsMargins(2, 2, 2, 2)
        extguilayout.addLayout(layout)
        extguilayout.addLayout(quickbuttonslayout)

        class ExtGuiFrame(QtWidgets.QFrame):
            def mouseDoubleClickEvent(_, event):
                self.toggle_ext_window_dockable(True)

            _size_hint = QtCore.QSize(options.win_x, options.ext_y)

            def sizeHint(self):
                return self._size_hint

        dockWidgetContents = ExtGuiFrame(self)
        dockWidgetContents.setLayout(extguilayout)
        dockWidgetContents.setObjectName("extgui")

        self.ext_window = \
            dockWidget = QtWidgets.QDockWidget(self)
        dockWidget.setWindowTitle(zh("External GUI"))
        dockWidget.setWidget(dockWidgetContents)
        if options.external_gui:
            dockWidget.setTitleBarWidget(QtWidgets.QWidget())
        else:
            dockWidget.hide()

        self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, dockWidget)

        # rearrange vertically if docking left or right
        @dockWidget.dockLocationChanged.connect
        def _(area):
            if area == Qt.DockWidgetArea.LeftDockWidgetArea or area == Qt.DockWidgetArea.RightDockWidgetArea:
                extguilayout.setDirection(QtWidgets.QBoxLayout.Direction.BottomToTop)
                quickbuttonslayout.takeAt(quickbuttons_stretch_index)
            else:
                extguilayout.setDirection(QtWidgets.QBoxLayout.Direction.LeftToRight)
                if quickbuttons_stretch_index >= quickbuttonslayout.count():
                    quickbuttonslayout.addStretch()

        # OpenGL Widget
        show_welcome = options.show_splash
        native_splash = '_do__ cmd.splash(1)'
        if options.deferred[:1] == [native_splash]:
            options.deferred.pop(0)
        # The Qt welcome page takes the place of the native startup banner.
        options.show_splash = 0
        try:
            self.pymolwidget = PyMOLGLWidget(self)
        finally:
            options.show_splash = show_welcome
        self.setCentralWidget(self.pymolwidget)

        cmd = self.cmd = self.pymolwidget.cmd

        '''
        # command completion
        completer = QtWidgets.QCompleter(cmd.kwhash.keywords, self)
        self.lineedit.setCompleter(completer)
        '''

        # overload <Tab> action
        self.lineedit.installEventFilter(self)
        self.pymolwidget.installEventFilter(self)

        # Quick Buttons（文案经 zh()；符号播放键保留符号并配中文提示）
        _qb_symbols = {'|<', '<', '>', '>|'}
        _qb_tips = {
            '|<': '回到开头',
            '<': '上一帧',
            'Stop': '暂停播放，保留当前帧',
            'Play': '开始播放',
            '>': '下一帧',
            '>|': '跳到结尾',
            'MClear': '清除动画帧',
            'Reset': '重置视角',
            'Zoom': '缩放至对象',
            'Orient': '定向对齐',
            'Draw/Ray': '绘制 / 光线追踪',
            'Unpick': '取消原子拾取',
            'Deselect': '取消选择',
            'Rock': '摇摆动画',
            'Get View': '复制当前视角矩阵',
            'Builder': '打开分子构建器',
            'Properties': '对象属性',
            'Rebuild': '重建几何',
        }
        for row in [
            [
                ('Reset', cmd.reset),
                ('Zoom', lambda: cmd.zoom(animate=1.0)),
                ('Orient', lambda: cmd.orient(animate=1.0)),

                # render dialog will be constructed when the menu is shown
                # for the first time. This way it's populated with the current
                # viewport and settings. Also defers parsing of the ui file.
                ('Draw/Ray', WidgetMenu(self).setSetupUi(self.render_dialog)),
            ],
            [
                ('Unpick', cmd.unpick),
                ('Deselect', cmd.deselect),
                ('Rock', cmd.rock),
                ('Get View', self.get_view),
            ],
            [
                ('|<', cmd.rewind),
                ('<', cmd.backward),
                ('Stop', cmd.mstop),
                ('Play', cmd.mplay),
                ('>', cmd.forward),
                ('>|', cmd.ending),
                ('MClear', cmd.mclear),
            ],
            [
                ('Builder', self.open_builder_panel),
                ('Properties', self.open_props_dialog),
                ('Rebuild', cmd.rebuild),
            ],
        ]:
            hbox = QtWidgets.QHBoxLayout()
            hbox.setSpacing(2)

            for name, callback in row:
                label = name if name in _qb_symbols else zh(name)
                btn = QtWidgets.QPushButton('暂停' if name == 'Stop' else label)
                if name in ('Stop', 'Play'):
                    btn.setObjectName('movie_pause' if name == 'Stop' else 'movie_play')
                btn.setProperty("quickbutton", True)
                btn.setAttribute(Qt.WidgetAttribute.WA_LayoutUsesWidgetRect) # OS X workaround
                if name in _qb_tips:
                    btn.setToolTip(_qb_tips[name])
                hbox.addWidget(btn)

                if callback is None:
                    btn.setEnabled(False)
                elif isinstance(callback, QtWidgets.QMenu):
                    btn.setMenu(callback)
                else:
                    btn.released.connect(callback)

            quickbuttonslayout.addLayout(hbox)
            if row[0][0] == '|<':
                from .movie_speed import MovieSpeedControl, PeriodicBoxButton
                self.movie_speed = MovieSpeedControl(cmd)
                quickbuttonslayout.addWidget(self.movie_speed)
                for setting in ('movie_fps', 'movie_delay'):
                    self.setting_callbacks[cmd.setting._get_index(setting)].append(self.movie_speed.sync)
                analysis_button = QtWidgets.QPushButton('轨迹分析…')
                analysis_button.setObjectName('trajectory_analysis')
                analysis_button.clicked.connect(self.open_trajectory_analysis)
                self.periodic_box_button = PeriodicBoxButton(
                    cmd, lambda message: self.statusBar().showMessage(message, 5000))
                tools = QtWidgets.QHBoxLayout()
                tools.setSpacing(2)
                tools.addWidget(analysis_button, 1)
                tools.addWidget(self.periodic_box_button)
                quickbuttonslayout.addLayout(tools)

        # progress bar
        hbox = QtWidgets.QHBoxLayout()
        self.progressbar = QtWidgets.QProgressBar()
        self.progressbar.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Minimum,
                QtWidgets.QSizePolicy.Policy.Minimum)
        hbox.addWidget(self.progressbar)
        self.abortbutton = QtWidgets.QPushButton(zh('Abort'))
        self.abortbutton.setStyleSheet("background: #FF0000; color: #FFFFFF")
        self.abortbutton.released.connect(cmd.interrupt)
        hbox.addWidget(self.abortbutton)
        quickbuttonslayout.addLayout(hbox)

        quickbuttonslayout.addStretch()
        quickbuttons_stretch_index = quickbuttonslayout.count() - 1

        # menu top level
        self.menubar = menubar = self.menuBar()

        # action groups
        actiongroups = {}

        def _addmenu(data, menu):
            '''Fill a menu from "data"'''
            menu.setTearOffEnabled(True)
            menu.setWindowTitle(menu.title())  # needed for Windows
            for item in data:
                if item[0] == 'separator':
                    menu.addSeparator()
                elif item[0] == 'menu':
                    _addmenu(item[2], menu.addMenu(zh(item[1]).replace('&', '&&')))
                elif item[0] == 'command':
                    command = item[2]
                    if command is None:
                        print('warning: skipping', item)
                    else:
                        if isinstance(command, str):
                            command = lambda c=command: cmd.do(c)
                        menu.addAction(zh(item[1]), command)
                elif item[0] == 'check':
                    if len(item) > 4:
                        menu.addAction(
                            SettingAction(self, cmd, item[2], zh(item[1]),
                                          item[3], item[4]))
                    else:
                        menu.addAction(
                            SettingAction(self, cmd, item[2], zh(item[1])))
                elif item[0] == 'radio':
                    label, name, value = item[1:4]
                    label = zh(label)
                    try:
                        group, type_, values = actiongroups[item[2]]
                    except KeyError:
                        group = QtWidgets.QActionGroup(self)
                        type_, values = cmd.get_setting_tuple(name)
                        actiongroups[item[2]] = group, type_, values
                    action = QtWidgets.QAction(label, self)
                    action.triggered.connect(lambda _=0, args=(name, value):
                                             cmd.set(*args, log=1, quiet=0))

                    self.setting_callbacks[cmd.setting._get_index(
                        name)].append(
                            lambda v, V=value, a=action: a.setChecked(v == V))

                    group.addAction(action)
                    menu.addAction(action)
                    action.setCheckable(True)
                    if values[0] == value:
                        action.setChecked(True)
                elif item[0] == 'open_recent_menu':
                    self.open_recent_menu = menu.addMenu(zh('Open Recent...'))
                else:
                    print('error:', item)

        # recent files menu
        self.open_recent_menu = None

        # for plugins
        self.menudict = {'': menubar}

        # menu
        for _, label, data in self.get_menudata(cmd):
            assert _ == 'menu'
            menu = menubar.addMenu(zh(label))
            self.menudict[label] = menu
            _addmenu(data, menu)

        self._plugin_menu = self.menudict['Plugin']
        self._add_builtin_plugin_actions()

        # hack for macOS to hide "Edit > Start Dictation"
        # https://bugreports.qt.io/browse/QTBUG-43217
        if pymol.IS_MACOS:
            self.menudict['Edit'].setTitle('Edit_')
            QtCore.QTimer.singleShot(10, lambda:
                    self.menudict['Edit'].setTitle(zh('Edit')))

        # recent files menu
        if self.open_recent_menu:
            @self.open_recent_menu.aboutToShow.connect
            def _():
                self.open_recent_menu.clear()
                for fname in self.recent_filenames:
                    self.open_recent_menu.addAction(
                            fname if len(fname) < 128 else '...' + fname[-120:],
                            lambda fname=fname: self.load_dialog(fname))
                if not self.open_recent_menu.actions():
                    self.open_recent_menu.addAction('没有最近文件').setEnabled(False)

        self.welcome_panel = None
        if show_welcome:
            from .welcome import WelcomePanel
            self.welcome_panel = WelcomePanel(self)

        # 显示菜单：选择与表示备用入口
        self.menudict['Display'].addSeparator()
        self.menudict['Display'].addAction(
            '选择与表示（V）…',
            lambda: self.open_visual_representations('all'))

        # some experimental window control
        menu = self.menudict['Display'].addSeparator()
        menu = self.menudict['Display'].addMenu(zh('External GUI'))
        menu.addAction(zh('Toggle dockable'), self.toggle_ext_window_dockable).setShortcut(
            QtGui.QKeySequence('Ctrl+E'))

        ext_vis_action = self.ext_window.toggleViewAction()
        ext_vis_action.setText(zh('Visible'))
        menu.addAction(ext_vis_action)

        # extra key mappings (MacPyMOL compatible)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+O'), self).activated.connect(self.file_open)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+S'), self).activated.connect(self.session_save)

        # feedback
        self.feedback_timer = QtCore.QTimer()
        self.feedback_timer.setSingleShot(True)
        self.feedback_timer.timeout.connect(self.update_feedback)
        self.feedback_timer.start(100)

        # 顶部出图工具栏（挂 feedback_timer 同步勾选；不自建 QTimer）
        from pmg_qt.figure_toolbar import install_figure_toolbar
        self.figure_toolbar = install_figure_toolbar(self, cmd)

        # legacy plugin system
        self.menudict['Plugin'].addAction(
            zh('Initialize Plugin System'), self.initializePlugins)

        # focus in command line
        if options.external_gui:
            self.lineedit.setFocus()
        else:
            self.pymolwidget.setFocus()

        # Apply PyMOL stylesheet
        try:
            with open(cmd.exp_path('$PYMOL_DATA/pmg_qt/styles/pymol.sty')) as f:
                style = f.read()
        except IOError:
            print('Could not read PyMOL stylesheet.')
            print('DEBUG: PYMOL_DATA=' + repr(os.getenv('PYMOL_DATA')))
            style = ""

        if style:
            self.setStyleSheet(style)

        # Load saved shortcuts on launch
        self.saved_shortcuts = pymol.save_shortcut.load_and_set(self.cmd)

    def lineeditKeyPressEventFilter(self, watched, event):
        key = event.key()
        if key == Qt.Key.Key_Tab:
            self.complete()
        elif key == Qt.Key.Key_Up:
            if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
                self.back_search()
            else:
                self.back()
        elif key == Qt.Key.Key_Down:
            self.forward()
        elif key == Qt.Key.Key_Return or key == Qt.Key.Key_Enter:
            # filter out "Return" instead of binding lineedit.returnPressed,
            # because otherwise OrthoKey would capture it as well.
            self.doPrompt()
        else:
            return False
        return True

    def eventFilter(self, watched, event):
        '''
        Filter out <Tab> event to do tab-completion instead of move focus
        '''
        type_ = event.type()
        if type_ == QtCore.QEvent.Type.KeyRelease:
            if event.key() == Qt.Key.Key_Tab:
                # silently skip tab release
                return True
        elif type_ == QtCore.QEvent.Type.KeyPress:
            if watched is self.lineedit:
                return self.lineeditKeyPressEventFilter(watched, event)
            elif event.key() == Qt.Key.Key_Tab:
                self.keyPressEvent(event)
                return True
        return False

    def toggle_ext_window_dockable(self, neverfloat=False):
        '''
        Toggle whether the "external" GUI is dockable
        '''
        dockWidget = self.ext_window

        if dockWidget.titleBarWidget() is None:
            tbw = QtWidgets.QWidget()
        else:
            tbw = None

        dockWidget.setFloating(tbw is None and not neverfloat)
        dockWidget.setTitleBarWidget(tbw)
        dockWidget.show()

    def toggle_fullscreen(self, toggle=-1):
        '''
        Full screen
        '''
        is_fullscreen = self.windowState() == Qt.WindowState.WindowFullScreen

        if toggle == -1:
            toggle = not is_fullscreen

        if not is_fullscreen:
            self._ext_window_visible = self.ext_window.isVisible()

        if toggle:
            self.menubar.hide()
            if not self.ext_window.isFloating():
                self.ext_window.hide()
            self.showFullScreen()
            self.pymolwidget.setFocus()
        else:
            self.menubar.show()
            if self._ext_window_visible:
                self.ext_window.show()
            self.showNormal()

    @property
    def initialdir(self):
        '''
        Be in sync with cd/pwd on the console until the first file has been
        browsed, then remember the last directory.
        '''
        return self._initialdir or os.getcwd()

    @initialdir.setter
    def initialdir(self, value):
        self._initialdir = value

    ##################
    # UI Forms
    ##################

    def load_form(self, name, dialog=None):
        '''Load a form from pmg_qt/forms/{name}.py'''
        import importlib
        if dialog is None:
            dialog = QtWidgets.QDialog(self)
            widget = dialog
        elif dialog == 'floating':
            widget = QtWidgets.QWidget(self)
        else:
            widget = dialog

        try:
            m = importlib.import_module('.forms.' + name, 'pmg_qt')
        except ImportError as e:
            if pymol.Qt.DEBUG:
                print('load_form import failed (%s)' % (e,))
            uifile = os.path.join(os.path.dirname(__file__), 'forms', '%s.ui' % name)
            form = pymol.Qt.utils.loadUi(uifile, widget)
        else:
            if hasattr(m, 'Ui_Form'):
                form = m.Ui_Form()
            else:
                form = m.Ui_Dialog()

            form.setupUi(widget)

        if dialog == 'floating':
            dialog = QtWidgets.QDockWidget(widget.windowTitle(), self)
            dialog.setFloating(True)
            dialog.setWidget(widget)
            dialog.resize(widget.size())

        form._dialog = dialog
        return form

    def edit_colors_dialog(self):
        form = self.load_form('colors')
        form.list_colors.setSortingEnabled(True)

        # populate list with named colors
        for color_index in self.cmd.get_color_indices():
            form.list_colors.addItem(color_index[0])

        # update spinboxes for given color
        def load_color(name):
            index = self.cmd.get_color_index(name)
            if index == -1:
                return
            rgb = self.cmd.get_color_tuple(index)
            form.input_R.setValue(rgb[0])
            form.input_G.setValue(rgb[1])
            form.input_B.setValue(rgb[2])

        # update spinbox from slider
        spinbox_lock = [False]
        def update_spinbox(spinbox, value):
            if not spinbox_lock[0]:
                spinbox.setValue(value / 100.)

        # update sliders and colored frame
        def update_gui(*args):
            spinbox_lock[0] = True
            R = form.input_R.value()
            G = form.input_G.value()
            B = form.input_B.value()
            form.slider_R.setValue(round(R * 100))
            form.slider_G.setValue(round(G * 100))
            form.slider_B.setValue(round(B * 100))
            form.frame_color.setStyleSheet(
                "background-color: rgb(%d,%d,%d)" % (
                    R * 0xFF, G * 0xFF, B * 0xFF))
            spinbox_lock[0] = False

        def run():
            name  = form.input_name.text()
            R = form.input_R.value()
            G = form.input_G.value()
            B = form.input_B.value()

            self.cmd.do('set_color %s, [%.2f, %.2f, %.2f]\nrecolor' %
                        (name, R, G, B))

            # if new color, insert and make current row
            if not form.list_colors.findItems(name, Qt.MatchFlag.MatchExactly):
                form.list_colors.addItem(name)
                form.list_colors.setCurrentItem(
                    form.list_colors.findItems(name, Qt.MatchFlag.MatchExactly)[0])

        # hook up events
        form.slider_R.valueChanged.connect(lambda v: update_spinbox(form.input_R, v))
        form.slider_G.valueChanged.connect(lambda v: update_spinbox(form.input_G, v))
        form.slider_B.valueChanged.connect(lambda v: update_spinbox(form.input_B, v))
        form.input_R.valueChanged.connect(update_gui)
        form.input_G.valueChanged.connect(update_gui)
        form.input_B.valueChanged.connect(update_gui)
        form.input_name.textChanged.connect(load_color)
        form.list_colors.currentTextChanged.connect(form.input_name.setText)
        form.button_apply.clicked.connect(run)

        form._dialog.show()

    def open_builder_panel(self):
        from pmg_qt.builder import BuilderPanelDocked
        from pymol import plugins

        app = plugins.get_pmgapp()
        if not self.builder:
            self.builder = BuilderPanelDocked(self, app)
            self.addDockWidget(Qt.DockWidgetArea.TopDockWidgetArea, self.builder)

        self.builder.show()
        self.builder.raise_()

    def open_props_dialog(self):
        from .properties_dialog import PropsDialog

        if not self.props_panel:
            self.props_panel = PropsDialog(self)

        self.props_panel.get_dialog().show()
        self.props_panel.get_dialog().raise_()

    def open_publication_dialog(self):
        """打开中文科研出图预览与导出对话框。"""
        from pmg_qt.publication_gui import open_publication_dialog
        open_publication_dialog(self)

    def open_ptm_dialog(self):
        if self._quit_requested or self._shutting_down:
            return
        from .ptm_gui import PTMDialog
        dialog = getattr(self, 'ptm_dialog', None)
        if dialog is None:
            dialog = self.ptm_dialog = PTMDialog(self)
        elif not dialog.isVisible():
            dialog.refresh()
        dialog.show()
        dialog.raise_()

    def open_md_setup_dialog(self):
        if self._quit_requested or self._shutting_down:
            return
        from .md_setup_gui import MDSetupDialog
        dialog = getattr(self, 'md_setup_dialog', None)
        if dialog is None:
            dialog = self.md_setup_dialog = MDSetupDialog(self)
        elif not dialog.is_running():
            dialog.refresh()
        dialog.show()
        dialog.raise_()

    def apply_publication_smart_theme(self, theme):
        """菜单入口：应用智能风格预览（与工具栏共用控制器）。"""
        from pmg_qt.publication_gui import apply_smart_theme
        apply_smart_theme(self, theme)

    def open_structure_compare(self):
        from .structure_compare_gui import StructureCompareDialog
        dialog = getattr(self, 'structure_compare_dialog', None)
        if dialog is None:
            dialog = self.structure_compare_dialog = StructureCompareDialog(self)
        else:
            dialog.refresh()
        dialog.show()
        dialog.raise_()

    def open_trajectory_analysis(self):
        from .trajectory_analysis_gui import TrajectoryAnalysisDialog
        dialog = getattr(self, 'trajectory_analysis_dialog', None)
        if dialog is None:
            dialog = self.trajectory_analysis_dialog = TrajectoryAnalysisDialog(self)
        dialog.show()
        dialog.raise_()

    def open_visual_representations(self, target='all'):
        """打开选择与表示窗口（Queued 钩子 / 命令栏 V / 显示菜单）。"""
        if self._quit_requested or self._shutting_down or getattr(self.cmd, '_COb', None) is None:
            return
        from pmg_qt.visual_representations_gui import open_visual_representations
        open_visual_representations(self, target=target or 'all')

    def open_color_palette(self, target='all'):
        """打开原生 C 对象色板（Queued 钩子）。"""
        if self._quit_requested or self._shutting_down or getattr(self.cmd, '_COb', None) is None:
            return
        from pmg_qt.color_palette import open_object_palette
        open_object_palette(self, target=target or 'all')

    def edit_pymolrc(self):
        from . import TextEditor
        from pymol import plugins
        TextEditor.edit_pymolrc(plugins.get_pmgapp())

    ##################
    # Menu callbacks
    ##################

    def file_open(self):
        fnames = getOpenFileNames(self, '打开结构或会话', self.initialdir)[0]
        partial = 0
        for fname in fnames:
            if not self.load_dialog(fname, partial=partial):
                break
            partial = 1

    def session_save(self):
        fname = self.cmd.get('session_file')
        fname = self.cmd.as_pathstr(fname)
        return self.session_save_as(fname)

    @PopupOnException.decorator
    def session_save_as(self, fname=''):
        formats = [
            f'{APP_NAME} Session File (*.pse *.pze *.pse.gz)',
            f'{APP_NAME} Show File (*.psw *.pzw *.psw.gz)',
        ]
        if not fname:
            fname = getSaveFileNameWithExt(
                self,
                'Save Session As...',
                self.initialdir,
                filter=';;'.join(formats))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.save(fname, format='pse', quiet=0)
            self.recent_filenames_add(fname)

    def render_dialog(self, widget=None):
        form = self.load_form('render', widget)
        lock = UpdateLock([ZeroDivisionError])

        def get_factor():
            units = form.input_units.currentText()
            factor = 1.0 if units == 'inch' else 2.54
            return factor / float(form.input_dpi.currentText())

        @lock.skipIfCircular
        def update_units(*args):
            width = form.input_width.value()
            height = form.input_height.value()
            factor = get_factor()
            form.input_width_units.setValue(width * factor)
            form.input_height_units.setValue(height * factor)

        @lock.skipIfCircular
        def update_pixels(*args):
            width = form.input_width_units.value()
            height = form.input_height_units.value()
            factor = get_factor()
            form.input_width.setValue(int(width / factor))
            form.input_height.setValue(int(height / factor))

        @lock.skipIfCircular
        def update_width(*args):
            if form.aspectratio > 0:
                width = form.input_height.value() * form.aspectratio
                form.input_width.setValue(int(width))
                form.input_width_units.setValue(width * get_factor())

        @lock.skipIfCircular
        def update_height(*args):
            if form.aspectratio > 0:
                height = form.input_width.value() / form.aspectratio
                form.input_height.setValue(int(height))
                form.input_height_units.setValue(height * get_factor())

        def update_aspectratio(checked=True):
            if checked:
                try:
                    form.aspectratio = (
                            float(form.input_width.value()) /
                            float(form.input_height.value()))
                except ZeroDivisionError:
                    form.button_lock.setChecked(False)
            else:
                form.aspectratio = 0

        def update_from_viewport():
            w, h = self.cmd.get_viewport()
            form.aspectratio = 0
            form.input_width.setValue(w)
            form.input_height.setValue(h)
            update_aspectratio(form.button_lock.isChecked())

        def run_draw(ray=False):
            width = form.input_width.value()
            height = form.input_height.value()
            if ray:
                self.cmd.set('opaque_background',
                        not form.input_transparent.isChecked())
                self.cmd.do('ray %d, %d, async=1' % (width, height))
            else:
                self.cmd.do('draw %d, %d' % (width, height))
            form.stack.setCurrentIndex(1)

        def run_ray():
            run_draw(ray=True)

        def run_save():
            fname = getSaveFileNameWithExt(self, 'Save As...', self.initialdir,
                    filter='PNG File (*.png)')
            if not fname:
                return
            self.initialdir = os.path.dirname(fname)
            self.cmd.png(fname, prior=1, dpi=form.input_dpi.currentText())

        def run_copy_clipboard():
            with PopupOnException():
                _copy_image(self.cmd, False, form.input_dpi.currentText())

        dpi = self.cmd.get_setting_int('image_dots_per_inch')
        if dpi > 0:
            form.input_dpi.setEditText(str(dpi))
        form.input_dpi.setValidator(QtGui.QIntValidator())

        # This connection used to be in the .ui file, but that fails with Qt6
        form.input_units.currentTextChanged.connect(lambda s: form.input_height_units.setSuffix(s))
        form.input_units.currentTextChanged.connect(lambda s: form.input_width_units.setSuffix(s))

        form.input_units.currentIndexChanged.connect(update_units)
        form.input_dpi.editTextChanged.connect(update_pixels)
        form.input_width.valueChanged.connect(update_units)
        form.input_height.valueChanged.connect(update_units)
        form.input_width_units.valueChanged.connect(update_pixels)
        form.input_height_units.valueChanged.connect(update_pixels)

        # set values before connecting mutual width<->height updates
        update_from_viewport()

        form.input_width.valueChanged.connect(update_height)
        form.input_height.valueChanged.connect(update_width)
        form.input_width_units.valueChanged.connect(update_height)
        form.input_height_units.valueChanged.connect(update_width)
        form.button_lock.toggled.connect(update_aspectratio)

        form.button_draw.clicked.connect(run_draw)
        form.button_ray.clicked.connect(run_ray)
        form.button_current.clicked.connect(update_from_viewport)
        form.button_back.clicked.connect(lambda: form.stack.setCurrentIndex(0))
        form.button_clip.clicked.connect(run_copy_clipboard)
        form.button_save.clicked.connect(run_save)

        if widget is None:
            form._dialog.show()
        return form

    @PopupOnException.decorator
    def _file_save(self, filter, format):
        fname = getSaveFileNameWithExt(
            self,
            'Save As...',
            self.initialdir,
            filter=filter)
        if fname:
            self.cmd.save(fname, format=format, quiet=0)

    def file_save_wrl(self):
        self._file_save('VRML 2 WRL File (*.wrl)', 'wrl')

    def file_save_dae(self):
        self._file_save('COLLADA File (*.dae)', 'dae')

    def file_save_pov(self):
        self._file_save('POV File (*.pov)', 'pov')

    def file_save_mpng(self):
        self.file_save_mpeg('png')

    def file_save_mov(self):
        self.file_save_mpeg('mov')

    def file_save_stl(self):
        self._file_save('STL File (*.stl)', 'stl')

    def file_save_glb(self):
        self._file_save('GLB File (*.glb)', 'glb')

    def file_save_gltf(self):
        self._file_save('GLTF File (*.gltf)', 'gltf')

    def file_save_usda(self):
        self._file_save('OpenUSD ASCII File (*.usda)', 'usda')

    def file_save_usdz(self):
        self._file_save('OpenUSD Package (*.usdz)', 'usdz')

    LOG_FORMATS = [
        f'{APP_NAME} Script (*.pml)',
        'Python Script (*.py *.pym)',
        'All (*)',
    ]

    def log_open(self, fname='', mode='w'):
        if not fname:
            fname = getSaveFileNameWithExt(self, 'Open Logfile...', self.initialdir,
                                    filter=';;'.join(self.LOG_FORMATS))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.log_open(fname, mode)

    def log_append(self):
        return self.log_open(mode='a')

    def log_resume(self):
        fname = getSaveFileNameWithExt(self, 'Open Logfile...', self.initialdir,
                                filter=';;'.join(self.LOG_FORMATS))
        if fname:
            self.initialdir = os.path.dirname(fname)
            self.cmd.resume(fname)

    def file_run(self):
        formats = [
            'All Runnable (*.pml *.py *.pym)',
            f'{APP_NAME} Command Script (*.pml)',
            f'{APP_NAME} Command Script (*.txt)',
            'Python Script (*.py *.pym)',
            'Python Script (*.txt)',
            'All Files(*)',
        ]
        fnames, selectedfilter = getOpenFileNames(
            self, 'Open file', self.initialdir, filter=';;'.join(formats))
        is_py = selectedfilter.startswith('Python')

        with PopupOnException():
            for fname in fnames:
                self.initialdir = os.path.dirname(fname)
                self.cmd.cd(self.initialdir, quiet=0)
                # detect: .py, .pym, .pyc, .pyo, .py.txt
                if is_py or re.search(r'\.py(|m|c|o|\.txt)$', fname, re.I):
                    self.cmd.run(fname)
                else:
                    self.cmd.do("@" + fname)

    def cd_dialog(self):
        dname = QFileDialog.getExistingDirectory(
            self, "Change Working Directory", self.initialdir)
        self.cmd.cd(dname or '.', quiet=0)

    def confirm_quit(self):
        # File/Quit 与 Cmd+Q：与窗口关闭同一优雅路径
        self.request_quit(0)

    def request_quit(self, code=0):
        """
        请求优雅退出（幂等，任意线程安全）。

        通过 QueuedConnection 投递到 GUI 线程；不在调用栈内直接 exit。
        """
        self.quitsignal.emit(int(code))

    def _handle_quit_request(self, code):
        """GUI 线程处理退出请求：打断渲染，必要时等待出图线程结束后再关机。"""
        if self._quit_requested:
            return
        self._quit_requested = True
        self._exit_code = int(code)
        self.feedback_timer.stop()
        self.setEnabled(False)
        try:
            self.cmd.interrupt()
        except Exception:
            pass
        pw = getattr(self, 'pymolwidget', None)
        if pw is not None:
            pw.stop_drawing()
        if self._background_export_busy():
            self._schedule_quit_when_export_idle()
            return
        self._begin_shutdown()

    def _background_export_busy(self):
        """Wait for figure and GROMACS export workers before destroying Qt."""
        md_dialog = getattr(self, 'md_setup_dialog', None)
        if md_dialog is not None and md_dialog.is_running():
            return True
        try:
            from pmg_qt.publication_gui import publication_export_running
            return publication_export_running()
        except Exception:
            return False

    def _schedule_quit_when_export_idle(self):
        """出图未结束时轮询等待；不销毁运行中的 QThread，不阻塞主线程等 ray。"""
        def _poll():
            if self._background_export_busy():
                QtCore.QTimer.singleShot(100, _poll)
                return
            self._begin_shutdown()
        QtCore.QTimer.singleShot(100, _poll)

    def _begin_shutdown(self):
        """将最终清理推迟到下一事件循环节拍（绝不在 closeEvent 内退出）。"""
        QtCore.QTimer.singleShot(0, self._finish_shutdown)

    def _close_command_log(self):
        """关闭命令日志文件（若有）。"""
        log = getattr(pymol, '_log_file', None)
        if log is None:
            return
        try:
            try:
                log.flush()
            except Exception:
                pass
            try:
                log.close()
            except Exception:
                pass
        finally:
            pymol._log_file = None

    def _finish_shutdown(self):
        """
        正式关机：停 GL → 关日志 → 带 GL 上下文停引擎 → 释放控件 → 退出事件循环。
        """
        if self._shutdown_done:
            return
        self._shutdown_done = True
        self._shutting_down = True
        self.feedback_timer.stop()

        pw = getattr(self, 'pymolwidget', None)
        if pw is not None:
            pw.stop_drawing()

        try:
            from pmg_qt import publication_gui as pg
            dlg = pg._dialog_ref
            if dlg is not None:
                try:
                    dlg.hide()
                except Exception:
                    pass
            pg.clear_publication_dialog_ref()
        except Exception:
            pass

        self._close_command_log()

        # pcatch 写入引擎反馈区；停引擎前解除重定向，避免清理异常写入已释放内存。
        import pcatch
        if sys.stdout is pcatch:
            sys.stdout = sys.__stdout__
        if sys.stderr is pcatch:
            sys.stderr = sys.__stderr__

        # 只保持 Qt 的 GL 上下文；不能使用 with pw，后者会在 stop 后再次访问引擎。
        if pw is not None and getattr(pw, 'pymol', None) is not None:
            pw.makeCurrent()
            try:
                pw.pymol.stop()
            finally:
                pw.doneCurrent()

        app = QtWidgets.QApplication.instance()
        try:
            # 仅隐藏；真正 deleteLater 在事件循环退出后、app 仍存活时执行
            self.hide()
        except Exception:
            pass

        if app is not None:
            app.exit(int(self._exit_code))

    def settings_edit_all_dialog(self):
        from .advanced_settings_gui import PyMOLAdvancedSettings
        if self.advanced_settings_dialog is None:
            self.advanced_settings_dialog = PyMOLAdvancedSettings(self,
                                                                  self.cmd)
        self.advanced_settings_dialog.show()

    def shortcut_menu_edit_dialog(self):
        from .shortcut_menu_gui import PyMOLShortcutMenu
        if self.shortcut_menu_filter_dialog is None:
            self.shortcut_menu_filter_dialog = PyMOLShortcutMenu(self, self.saved_shortcuts, self.cmd)
        self.shortcut_menu_filter_dialog.show()

    def scene_panel_menu_dialog(self):
        from .scene_bin_gui import ScenePanel

        if self.scene_panel_dialog is None:
            self.scene_panel_dialog = ScenePanel(self)

        self.scene_panel_dialog.show()

    def show_about(self):
        build = local_build_info()
        msg = [
            f'{APP_NAME} {self.cmd.get_version()[0]}',
            '基于 Open-Source PyMOL 的第三方分子可视化工作台。',
            '选择与表示参考 VMD，快捷工具栏参考 ChimeraX。',
            '',
            'Open-Source PyMOL is Copyright (C) Schrodinger, LLC.',
            'All Rights Reserved.',
            'PyMOL is a trademark of Schrodinger, LLC.',
            'https://pymol.org',
        ]
        if build.get('commit') and build.get('built_at'):
            source = '（含工作区改动）' if build.get('dirty') else ''
            msg[1:1] = [f"本机构建 {build['commit']}{source}", f"更新于 {build['built_at']}"]
        dialog = QtWidgets.QMessageBox(self)
        dialog.setWindowTitle(f'关于 {APP_NAME}')
        dialog.setText('\n'.join(msg))
        dialog.setIcon(QtWidgets.QMessageBox.Icon.Information)
        try:
            with open(self.cmd.exp_path('$PYMOL_PATH/LICENSE'), encoding='utf-8') as license_file:
                dialog.setDetailedText(license_file.read())
        except OSError:
            pass
        dialog.exec()

    #################
    # GUI callbacks
    #################

    def command_get(self):
        return self.lineedit.text()

    def command_set(self, v):
        return self.lineedit.setText(v)

    def command_set_cursor(self, i):
        return self.lineedit.setCursorPosition(i)

    def update_progress(self):
        if self._shutting_down:
            return
        progress = int(self.cmd.get_progress() * 100)
        if progress >= 0:
            self.progressbar.setValue(progress)
            self.progressbar.show()
            self.abortbutton.show()
        else:
            self.progressbar.hide()
            self.abortbutton.hide()

    def update_feedback(self):
        if self._shutting_down:
            return
        self.update_progress()
        if self.welcome_panel is not None:
            self.welcome_panel.refresh()

        self.periodic_box_button.sync()
        feedback = self.cmd._get_feedback()
        if feedback:
            html = colorprinting.text2html('\n'.join(feedback))
            self.browser.appendHtml(html)

            scrollbar = self.browser.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

        for setting in self.cmd.get_setting_updates() or ():
            if setting in self.setting_callbacks:
                current_value = self.cmd.get_setting_tuple(setting)[1][0]
                for callback in self.setting_callbacks[setting]:
                    callback(current_value)

        self.feedback_timer.start(500)

    def doPrompt(self):
        command = self.command_get()
        if self.welcome_panel is not None and command.strip():
            self.welcome_panel.dismiss()
        self.doTypedCommand(command)
        self.pymolwidget._pymolProcess()
        self.lineedit.clear()
        self.feedback_timer.start(0)

    ##########################
    # legacy plugin system
    ##########################

    def _add_builtin_plugin_actions(self):
        from pymol import publication

        menu = self._plugin_menu
        pub_menu = menu.addMenu('科研出图')
        for theme in publication.list_themes():
            pub_menu.addAction(
                theme['label'],
                lambda checked=False, theme_id=theme['id']: self.apply_publication_smart_theme(theme_id),
            )
        pub_menu.addSeparator()
        pub_menu.addAction('预览与导出…', self.open_publication_dialog)
        menu.addAction('结构对照…', self.open_structure_compare)
        menu.addAction('轨迹分析…', self.open_trajectory_analysis)
        menu.addAction('蛋白质翻译后修饰（PTM）…', self.open_ptm_dialog)
        menu.addAction('生成 GROMACS 模拟包…', self.open_md_setup_dialog)
        menu.addSeparator()

    @PopupOnException.decorator
    def initializePlugins(self):
        from pymol import plugins
        from . import mimic_tk

        if 'PluginQt' in self.menudict:
            return
        self._plugin_menu.clear()
        self._add_builtin_plugin_actions()

        app = plugins.get_pmgapp()

        plugins.legacysupport.addPluginManagerMenuItem()

        # Redirect to Legacy submenu
        self.menudict['PluginQt'] = self.menudict['Plugin']
        self.menudict['Plugin'] = self.menudict['PluginQt'].addMenu(zh('Legacy Plugins'))
        self.menudict['Plugin'].setTearOffEnabled(True)
        self.menudict['PluginQt'].addSeparator()

        plugins.HAVE_QT = True
        plugins.initialize(app)

    def createlegacypmgapp(self):
        from . import mimic_pmg_tk as mimic
        pmgapp = mimic.PMGApp()
        pmgapp.menuBar = mimic.PmwMenuBar(self.menudict)
        return pmgapp

    def window_cmd(self, action, x, y, w, h):
        if action == 0: # hide
            self.hide()
        elif action == 1: # show
            self.show()
        elif action == 2: # position
            self.move(x, y)
        elif action == 3: # size (first two arguments)
            self.resize(x, y)
        elif action == 4: # box
            self.move(x, y)
            self.resize(w, h)
        elif action == 5: # maximize
            self.showMaximized()
        elif action == 6: # fit
            if hasattr(QtGui, 'QWindow') and self.windowHandle().visibility() in (
                    QtGui.QWindow.Visibility.Maximized, QtGui.QWindow.Visibility.FullScreen):
                return
            a = self.screen().availableGeometry()
            g = self.geometry()
            f = self.frameGeometry()
            w = min(f.width(), a.width())
            h = min(f.height(), a.height())
            x = max(min(f.x(), a.right() - w), a.x())
            y = max(min(f.y(), a.bottom() - h), a.y())
            self.setGeometry(
                x - f.x() + g.x(),
                y - f.y() + g.y(),
                w - f.width() + g.width(),
                h - f.height() + g.height(),
            )
        elif action == 7: # focus
            self.setFocus(Qt.FocusReason.OtherFocusReason)
        elif action == 8: # defocus
            self.clearFocus()


def commandoverloaddecorator(func):
    name = func.__name__
    func.__doc__ = getattr(pymol.cmd, name).__doc__
    setattr(pymol.cmd, name, func)
    pymol.cmd.extend(func)
    return func


def SettingAction(parent, cmd, name, label='', true_value=1, false_value=0,
                  command=None):
    '''
    Menu toggle action for a PyMOL setting

    parent: parent QObject
    cmd: PyMOL instance
    name: setting name
    label: menu item text
    '''
    if not label:
        label = name

    index = cmd.setting._get_index(name)
    type_, values = cmd.get_setting_tuple(index)
    action = QtWidgets.QAction(label, parent)

    if not command:
        command = lambda: cmd.set(
            index,
            true_value if action.isChecked() else false_value,
            log=1,
            quiet=0)

    parent.setting_callbacks[index].append(
        lambda v: action.setChecked(v != false_value))

    if type_ in (
            1,  # bool
            2,  # int
            3,  # float
            5,  # color
            6,  # str
    ):
        action.setCheckable(True)
        if values[0] == true_value:
            action.setChecked(True)
    else:
        print('TODO', type_, name)

    action.triggered.connect(command)
    return action

window = None


class CommandLineEdit(QtWidgets.QLineEdit):
    '''
    Line edit widget with instant text insert on drag-enter
    '''
    _saved_pos = -1

    def dragMoveEvent(self, event):
        pass

    def dropEvent(self, event):
        if event.mimeData().hasText():
            event.acceptProposedAction()

    def dragEnterEvent(self, event):
        if not event.mimeData().hasText():
            self._saved_pos = -1
            return

        event.acceptProposedAction()

        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            droppedtext = urls[0].toLocalFile()
        else:
            droppedtext = event.mimeData().text()

        pos = self.cursorPosition()
        text = self.text()
        self._saved_pos = pos
        self._saved_text = text

        self.setText(text[:pos] + droppedtext + text[pos:])
        self.setSelection(pos, len(droppedtext))

    def dragLeaveEvent(self, event):
        if self._saved_pos != -1:
            self.setText(self._saved_text)
            self.setCursorPosition(self._saved_pos)


class PyMOLApplication(QtWidgets.QApplication):
    '''
    Catch drop events on app icon
    '''
    # FileOpen event is only activated after the first
    # application state change, otherwise sys.argv would be
    # handled by Qt, we don't want that.

    def handle_file_open(self, ev):
        if ev.type() == QtCore.QEvent.Type.ApplicationActivate:
            self.handle_file_open = self.handle_file_open_active
        return False

    def handle_file_open_active(self, ev):
        if ev.type() != QtCore.QEvent.Type.FileOpen:
            return False

        # When double clicking a file in Finder, open it in a new instance
        if not pymol.invocation.options.reuse_helper and pymol.cmd.get_names():
            window.new_window([ev.file()])
            return True

        # pymol -I -U
        if pymol.invocation.options.auto_reinitialize:
            pymol.cmd.reinitialize()

        # PyMOL Show
        if ev.file().endswith('.psw'):
            pymol.cmd.set('presentation')
            pymol.cmd.set('internal_gui', 0)
            pymol.cmd.set('internal_feedback', 0)
            pymol.cmd.full_screen('on')

        window.load_dialog(ev.file())
        return True

    def event(self, ev):
        if ev.type() == QtCore.QEvent.Type.Quit:
            w = globals().get('window')
            if w is not None and not w._shutting_down:
                w.request_quit(0)
                return True
        if self.handle_file_open(ev):
            return True
        return super(PyMOLApplication, self).event(ev)


# like pymol.internal._copy_image
def _copy_image(_self=pymol.cmd, quiet=1, dpi=-1):
    import tempfile
    fname = tempfile.mktemp('.png')

    if not _self.png(fname, prior=1, dpi=dpi):
        print("no prior image")
        return

    try:
        qim = QtGui.QImage(fname)
        QtWidgets.QApplication.clipboard().setImage(qim)
    finally:
        os.unlink(fname)

    if not quiet:
        print(" Image copied to clipboard")


def make_pymol_qicon():
    # 窗口图标使用 PymolX 图；数据目录里没有新图时回退到上游图标
    icons_dir = os.path.expandvars('$PYMOL_DATA/pymol/icons')
    fp = os.path.join(icons_dir, 'pymolx.png')
    if not os.path.isfile(fp):
        fp = os.path.join(icons_dir, 'icon2.svg')
    return QtGui.QIcon(fp)


def _show_qt_plugins() -> None:
    """去掉 Qt 插件上的隐藏旗标。

    macOS 有时会把 .venv 里的库标成隐藏。Qt 列目录时会跳过这些文件，
    于是找不到 cocoa 插件，窗口还没出现程序就被系统中止。
    """
    plug = QtCore.QLibraryInfo.path(QtCore.QLibraryInfo.LibraryPath.PluginsPath)
    cocoa_dir = os.path.join(plug, "platforms")
    cocoa = os.path.join(cocoa_dir, "libqcocoa.dylib")
    if not os.path.isfile(cocoa):
        return
    # 已经能被 Qt 看到就不用改旗标
    if "libqcocoa.dylib" in QtCore.QDir(cocoa_dir).entryList():
        return
    import subprocess
    subprocess.run(["/usr/bin/chflags", "-R", "nohidden", plug], check=False)


_MISSING_QT_HOOK = object()

def _teardown_qt_hooks(app=None, previous_hooks=None):
    """移除 Qt 注册的 cmd 钩子与模块级窗口引用，避免单例悬挂。"""
    global window
    w = window
    window = None
    cg_gui = sys.modules.get("pmg_qt.coarse_grained_gui")
    if cg_gui is not None:
        cg_gui.clear_coarse_grained_dialog_ref()
    vis_gui = sys.modules.get("pmg_qt.visual_representations_gui")
    if vis_gui is not None:
        vis_gui.clear_visual_representations_ref()
    pal_gui = sys.modules.get("pmg_qt.color_palette")
    if pal_gui is not None:
        pal_gui.clear_color_palette_refs()
    if w is not None:
        try:
            w.deleteLater()
        except Exception:
            pass
    for name, original in (previous_hooks or {}).items():
        if original is _MISSING_QT_HOOK:
            if hasattr(pymol.cmd, name):
                delattr(pymol.cmd, name)
        else:
            setattr(pymol.cmd, name, original)
    try:
        import pymol.gui as _gui
        if getattr(_gui, 'createlegacypmgapp', None) is not None:
            _gui.createlegacypmgapp = None
    except Exception:
        pass
    if app is not None:
        try:
            app.processEvents()
            app.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete)
        except Exception:
            pass


def execapp():
    '''
    Run PyMOL as a Qt application

    Returns the process exit code for launch/.app wrappers.
    '''
    global window
    global pymol

    previous_hooks = {
        name: getattr(pymol.cmd, name, _MISSING_QT_HOOK)
        for name in ('_qt_quit', '_copy_image', '_call_in_gui_thread',
                     '_call_with_opengl_context', 'viewport', 'full_screen',
                     '_qt_visual_panel', '_qt_color_panel')
    }
    # 先保证 cocoa 插件可见，否则 QApplication 会直接 abort
    _show_qt_plugins()

    # don't let exceptions stop PyMOL
    import traceback
    sys.excepthook = traceback.print_exception

    # use QT_OPENGL=desktop (auto-detection may fail on Windows)
    if hasattr(Qt, 'AA_UseDesktopOpenGL') and pymol.IS_WINDOWS:
        QtCore.QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_UseDesktopOpenGL)

    # enable 4K scaling on Windows and Linux
    if hasattr(Qt, 'AA_EnableHighDpiScaling') and not any(
            v in os.environ
            for v in ['QT_SCALE_FACTOR', 'QT_SCREEN_SCALE_FACTORS']):
        QtCore.QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_EnableHighDpiScaling)

    # fix Windows taskbar icon
    if pymol.IS_WINDOWS:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                u'org.pymolx.PymolX')

    # 标记 GUI 线程，供 is_gui_thread / 后台 quit 路径使用
    pymol.prime_pymol()

    app = PyMOLApplication([APP_NAME])
    app.setApplicationName(APP_NAME)
    app.setApplicationDisplayName(APP_NAME)
    app.setWindowIcon(make_pymol_qicon())

    window = PyMOLQtGUI()
    window.setWindowTitle(APP_NAME)

    # fix gnome/wayland dash icon/missing wmclass
    app.setDesktopFileName("org.pymolx.PymolX")

    @commandoverloaddecorator
    def viewport(w=-1, h=-1, _self=None):
        window.viewportsignal.emit(int(w), int(h))

    @commandoverloaddecorator
    def full_screen(toggle=-1, _self=None):
        from pymol import viewing as v
        toggle = v.toggle_dict[v.toggle_sc.auto_err(str(toggle), 'toggle')]
        window.toggle_fullscreen(toggle)

    import pymol.gui
    pymol.gui.createlegacypmgapp = window.createlegacypmgapp

    pymol.cmd._copy_image = _copy_image
    pymol.cmd._call_in_gui_thread = MainThreadCaller()

    # Assume GUI thread, make OpenGL context current before calling func().
    def _call_with_opengl_context_gui_thread(func):
        with window.pymolwidget:
            return func()

    # Dispatch to GUI thread and make OpenGL context current before calling func().
    pymol.cmd._call_with_opengl_context = lambda func: pymol.cmd._call_in_gui_thread(
        lambda: _call_with_opengl_context_gui_thread(func))

    # 仅 Qt 前端注册：cmd.quit / _quit 走优雅生命周期，而非 PExit
    def _qt_quit(code=0, _self=None):
        w = window
        if w is not None:
            w.request_quit(int(code))
        return None

    pymol.cmd._qt_quit = _qt_quit

    # 原生对象列表 V/C：接收 target 后 Queued emit，禁止在调用栈内直接建窗
    def _qt_visual_panel(sele='all', _self=None):
        w = window
        if w is not None:
            w.visualsignal.emit(str(sele if sele is not None else 'all'))
        return None

    def _qt_color_panel(sele='all', _self=None):
        w = window
        if w is not None:
            w.palettesignal.emit(str(sele if sele is not None else 'all'))
        return None

    pymol.cmd._qt_visual_panel = _qt_visual_panel
    pymol.cmd._qt_color_panel = _qt_color_panel

    window.show()
    window.raise_()

    # window size according to -W -H options
    options = pymol.invocation.options
    if options.win_xy_set:
        scale = window.pymolwidget.fb_scale
        viewport(scale * options.win_x, scale * options.win_y)

    # load plugins
    if options.plugins:
        window.initializePlugins()

    exit_code = 0
    try:
        exit_code = app.exec()
    finally:
        # 事件循环结束后、QApplication 析构前：补齐清理并清钩子
        try:
            if window is not None and not getattr(window, '_shutdown_done', False):
                window._exit_code = int(exit_code or 0)
                window._finish_shutdown()
        except Exception:
            pass
        _teardown_qt_hooks(app, previous_hooks)
        # PySide 在解释器结束前主动释放 QApplication 及其绑定对象。
        shutdown = getattr(app, 'shutdown', None)
        if callable(shutdown):
            shutdown()

    return int(exit_code) if exit_code is not None else 0
