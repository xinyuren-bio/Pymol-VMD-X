import os
import sys
import math
from pymol2 import SingletonPyMOL as PyMOL

import pymol

from pymol.Qt import QtCore
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets
Gesture = QtCore.QEvent.Type.Gesture
Qt = QtCore.Qt

from .keymapping import get_modifiers
from .keymapping import get_wheel_button

# don't import the heavy OpenGL (PyOpenGL) module
from pymol._cmd import glViewport

# QOpenGLWidget is supposed to supersede QGLWidget, but has issues (e.g.
# no stereo support)
USE_QOPENGLWIDGET = int(
    os.getenv("PYMOL_USE_QOPENGLWIDGET") or
    (QtCore.QT_VERSION >= 0x50400 and pymol.IS_MACOS or
     QtCore.QT_VERSION >= 0x60000))

if USE_QOPENGLWIDGET:
    BaseGLWidget = QtWidgets.QOpenGLWidget
    AUTO_DETECT_STEREO = False
else:
    from pymol.Qt import QtOpenGL
    BaseGLWidget = QtOpenGL.QGLWidget
    # only attempt stereo detection in Qt <= 5.6 (with 5.9+ on Linux I
    # get GL_DOUBLEBUFFER=0 with flickering when requesting stereo)
    AUTO_DETECT_STEREO = pymol.IS_WINDOWS or QtCore.QT_VERSION < 0x50700


class PyMOLGLWidget(BaseGLWidget):
    '''
    PyMOL OpenGL Widget

    Can be used as a context manager to make OpenGL context current.
    '''

    # mouse button map
    _buttonMap = {
        Qt.MouseButton.LeftButton: 0,
        Qt.MouseButton.MiddleButton: 1,
        Qt.MouseButton.RightButton: 2,
    }

    def __enter__(self):
        '''
        Context manager to make the OpenGL context "current".

        Fixes depth-buffer issue with QOpenGLWidget
        https://github.com/schrodinger/pymol-open-source/issues/25
        '''
        self.makeCurrent()
        pymol._cmd._pushValidContext(self.cmd._COb)

    def __exit__(self, exc_type, exc_value, traceback):
        pymol._cmd._popValidContext(self.cmd._COb)

    def __init__(self, parent):
        self.gui = parent
        self.fb_scale = 1.0

        # OpenGL context setup
        if USE_QOPENGLWIDGET:
            f = QtGui.QSurfaceFormat()
        else:
            f = QtOpenGL.QGLFormat()

        from pymol.invocation import options

        # logic equivalent to layer5/main.cpp:launch

        if options.multisample:
            f.setSamples(4)

        if options.force_stereo != -1:
            # See layer1/Setting.h for stereo modes

            if options.stereo_mode in (1, 12) or (
                    options.stereo_mode == 0 and AUTO_DETECT_STEREO):
                f.setStereo(True)

            if options.stereo_mode in (11, 12) and not USE_QOPENGLWIDGET:
                f.setAccum(True)

        if USE_QOPENGLWIDGET:
            super(PyMOLGLWidget, self).__init__(parent=parent)
            self.setFormat(f)
            self.setUpdateBehavior(QtWidgets.QOpenGLWidget.UpdateBehavior.PartialUpdate)
        else:
            super(PyMOLGLWidget, self).__init__(f, parent=parent)

        # pymol instance
        self.pymol = PyMOL()
        self.pymol.start()
        self.cmd = self.pymol.cmd

        # capture python output for feedback
        import pcatch
        pcatch._install()

        # for passive move drag
        self.setMouseTracking(True)

        # for accepting keyboard input (command line, shortcuts)
        self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)

        # for idle rendering
        self._timer = QtCore.QTimer()
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._pymolProcess)
        # 退出流程中置位，阻止继续调度绘制
        self._drawing_stopped = False

        # drag n drop
        self.setAcceptDrops(True)

        # pinch-zoom
        self.grabGesture(Qt.GestureType.PinchGesture)

    def sizeHint(self):
        # default 640 + internal_gui, 480 + internal_feedback
        return QtCore.QSize(860, 498)

    ##########################
    # Input Events
    ##########################

    def event(self, ev):
        if ev.type() == Gesture:
            return self.gestureEvent(ev)

        result = super(PyMOLGLWidget, self).event(ev)
        scale_events = tuple(getattr(QtCore.QEvent.Type, name, None) for name in
                             ('DevicePixelRatioChange', 'ScreenChangeInternal'))
        if ev.type() in scale_events and hasattr(self, 'cmd'):
            self.updateFbScale(self.window().windowHandle())
        return result

    def gestureEvent(self, ev):
        gesture = ev.gesture(Qt.GestureType.PinchGesture)

        if gesture is None:
            return False

        if gesture.state() == Qt.GestureState.GestureStarted:
            self.pinch_start_z = self.cmd.get_view()[11]

        changeFlags = gesture.changeFlags()

        if changeFlags & QtWidgets.QPinchGesture.ChangeFlag.RotationAngleChanged:
            delta = gesture.lastRotationAngle() - gesture.rotationAngle()
            self.cmd.turn('z', delta)

        if changeFlags & QtWidgets.QPinchGesture.ChangeFlag.ScaleFactorChanged:
            view = list(self.cmd.get_view())

            # best guess for https://bugreports.qt.io/browse/QTBUG-48138
            totalscalefactor = gesture.totalScaleFactor()
            if totalscalefactor == 1.0:
                totalscalefactor = gesture.scaleFactor()

            z = self.pinch_start_z / totalscalefactor
            delta = z - view[11]
            view[11] = z
            view[15] -= delta
            view[16] -= delta
            self.cmd.set_view(view)

        return True

    def _event_x_y_mod(self, ev):
        # Pointer positions and reshape use the same physical framebuffer size.
        # Qt's origin is top-left; OpenGL pixel indices end at height - 1.
        self.updateFbScale(self.window().windowHandle())
        pos = ev.position() if hasattr(ev, "position") else ev.pos()
        width, height = self._framebuffer_size
        return (
            math.floor(pos.x() * width / max(1, self.width())),
            height - 1 - math.floor(pos.y() * height / max(1, self.height())),
            get_modifiers(ev),
        )

    def mouseMoveEvent(self, ev):
        if getattr(self, '_drawing_stopped', False):
            return
        self.pymol.drag(*self._event_x_y_mod(ev))

    def mousePressEvent(self, ev, state=0):
        if getattr(self, '_drawing_stopped', False):
            return
        if ev.button() not in self._buttonMap:
            return
        self.pymol.button(self._buttonMap[ev.button()], state,
                          *self._event_x_y_mod(ev))

    def mouseReleaseEvent(self, ev):
        self.mousePressEvent(ev, 1)

    def wheelEvent(self, ev):
        if getattr(self, '_drawing_stopped', False):
            return
        button = get_wheel_button(ev)
        if not button:
            return
        args = self._event_x_y_mod(ev)
        self.pymol.button(button, 0, *args)
        self.pymol.button(button, 1, *args)

    ##########################
    # OpenGL
    ##########################

    def stop_drawing(self):
        """停止空闲绘制定时器，避免退出过程中继续 paint/idle。"""
        self._drawing_stopped = True
        try:
            self._timer.stop()
        except Exception:
            pass

    def paintGL(self):
        if getattr(self, '_drawing_stopped', False):
            return
        self.updateFbScale(self.window().windowHandle())
        if not USE_QOPENGLWIDGET:
            glViewport(0, 0, int(self.fb_scale * self.width()),
                         int(self.fb_scale * self.height()))
        self.pymol.draw()
        self._timer.start(0)

    def resizeGL(self, w, h):
        if getattr(self, '_drawing_stopped', False):
            return
        self.updateFbScale(self.window().windowHandle(), reshape=False)
        if USE_QOPENGLWIDGET:
            w, h = self._physical_size(w, h)

        self._framebuffer_size = (w, h)
        self.pymol.reshape(w, h, True)

    def _physical_size(self, w, h):
        """Use Qt's nearest-pixel rounding for fractional display scales."""
        return (max(1, math.floor(w * self.fb_scale + 0.5)),
                max(1, math.floor(h * self.fb_scale + 0.5)))

    def updateFbScale(self, context=None, reshape=True):
        '''Update PyMOL's display scale factor from the window or screen context
        @type context: QWindow or QScreen
        '''
        if getattr(self, '_drawing_stopped', False):
            return
        try:
            scale = float(self.devicePixelRatioF())
        except AttributeError:
            scale = float(context.devicePixelRatio()) if context is not None else 1.0
        if not math.isfinite(scale) or scale <= 0:
            scale = 1.0
        changed = scale != self.fb_scale
        self.fb_scale = scale
        size = self._physical_size(self.width(), self.height())
        if changed:
            try:
                self.cmd.set('display_scale_factor', max(1, int(scale)))
            except Exception as exc:
                print(exc)
        if reshape and size != getattr(self, '_framebuffer_size', None):
            self._framebuffer_size = size
            self.pymol.reshape(*size, True)
            self.update()

    def initializeGL(self):
        # Scale framebuffer for Retina displays
        try:
            window = self.windowHandle()

            # QOpenGLWidget workaround
            if window is None:
                window = self.parent().windowHandle()

            self.updateFbScale(window)
            window.screenChanged.connect(self.updateFbScale)
            window.screen().physicalDotsPerInchChanged.connect(
                    lambda dpi: self.updateFbScale(window))

        except AttributeError:
            # Fallback for Qt4
            pass

    def _pymolProcess(self):
        if getattr(self, '_drawing_stopped', False):
            return
        idle = self.pymol.idle()
        if idle or self.pymol.getRedisplay():
            self.update()

        self._timer.start(20)

    ##########################
    # drag n drop
    ##########################

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    url = url.toLocalFile()
                else:
                    url = url.toString()
                self.gui.load_dialog(url)
            event.accept()
