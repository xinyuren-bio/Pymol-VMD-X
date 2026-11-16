"""The input adapter and GL resize must share one physical pixel space."""
import unittest
from types import SimpleNamespace
from pymol.Qt import QtCore
from pmg_qt.pymol_gl_widget import PyMOLGLWidget


class FakeViewport:
    _physical_size = PyMOLGLWidget._physical_size
    updateFbScale = PyMOLGLWidget.updateFbScale
    resizeGL = PyMOLGLWidget.resizeGL
    _event_x_y_mod = PyMOLGLWidget._event_x_y_mod

    def __init__(self, scale):
        self.scale=scale; self.fb_scale=1.; self.shapes=[]; self.settings=[]
        self.cmd=SimpleNamespace(set=lambda *a:self.settings.append(a))
        self.pymol=SimpleNamespace(reshape=lambda *a:self.shapes.append(a))

    def devicePixelRatioF(self): return self.scale
    def width(self): return 101
    def height(self): return 81
    def window(self): return self
    def windowHandle(self): return self
    def update(self): pass


class CoordinateTests(unittest.TestCase):
    def test_pixel_edges_at_integer_and_fractional_scale(self):
        for scale in (1.,1.25,1.5,2.,2.5):
            view=FakeViewport(scale);view.resizeGL(view.width(),view.height())
            width,height=view._framebuffer_size
            def convert(x,y):
                event=SimpleNamespace(position=lambda:QtCore.QPointF(x,y),
                    modifiers=lambda:QtCore.Qt.KeyboardModifier.NoModifier)
                return view._event_x_y_mod(event)[:2]
            self.assertEqual(convert(0,0),(0,height-1))
            self.assertEqual(convert(100.999,80.999),(width-1,0))
            self.assertEqual(view.shapes[-1],(width,height,True))
            self.assertLess(convert(-.1,0)[0],0)

    def test_scale_change_reshapes_before_next_click(self):
        view=FakeViewport(2.);view.resizeGL(101,81)
        before=view._framebuffer_size
        view.scale=1.25
        event=SimpleNamespace(position=lambda:QtCore.QPointF(0,0),
            modifiers=lambda:QtCore.Qt.KeyboardModifier.NoModifier)
        view._event_x_y_mod(event)
        self.assertNotEqual(before,view._framebuffer_size)
        self.assertEqual(view.shapes[-1],(*view._framebuffer_size,True))


if __name__=='__main__': unittest.main()
