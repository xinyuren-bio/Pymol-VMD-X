"""Run in a fresh PymolX process with -r to verify native panel hit regions.
Creates test objects, checks rendered pixel regions and exits the test app.
PYMOLX_TEST_ARTIFACTS sets the result directory; QT_SCALE_FACTOR tests fractional DPR.
"""
import json, math, traceback
from pathlib import Path
from pymol import cmd, menu
from pymol.Qt import QtCore, QtGui, QtWidgets
from pmg_qt import pymol_qt_gui as gui
from PySide6.QtTest import QTest

import os, tempfile
art=Path(os.environ.get('PYMOLX_TEST_ARTIFACTS') or tempfile.mkdtemp(prefix='pymolx_hits_'))
art.mkdir(parents=True, exist_ok=True)
tag=__import__('os').environ.get('PYMOLX_HIT_TAG','validation')
calls=[]
for name,col in [('all_action','A'),('mol_action','A'),('sele_action','A'),
                 ('mol_show','S'),('mol_hide','H'),('mol_labels','L'),
                 ('color_panel','C'),('visual_panel','V'),('obj_motion','M')]:
    def record(_cmd,target,_col=col):
        calls.append((_col,target.strip('()')));return []
    setattr(menu,name,record)
menu.camera_motion=lambda _cmd: calls.append(('M','all')) or []
for i in range(60): cmd.pseudoatom('obj%02d'%i,pos=[i,0,0])
cmd.order('obj*',sort=1)


def runs(values):
    out=[]
    for v in values:
        if out and v==out[-1][-1]+1: out[-1].append(v)
        else: out.append([v])
    return [(x[0],x[-1])for x in out if len(x)>4]


def scan(pw, columns):
    image=pw.grabFramebuffer()
    scale=pw.fb_scale
    display_scale=cmd.get_setting_int('display_scale_factor')
    x=image.width()-int((17 if columns==6 else 34)*display_scale)+2
    ys=[]
    for y in range(image.height()):
        c=image.pixelColor(x,y)
        if abs(c.red()-115)<3 and abs(c.green()-196)<3 and abs(c.blue()-209)<3:
            ys.append(y)
    bounds=runs(ys)
    assert bounds,('No V buttons',image.width(),x,scale,display_scale)
    y=bounds[0][0]+2
    xs=[]
    for xx in range(max(0,x-int(18*display_scale)),min(image.width(),x+int(18*display_scale))):
        c=image.pixelColor(xx,y)
        if abs(c.red()-115)<3 and abs(c.green()-196)<3 and abs(c.blue()-209)<3:
            xs.append(xx)
    xr=max(runs(xs),key=lambda r:r[1]-r[0])
    # Interior flat fill plus its one framebuffer pixel outline.
    vleft=xr[0]-1
    return image, bounds, vleft, display_scale


def click_fb(pw,x,y,image):
    # Map a rendered framebuffer pixel center to Qt's logical coordinates.
    local=QtCore.QPointF((x+.5)*pw.width()/image.width(),(y+.5)*pw.height()/image.height())
    global_pos=QtCore.QPointF(pw.mapToGlobal(local.toPoint()))
    cover=QtWidgets.QApplication.widgetAt(global_pos.toPoint())
    if os.environ.get('PYMOLX_INVALIDATE_BEFORE_CLICK') == '1':
        # Coloring clears Executive::Panel. Do not repaint or move first.
        click_fb.counter = getattr(click_fb, 'counter', 0) + 1
        cmd.color('0x88aabb' if click_fb.counter%2 else '0x99bbcc', 'obj00')
    for kind,buttons in [(QtCore.QEvent.Type.MouseButtonPress,QtCore.Qt.MouseButton.LeftButton),
                         (QtCore.QEvent.Type.MouseButtonRelease,QtCore.Qt.MouseButton.NoButton)]:
        event=QtGui.QMouseEvent(kind,local,global_pos,QtCore.Qt.MouseButton.LeftButton,buttons,QtCore.Qt.KeyboardModifier.NoModifier)
        QtWidgets.QApplication.sendEvent(pw,event)
    return None if cover is None else type(cover).__name__


def check():
    results=[]; failures=[]
    try:
        w=gui.window;pw=w.pymolwidget
        fixtures=[('three_button_viewing',1250,900,0,None),
                  ('three_button_motions',1250,900,0,None),
                  ('three_button_viewing',900,700,0,None),
                  ('three_button_viewing',900,700,7,None)]
        if os.environ.get('PYMOLX_MOVE_SCREENS') == '1':
            fixtures=[('three_button_viewing',900,700,0,screen)
                      for screen in QtGui.QGuiApplication.screens()]
        for mode,width,height,scroll,screen in fixtures:
            cmd.config_mouse(mode)
            if screen is not None:
                w.windowHandle().setScreen(screen)
                rect=screen.availableGeometry()
                width=min(width,rect.width()-30);height=min(height,rect.height()-70)
                w.move(rect.x()+10,rect.y()+30)
            else:
                w.move(10,30)
            w.resize(width,height);QTest.qWait(250)
            if scroll:
                # Native list scrollbar: wheel scrolls the executive list, not camera.
                for _ in range(scroll):
                    pw.pymol.button(4,0,int(pw.width()*pw.fb_scale-4),int(pw.height()*pw.fb_scale-40),0)
                    pw.pymol.button(4,1,int(pw.width()*pw.fb_scale-4),int(pw.height()*pw.fb_scale-40),0)
                QTest.qWait(120)
            cols=7 if 'motions' in mode else 6
            image,bounds,vleft,ds=scan(pw,cols)
            rows=sorted(set([0,len(bounds)//2,len(bounds)-1]))
            for row in rows:
                top,bottom=bounds[row]
                # Derive the expected target from the center before testing edges.
                calls.clear();click_fb(pw,vleft+8*ds,(top+bottom)//2,image)
                panel_index=row+scroll
                target='all' if panel_index==0 else 'obj%02d'%(panel_index-1)
                for col,label in enumerate('ASHLCV'+('M'if cols==7 else '')):
                    left=vleft-(5-col)*17*ds
                    right=left+16*ds
                    points={'center':((left+right)//2,(top+bottom)//2),
                            'left':(left+1,(top+bottom)//2),
                            'right':(right-1,(top+bottom)//2),
                            'top':((left+right)//2,top),
                            'bottom':((left+right)//2,bottom),
                            'top_left':(left+1,top),'top_right':(right-1,top),
                            'bottom_left':(left+1,bottom),'bottom_right':(right-1,bottom)}
                    for point,(x,y) in points.items():
                        calls.clear();cover=click_fb(pw,x,y,image)
                        got=list(calls);expected=[(label,target)]
                        item={'mode':mode,'size':[width,height],'scroll':scroll,'row':row,'col':label,'point':point,'expected':expected,'got':got,'cover':cover,'fb':[image.width(),image.height()],'fb_scale':pw.fb_scale,'ds':ds,'screen':w.windowHandle().screen().name()}
                        results.append(item)
                        if got!=expected:failures.append(item)
            image.save(str(art/('pymolx_hit_'+tag+'_'+mode+'_'+str(scroll)+'.png')))
        payload={'total':len(results),'failures':failures,'screen_count':len(QtGui.QGuiApplication.screens()),'coverage_widgets':sorted(set(str(x['cover'])for x in results)),'results':results}
        path=art/('pymolx_hit_matrix_'+tag+'.json');path.write_text(json.dumps(payload,indent=2))
        print('HIT_MATRIX',tag,'total',len(results),'failures',len(failures),'cover',payload['coverage_widgets'],flush=True)
        print('FIRST_FAILURES',failures[:5],flush=True)
        cmd.quit(1 if failures else 0)
    except BaseException:
        traceback.print_exc();cmd.quit(1)

QtCore.QTimer.singleShot(1000,check)
