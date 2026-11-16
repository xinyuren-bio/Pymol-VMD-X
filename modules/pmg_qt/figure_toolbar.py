# ==================================================
# 功能说明：ChimeraX 风格顶部出图工具栏（智能风格/背景/光照/效果/结构显示）。
# 使用方法：由 pymol_qt_gui.install_figure_toolbar 挂载；测试可直接调用本模块纯函数。
# 依赖环境：PyMOL + PySide6（pymol.Qt）；不下载外部图标资源。
# 生成时间：2026-11-14 13:30
# ==================================================

"""顶部分组出图工具栏：真实 PyMOL 设置与表示，不含自动 ray。"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from pymol.Qt import QtCore, QtGui, QtWidgets

Qt = QtCore.Qt

# 会话私有标志：工具栏是否改过背景/光照/轮廓阴影等外观
_APPEARANCE_FLAG = "_pymolx_figure_appearance_customized"

# 结构类表示（互斥切换；不触碰 labels/测量等）
_STRUCTURE_REPS: Tuple[str, ...] = ("cartoon", "sticks", "spheres", "surface")

# 背景色名 -> 近似 RGB（用于勾选同步）
_BG_RGB: Dict[str, Tuple[float, float, float]] = {
    "white": (1.0, 1.0, 1.0),
    "grey": (0.5, 0.5, 0.5),
    "gray": (0.5, 0.5, 0.5),
    "black": (0.0, 0.0, 0.0),
}

# 光照预设：每套完整、与前一预设无关的确定性参数包（不含 ray_shadows / ray_trace_mode）
LIGHT_PRESETS: Dict[str, Dict[str, float]] = {
    "简洁": {
        "ambient": 0.35,
        "direct": 0.55,
        "reflect": 0.35,
        "specular": 0.25,
        "shininess": 40.0,
        "light_count": 2.0,
        "ambient_occlusion_mode": 0.0,
    },
    "柔和": {
        "ambient": 0.48,
        "direct": 0.38,
        "reflect": 0.42,
        "specular": 0.12,
        "shininess": 18.0,
        "light_count": 4.0,
        "ambient_occlusion_mode": 0.0,
    },
    "立体": {
        "ambient": 0.18,
        "direct": 0.32,
        "reflect": 0.72,
        "specular": 0.48,
        "shininess": 60.0,
        "light_count": 2.0,
        "ambient_occlusion_mode": 1.0,
    },
    "平光": {
        "ambient": 0.72,
        "direct": 0.28,
        "reflect": 0.08,
        "specular": 0.04,
        "shininess": 8.0,
        "light_count": 1.0,
        "ambient_occlusion_mode": 0.0,
    },
}


def is_appearance_customized(cmd: Any) -> bool:
    """工具栏是否已定制过背景/光照/效果外观。"""
    return bool(getattr(cmd, _APPEARANCE_FLAG, False))


def mark_appearance_customized(cmd: Any) -> None:
    """标记当前会话外观已被工具栏改过（供出图对话框默认勾选）。"""
    setattr(cmd, _APPEARANCE_FLAG, True)


def clear_appearance_customized(cmd: Any) -> None:
    """清除外观定制标志（测试或重新初始化用）。"""
    if hasattr(cmd, _APPEARANCE_FLAG):
        try:
            delattr(cmd, _APPEARANCE_FLAG)
        except Exception:
            setattr(cmd, _APPEARANCE_FLAG, False)


def _cmd_alive(cmd: Any) -> bool:
    """引擎是否仍可用（关机后 _COb 为 None）。"""
    return cmd is not None and getattr(cmd, "_COb", None) is not None


def enabled_mol_objects(cmd: Any) -> List[str]:
    """返回当前启用的分子对象名（排除非分子类型）。"""
    if not _cmd_alive(cmd):
        return []
    try:
        names = cmd.get_names("public_objects", enabled_only=1) or []
    except Exception:
        return []
    out: List[str] = []
    for n in names:
        try:
            if cmd.get_type(n) == "object:molecule":
                out.append(n)
        except Exception:
            continue
    return out


def _scope_sele(cmd: Any, objs: Optional[Sequence[str]] = None) -> Optional[str]:
    """
    构造启用分子对象并排除溶剂的选择串。

    空场景返回 None，调用方应优雅跳过表示操作。
    """
    names = list(objs) if objs is not None else enabled_mol_objects(cmd)
    if not names:
        return None
    parts = []
    for o in names:
        safe = o.replace("\\", "\\\\").replace('"', '\\"')
        parts.append('(model "%s")' % safe)
    return "((%s) and not solvent)" % (" or ".join(parts))


def _get_float(cmd: Any, name: str, default: float = 0.0) -> float:
    """安全读取浮点设置。"""
    try:
        return float(cmd.get_setting_float(name))
    except Exception:
        return default


def _get_int(cmd: Any, name: str, default: int = 0) -> int:
    """安全读取整数设置。"""
    try:
        return int(round(float(cmd.get_setting_float(name))))
    except Exception:
        return default


def read_bg_rgb(cmd: Any) -> Optional[Tuple[float, float, float]]:
    """读取当前背景 RGB（0–1）；失败返回 None。"""
    if not _cmd_alive(cmd):
        return None
    try:
        t = cmd.get_setting_tuple("bg_rgb")
        v = t[1]
        if isinstance(v, (list, tuple)) and len(v) >= 3:
            return float(v[0]), float(v[1]), float(v[2])
        idx = int(v[0]) if isinstance(v, (list, tuple)) else int(v)
        rgb = cmd.get_color_tuple(idx)
        if rgb and len(rgb) >= 3:
            return float(rgb[0]), float(rgb[1]), float(rgb[2])
    except Exception:
        return None
    return None


def match_background(cmd: Any) -> Optional[str]:
    """将当前背景匹配为 white/grey/black；无法判断则 None。"""
    if _get_int(cmd, "bg_gradient", 0):
        return None
    rgb = read_bg_rgb(cmd)
    if rgb is None:
        return None
    best = None
    best_d = 0.12
    for key, ref in (("white", _BG_RGB["white"]), ("grey", _BG_RGB["grey"]), ("black", _BG_RGB["black"])):
        d = abs(rgb[0] - ref[0]) + abs(rgb[1] - ref[1]) + abs(rgb[2] - ref[2])
        if d < best_d:
            best_d = d
            best = key
    return best


def set_background(cmd: Any, color: str) -> None:
    """
    设置纯色背景并关闭 bg_gradient，使色块与实际背景一致。

    color: white | grey | gray | black
    """
    if not _cmd_alive(cmd):
        return
    name = "grey" if color == "gray" else color
    if name not in ("white", "grey", "black"):
        name = "white"
    try:
        cmd.set("bg_gradient", 0)
    except Exception:
        pass
    cmd.bg_color(name)
    mark_appearance_customized(cmd)


def apply_light_preset(cmd: Any, key: str) -> None:
    """
    应用一套完整光照预设（确定性，与先前预设无关）。

    不修改 ray_shadows / ray_trace_mode，效果开关保持独立。
    """
    if not _cmd_alive(cmd):
        return
    preset = LIGHT_PRESETS.get(key)
    if not preset:
        return
    for name, val in preset.items():
        try:
            cmd.set(name, val)
        except Exception:
            # ambient_occlusion_mode 等在旧内核可能不可用，跳过单项
            continue
    mark_appearance_customized(cmd)


def match_light_preset(cmd: Any) -> Optional[str]:
    """若当前光照与某预设一致（容差内）则返回键名，否则 None（避免陈旧勾选）。"""
    if not _cmd_alive(cmd):
        return None
    for key, preset in LIGHT_PRESETS.items():
        ok = True
        for name, val in preset.items():
            cur = _get_float(cmd, name, default=-999.0)
            if abs(cur - float(val)) > 0.06:
                ok = False
                break
        if ok:
            return key
    return None


def get_effect_state(cmd: Any) -> Dict[str, bool]:
    """读取阴影/轮廓/雾化当前开关（阴影与轮廓主要在 ray 时生效）。"""
    return {
        "shadow": _get_int(cmd, "ray_shadows", 0) != 0,
        "outline": _get_int(cmd, "ray_trace_mode", 0) == 1,
        "fog": _get_int(cmd, "depth_cue", 0) != 0,
    }


def set_effect_shadow(cmd: Any, on: bool) -> None:
    """切换 ray_shadows（光线追踪/导出时阴影，非实时预览保证）。"""
    if not _cmd_alive(cmd):
        return
    cmd.set("ray_shadows", 1 if on else 0)
    mark_appearance_customized(cmd)


def set_effect_outline(cmd: Any, on: bool) -> None:
    """切换轮廓：ray_trace_mode 1/0（主要在 ray/导出时可见）。"""
    if not _cmd_alive(cmd):
        return
    cmd.set("ray_trace_mode", 1 if on else 0)
    mark_appearance_customized(cmd)


def set_effect_fog(cmd: Any, on: bool) -> None:
    """切换雾化 depth_cue。"""
    if not _cmd_alive(cmd):
        return
    cmd.set("depth_cue", 1 if on else 0)
    mark_appearance_customized(cmd)


def apply_representation(cmd: Any, rep: str) -> None:
    """
    在启用分子对象上切换结构表示（排除溶剂）。

    仅隐藏/显示 cartoon/sticks/spheres/surface，保留标签与测量对象；
    不删除任何数据；禁用对象与空场景安全跳过。
    """
    if not _cmd_alive(cmd):
        return
    if rep not in _STRUCTURE_REPS:
        return
    sele = _scope_sele(cmd)
    if not sele:
        return
    for r in _STRUCTURE_REPS + ("lines", "ribbon", "nonbonded", "nb_spheres"):
        cmd.hide(r, sele)
    if rep == "cartoon":
        cmd.show("cartoon", "(%s) and polymer" % sele)
        cmd.show("sticks", "(%s) and organic and not polymer" % sele)
    else:
        cmd.show(rep, sele)


def dominant_structure_rep(cmd: Any) -> Optional[str]:
    """根据启用分子上可见原子数推断当前主结构表示。"""
    sele = _scope_sele(cmd)
    if not sele:
        return None
    best = None
    best_n = 0
    for r in _STRUCTURE_REPS:
        try:
            n = int(cmd.count_atoms("(%s) and rep %s" % (sele, r)))
        except Exception:
            n = 0
        if n > best_n:
            best_n = n
            best = r
    return best if best_n > 0 else None


# ---------------------------------------------------------------------------
# 图标绘制（原生 QPainter，约 26px，无外部资源/emoji）
# ---------------------------------------------------------------------------

def _mk_icon(paint: Callable[[QtGui.QPainter, int], None], size: int = 26) -> QtGui.QIcon:
    """用绘制回调生成 QIcon。"""
    pm = QtGui.QPixmap(size * 2, size * 2)
    pm.setDevicePixelRatio(2.0)
    pm.fill(QtCore.Qt.GlobalColor.transparent)
    p = QtGui.QPainter(pm)
    p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    paint(p, size)
    p.end()
    return QtGui.QIcon(pm)


def _icon_swatch(rgb: Tuple[int, int, int]) -> QtGui.QIcon:
    """背景色块图标。"""
    def _p(p: QtGui.QPainter, s: int) -> None:
        p.setPen(QtGui.QPen(QtGui.QColor(120, 120, 120), 1))
        p.setBrush(QtGui.QColor(*rgb))
        m = 3
        p.drawRoundedRect(m, m, s - 2 * m, s - 2 * m, 3, 3)
    return _mk_icon(_p)


def _icon_light(kind: str) -> QtGui.QIcon:
    """光照预设小图标。"""
    def _p(p: QtGui.QPainter, s: int) -> None:
        cx, cy = s / 2, s / 2
        if kind == "简洁":
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(90, 90, 95))
            p.drawEllipse(QtCore.QPointF(cx, cy), 7, 7)
            p.setBrush(QtGui.QColor(210, 210, 215))
            p.drawEllipse(QtCore.QPointF(cx - 2, cy - 2), 3, 3)
        elif kind == "柔和":
            for i, a in enumerate((40, 70, 110)):
                p.setPen(QtCore.Qt.PenStyle.NoPen)
                p.setBrush(QtGui.QColor(160, 160, 165, a))
                p.drawEllipse(QtCore.QPointF(cx, cy), 10 - i * 2, 10 - i * 2)
        elif kind == "立体":
            path = QtGui.QPainterPath()
            path.moveTo(cx, 5)
            path.lineTo(s - 5, s - 6)
            path.lineTo(5, s - 6)
            path.closeSubpath()
            p.setPen(QtGui.QPen(QtGui.QColor(70, 70, 75), 1))
            grad = QtGui.QLinearGradient(5, 5, s - 5, s - 5)
            grad.setColorAt(0.0, QtGui.QColor(230, 230, 235))
            grad.setColorAt(1.0, QtGui.QColor(80, 80, 85))
            p.setBrush(grad)
            p.drawPath(path)
        else:  # 平光
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(190, 190, 195))
            p.drawRect(4, 4, s - 8, s - 8)
            p.setBrush(QtGui.QColor(150, 150, 155))
            p.drawEllipse(QtCore.QPointF(cx, cy), 5, 5)
    return _mk_icon(_p)


def _icon_effect(kind: str) -> QtGui.QIcon:
    """效果开关图标。"""
    def _p(p: QtGui.QPainter, s: int) -> None:
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        if kind == "shadow":
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(180, 180, 185))
            p.drawEllipse(8, 5, 12, 12)
            p.setBrush(QtGui.QColor(120, 120, 125, 120))
            p.drawEllipse(10, 14, 12, 7)
        elif kind == "outline":
            p.setBrush(QtGui.QColor(210, 210, 215))
            p.setPen(QtGui.QPen(QtGui.QColor(40, 40, 45), 2.2))
            p.drawEllipse(6, 6, s - 12, s - 12)
        else:  # fog
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            for i, y in enumerate((8, 13, 18)):
                p.setBrush(QtGui.QColor(140, 140, 145, 60 + i * 40))
                p.drawRoundedRect(4, y, s - 8, 5, 2, 2)
    return _mk_icon(_p)


def _icon_rep(kind: str) -> QtGui.QIcon:
    """结构表示图标。"""
    def _p(p: QtGui.QPainter, s: int) -> None:
        p.setPen(QtGui.QPen(QtGui.QColor(70, 70, 75), 1.4))
        if kind == "cartoon":
            path = QtGui.QPainterPath()
            path.moveTo(4, s - 6)
            path.cubicTo(8, 4, 12, s - 2, 16, 8)
            path.cubicTo(18, 4, 20, 14, s - 4, 10)
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.drawPath(path)
        elif kind == "sticks":
            p.drawLine(6, s - 7, s - 7, 7)
            p.setBrush(QtGui.QColor(90, 90, 95))
            p.drawEllipse(4, s - 11, 7, 7)
            p.drawEllipse(s - 11, 4, 7, 7)
        elif kind == "spheres":
            p.setBrush(QtGui.QColor(150, 150, 155))
            p.drawEllipse(4, 8, 10, 10)
            p.drawEllipse(12, 6, 11, 11)
        else:  # surface
            path = QtGui.QPainterPath()
            path.moveTo(4, s - 5)
            path.cubicTo(8, 4, 14, 4, s - 4, s - 6)
            path.lineTo(s - 4, s - 4)
            path.lineTo(4, s - 4)
            path.closeSubpath()
            p.setBrush(QtGui.QColor(170, 170, 175, 160))
            p.drawPath(path)
    return _mk_icon(_p)


def _tool_btn(
    text: str,
    icon: QtGui.QIcon,
    tip: str,
    checkable: bool = True,
    width: int = 52,
) -> QtWidgets.QToolButton:
    """创建图标在上、文字在下的紧凑工具按钮。"""
    b = QtWidgets.QToolButton()
    b.setText(text)
    b.setIcon(icon)
    b.setToolTip(tip)
    b.setCheckable(checkable)
    b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
    b.setIconSize(QtCore.QSize(26, 26))
    b.setAutoRaise(True)
    b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    b.setFixedWidth(int(width))
    return b


def _icon_theme(kind: str) -> QtGui.QIcon:
    """
    三种智能风格原生图标（期刊简约 / 卡通插画 / 通透表面）。

    设计思路：用色块与简单形状区分主题，避免 emoji 与外部资源。
    """
    def _p(p: QtGui.QPainter, s: int) -> None:
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        if kind == "journal":
            # 白底 + 蓝蛋白条带 + 橙配体点
            p.setPen(QtGui.QPen(QtGui.QColor(160, 160, 165), 1))
            p.setBrush(QtGui.QColor(255, 255, 255))
            p.drawRoundedRect(3, 3, s - 6, s - 6, 3, 3)
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(90, 130, 180))
            path = QtGui.QPainterPath()
            path.moveTo(6, s - 8)
            path.cubicTo(10, 6, 14, s - 4, 18, 10)
            p.setPen(QtGui.QPen(QtGui.QColor(70, 110, 160), 2.2))
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.drawPath(path)
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(230, 130, 40))
            p.drawEllipse(QtCore.QPointF(s - 8, s - 9), 3.2, 3.2)
        elif kind == "illustration":
            # 青绿卡通弧 + 珊瑚配体
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(235, 248, 248))
            p.drawRoundedRect(3, 3, s - 6, s - 6, 3, 3)
            p.setPen(QtGui.QPen(QtGui.QColor(40, 140, 150), 3.0))
            path = QtGui.QPainterPath()
            path.moveTo(5, s - 7)
            path.cubicTo(9, 5, 13, s - 3, 17, 9)
            path.cubicTo(19, 5, 21, 14, s - 4, 11)
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.drawPath(path)
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(240, 120, 110))
            p.drawEllipse(QtCore.QPointF(s - 8, s - 8), 3.5, 3.5)
        else:  # surface
            # 半透明冷色表面块 + 不透明配体棒
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(210, 225, 235, 200))
            path = QtGui.QPainterPath()
            path.moveTo(4, s - 5)
            path.cubicTo(8, 5, 14, 5, s - 5, s - 6)
            path.lineTo(s - 5, s - 4)
            path.lineTo(4, s - 4)
            path.closeSubpath()
            p.drawPath(path)
            p.setPen(QtGui.QPen(QtGui.QColor(220, 120, 40), 2.2))
            p.drawLine(8, s - 10, s - 8, 8)
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(220, 120, 40))
            p.drawEllipse(QtCore.QPointF(8, s - 10), 2.5, 2.5)
            p.drawEllipse(QtCore.QPointF(s - 8, 8), 2.5, 2.5)
    return _mk_icon(_p)


def _group_widget(caption: str, buttons: Sequence[QtWidgets.QToolButton]) -> QtWidgets.QWidget:
    """带分组标题的一列按钮组。"""
    w = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(w)
    lay.setContentsMargins(6, 4, 6, 4)
    lay.setSpacing(2)
    cap = QtWidgets.QLabel(caption)
    cap.setAlignment(Qt.AlignmentFlag.AlignHCenter)
    cap.setStyleSheet("color: #555; font-size: 11px;")
    lay.addWidget(cap)
    row = QtWidgets.QHBoxLayout()
    row.setSpacing(2)
    row.setContentsMargins(0, 0, 0, 0)
    for b in buttons:
        row.addWidget(b)
    lay.addLayout(row)
    return w


_TOOLBAR_STYLE = """
QToolBar#figure_toolbar {
    background: #F2F2F2;
    border: none;
    border-bottom: 1px solid #D0D0D0;
    spacing: 0px;
    padding: 0px;
}
QToolBar#figure_toolbar QToolButton {
    color: #333;
    background: transparent;
    border: 1px solid transparent;
    border-radius: 4px;
    padding: 2px 1px 1px 1px;
    font-size: 11px;
}
QToolBar#figure_toolbar QToolButton:hover {
    background: #E6E6E6;
    border-color: #C8C8C8;
}
QToolBar#figure_toolbar QToolButton:checked {
    background: #DDDDDD;
    border-color: #B0B0B0;
}
QToolBar#figure_toolbar QFrame[figureSep="true"] {
    background: #C8C8C8;
    max-width: 1px;
    min-width: 1px;
    margin: 10px 4px;
}
"""


class FigureToolBar(QtWidgets.QToolBar):
    """
    顶部横向分组出图工具栏（约 100px 高，窄窗可横向滚动）。

    状态刷新挂在主窗口 feedback_timer / showEvent 上，不自建 QTimer。
    """

    def __init__(self, window: QtWidgets.QMainWindow, cmd: Any) -> None:
        super().__init__("出图工具栏", window)
        self._window = window
        self.cmd = cmd
        self._refreshing = False
        self.setObjectName("figure_toolbar")
        self.setMovable(False)
        self.setFloatable(False)
        self.setAllowedAreas(Qt.ToolBarArea.TopToolBarArea)
        self.setIconSize(QtCore.QSize(26, 26))
        self.setStyleSheet(_TOOLBAR_STYLE)
        self.setMinimumHeight(96)
        self.setMaximumHeight(104)

        # 可横向滚动内容，避免逼迫主窗口变宽
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFixedHeight(100)

        inner = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(inner)
        row.setContentsMargins(4, 0, 8, 0)
        row.setSpacing(0)
        row.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)

        # --- 智能风格（首组，约 72px 宽中文标签）---
        from pymol import publication as pub

        self._theme_btns: Dict[str, QtWidgets.QToolButton] = {}
        theme_specs = (
            (pub.THEME_JOURNAL, "期刊简约", "期刊简约（Nature 风格）：蓝色蛋白、橙色配体、白底柔和光照"),
            (pub.THEME_ILLUSTRATION, "卡通插画", "卡通插画：青绿色蛋白、珊瑚色配体；渲染时带轮廓描边"),
            (pub.THEME_SURFACE, "通透表面", "通透表面：蛋白半透明冷色表面 + 不透明配体球棍"),
        )
        theme_list = []
        for tid, lab, tip in theme_specs:
            b = _tool_btn(lab, _icon_theme(tid), tip, checkable=False, width=72)
            self._theme_btns[tid] = b
            b.clicked.connect(lambda _=False, t=tid: self._on_theme(t))
            theme_list.append(b)
        row.addWidget(_group_widget("智能风格", theme_list))
        row.addWidget(self._sep())

        # 粗粒化显示与三个智能风格独立，保留粗粒化坐标。
        self.btn_cg = _tool_btn("粗粒化", _icon_rep("sticks"),
                                "连接粗粒化珠，显示脂质、蛋白和核酸的连续形态", checkable=False)
        self.btn_cg.setFixedWidth(64)
        self.btn_cg.clicked.connect(self._on_cg)
        row.addWidget(_group_widget("粗粒化", (self.btn_cg,)))
        row.addWidget(self._sep())

        # --- 背景 ---
        self.btn_bg_white = _tool_btn(
            "白", _icon_swatch((255, 255, 255)), "白色背景（关闭渐变）")
        self.btn_bg_grey = _tool_btn(
            "灰", _icon_swatch((128, 128, 128)), "灰色背景（关闭渐变）")
        self.btn_bg_black = _tool_btn(
            "黑", _icon_swatch((20, 20, 20)), "黑色背景（关闭渐变）")
        # 非强制互斥组：允许“无一匹配”时全部不勾，避免陈旧预设勾选
        self._bg_group = QtWidgets.QButtonGroup(self)
        self._bg_group.setExclusive(False)
        for b in (self.btn_bg_white, self.btn_bg_grey, self.btn_bg_black):
            self._bg_group.addButton(b)
        self.btn_bg_white.clicked.connect(lambda: self._on_bg("white"))
        self.btn_bg_grey.clicked.connect(lambda: self._on_bg("grey"))
        self.btn_bg_black.clicked.connect(lambda: self._on_bg("black"))
        row.addWidget(_group_widget("背景", (self.btn_bg_white, self.btn_bg_grey, self.btn_bg_black)))
        row.addWidget(self._sep())

        # --- 光照 ---
        self._light_btns: Dict[str, QtWidgets.QToolButton] = {}
        light_list = []
        tips = {
            "简洁": "简洁光照：中等环境光与少量光源",
            "柔和": "柔和光照：偏高环境光、低高光",
            "立体": "立体光照：对比更强，利于形体感",
            "平光": "平光：高环境光、弱反射，明暗差小",
        }
        self._light_group = QtWidgets.QButtonGroup(self)
        self._light_group.setExclusive(False)
        for key in ("简洁", "柔和", "立体", "平光"):
            b = _tool_btn(key, _icon_light(key), tips[key])
            self._light_btns[key] = b
            self._light_group.addButton(b)
            b.clicked.connect(lambda _=False, k=key: self._on_light(k))
            light_list.append(b)
        row.addWidget(_group_widget("光照", light_list))
        row.addWidget(self._sep())

        # --- 效果 ---
        self.btn_shadow = _tool_btn(
            "阴影",
            _icon_effect("shadow"),
            "开启或关闭阴影；光线追踪渲染或导出时生效",
        )
        self.btn_outline = _tool_btn(
            "轮廓",
            _icon_effect("outline"),
            "开启或关闭轮廓描边；光线追踪渲染或导出时生效",
        )
        self.btn_fog = _tool_btn(
            "雾化",
            _icon_effect("fog"),
            "开启或关闭远处结构的雾化效果，实时预览可见",
        )
        self.btn_shadow.clicked.connect(self._on_shadow)
        self.btn_outline.clicked.connect(self._on_outline)
        self.btn_fog.clicked.connect(self._on_fog)
        row.addWidget(_group_widget("效果", (self.btn_shadow, self.btn_outline, self.btn_fog)))
        row.addWidget(self._sep())

        # --- 结构显示 ---
        self._rep_btns: Dict[str, QtWidgets.QToolButton] = {}
        rep_labels = (
            ("cartoon", "卡通", "蛋白和核酸显示为卡通，小分子显示为棒状；保留标签与测距"),
            ("sticks", "棒状", "棒状表示（仅启用分子，排除溶剂；不改标签/测量）"),
            ("spheres", "球体", "球体表示（仅启用分子，排除溶剂；不改标签/测量）"),
            ("surface", "表面", "表面表示（仅启用分子，排除溶剂；不改标签/测量）"),
        )
        self._rep_group = QtWidgets.QButtonGroup(self)
        self._rep_group.setExclusive(False)
        rep_list = []
        for rep, lab, tip in rep_labels:
            b = _tool_btn(lab, _icon_rep(rep), tip)
            self._rep_btns[rep] = b
            self._rep_group.addButton(b)
            b.clicked.connect(lambda _=False, r=rep: self._on_rep(r))
            rep_list.append(b)
        row.addWidget(_group_widget("结构显示", rep_list))
        row.addStretch(1)

        scroll.setWidget(inner)
        self.addWidget(scroll)

        # 挂到主反馈定时器，避免自建 QTimer
        timer = getattr(window, "feedback_timer", None)
        if timer is not None:
            timer.timeout.connect(self._safe_refresh)

        self._safe_refresh()

    def _sep(self) -> QtWidgets.QFrame:
        """细分隔线。"""
        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        line.setProperty("figureSep", True)
        line.setFixedWidth(1)
        return line

    def _guard(self) -> bool:
        """退出/关机/引擎已释放时禁止刷新或写设置。"""
        w = self._window
        if getattr(w, "_quit_requested", False) or getattr(w, "_shutting_down", False):
            return False
        return _cmd_alive(self.cmd)

    def _safe_refresh(self) -> None:
        """带守卫的勾选状态同步。"""
        if not self._guard():
            return
        self.refresh_checked_states()

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        """显示时同步一次状态。"""
        super().showEvent(event)
        self._safe_refresh()

    def refresh_checked_states(self) -> None:
        """根据真实设置/表示刷新勾选，避免外部命令或会话恢复后陈旧状态。"""
        if self._refreshing or not self._guard():
            return
        self._refreshing = True
        # 允许真实状态不匹配任何预设（例如自定义背景），避免独占组留下旧勾选。
        groups = (self._bg_group, self._light_group, self._rep_group)
        for group in groups:
            group.setExclusive(False)
        try:
            # 导出中禁用风格按钮
            busy = False
            try:
                from pmg_qt.publication_gui import publication_export_running
                busy = publication_export_running()
            except Exception:
                busy = False
            for b in self._theme_btns.values():
                b.setEnabled(not busy)

            bg = match_background(self.cmd)
            self.btn_bg_white.setChecked(bg == "white")
            self.btn_bg_grey.setChecked(bg == "grey")
            self.btn_bg_black.setChecked(bg == "black")

            lit = match_light_preset(self.cmd)
            for key, b in self._light_btns.items():
                b.setChecked(lit == key)

            fx = get_effect_state(self.cmd)
            self.btn_shadow.setChecked(fx["shadow"])
            self.btn_outline.setChecked(fx["outline"])
            self.btn_fog.setChecked(fx["fog"])

            dom = dominant_structure_rep(self.cmd)
            for rep, b in self._rep_btns.items():
                b.setChecked(dom == rep)
        finally:
            for group in groups:
                group.setExclusive(True)
            self._refreshing = False

    def _on_theme(self, theme: str) -> None:
        """智能风格是一次性预览操作，不把上次点击误示为当前精确外观。"""
        if not self._guard():
            return
        from pmg_qt.publication_gui import apply_smart_theme
        apply_smart_theme(self._window, theme)
        self.refresh_checked_states()

    def _on_cg(self):
        if self._guard():
            from pmg_qt.coarse_grained_gui import open_coarse_grained_dialog
            open_coarse_grained_dialog(self._window)

    def _on_bg(self, color: str) -> None:
        if not self._guard():
            return
        set_background(self.cmd, color)
        self.refresh_checked_states()

    def _on_light(self, key: str) -> None:
        if not self._guard():
            return
        apply_light_preset(self.cmd, key)
        self.refresh_checked_states()

    def _on_shadow(self) -> None:
        if not self._guard():
            return
        set_effect_shadow(self.cmd, self.btn_shadow.isChecked())
        self.refresh_checked_states()

    def _on_outline(self) -> None:
        if not self._guard():
            return
        set_effect_outline(self.cmd, self.btn_outline.isChecked())
        self.refresh_checked_states()

    def _on_fog(self) -> None:
        if not self._guard():
            return
        set_effect_fog(self.cmd, self.btn_fog.isChecked())
        self.refresh_checked_states()

    def _on_rep(self, rep: str) -> None:
        if not self._guard():
            return
        apply_representation(self.cmd, rep)
        self.refresh_checked_states()


def install_figure_toolbar(window: QtWidgets.QMainWindow, cmd: Any) -> FigureToolBar:
    """
    将出图工具栏加入主窗口顶部，并在“显示”菜单提供显隐。

    放在外置命令坞上方的 ToolBar 区域；窄窗口内容可横向滚动。
    """
    bar = FigureToolBar(window, cmd)
    window.addToolBar(Qt.ToolBarArea.TopToolBarArea, bar)
    window.figure_toolbar = bar  # type: ignore[attr-defined]

    menudict = getattr(window, "menudict", None) or {}
    disp = menudict.get("Display")
    if disp is not None:
        act = bar.toggleViewAction()
        act.setText("出图工具栏")
        act.setToolTip("显示或隐藏顶部出图工具栏")
        disp.addSeparator()
        disp.addAction(act)
    if disp is not None:
        disp.addAction("粗粒化显示…", bar._on_cg)
    return bar
