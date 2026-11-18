# ==================================================
# 功能说明：VMD 常用选择语法安全解析与多层独立图形表示后端（无 Qt）。
# 使用方法：由 pmg_qt.visual_representations_gui / 测试导入；通过 RepresentationController 预览与 Apply/Restore。
# 依赖环境：已安装的 PyMOL（pymol.cmd、chempy）；可选 pymol.coarse_grained / cg_templates 做 CG 分类与连接。
# 生成时间：2026-11-14 16:50
# ==================================================

"""VMD 风格选择与独立表示层纯后端（无 Qt 依赖）。

科学与安全约定：
- 选择表达式解析为 AST 后在 get_model 原子集上求值；不执行输入为 Python / cmd.do。
- 未经完整校验前不改场景；Apply 事务式：新显示全部建好再替换旧副本。
- 源对象坐标 / 颜色 / 键 / 显示设置不变；仅禁用源并管理控制器拥有的副本。
- 本模块支持的是 VMD 常用子集，非完整 VMD/Tcl。
"""

from __future__ import annotations

import fnmatch
import math
import operator
import re
import uuid
from dataclasses import dataclass, field
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

from pymol import cmd as _global_cmd

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

OWNED_PREFIX = "visrep_"

STYLE_ORIGINAL = "original"
STYLE_LINES = "lines"
STYLE_STICKS = "sticks"
STYLE_BALLSTICK = "ballstick"
STYLE_SPHERES = "spheres"
STYLE_CARTOON = "cartoon"
STYLE_SURFACE = "surface"

STYLE_ORDER: Tuple[str, ...] = (
    STYLE_ORIGINAL,
    STYLE_LINES,
    STYLE_STICKS,
    STYLE_BALLSTICK,
    STYLE_SPHERES,
    STYLE_CARTOON,
    STYLE_SURFACE,
)

STYLE_LABELS: Dict[str, str] = {
    STYLE_ORIGINAL: "当前显示（原始）",
    STYLE_LINES: "线框",
    STYLE_STICKS: "棍状（sticks）",
    STYLE_BALLSTICK: "球棍（ball-and-stick）",
    STYLE_SPHERES: "球体",
    STYLE_CARTOON: "卡通",
    STYLE_SURFACE: "表面",
}

COLOR_ORIGINAL = "original"
COLOR_SOLID = "solid"
COLOR_BYRESN = "byresn"
COLOR_BYELEMENT = "byelement"

COLOR_MODE_ORDER: Tuple[str, ...] = (
    COLOR_ORIGINAL,
    COLOR_SOLID,
    COLOR_BYRESN,
    COLOR_BYELEMENT,
)

COLOR_MODE_LABELS: Dict[str, str] = {
    COLOR_ORIGINAL: "原色",
    COLOR_SOLID: "单色",
    COLOR_BYRESN: "按残基类型",
    COLOR_BYELEMENT: "按元素",
}

# 按残基类型固定调色（同 resn 同色；核酸统一橙色）
RESN_PALETTE: Dict[str, str] = {
    "AL31": "0xcc6677",
    "AL3": "0xcc6677",
    "AL15": "0x8856a7",
    "AL1": "0x8856a7",
    "DPPC": "0x4477aa",
    "DPP": "0x4477aa",
    "CHOL": "0xb89939",
    "CHO": "0xb89939",
    "POPC": "0x44aa99",
    "POPE": "0x117733",
    "DOPC": "0x332288",
}
NUCLEIC_COLOR = "0xe69f00"
PROTEIN_COLOR = "0x999999"
LIPID_FALLBACK = "0x56b4e9"
WATER_COLOR = "0x88ccee"
DEFAULT_SOLID = (0.80, 0.20, 0.20)

# 关键词宏
_KW_SIMPLE = frozenset({
    "all", "none", "protein", "nucleic", "lipid", "water",
    "hydrogen", "backbone", "sidechain",
})
_FIELD_ALIASES = {
    "resname": "resname",
    "resn": "resname",
    "name": "name",
    "chain": "chain",
    "segname": "segname",
    "segi": "segname",
    "resid": "resid",
    "resi": "resid",
    "index": "index",
    "serial": "serial",
    "element": "element",
}
_INJECT_RE = re.compile(r"[;\n\r`]|cmd\.|python|import\s|exec\s|eval\s", re.I)

# 标准氨基酸 / 核酸（全原子回退）
_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP",
})
_NA = frozenset({
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
    "RA", "RC", "RG", "RU", "ADE", "CYT", "GUA", "URA", "THY",
})
_WATER = frozenset({
    "HOH", "WAT", "H2O", "TIP", "TIP3", "TIP4", "T3P", "T4P", "SPC", "SOL",
    "W", "PW", "WN", "WF",
})
_BB_AA = frozenset({"N", "CA", "C", "O", "OXT", "HN", "H", "HA", "BB"})
_BB_NA = frozenset({
    "P", "O1P", "O2P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'",
    "O3'", "C2'", "C1'", "BB", "BB1", "BB2", "BB3",
})


# ---------------------------------------------------------------------------
# 异常与数据结构
# ---------------------------------------------------------------------------

class SelectionError(ValueError):
    """选择表达式错误（友好拒绝，不改场景）。"""


@dataclass
class AtomRec:
    """源对象中的单个原子/珠快照。"""

    model_idx: int
    index: int
    serial: int
    name: str
    resn: str
    resi: int
    resi_str: str
    chain: str
    segi: str
    element: str
    coord: Tuple[float, float, float]
    kind: str = "other"
    component: str = "other"


@dataclass
class RepresentationLayer:
    """单行独立表示规则。"""

    layer_id: str
    enabled: bool = True
    style: str = STYLE_ORIGINAL
    color_mode: str = COLOR_ORIGINAL
    color_rgb: Tuple[float, float, float] = DEFAULT_SOLID
    selection: str = "all"
    owned_object: str = ""
    label: str = ""


@dataclass
class PreviewResult:
    """预览：匹配原子数与错误信息。"""

    ok: bool
    count: int = 0
    message: str = ""
    indices: List[int] = field(default_factory=list)
    bond_count: int = 0
    external_bond_count: int = 0


@dataclass
class ApplyReport:
    """一次 Apply 的结果说明。"""

    ok: bool = False
    message: str = ""
    owned_objects: List[str] = field(default_factory=list)
    layer_counts: Dict[str, int] = field(default_factory=dict)
    bond_note: str = ""
    warnings: List[str] = field(default_factory=list)
    layer_bond_counts: Dict[str, int] = field(default_factory=dict)


# AST 节点类型
AstNode = Union[
    "KwNode",
    "FieldNode",
    "ResidComparisonNode",
    "NotNode",
    "BinNode",
    "WithinNode",
    "SameResidueNode",
]


@dataclass
class KwNode:
    """关键词宏节点。"""

    kw: str


@dataclass
class FieldNode:
    """字段匹配节点（可含多值与 resid 范围）。"""

    field: str
    values: List[str]
    range_lo: Optional[int] = None
    range_hi: Optional[int] = None


_RESID_COMPARISONS = {
    '<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge,
    '=': operator.eq, '==': operator.eq, '!=': operator.ne,
}


@dataclass
class ResidComparisonNode:
    comparison: str
    value: int


@dataclass
class NotNode:
    """逻辑非。"""

    child: AstNode


@dataclass
class BinNode:
    """逻辑与/或。"""

    op: str  # and / or
    left: AstNode
    right: AstNode


@dataclass
class WithinNode:
    """within X of sel。"""

    distance: float
    child: AstNode


@dataclass
class SameResidueNode:
    """same residue as sel。"""

    child: AstNode


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------

def _quote_sel(name: str) -> str:
    """选择串中的对象名加引号。"""
    return "'" + str(name).replace("'", r"\'") + "'"


def is_owned_name(name: str) -> bool:
    """判断是否为本模块生成的显示副本。"""
    return str(name).startswith(OWNED_PREFIX)


def _norm(s: str) -> str:
    """规范化大小写。"""
    return (s or "").strip().upper()


def validate_style(style: str) -> str:
    """校验显示样式。"""
    key = (style or "").strip().lower()
    aliases = {
        "当前显示": STYLE_ORIGINAL,
        "原始": STYLE_ORIGINAL,
        "线框": STYLE_LINES,
        "棒状": STYLE_STICKS,
        "棍状": STYLE_STICKS,
        "棍状（sticks）": STYLE_STICKS,
        "球棍": STYLE_BALLSTICK,
        "球棍（ball-and-stick）": STYLE_BALLSTICK,
        "球体": STYLE_SPHERES,
        "卡通": STYLE_CARTOON,
        "表面": STYLE_SURFACE,
        "ball_stick": STYLE_BALLSTICK,
    }
    key = aliases.get(key, key)
    if key not in STYLE_LABELS:
        raise ValueError("未知显示方式“%s”。" % style)
    return key


def validate_color_mode(mode: str) -> str:
    """校验颜色模式。"""
    key = (mode or "").strip().lower()
    aliases = {
        "原色": COLOR_ORIGINAL,
        "单色": COLOR_SOLID,
        "按残基类型": COLOR_BYRESN,
        "按元素": COLOR_BYELEMENT,
        "element": COLOR_BYELEMENT,
        "resn": COLOR_BYRESN,
    }
    key = aliases.get(key, key)
    if key not in COLOR_MODE_LABELS:
        raise ValueError("未知颜色模式“%s”。" % mode)
    return key


def validate_rgb(rgb: Sequence[float]) -> Tuple[float, float, float]:
    """校验 RGB 分量在 0..1。"""
    if rgb is None or len(rgb) < 3:
        raise ValueError("RGB 需要三个 0..1 分量。")
    out = []
    for v in rgb[:3]:
        f = float(v)
        if f < 0.0 or f > 1.0:
            raise ValueError("RGB 分量须在 0..1，收到 %s。" % f)
        out.append(f)
    return (out[0], out[1], out[2])


def rgb_to_hex_color(rgb: Sequence[float]) -> str:
    """将 0..1 RGB 转为 PyMOL truecolor 字符串 0xRRGGBB。"""
    r, g, b = validate_rgb(rgb)
    return "0x%02x%02x%02x" % (
        int(round(r * 255)),
        int(round(g * 255)),
        int(round(b * 255)),
    )


def new_layer_id() -> str:
    """生成层唯一 ID。"""
    return uuid.uuid4().hex[:12]


def list_styles() -> List[Dict[str, str]]:
    """样式元数据列表。"""
    return [{"id": k, "label": STYLE_LABELS[k]} for k in STYLE_ORDER]


def list_color_modes() -> List[Dict[str, str]]:
    """颜色模式元数据列表。"""
    return [{"id": k, "label": COLOR_MODE_LABELS[k]} for k in COLOR_MODE_ORDER]


# ---------------------------------------------------------------------------
# 词法 / 语法分析（VMD 子集）
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(
    r"""
    \s+
    |(?P<quoted>"[^"\n\r]*"|'[^'\n\r]*')
    |(?P<cmp><=|>=|==|!=|[<>=])
    |(?P<op>\(|\)|and\b|or\b|not\b|to\b|of\b|within\b|same\b|residue\b|as\b)
    |(?P<field>(?:resname|resn|name|chain|segname|segi|resid|resi|index|serial|element)\b)
    |(?P<kw>(?:all|none|protein|nucleic|lipid|water|hydrogen|backbone|sidechain)\b)
    |(?P<num>-?\d+(?:\.\d+)?)(?![A-Za-z0-9_.*?'])
    |(?P<word>[A-Za-z0-9_.*?+'\-]+)
    |(?P<bad>.)
    """,
    re.IGNORECASE | re.VERBOSE,
)


@dataclass
class _Tok:
    """词法单元。"""

    kind: str
    value: str
    pos: int


def tokenize(expr: str) -> List[_Tok]:
    """
    将选择表达式拆为词法单元。

    设计思路：拒绝分号/注入字符与未知符号；不把输入交给 Tcl/Python。
    """
    if expr is None:
        raise SelectionError("选择表达式为空。")
    text = str(expr).strip()
    if not text:
        raise SelectionError("选择表达式为空。")
    if _INJECT_RE.search(text):
        raise SelectionError("表达式含非法字符或疑似命令注入，已拒绝。")
    out: List[_Tok] = []
    for m in _TOKEN_RE.finditer(text):
        if m.group(0).isspace():
            continue
        if m.group("bad"):
            raise SelectionError(
                "未知符号“%s”（位置 %d）。本功能仅支持 VMD 常用子集，不会当作 PyMOL 命令执行。"
                % (m.group("bad"), m.start())
            )
        if m.group("quoted"):
            out.append(_Tok("word", m.group("quoted")[1:-1], m.start()))
        elif m.group("cmp"):
            out.append(_Tok("cmp", m.group("cmp"), m.start()))
        elif m.group("op"):
            out.append(_Tok("op", m.group("op").lower(), m.start()))
        elif m.group("field"):
            out.append(_Tok("field", m.group("field").lower(), m.start()))
        elif m.group("kw"):
            out.append(_Tok("kw", m.group("kw").lower(), m.start()))
        elif m.group("num"):
            out.append(_Tok("num", m.group("num"), m.start()))
        elif m.group("word"):
            out.append(_Tok("word", m.group("word"), m.start()))
    if not out:
        raise SelectionError("选择表达式为空。")
    return out


class _Parser:
    """递归下降：not > and > or。"""

    def __init__(self, tokens: List[_Tok]) -> None:
        self.tokens = tokens
        self.i = 0

    def _peek(self) -> Optional[_Tok]:
        if self.i < len(self.tokens):
            return self.tokens[self.i]
        return None

    def _take(self, *ok: str) -> _Tok:
        t = self._peek()
        if t is None:
            raise SelectionError("表达式不完整。")
        if ok and t.value not in ok and t.kind not in ok:
            raise SelectionError("期望 %s，却遇到“%s”。" % ("/".join(ok), t.value))
        self.i += 1
        return t

    def parse(self) -> AstNode:
        """解析完整表达式。"""
        node = self._parse_or()
        if self._peek() is not None:
            raise SelectionError("表达式末尾有多余内容“%s”。" % self._peek().value)
        return node

    def _parse_or(self) -> AstNode:
        left = self._parse_and()
        while True:
            t = self._peek()
            if t and t.kind == "op" and t.value == "or":
                self.i += 1
                right = self._parse_and()
                left = BinNode("or", left, right)
            else:
                break
        return left

    def _parse_and(self) -> AstNode:
        left = self._parse_not()
        while True:
            t = self._peek()
            # 隐式 and：两个原子项相邻
            if t and t.kind == "op" and t.value == "and":
                self.i += 1
                right = self._parse_not()
                left = BinNode("and", left, right)
            elif t and self._starts_atom(t):
                right = self._parse_not()
                left = BinNode("and", left, right)
            else:
                break
        return left

    def _starts_atom(self, t: _Tok) -> bool:
        if t.kind in ("kw", "field"):
            return True
        if t.kind == "op" and t.value in ("(", "not", "within", "same"):
            return True
        return False

    def _parse_not(self) -> AstNode:
        t = self._peek()
        if t and t.kind == "op" and t.value == "not":
            self.i += 1
            return NotNode(self._parse_not())
        return self._parse_atom()

    def _parse_atom(self) -> AstNode:
        t = self._peek()
        if t is None:
            raise SelectionError("表达式不完整。")
        if t.kind == "op" and t.value == "(":
            self.i += 1
            node = self._parse_or()
            self._take(")")
            return node
        if t.kind == "op" and t.value == "within":
            return self._parse_within()
        if t.kind == "op" and t.value == "same":
            return self._parse_same()
        if t.kind == "kw":
            self.i += 1
            return KwNode(t.value)
        if t.kind == "field":
            return self._parse_field()
        raise SelectionError("无法解析“%s”。" % t.value)

    def _parse_within(self) -> AstNode:
        self._take("within")
        n = self._take()
        if n.kind != "num":
            raise SelectionError("within 后需要距离（Å）。")
        dist = float(n.value)
        if dist < 0:
            raise SelectionError("within 距离不能为负。")
        self._take("of")
        child = self._parse_not()
        return WithinNode(dist, child)

    def _parse_same(self) -> AstNode:
        self._take("same")
        self._take("residue")
        self._take("as")
        child = self._parse_not()
        return SameResidueNode(child)

    def _parse_field(self) -> AstNode:
        ft = self._take("field")
        field = _FIELD_ALIASES[ft.value]
        term = self._peek()
        if term and term.kind == "cmp":
            if field != "resid":
                raise SelectionError("比较运算符用于 resid 编号，例如 resid < 90。%s 请直接写值，如 %s A。"
                                     % (ft.value, ft.value))
            self._take("cmp")
            value = self._peek()
            if not value or value.kind not in ("num", "word") or not re.fullmatch(r"[+-]?\d+", value.value):
                raise SelectionError("resid %s 后需要一个整数，例如 resid %s 90。" % (term.value, term.value))
            self._take()
            return ResidComparisonNode(term.value, int(value.value))
        values: List[str] = []
        range_lo: Optional[int] = None
        range_hi: Optional[int] = None

        # resid/index/serial: 支持 "1 to 10" 或一列数值
        if field in ("resid", "index", "serial"):
            first = self._peek()
            if first is None or first.kind not in ("num", "word"):
                raise SelectionError("%s 后需要数值。" % ft.value)
            self.i += 1
            lo_tok = first.value
            nxt = self._peek()
            if nxt and nxt.kind == "op" and nxt.value == "to":
                self.i += 1
                hi = self._take()
                if hi.kind not in ("num", "word"):
                    raise SelectionError("范围终点无效。")
                try:
                    lo_i = int(lo_tok)
                    hi_i = int(hi.value)
                except ValueError as exc:
                    raise SelectionError("范围必须是整数。") from exc
                if lo_i > hi_i:
                    raise SelectionError("范围颠倒：%s to %s。" % (lo_i, hi_i))
                range_lo, range_hi = lo_i, hi_i
            else:
                values.append(lo_tok)
                while True:
                    p = self._peek()
                    if p and p.kind in ("num", "word") and not self._is_keywordish(p):
                        values.append(p.value)
                        self.i += 1
                    else:
                        break
                try:
                    for value in values:
                        number = int(value)
                        if field in ("index", "serial") and number < 0:
                            raise ValueError(value)
                except ValueError as exc:
                    raise SelectionError("%s 需要整数编号（index/serial 不接受负值）。" % field) from exc
            if field in ("index", "serial") and range_lo is not None and range_lo < 0:
                raise SelectionError("%s 不接受负值范围。" % field)
            return FieldNode(field, values, range_lo, range_hi)

        # 多值：resname AL31 DPPC；name BB*
        # 关键词宏（protein 等）留给上层隐式 and，不当作字段值
        while True:
            p = self._peek()
            if p is None:
                break
            if p.kind == "kw":
                break
            if p.kind in ("word", "num"):
                values.append(p.value)
                self.i += 1
                continue
            break
        if not values:
            raise SelectionError("%s 后需要至少一个值。" % ft.value)
        return FieldNode(field, values)

    def _is_keywordish(self, p: _Tok) -> bool:
        if p.kind == "op":
            return True
        if p.kind == "field":
            return True
        if p.kind == "kw" and p.value in _KW_SIMPLE:
            return True
        return False


def parse_selection(expr: str) -> AstNode:
    """解析选择表达式为 AST；失败抛 SelectionError。"""
    return _Parser(tokenize(expr)).parse()


# ---------------------------------------------------------------------------
# 原子快照与分类
# ---------------------------------------------------------------------------

def _classify_atom(resn: str, name: str, element: str, names_in_resi: Sequence[str]) -> Tuple[str, str]:
    """
    分类单个原子所属种类。

    结合 coarse_grained.classify_residue（CG）与标准 AA/NA；
    不对任意未知残基做模糊前缀匹配。
    """
    rn = _norm(resn)
    nm = _norm(name)
    el = _norm(element) or (nm[:1] if nm else "")

    if rn in _WATER or nm in ("OW", "OH2", "W"):
        return "water", "water"

    try:
        from pymol import coarse_grained as cg

        kind, comp = cg.classify_residue(rn, list(names_in_resi) or [name])
        if kind == cg.KIND_NA:
            return "nucleic", comp
        if kind == cg.KIND_PROTEIN:
            return "protein", comp
        if kind in (cg.KIND_LIPID, cg.KIND_STEROL):
            return "lipid", comp
    except Exception:
        pass

    # 全原子回退
    if rn in _NA:
        return "nucleic", "rna"
    if rn in _AA:
        return "protein", "protein"
    if el == "H":
        return "other", "other"
    return "other", "other"


def load_atom_records(obj: str, cmd: Any = None, *, state: int = 1) -> List[AtomRec]:
    """从分子对象读取原子记录（严格限定该对象）。"""
    c = cmd or _global_cmd
    names = c.get_names("public_objects", enabled_only=0) or []
    if obj not in names:
        # 保留字对象名：仍按精确名字取
        if obj not in (c.get_names("all") or []):
            raise SelectionError("对象“%s”不存在。" % obj)
    try:
        if c.get_type(obj) != "object:molecule":
            raise SelectionError("“%s”不是分子对象。" % obj)
    except SelectionError:
        raise
    except Exception as exc:
        raise SelectionError("无法读取对象类型：%s" % exc) from exc

    model = c.get_model(_quote_sel(obj), state=state)
    # 先按残基收集珠名，供 classify_residue
    res_names: Dict[Tuple[str, str, str], List[str]] = {}
    raw: List[Tuple[int, Any]] = []
    for i, at in enumerate(model.atom):
        raw.append((i, at))
        key = (
            (getattr(at, "segi", "") or "").strip(),
            (getattr(at, "chain", "") or "").strip(),
            str(getattr(at, "resi", "") or "").strip(),
        )
        res_names.setdefault(key, []).append((getattr(at, "name", "") or "").strip())

    out: List[AtomRec] = []
    for i, at in raw:
        segi = (getattr(at, "segi", "") or "").strip()
        chain = (getattr(at, "chain", "") or "").strip()
        resi_str = str(getattr(at, "resi", "") or "").strip()
        resn = (getattr(at, "resn", "") or "").strip()
        name = (getattr(at, "name", "") or "").strip()
        element = (getattr(at, "symbol", None) or getattr(at, "elem", None) or "")
        element = str(element or "").strip()
        try:
            resi_i = int(str(getattr(at, "resi", "0")).strip())
        except Exception:
            m = re.match(r"^(-?\d+)", resi_str)
            resi_i = int(m.group(1)) if m else 0
        idx = int(getattr(at, "index", i + 1) or (i + 1))
        serial = int(getattr(at, "id", idx) or idx)
        coord = tuple(float(x) for x in at.coord[:3])
        kind, comp = _classify_atom(resn, name, element, res_names[(segi, chain, resi_str)])
        out.append(
            AtomRec(
                model_idx=i,
                index=idx,
                serial=serial,
                name=name,
                resn=resn,
                resi=resi_i,
                resi_str=resi_str,
                chain=chain,
                segi=segi,
                element=element,
                coord=coord,  # type: ignore[arg-type]
                kind=kind,
                component=comp,
            )
        )
    return out


# ---------------------------------------------------------------------------
# 求值
# ---------------------------------------------------------------------------

def _match_wild(pat: str, text: str) -> bool:
    """大小写不敏感通配（* ?）。"""
    return fnmatch.fnmatchcase(_norm(text), _norm(pat))


def _resn_aliases_for_query(q: str) -> Set[str]:
    """
    查询残基名的可接受别名集合。

    AL31 同时匹配 PDB 截断 AL3，但必须结合珠特征分类，不对其它未知残基前缀模糊匹配。
    """
    u = _norm(q)
    out = {u}
    # 用户写完整名时，允许截断三字符 PDB 名，稍后用 component 校验
    mapping = {
        "AL31": {"AL31", "AL3"},
        "AL15": {"AL15", "AL1", "AL15M"},
        "DPPC": {"DPPC", "DPP"},
        "CHOL": {"CHOL", "CHO", "CHL1"},
    }
    if u in mapping:
        out |= mapping[u]
    return out


def _resn_query_ok(atom: AtomRec, query: str) -> bool:
    """resname 匹配：精确/通配，并对截断脂名结合 classify 组件。"""
    q = query.strip()
    if any(ch in q for ch in "*?[]"):
        return _match_wild(q, atom.resn)
    aliases = _resn_aliases_for_query(q)
    rn = _norm(atom.resn)
    if rn not in {_norm(a) for a in aliases}:
        return False
    uq = _norm(q)
    # 截断名必须结合组件，避免误伤其它 AL3*
    if uq == "AL31" and rn in ("AL3", "AL31"):
        return atom.component == "al31" or rn == "AL31"
    if uq == "AL15" and rn in ("AL1", "AL15", "AL15M"):
        return atom.component == "al15" or rn.startswith("AL15")
    if uq == "DPPC" and rn in ("DPP", "DPPC"):
        return atom.component == "dppc" or rn == "DPPC"
    if uq == "CHOL" and rn in ("CHO", "CHOL", "CHL1"):
        return atom.component == "chol" or rn in ("CHOL", "CHL1")
    return True


def _kw_match(atom: AtomRec, kw: str) -> bool:
    """关键词宏。"""
    if kw == "all":
        return True
    if kw == "none":
        return False
    if kw == "protein":
        return atom.kind == "protein"
    if kw == "nucleic":
        return atom.kind == "nucleic"
    if kw == "lipid":
        return atom.kind == "lipid"
    if kw == "water":
        return atom.kind == "water"
    if kw == "hydrogen":
        return _norm(atom.element) == "H" or _norm(atom.name).startswith("H")
    if kw == "backbone":
        nm = _norm(atom.name)
        if atom.kind == "nucleic":
            return nm in _BB_NA or nm.startswith("BB")
        if atom.kind == "protein":
            return nm in _BB_AA or nm == "BB"
        return nm in _BB_AA or nm in _BB_NA or nm.startswith("BB")
    if kw == "sidechain":
        if atom.kind not in ("protein", "nucleic"):
            return False
        return not _kw_match(atom, "backbone") and not _kw_match(atom, "hydrogen")
    return False


def _field_match(atom: AtomRec, node: FieldNode) -> bool:
    """字段节点匹配。"""
    f = node.field
    if node.range_lo is not None and node.range_hi is not None:
        if f == "resid":
            return node.range_lo <= atom.resi <= node.range_hi
        if f == "index":
            # index 0 = get_model 顺序第 1 个原子 → 映射到该原子的 source.index
            # 用户写的是 0-based 模型序，不是 source.index 本身
            # 约定：index N 指 model 顺序下标 N（0 起），匹配 atoms[N].index
            return node.range_lo <= atom.model_idx <= node.range_hi
        if f == "serial":
            return node.range_lo <= atom.serial <= node.range_hi

    vals = node.values
    if f == "resname":
        return any(_resn_query_ok(atom, v) for v in vals)
    if f == "name":
        return any(_match_wild(v, atom.name) for v in vals)
    if f == "chain":
        return any(_match_wild(v, atom.chain) for v in vals)
    if f == "segname":
        return any(_match_wild(v, atom.segi) for v in vals)
    if f == "element":
        return any(_match_wild(v, atom.element) for v in vals)
    if f == "resid":
        want = set()
        for v in vals:
            try:
                want.add(int(v))
            except ValueError as exc:
                raise SelectionError("resid 值无效：%s" % v) from exc
        return atom.resi in want
    if f == "index":
        # 用户 index 为 get_model 0-based 序
        want = set()
        for v in vals:
            try:
                want.add(int(v))
            except ValueError as exc:
                raise SelectionError("index 值无效：%s" % v) from exc
        return atom.model_idx in want
    if f == "serial":
        want = set()
        for v in vals:
            try:
                want.add(int(v))
            except ValueError as exc:
                raise SelectionError("serial 值无效：%s" % v) from exc
        return atom.serial in want
    raise SelectionError("未知字段 %s。" % f)


def _eval_set(node: AstNode, atoms: Sequence[AtomRec]) -> Set[int]:
    """
    求值 AST，返回匹配的 model_idx 集合。

    within 使用网格加速，避免对大规模体系做完整 NxN。
    """
    if isinstance(node, KwNode):
        return {a.model_idx for a in atoms if _kw_match(a, node.kw)}
    if isinstance(node, FieldNode):
        return {a.model_idx for a in atoms if _field_match(a, node)}
    if isinstance(node, ResidComparisonNode):
        compare = _RESID_COMPARISONS[node.comparison]
        return {a.model_idx for a in atoms if compare(a.resi, node.value)}
    if isinstance(node, NotNode):
        all_i = {a.model_idx for a in atoms}
        return all_i - _eval_set(node.child, atoms)
    if isinstance(node, BinNode):
        L = _eval_set(node.left, atoms)
        R = _eval_set(node.right, atoms)
        if node.op == "and":
            return L & R
        return L | R
    if isinstance(node, SameResidueNode):
        base = _eval_set(node.child, atoms)
        by_model = {a.model_idx: a for a in atoms}
        keys = {
            (by_model[i].segi, by_model[i].chain, by_model[i].resi_str, _norm(by_model[i].resn))
            for i in base
            if i in by_model
        }
        return {
            a.model_idx
            for a in atoms
            if (a.segi, a.chain, a.resi_str, _norm(a.resn)) in keys
        }
    if isinstance(node, WithinNode):
        base = _eval_set(node.child, atoms)
        return _within_indices(atoms, base, node.distance)
    raise SelectionError("内部 AST 类型错误。")


def _within_indices(
    atoms: Sequence[AtomRec], base: Set[int], distance: float
) -> Set[int]:
    """网格分块 within 查询（Å），严格在当前 atoms 范围内。"""
    if not base or distance < 0:
        return set()
    cell = max(distance, 1e-6)
    grid: Dict[Tuple[int, int, int], List[int]] = {}
    by_model = {a.model_idx: a for a in atoms}
    for i in base:
        if i not in by_model:
            continue
        x, y, z = by_model[i].coord
        key = (int(math.floor(x / cell)), int(math.floor(y / cell)), int(math.floor(z / cell)))
        grid.setdefault(key, []).append(i)

    d2 = distance * distance
    hit: Set[int] = set()
    # 自身也算 within
    hit |= base & by_model.keys()
    for a in atoms:
        x, y, z = a.coord
        ix, iy, iz = (
            int(math.floor(x / cell)),
            int(math.floor(y / cell)),
            int(math.floor(z / cell)),
        )
        found = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((ix + dx, iy + dy, iz + dz), ()):
                        bx, by, bz = by_model[j].coord
                        dd = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
                        if dd <= d2:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if found:
            hit.add(a.model_idx)
    return hit


def evaluate_selection(
    expr: str, atoms: Sequence[AtomRec]
) -> Tuple[List[int], List[int]]:
    """
    求值表达式。

    返回 (source_index 列表, model_idx 列表)。
    index 用户语法按 model 序；返回给 select_list 的是 atom.index。
    """
    ast = parse_selection(expr)
    mids = sorted(_eval_set(ast, atoms))
    by_model = {a.model_idx: a for a in atoms}
    indices = [by_model[i].index for i in mids]
    return indices, mids


def preview_selection(obj: str, expr: str, cmd: Any = None) -> PreviewResult:
    """预览选择：不改场景。"""
    c = cmd or _global_cmd
    try:
        atoms = load_atom_records(obj, c)
        indices, mids = evaluate_selection(expr, atoms)
        model = c.get_model(_quote_sel(obj))
        internal, external = selection_bond_counts(mids, [b.index for b in model.bond])
        return PreviewResult(True, len(indices),
            "匹配 %d 个原子/珠；源内键数 %d。" % (len(indices), internal),
            indices, internal, external)
    except SelectionError as exc:
        return PreviewResult(False, 0, str(exc))
    except Exception as exc:
        return PreviewResult(False, 0, "预览失败：%s" % exc)


# ---------------------------------------------------------------------------
# 着色与样式
# ---------------------------------------------------------------------------

def _resn_color(resn: str, kind: str, component: str) -> str:
    """按残基类型取色。"""
    rn = _norm(resn)
    if kind == "nucleic":
        return NUCLEIC_COLOR
    if kind == "protein":
        return PROTEIN_COLOR
    if kind == "water":
        return WATER_COLOR
    if rn in RESN_PALETTE:
        return RESN_PALETTE[rn]
    if component in ("al31", "al15", "dppc", "chol"):
        return {
            "al31": RESN_PALETTE["AL31"],
            "al15": RESN_PALETTE["AL15"],
            "dppc": RESN_PALETTE["DPPC"],
            "chol": RESN_PALETTE["CHOL"],
        }[component]
    if kind == "lipid":
        return LIPID_FALLBACK
    # 确定性哈希色（不依赖 PYTHONHASHSEED），保证同 resn 同色
    h = 0
    for i, ch in enumerate(rn):
        h = (h * 131 + ord(ch) * (i + 1)) & 0xFFFFFF
    return "0x%06x" % h


def selection_bond_counts(model_indices, pairs):
    """Count only bonds wholly inside this rule; report boundary bonds separately."""
    keep = set(model_indices)
    internal = external = 0
    for i, j in pairs:
        if i in keep and j in keep:
            internal += 1
        elif i in keep or j in keep:
            external += 1
    return internal, external


def _cg_nucleic_indices(atoms: Sequence[AtomRec]) -> Set[int]:
    from pymol.cg_cartoon import BACKBONE_NAMES

    def residue(atom):
        return atom.segi, atom.chain, atom.resi_str, atom.resn

    cg_residues = {residue(atom) for atom in atoms
                   if atom.kind == "nucleic" and _norm(atom.name) in BACKBONE_NAMES}
    return {atom.index for atom in atoms if residue(atom) in cg_residues}


def _apply_cartoon(cmd: Any, name: str, warning_sink: List[str]) -> List[List[float]]:
    from pymol.cg_cartoon import display_states

    sel = _quote_sel(name)
    atoms = load_atom_records(name, cmd)
    cg_indices = _cg_nucleic_indices(atoms)
    remaining = [atom for atom in atoms if atom.index not in cg_indices]
    if remaining:
        has_guide = any(
            (atom.kind == "protein" and _norm(atom.name) == "CA") or
            (atom.kind == "nucleic" and _norm(atom.name) in ("C4'", "C4*", "C3'", "C3*"))
            for atom in remaining
        )
        cmd.show("cartoon" if has_guide else "lines", sel)
        if not has_guide:
            warning_sink.append("其余选中原子缺少卡通所需的全原子骨架，使用连接线框显示。")
    if not cg_indices:
        return []

    # A mixed layer may show native cartoon and CG geometry side by side.
    if remaining:
        temporary = cmd.get_unused_name("_visrep_cg_nucleic")
        try:
            cmd.select_list(temporary, name, sorted(cg_indices), mode="index")
            cmd.hide("everything", temporary)
        finally:
            cmd.delete(temporary)
    states, has_backbone_bond = display_states(cmd, sel, cg_indices)
    if not has_backbone_bond:
        warning_sink.append("粗粒化核酸缺少所选范围内的连续骨架键，仅标出单独的骨架珠。")
    return states


def _apply_style(cmd: Any, name: str, style: str, warning_sink: List[str]) -> List[List[float]]:
    """在显示副本上设置表示。"""
    sel = _quote_sel(name)
    st = validate_style(style)
    if st == STYLE_ORIGINAL:
        return []
    if st in (STYLE_STICKS, STYLE_BALLSTICK):
        # Clear inherited invisible-radius/transparency overrides on this copy.
        cmd.unset_deep("stick_radius stick_transparency stick_ball", name, updates=0)
        cmd.unset_bond("stick_radius", sel, sel)
        cmd.unset_bond("stick_transparency", sel, sel)
        cmd.set("stick_radius", 0.18, name)
        cmd.set("stick_transparency", 0.0, name)
        cmd.set("stick_ball", 0, name)
    cmd.hide("everything", sel)
    if st == STYLE_LINES:
        cmd.show("lines", sel)
    elif st == STYLE_STICKS:
        cmd.show("sticks", sel)
    elif st == STYLE_BALLSTICK:
        cmd.unset_deep("sphere_scale sphere_transparency", name, updates=0)
        cmd.show("sticks", sel)
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 0.25, name)
        cmd.set("sphere_transparency", 0.0, name)
    elif st == STYLE_SPHERES:
        cmd.show("spheres", sel)
    elif st == STYLE_CARTOON:
        return _apply_cartoon(cmd, name, warning_sink)
    elif st == STYLE_SURFACE:
        cmd.show("surface", sel)
    return []


def _apply_color(
    cmd: Any,
    name: str,
    mode: str,
    rgb: Tuple[float, float, float],
    atoms: Optional[Sequence[AtomRec]] = None,
    *,
    selection: Optional[str] = None,
) -> None:
    """给整个副本或来源对象内的一个选区着色。"""
    sel = selection if selection is not None else _quote_sel(name)
    m = validate_color_mode(mode)
    if m == COLOR_ORIGINAL:
        return
    if m == COLOR_SOLID:
        cmd.color(rgb_to_hex_color(rgb), sel)
        return
    if m == COLOR_BYELEMENT:
        try:
            cmd.util.cbag(sel)
        except Exception:
            try:
                from pymol import util as _util

                _util.cbag(sel)
            except Exception:
                cmd.color("atomic", sel)
        return
    # byresn
    recs = list(atoms) if atoms is not None else load_atom_records(name, cmd)
    groups: Dict[str, List[int]] = {}
    for a in recs:
        col = _resn_color(a.resn, a.kind, a.component)
        groups.setdefault(col, []).append(a.index)
    for col, idxs in groups.items():
        if not idxs:
            continue
        # 分块避免超长选择串
        for k in range(0, len(idxs), 800):
            chunk = idxs[k : k + 800]
            cmd.color(
                col,
                "(%s and index %s)" % (sel, "+".join(str(i) for i in chunk)),
            )


_SOURCE_REPS = {
    STYLE_LINES: ("lines",),
    STYLE_STICKS: ("sticks",),
    STYLE_BALLSTICK: ("sticks", "spheres"),
    STYLE_SPHERES: ("spheres",),
    STYLE_CARTOON: ("cartoon",),
    STYLE_SURFACE: ("surface",),
}


def _apply_source_style(cmd, selection, style):
    if style == STYLE_ORIGINAL:
        return
    cmd.hide("everything", selection)
    for rep in _SOURCE_REPS[style]:
        cmd.show(rep, selection)
    if style == STYLE_BALLSTICK:
        cmd.set("stick_ball", 0, selection)
        cmd.set("sphere_scale", 0.25, selection)
    elif style == STYLE_SPHERES:
        cmd.set("sphere_scale", 1.0, selection)


def _use_atom_colors(cmd, selection):
    """Let the new atom colors show through any fixed representation colors."""
    for rep, setting in (("cartoon", "cartoon_color"), ("surface", "surface_color"),
                         ("spheres", "sphere_color"), ("ribbon", "ribbon_color"),
                         ("mesh", "mesh_color"), ("dots", "dot_color")):
        visible = f"({selection} and rep {rep})"
        if cmd.count_atoms(visible):
            cmd.set(setting, -1, visible)
    for rep, setting in (("sticks", "stick_color"), ("lines", "line_color")):
        visible = f"({selection} and rep {rep})"
        if cmd.count_atoms(visible):
            cmd.set_bond(setting, -1, visible, visible)


def _style_needs_bonds(style: str) -> bool:
    """该样式是否依赖化学键。"""
    return validate_style(style) in (STYLE_STICKS, STYLE_BALLSTICK, STYLE_LINES, STYLE_CARTOON)


def _object_bond_count(cmd: Any, obj: str) -> int:
    """源对象当前键数。"""
    try:
        model = cmd.get_model(_quote_sel(obj))
        return len(getattr(model, "bond", []) or [])
    except Exception:
        return 0


def _build_cg_pairs_for_object(
    cmd: Any, obj: str
) -> Tuple[List[Tuple[int, int]], str]:
    """
    为无键 CG 源构建显示用连接（0-based model 下标对）。

    只复用已保存且完整匹配的模板；其余保留源连接，不推断新键。
    """
    try:
        from pymol import coarse_grained as cg
        from pymol import cg_templates as cgt

        det = cg.detect_cg(obj, cmd)
        if det.empty:
            return [], ""
        match = cgt.match_templates(det.residues)
        # PDB distance guessing can create a few incidental bonds in a CG model.
        # Complete saved templates still take precedence in the display copy.
        if not match.matched_beads:
            return [], ""
        matched_indices = {bead.model_idx for res in det.residues
                           if res.key in match.matched_keys for bead in res.beads}
        model = cmd.get_model(_quote_sel(obj))
        template_pairs, skipped = cg.filter_bonds_by_distance(
            match.pairs, [at.coord for at in model.atom], cg.PBC_REJECT)
        existing = {tuple(sorted(b.index)) for b in model.bond
                    if not set(b.index).issubset(matched_indices)}
        pairs = sorted(existing | set(template_pairs))
        note = "连接来源：已保存 CG 模板；其余保留源连接，未猜测新键。"
        if skipped:
            note += " 已跳过 %d 处超长连接。" % skipped
        return pairs, note
    except Exception as exc:
        return [], "CG 补键失败：%s" % exc


def _copy_bonded_subset(
    cmd: Any,
    obj: str,
    model_indices: Sequence[int],
    pairs: Sequence[Tuple[int, int]],
    new_name: str,
) -> None:
    """复制选中珠及内部连接，保留所有坐标状态和原子颜色。"""
    from chempy import Bond

    model = cmd.get_model(_quote_sel(obj), state=1)
    keep = set(int(i) for i in model_indices)
    old_to_new = {}
    from chempy.models import Indexed

    out = Indexed()
    out.connect_mode = 1  # Explicit connectivity only; never infer new bonds.
    source_orders = {tuple(sorted(b.index)): b.order for b in model.bond}
    for i, at in enumerate(model.atom):
        if i in keep:
            old_to_new[i] = len(out.atom)
            out.add_atom(at)
    out.bond = []
    for i, j in pairs:
        if i in old_to_new and j in old_to_new:
            b = Bond()
            b.index = [old_to_new[i], old_to_new[j]]
            b.order = source_orders.get(tuple(sorted((i, j))), 1)
            out.add_bond(b)
    if new_name in (cmd.get_names("all") or []):
        cmd.delete(new_name)
    cmd.load_model(out, new_name, zoom=0, finish=0)
    for state in range(2, cmd.count_states(_quote_sel(obj)) + 1):
        coords = cmd.get_coords(_quote_sel(obj), state=state)
        if coords is None or len(coords) != len(model.atom):
            raise SelectionError("多帧 CG 对象的原子数不一致，已取消显示复制。")
        cmd.load_coordset(coords[sorted(keep)], new_name, state=state)
    colors = {}
    cmd.iterate(_quote_sel(obj), "colors[index] = color", space={"colors": colors})
    groups = {}
    for old_idx, new_idx in old_to_new.items():
        color = colors.get(model.atom[old_idx].index)
        if color is not None:
            groups.setdefault(color, []).append(str(new_idx + 1))
    for color, indices in groups.items():
        cmd.color(str(color), "(%s and index %s)" % (
            _quote_sel(new_name), "+".join(indices)))


# ---------------------------------------------------------------------------
# 控制器
# ---------------------------------------------------------------------------

def list_candidate_objects(cmd: Any = None, exclude_owned: bool = True) -> List[str]:
    """
    对象下拉候选项。

    含当前可见分子、禁用的原始源、CG 连接副本（cgvis_*）；
    默认排除本模块 visrep_* 副本。
    """
    c = cmd or _global_cmd
    out: List[str] = []
    try:
        names = c.get_names("public_objects", enabled_only=0) or []
    except Exception:
        return out
    for n in names:
        if exclude_owned and is_owned_name(n):
            continue
        try:
            if c.get_type(n) == "object:molecule":
                out.append(n)
        except Exception:
            continue
    return out


class RepresentationController:
    """
    多层独立图形表示控制器。

    每层真正独立显示副本；源名字不覆写/删除；仅删除 OWNED_PREFIX 对象。
    """

    def __init__(self, cmd: Any = None, instance_prefix: Optional[str] = None) -> None:
        self.cmd = cmd or _global_cmd
        self.instance_prefix = instance_prefix or ("vr" + uuid.uuid4().hex[:6])
        self.source: Optional[str] = None
        self._source_was_enabled: bool = True
        self.layers: List[RepresentationLayer] = []
        self.owned: List[str] = []
        self.related_hide: List[str] = []  # 接管时暂时隐藏的相关对象（如 CG guide）
        self._related_was_enabled: Dict[str, bool] = {}
        self.last_report: Optional[ApplyReport] = None
        self._sources_were_enabled: Dict[str, bool] = {}
        self.object_layers: Dict[str, str] = {}
        self._ensure_default_layer()

    def _ensure_default_layer(self) -> None:
        """窗口初始一行：all / 原色 / 当前显示。"""
        if not self.layers:
            self.layers.append(
                RepresentationLayer(
                    layer_id=new_layer_id(),
                    enabled=True,
                    style=STYLE_ORIGINAL,
                    color_mode=COLOR_ORIGINAL,
                    selection="all",
                    label="原始表示",
                )
            )

    def _owned_basename(self, src: str, layer_id: str) -> str:
        """生成自有对象名候选。"""
        safe = re.sub(r"[^\w\-]+", "_", src)[:28]
        return "%s%s_%s_%s" % (OWNED_PREFIX, self.instance_prefix, safe, layer_id[:8])

    def set_related_objects(self, names: Optional[Sequence[str]]) -> None:
        """设置接管源时需暂隐的相关对象（如 cgvis_*_guide）。"""
        self.related_hide = [str(x) for x in (names or []) if x]

    def set_layers(self, layers: Sequence[RepresentationLayer]) -> None:
        """替换层定义（不立即改场景）。"""
        self.layers = list(layers)

    def add_layer(
        self,
        selection: str = "all",
        style: str = STYLE_STICKS,
        color_mode: str = COLOR_BYRESN,
        color_rgb: Tuple[float, float, float] = DEFAULT_SOLID,
        enabled: bool = True,
        label: str = "",
    ) -> RepresentationLayer:
        """新建一层规则。"""
        layer = RepresentationLayer(
            layer_id=new_layer_id(),
            enabled=enabled,
            style=validate_style(style),
            color_mode=validate_color_mode(color_mode),
            color_rgb=validate_rgb(color_rgb),
            selection=selection or "all",
            label=label,
        )
        self.layers.append(layer)
        return layer

    def remove_layer(self, layer_id: str) -> None:
        """删除一层定义（不立即改场景；Apply 后生效）。"""
        self.layers = [L for L in self.layers if L.layer_id != layer_id]

    def suggest_component_layers(self, obj: str, *, include_disabled=False) -> List[RepresentationLayer]:
        """
        一键按组分分层：各脂质残基类型、核酸整体、蛋白整体；排除水/离子。
        """
        atoms = [a for source in self._resolve_sources(obj, include_disabled=include_disabled)
                 for a in load_atom_records(source, self.cmd)]
        lipids: Dict[str, str] = {}
        has_na = False
        has_prot = False
        for a in atoms:
            if a.kind == "lipid":
                # 优先完整常见名
                key = _norm(a.resn)
                if a.component == "al31":
                    key = "AL31"
                elif a.component == "al15":
                    key = "AL15"
                elif a.component == "dppc":
                    key = "DPPC"
                elif a.component == "chol":
                    key = "CHOL"
                lipids[key] = a.resn
            elif a.kind == "nucleic":
                has_na = True
            elif a.kind == "protein":
                has_prot = True
        layers: List[RepresentationLayer] = []
        palette_cycle = [
            (0.80, 0.40, 0.47),
            (0.27, 0.47, 0.67),
            (0.53, 0.34, 0.66),
            (0.72, 0.60, 0.22),
            (0.27, 0.67, 0.60),
        ]
        i = 0
        for key in sorted(lipids.keys()):
            rgb = palette_cycle[i % len(palette_cycle)]
            i += 1
            layers.append(
                RepresentationLayer(
                    layer_id=new_layer_id(),
                    enabled=True,
                    style=STYLE_STICKS,
                    color_mode=COLOR_SOLID,
                    color_rgb=rgb,
                    selection="resname %s" % key,
                    label=key,
                )
            )
        if has_na:
            layers.append(
                RepresentationLayer(
                    layer_id=new_layer_id(),
                    enabled=True,
                    style=STYLE_STICKS,
                    color_mode=COLOR_SOLID,
                    color_rgb=(0.90, 0.62, 0.00),
                    selection="nucleic",
                    label="核酸",
                )
            )
        if has_prot:
            layers.append(
                RepresentationLayer(
                    layer_id=new_layer_id(),
                    enabled=True,
                    style=STYLE_STICKS,
                    color_mode=COLOR_SOLID,
                    color_rgb=(0.60, 0.60, 0.60),
                    selection="protein",
                    label="蛋白",
                )
            )
        if not layers:
            layers.append(
                RepresentationLayer(
                    layer_id=new_layer_id(),
                    selection="all",
                    style=STYLE_ORIGINAL,
                    color_mode=COLOR_ORIGINAL,
                    label="原始表示",
                )
            )
        return layers

    def preview(self, obj: str, expr: str, *, include_disabled=False) -> PreviewResult:
        """预览表达式匹配数。"""
        try:
            sources = self._resolve_sources(obj, include_disabled=include_disabled)
            mask = None
            if obj in self.cmd.get_names("selections"):
                mask = set(self.cmd.index(_quote_sel(obj)))
            count = 0
            bonds = external_bonds = 0
            for source in sources:
                atoms = load_atom_records(source, self.cmd)
                if mask is not None:
                    atoms = [a for a in atoms if (source, a.index) in mask]
                indices, mids = evaluate_selection(expr, atoms)
                count += len(indices)
                model = self.cmd.get_model(_quote_sel(source))
                internal, external = selection_bond_counts(mids, [b.index for b in model.bond])
                bonds += internal
                external_bonds += external
            return PreviewResult(True, count,
                "匹配 %s 个原子/珠；源内键数 %s。" % (format(count, ","), format(bonds, ",")),
                bond_count=bonds, external_bond_count=external_bonds)
        except Exception as exc:
            return PreviewResult(False, 0, str(exc))

    def _resolve_sources(self, target: str, *, include_disabled=False) -> List[str]:
        """将 all / 对象名解析为分子对象列表。"""
        c = self.cmd
        t = (target or "all").strip()
        if t.lower() == "all":
            if include_disabled:
                return list_candidate_objects(c)
            names = c.get_names("public_objects", enabled_only=1) or []
            out = []
            for n in names:
                if is_owned_name(n):
                    continue
                try:
                    if c.get_type(n) == "object:molecule":
                        out.append(n)
                except Exception:
                    continue
            if self.source == "all":
                for n in self._sources_were_enabled:
                    if n not in out and n in c.get_names("public_objects", enabled_only=0):
                        out.append(n)
            return out
        # 精确名字（即使叫 all/protein）
        all_names = c.get_names("public_objects", enabled_only=0) or []
        if t in all_names:
            try:
                if c.get_type(t) == "object:molecule":
                    return [t]
            except Exception:
                pass
        # 选择名：收集涉及的对象
        try:
            models = c.get_object_list("(" + t + ")") or []
            return [m for m in models if not is_owned_name(m)]
        except Exception:
            pass
        raise SelectionError("无法解析目标“%s”。" % target)

    def _tmp_sel_name(self) -> str:
        return self.cmd.get_unused_name("_%s_sel" % self.instance_prefix)

    def _hide_related(self) -> None:
        c = self.cmd
        enabled = set(c.get_names("public_objects", enabled_only=1) or [])
        for n in self.related_hide:
            if n in (c.get_names("all") or []):
                self._related_was_enabled.setdefault(n, n in enabled)
                try:
                    c.disable(n)
                except Exception:
                    pass

    def _restore_related(self) -> None:
        c = self.cmd
        for n, was in list(self._related_was_enabled.items()):
            if n not in (c.get_names("all") or []):
                continue
            try:
                (c.enable if was else c.disable)(n)
            except Exception:
                pass
        self._related_was_enabled = {}

    def _clean_owned(self, names: Optional[Sequence[str]] = None) -> None:
        targets = list(names) if names is not None else list(self.owned)
        for n in targets:
            if is_owned_name(n):
                try:
                    if n in (self.cmd.get_names("all") or []):
                        self.cmd.delete(n)
                except Exception:
                    pass
        if names is None:
            self.owned = []
        else:
            drop = set(targets)
            self.owned = [x for x in self.owned if x not in drop]

    def restore(self) -> str:
        """清理自有副本并恢复源与相关对象的启用状态。"""
        self._clean_owned()
        self._restore_related()
        msg = "已恢复源对象显示。"
        for obj, was_enabled in self._sources_were_enabled.items():
            if obj in self.cmd.get_names("all"):
                (self.cmd.enable if was_enabled else self.cmd.disable)(obj)
        self._sources_were_enabled = {}
        self.object_layers = {}
        for L in self.layers:
            L.owned_object = ""
        self.source = None
        self.last_report = ApplyReport(ok=True, message=msg)
        return msg

    def apply_to_sources(self, source: str, *, colors_only=False) -> ApplyReport:
        """Apply rules to matched atoms in place without creating objects."""
        c = self.cmd
        report = ApplyReport()
        jobs = []
        try:
            rules = []
            for layer in self.layers:
                if not layer.enabled:
                    continue
                mode = validate_color_mode(layer.color_mode)
                style = STYLE_ORIGINAL if colors_only else validate_style(layer.style)
                if style == STYLE_ORIGINAL and mode == COLOR_ORIGINAL:
                    continue
                rgb = validate_rgb(layer.color_rgb) if mode == COLOR_SOLID else layer.color_rgb
                rules.append((layer.selection, parse_selection(layer.selection), style, mode, rgb))
            if not rules:
                raise SelectionError("没有启用的着色规则。" if colors_only else "没有需要应用的表示或颜色改动。")

            sources = self._resolve_sources(source, include_disabled=True)
            if not sources:
                raise SelectionError("没有可用的分子对象。")
            mask = set(c.index(_quote_sel(source))) if source in c.get_names("selections") else None
            has_sticks = any(rule[2] in (STYLE_STICKS, STYLE_BALLSTICK) for rule in rules)
            matched_rules = set()
            for obj in sources:
                atoms = load_atom_records(obj, c)
                if mask is not None:
                    atoms = [atom for atom in atoms if (obj, atom.index) in mask]
                source_bonds = [bond.index for bond in c.get_model(_quote_sel(obj)).bond] if has_sticks else []
                for i, (expr, ast, style, mode, rgb) in enumerate(rules):
                    mids = _eval_set(ast, atoms)
                    if not mids:
                        continue
                    if style in (STYLE_STICKS, STYLE_BALLSTICK):
                        internal, _ = selection_bond_counts(mids, source_bonds)
                        if not internal:
                            raise SelectionError(
                                f"{obj} 的「{expr}」有 0 条内部键，无法原位显示{STYLE_LABELS[style]}。"
                                "可选球体，或在单对象的独立副本中使用已有的 CG 拓扑模板。")
                    selected = [atom for atom in atoms if atom.model_idx in mids]
                    jobs.append((obj, selected, style, mode, rgb))
                    matched_rules.add(i)
            if len(matched_rules) != len(rules):
                missing = [rule[0] for i, rule in enumerate(rules) if i not in matched_rules]
                raise SelectionError("这些规则没有匹配原子或珠：" + "；".join(missing))
        except Exception as exc:
            report.message = f"应用失败，已保留前次显示：{exc}"
            self.last_report = report
            return report

        # All rules and their bond requirements are checked before any changes.
        counts = {}
        try:
            for obj, atoms, style, mode, rgb in jobs:
                selection = self._tmp_sel_name()
                try:
                    c.select_list(selection, obj, [atom.index for atom in atoms], mode="index")
                    _apply_source_style(c, selection, style)
                    if mode != COLOR_ORIGINAL:
                        _apply_color(c, obj, mode, rgb, atoms, selection=selection)
                        _use_atom_colors(c, selection)
                finally:
                    c.delete(selection)
                counts[obj] = counts.get(obj, 0) + len(atoms)
            covered_by_copies = bool(set(counts) & self._sources_were_enabled.keys())
            report.ok = True
            action = "仅更新颜色" if colors_only else "原位更新显示与颜色"
            report.message = f"已对 {len(counts)} 个分子对象{action}，共应用 {len(matched_rules)} 条规则。"
            if covered_by_copies and not colors_only and source == self.source and source.lower() != "all":
                self.restore()
                report.message += " 已切回原对象显示。"
            elif covered_by_copies:
                report.message += " 分层副本仍在显示，可点「恢复源对象」查看原对象的更新。"
            report.layer_counts = counts
        except Exception as exc:
            report.message = f"更新对象时出错：{exc}"
        self.last_report = report
        return report

    def apply(self, source: str, layers=None) -> ApplyReport:
        """Build every requested layer before replacing the previous display."""
        c = self.cmd
        report = ApplyReport()
        if c is None or getattr(c, "_COb", None) is None:
            report.message = "PyMOL 引擎不可用。"
            self.last_report = report
            return report
        if layers is not None:
            self.set_layers(layers)
        active = [layer for layer in self.layers if layer.enabled]
        if not active:
            report.ok = True
            report.message = self.restore()
            self.last_report = report
            return report
        new_owned = []
        selections = []
        outputs = {}
        notes = []
        warnings = []
        counts = {}
        bond_counts = {}
        layer_outputs = {layer.layer_id: [] for layer in active}
        try:
            sources = self._resolve_sources(source)
            if not sources:
                raise SelectionError("没有可用的分子对象。")
            for layer in active:
                validate_style(layer.style)
                validate_color_mode(layer.color_mode)
                validate_rgb(layer.color_rgb)
                parse_selection(layer.selection)
            enabled = set(c.get_names("public_objects", enabled_only=1) or [])
            old_visibility = dict(self._sources_were_enabled)
            source_visibility = {obj: old_visibility.get(obj, obj in enabled) for obj in sources}
            mask = None
            if source in (c.get_names("selections") or []):
                mask = {}
                for obj, index in c.index(_quote_sel(source)):
                    mask.setdefault(obj, set()).add(index)
            jobs = []
            has_cartoon = any(validate_style(layer.style) == STYLE_CARTOON for layer in active)
            for obj in sources:
                atoms = load_atom_records(obj, c)
                cg_nucleic_indices = _cg_nucleic_indices(atoms) if has_cartoon else set()
                source_model = c.get_model(_quote_sel(obj))
                source_pairs = [b.index for b in source_model.bond]
                if mask is not None:
                    atoms = [atom for atom in atoms if atom.index in mask.get(obj, set())]
                pairs, note = ([], "")
                if any(_style_needs_bonds(layer.style) for layer in active):
                    pairs, note = _build_cg_pairs_for_object(c, obj)
                used_templates = False
                for layer in active:
                    indices, mids = evaluate_selection(layer.selection, atoms)
                    if indices:
                        style = validate_style(layer.style)
                        # Preserve explicit CG nucleic-acid breaks. Template filling
                        # for other styles must not reconnect a cartoon backbone.
                        display_pairs = ([] if style == STYLE_CARTOON and
                                         cg_nucleic_indices.intersection(indices) else pairs)
                        used_templates |= bool(display_pairs and _style_needs_bonds(style))
                        internal, external = selection_bond_counts(mids, display_pairs or source_pairs)
                        if style in (STYLE_STICKS, STYLE_BALLSTICK) and not internal:
                            boundary = ("另有 %d 条连接指向选择外，未扩大选择范围。" % external) if external else ""
                            raise SelectionError(
                                "对象“%s”的选择“%s”匹配 %d 个原子/珠，但有 0 条内部键，无法显示%s。"
                                "%s可改用球体，或导入含配体连接的 MOL2/SDF/已确认拓扑。"
                                "未自动添加化学键。" %
                                (obj, layer.selection, len(indices), STYLE_LABELS[style], boundary))
                        jobs.append((obj, layer, indices, mids, display_pairs, internal))
                        layer_outputs[layer.layer_id].append(True)
                if used_templates and note:
                    notes.append("%s：%s" % (obj, note))
            for layer in active:
                if not layer_outputs[layer.layer_id]:
                    raise SelectionError("表示“%s”无匹配原子：%s" % (
                        layer.label or layer.layer_id, layer.selection))
                layer_outputs[layer.layer_id] = []
            # Every expression has now been validated without changing the scene.
            for obj, layer, indices, mids, pairs, internal in jobs:
                name = c.get_unused_name(self._owned_basename(obj, layer.layer_id))
                new_owned.append(name)
                if pairs and _style_needs_bonds(layer.style):
                    _copy_bonded_subset(c, obj, mids, pairs, name)
                else:
                    temporary = self._tmp_sel_name()
                    selections.append(temporary)
                    c.select_list(temporary, obj, indices, mode="index")
                    c.create(name, temporary, source_state=0, target_state=0,
                             zoom=0, copy_properties=True)
                    c.delete(temporary)
                    selections.remove(temporary)
                _apply_color(c, name, layer.color_mode, layer.color_rgb)
                cartoon_states = _apply_style(c, name, layer.style, warnings)
                if any(cartoon_states):
                    cartoon_name = c.get_unused_name(name + "_cartoon")
                    new_owned.append(cartoon_name)
                    for state, data in enumerate(cartoon_states, 1):
                        c.load_cgo(data, cartoon_name, state=state, zoom=0)
                    outputs[cartoon_name] = layer.layer_id
                layer_outputs[layer.layer_id].append(name)
                outputs[name] = layer.layer_id
                counts[name] = len(indices)
                bond_counts[name] = internal
            # Commit after all copies and styles succeed.
            old_owned = list(self.owned)
            self._restore_related()
            for obj, was_enabled in old_visibility.items():
                if obj not in source_visibility and obj in c.get_names("all"):
                    (c.enable if was_enabled else c.disable)(obj)
            self._clean_owned(old_owned)
            self.owned = list(new_owned)
            self.object_layers = outputs
            self._sources_were_enabled = source_visibility
            self.source = source
            self._source_was_enabled = source_visibility[sources[0]]
            for layer in self.layers:
                names = layer_outputs.get(layer.layer_id, [])
                layer.owned_object = names[0] if names else ""
            for obj in sources:
                c.disable(obj)
            self._hide_related()
            report.ok = True
            report.owned_objects = list(self.owned)
            report.layer_counts = counts
            report.layer_bond_counts = bond_counts
            report.bond_note = " ".join(dict.fromkeys(notes))
            report.warnings = list(dict.fromkeys(warnings))
            report.message = "已应用 %d 条表示规则，覆盖 %d 个分子对象。" % (len(active), len(sources))
            if report.bond_note:
                report.message += " " + report.bond_note
        except Exception as exc:
            self._clean_owned(new_owned)
            report.message = "应用失败，已保留前次显示：%s" % exc
            report.warnings = warnings
        finally:
            for temporary in selections:
                c.delete(temporary)
        self.last_report = report
        return report



def color_target_selection(cmd, target):
    """Accept an exact object/selection name or the explicit global all scope."""
    if target == "all":
        return "all"
    if target not in (cmd.get_names("all") or []):
        raise SelectionError("着色目标“%s”不存在。" % target)
    return _quote_sel(target)


def color_selection_truecolor(
    cmd: Any,
    target: str,
    rgb: Sequence[float],
) -> str:
    """
    对对象/选择/all 做 truecolor 着色（不 set_color 全局名）。

    供原生 C 色板使用。
    """
    if cmd is None or getattr(cmd, "_COb", None) is None:
        raise RuntimeError("PyMOL 引擎不可用。")
    hexcol = rgb_to_hex_color(rgb)
    t = (target or "all").strip() or "all"
    cmd.color(hexcol, color_target_selection(cmd, t))
    return hexcol
