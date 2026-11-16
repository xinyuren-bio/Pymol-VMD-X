# ==================================================
# 功能说明：Martini 风格粗粒化（CG）连通显示后端——检测分组、已保存模板/近似/外部键与独立显示副本。
# 使用方法：由 pmg_qt.coarse_grained_gui 或测试导入；通过 CGDisplayController 调用 detect/apply/restore。
# 依赖环境：已安装的 PyMOL（pymol.cmd、chempy）；可选 pymol.cg_templates 自动匹配。
# 生成时间：2026-11-14 16:30
# ==================================================

"""粗粒化连通可视化纯后端（无 Qt 依赖）。

科学诚实性约定：
- 无拓扑时键图为最近邻生成树近似，仅用于显示，非力场拓扑。
- 保留原始 CG 坐标与珠名，不伪造全原子回溯或二级结构。
- 源对象坐标/键/原子数在 apply 前后不变；仅生成控制器拥有的显示副本。
"""

from __future__ import annotations

import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from pymol import cmd as _global_cmd

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

# 自有显示对象前缀（禁止再作为检测源）
OWNED_PREFIX = "cgvis_"

# 显示模式
MODE_CONTINUOUS = "continuous"  # 连续形态
MODE_BALLSTICK = "ballstick"  # 球棍连接
MODE_ORDER: Tuple[str, ...] = (MODE_CONTINUOUS, MODE_BALLSTICK)

MODE_META: Dict[str, Dict[str, str]] = {
    MODE_CONTINUOUS: {
        "id": MODE_CONTINUOUS,
        "label": "连续形态",
        "description": "细管/棍连接脂类与固醇；蛋白/核酸主链光滑引导（仅插值显示）。",
    },
    MODE_BALLSTICK: {
        "id": MODE_BALLSTICK,
        "label": "球棍连接",
        "description": "按重建 CG 键绘制较小球体与棍。",
    },
}

# 距离上限（Å）——可视化近似，非力场参数
INTRA_CAP = 8.0
PROTEIN_BB_CAP = 6.0
NA_BB_CAP = 10.0
# 外部/拓扑键的周期跨越拒绝阈值
PBC_REJECT = 15.0

PROVENANCE_APPROX = "approximate_geometry"
PROVENANCE_BONDS_FILE = "user_bond_pairs"
PROVENANCE_ITP = "topology_itp"
PROVENANCE_SAVED = "saved_template"
PROVENANCE_MIXED = "mixed_template_geometry"

PROVENANCE_LABELS = {
    PROVENANCE_APPROX: "几何近似（仅用于显示，非力场拓扑）",
    PROVENANCE_BONDS_FILE: "用户提供的全局零基键对文件",
    PROVENANCE_ITP: "原子顺序匹配的单文件拓扑",
    PROVENANCE_SAVED: "已保存拓扑连接模板（键/约束；虚拟点另示）",
    PROVENANCE_MIXED: "混合：已保存连接模板匹配 + 几何近似（仅用于显示）",
}

# 标准氨基酸（Martini 蛋白残基名）
_AA_RESN = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP",
})

# 核酸残基
_NA_RESN = frozenset({
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
    "RA", "RC", "RG", "RU", "ADE", "CYT", "GUA", "URA", "THY",
})

# 已知脂类 / 固醇残基名前缀或全名
_LIPID_EXACT = frozenset({
    "DPPC", "DPPE", "DPPS", "DPPG", "POPC", "POPE", "POPS", "POPG",
    "DOPC", "DOPE", "DLPC", "DMPC", "DSPC", "AL31",
})
_STEROL_EXACT = frozenset({"CHOL", "CHL1", "CLR", "ERG"})

# 水 / 离子 / 溶剂排除
_EXCLUDED_RESN = frozenset({
    "HOH", "WAT", "H2O", "TIP", "TIP3", "TIP4", "T3P", "T4P", "SPC", "SOL",
    "W", "PW", "WN", "WF",
    "NA", "NA+", "NAT", "CL", "CL-", "CLA", "K", "K+", "CA", "CA2", "MG",
    "MG2", "ZN", "ZN2", "ION", "ION+",
})

# 头基珠名（用于着色/标记）
_HEAD_NAMES = frozenset({
    "NC3", "PO4", "NH3", "NC4", "POH", "ROH", "LNK", "GL1", "GL2",
})

# 组件调色板（确定性）
COMPONENT_COLORS: Dict[str, Dict[str, str]] = {
    "rna": {"backbone": "tv_orange", "base": "teal", "default": "tv_orange"},
    "dna": {"backbone": "tv_orange", "base": "teal", "default": "tv_orange"},
    "protein": {"backbone": "gray70", "sidechain": "gray50", "default": "gray60"},
    "dppc": {"head": "0x4477aa", "tail": "0x88aadd", "default": "0x4477aa"},
    "al31": {"head": "0xcc6677", "tail": "0xd99ba5", "default": "0xcc6677"},
    "al15": {"head": "0x8856a7", "tail": "0xab8ac4", "default": "0x8856a7"},
    "chol": {"head": "0xb89939", "tail": "0xddbb66", "default": "0xb89939"},
    "lipid": {"head": "cyan", "tail": "bluewhite", "default": "cyan"},
    "other": {"default": "white"},
}

# 组件类型
KIND_PROTEIN = "protein"
KIND_NA = "nucleic"
KIND_LIPID = "lipid"
KIND_STEROL = "sterol"
KIND_OTHER = "other"


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class CGBead:
    """单个 CG 珠（保留源原子索引与名称）。"""

    index: int  # PyMOL 1-based index（源对象内）
    name: str
    resn: str
    resi: str
    chain: str
    segi: str
    coord: Tuple[float, float, float]
    model_idx: int = 0  # get_model 中 0-based 下标


@dataclass
class CGResidue:
    """按 model/chain/segi/resi/resn 分组的 CG 残基。"""

    object: str
    segi: str
    chain: str
    resi: str
    resn: str
    kind: str
    beads: List[CGBead] = field(default_factory=list)
    component: str = "other"  # 调色板键

    @property
    def key(self) -> Tuple[str, str, str, str, str]:
        """残基唯一键。"""
        return (self.object, self.segi, self.chain, self.resi, self.resn)


@dataclass
class DetectionReport:
    """检测结果摘要（供 GUI 中文汇报）。"""

    object: str = ""
    bead_count: int = 0
    residue_count: int = 0
    by_kind: Dict[str, int] = field(default_factory=dict)
    by_component: Dict[str, int] = field(default_factory=dict)
    by_resn: Dict[str, int] = field(default_factory=dict)
    excluded_beads: int = 0
    note: str = ""
    residues: List[CGResidue] = field(default_factory=list)
    empty: bool = False


@dataclass
class ApplyReport:
    """一次 Apply 的结果与局限说明。"""

    ok: bool = False
    mode: str = MODE_CONTINUOUS
    provenance: str = PROVENANCE_APPROX
    provenance_label: str = PROVENANCE_LABELS[PROVENANCE_APPROX]
    owned_objects: List[str] = field(default_factory=list)
    bond_count: int = 0
    skipped_long: int = 0
    skipped_inter: int = 0
    message: str = ""
    detection: Optional[DetectionReport] = None
    matched_beads: int = 0
    approximate_beads: int = 0
    template_counts: Dict[str, int] = field(default_factory=dict)
    template_note: str = ""
    isolated_indices: List[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def validate_mode(mode: str) -> str:
    """校验显示模式标识。"""
    key = (mode or "").strip().lower()
    if key in ("连续形态", "continuous", "tube", "smooth"):
        return MODE_CONTINUOUS
    if key in ("球棍连接", "ballstick", "ball_stick", "sticks"):
        return MODE_BALLSTICK
    if key not in MODE_META:
        raise ValueError(
            "未知显示模式“%s”，可选：%s。"
            % (mode, "、".join(MODE_ORDER))
        )
    return key


def list_modes() -> List[Dict[str, str]]:
    """返回模式元数据列表。"""
    return [dict(MODE_META[k]) for k in MODE_ORDER]


def _quote_sel(name: str) -> str:
    """选择串中的对象名加引号。"""
    return "'" + str(name).replace("'", r"\'") + "'"


def is_owned_name(name: str) -> bool:
    """判断是否为本模块生成的显示对象。"""
    return str(name).startswith(OWNED_PREFIX)


def _norm_name(n: str) -> str:
    """规范化珠名（去空白、大写）。"""
    return (n or "").strip().upper()


def _norm_resn(r: str) -> str:
    """规范化残基名。"""
    return (r or "").strip().upper()


def _resi_int(r: str) -> Optional[int]:
    """残基号转整数；失败返回 None。"""
    try:
        return int(str(r).strip())
    except (TypeError, ValueError):
        m = re.match(r"^(-?\d+)", str(r).strip())
        if m:
            return int(m.group(1))
    return None


def _dist2(a: Sequence[float], b: Sequence[float]) -> float:
    """平方距离。"""
    dx = float(a[0]) - float(b[0])
    dy = float(a[1]) - float(b[1])
    dz = float(a[2]) - float(b[2])
    return dx * dx + dy * dy + dz * dz


def _dist(a: Sequence[float], b: Sequence[float]) -> float:
    """欧氏距离。"""
    return _dist2(a, b) ** 0.5


def classify_residue(resn: str, names: Sequence[str]) -> Tuple[str, str]:
    """
    分类残基种类与调色板组件键。

    返回 (kind, component)。排除溶剂时调用方应先过滤。
    设计思路：优先精确脂/固醇名，再核酸/蛋白标准残基，再用珠名启发式识别脂尾。
    """
    rn = _norm_resn(resn)
    nm = {_norm_name(x) for x in names}

    if rn in _EXCLUDED_RESN:
        return KIND_OTHER, "other"

    # 固醇（含 CHOL→CHO 截断）
    if rn in _STEROL_EXACT or rn.startswith("CHOL") or rn == "CHO":
        return KIND_STEROL, "chol"
    if "ROH" in nm and any(
        x in nm for x in ("R1", "R2", "R3", "R4", "R5")
    ):
        if not nm.intersection({"BB", "BB1", "BB2", "BB3"}):
            return KIND_STEROL, "chol"

    # 已知脂类（含 PDB 四字符残基名截断时的珠名回退）
    if rn.startswith("AL31") or (
        rn.startswith("AL3") and (("POH" in nm and "NC4" in nm) or "K2A" in nm)
    ):
        return KIND_LIPID, "al31"
    if rn.startswith("AL15") or (
        rn.startswith("AL1")
        and (
            "LNK" in nm
            or any(x.startswith("EC") or re.match(r"^\dEC", x) for x in nm)
        )
    ):
        return KIND_LIPID, "al15"
    if rn.startswith("DPPC") or (
        rn.startswith("DPP") and "NC3" in nm and "PO4" in nm
    ):
        return KIND_LIPID, "dppc"
    if rn in _LIPID_EXACT:
        return KIND_LIPID, "lipid"

    # 核酸
    if rn in _NA_RESN and nm & {"BB1", "BB2", "BB3"}:
        comp = "dna" if rn.startswith("D") or rn in ("DT", "DA", "DC", "DG", "THY") else "rna"
        return KIND_NA, comp
    if nm & {"BB1", "BB2", "BB3"} and any(x.startswith("SC") for x in nm):
        return KIND_NA, "rna"

    # 蛋白
    if rn in _AA_RESN and ("BB" in nm or any(x.startswith("SC") for x in nm)):
        return KIND_PROTEIN, "protein"
    if "BB" in nm and any(x.startswith("SC") for x in nm):
        return KIND_PROTEIN, "protein"

    # 脂类启发式：头基 + CnA/CnB 尾；并尽量细分调色板
    has_ec = any(x.startswith("EC") or re.match(r"^\dEC\d*$", x) for x in nm)
    has_tail = any(re.match(r"^C\d+[AB]$", x) for x in nm) or has_ec
    has_head = bool(nm & _HEAD_NAMES) or "LNK" in nm
    if has_tail and (has_head or any(x.startswith("GL") for x in nm)):
        if "LNK" in nm or has_ec:
            return KIND_LIPID, "al15"
        if "POH" in nm or "K2A" in nm:
            return KIND_LIPID, "al31"
        if "NC3" in nm and "PO4" in nm:
            return KIND_LIPID, "dppc"
        return KIND_LIPID, "lipid"

    return KIND_OTHER, "other"


def list_source_objects(cmd: Any = None) -> List[str]:
    """
    列出可选源分子对象（含禁用），排除本模块生成的显示副本。
    """
    c = cmd or _global_cmd
    names = c.get_names("public_objects", enabled_only=0) or []
    out: List[str] = []
    for n in names:
        if is_owned_name(n):
            continue
        try:
            if c.get_type(n) == "object:molecule":
                out.append(n)
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------------
# 检测
# ---------------------------------------------------------------------------

def detect_cg(obj: str, cmd: Any = None) -> DetectionReport:
    """
    在指定分子对象中检测并分组 CG 残基。

    使用安全的 model/chain/segi/resi/resn 键；保留原始原子索引。
    排除水/离子；不修改源对象。
    """
    c = cmd or _global_cmd
    report = DetectionReport(object=obj)
    if not obj or obj not in (c.get_names("public_objects", enabled_only=0) or []):
        report.empty = True
        report.note = "对象“%s”不存在。" % obj
        return report
    if is_owned_name(obj):
        report.empty = True
        report.note = "不能对粗粒化显示副本再次处理。"
        return report
    try:
        if c.get_type(obj) != "object:molecule":
            report.empty = True
            report.note = "“%s”不是分子对象。" % obj
            return report
    except Exception:
        report.empty = True
        report.note = "无法读取对象类型。"
        return report

    # 一次取出模型，保证 model_idx 与全局零基下标一致
    try:
        model = c.get_model(_quote_sel(obj))
    except Exception as exc:
        report.empty = True
        report.note = "读取模型失败：%s" % exc
        return report

    buckets: Dict[Tuple[str, str, str, str], List[CGBead]] = defaultdict(list)
    excluded = 0
    for i, at in enumerate(model.atom):
        rn = _norm_resn(getattr(at, "resn", "") or "")
        name = (getattr(at, "name", "") or "").strip()
        chain = (getattr(at, "chain", "") or "").strip()
        segi = (getattr(at, "segi", "") or "").strip()
        resi = str(getattr(at, "resi", "") or "").strip()
        idx = int(getattr(at, "index", i + 1) or (i + 1))
        coord = tuple(float(x) for x in at.coord[:3])
        bead = CGBead(
            index=idx,
            name=name,
            resn=rn or (getattr(at, "resn", "") or ""),
            resi=resi,
            chain=chain,
            segi=segi,
            coord=coord,  # type: ignore[arg-type]
            model_idx=i,
        )
        if rn in _EXCLUDED_RESN:
            excluded += 1
            continue
        key = (segi, chain, resi, rn or bead.resn)
        buckets[key].append(bead)

    residues: List[CGResidue] = []
    by_kind: Dict[str, int] = defaultdict(int)
    by_comp: Dict[str, int] = defaultdict(int)
    by_resn: Dict[str, int] = defaultdict(int)

    for (segi, chain, resi, resn), beads in sorted(
        buckets.items(), key=lambda x: (x[0][1], _resi_int(x[0][2]) or 0, x[0][2], x[0][3])
    ):
        names = [b.name for b in beads]
        kind, comp = classify_residue(resn, names)
        if kind == KIND_OTHER and resn in _EXCLUDED_RESN:
            excluded += len(beads)
            continue
        # 无明确 CG 特征的孤立珠：若珠名像溶剂已排除；其余记 other 仍纳入以便透明
        res = CGResidue(
            object=obj,
            segi=segi,
            chain=chain,
            resi=resi,
            resn=resn,
            kind=kind,
            beads=beads,
            component=comp,
        )
        residues.append(res)
        by_kind[kind] += 1
        by_comp[comp] += 1
        by_resn[resn] += 1

    # 仅保留可识别 CG 或脂/固醇/核酸/蛋白；other 若完全无特征可丢掉
    kept: List[CGResidue] = []
    for r in residues:
        names = {_norm_name(b.name) for b in r.beads}
        if r.kind != KIND_OTHER:
            kept.append(r)
            continue
        # 保留带 LNK/EC/尾链启发式但仍标 other 的极少情况已在 classify 处理
        if names & _HEAD_NAMES or any(x.startswith("BB") for x in names):
            kept.append(r)

    bead_count = sum(len(r.beads) for r in kept)
    report.residues = kept
    report.bead_count = bead_count
    report.residue_count = len(kept)
    report.by_kind = dict(by_kind)
    report.by_component = dict(by_comp)
    report.by_resn = dict(sorted(by_resn.items(), key=lambda x: (-x[1], x[0])))
    report.excluded_beads = excluded
    if bead_count == 0:
        report.empty = True
        report.note = "未检测到粗粒化珠（或仅含溶剂/离子）。"
    else:
        report.note = (
            "检测到 %d 个 CG 珠、%d 个残基组分；键图在无拓扑时为显示近似。"
            % (bead_count, len(kept))
        )
    return report


# ---------------------------------------------------------------------------
# 键图构建（可视化近似）
# ---------------------------------------------------------------------------

def _mst_edges(
    beads: Sequence[CGBead], cap: float
) -> Tuple[List[Tuple[int, int]], int]:
    """
    残基内最近邻生成树（Kruskal）；边长 > cap 丢弃。

    返回 (model_idx 对列表, 无法连通的缺口数)。
    组内规模通常 ≤56，O(n² log n) 可接受。
    """
    n = len(beads)
    if n < 2:
        return [], 0
    skipped = 0
    edges: List[Tuple[float, int, int]] = []
    cap2 = cap * cap
    for i in range(n):
        for j in range(i + 1, n):
            d2 = _dist2(beads[i].coord, beads[j].coord)
            if d2 > cap2:
                skipped += 1
                continue
            edges.append((d2, i, j))
    edges.sort(key=lambda x: x[0])
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        parent[rb] = ra
        return True

    out: List[Tuple[int, int]] = []
    for _, i, j in edges:
        if union(i, j):
            a, b = beads[i].model_idx, beads[j].model_idx
            out.append((a, b) if a < b else (b, a))
            if len(out) >= n - 1:
                break
    # 统计因距离限制无法连通的最少缺口，不把所有远距离点对当作“跳过的键”。
    return out, max(0, n - 1 - len(out))


def _bead_by_name(res: CGResidue, name: str) -> Optional[CGBead]:
    """按珠名取第一个匹配珠。"""
    target = _norm_name(name)
    for b in res.beads:
        if _norm_name(b.name) == target:
            return b
    return None


def _backbone_beads_protein(res: CGResidue) -> List[CGBead]:
    """蛋白主链珠（BB）。"""
    return [b for b in res.beads if _norm_name(b.name) == "BB"]


def _backbone_path_na(res: CGResidue) -> List[CGBead]:
    """核酸主链路径 BB1→BB2→BB3（缺则按名称序回退）。"""
    path: List[CGBead] = []
    for nm in ("BB1", "BB2", "BB3"):
        b = _bead_by_name(res, nm)
        if b is not None:
            path.append(b)
    if path:
        return path
    bbs = [b for b in res.beads if _norm_name(b.name).startswith("BB")]
    bbs.sort(key=lambda x: _norm_name(x.name))
    return bbs


def build_approximate_bonds(
    residues: Sequence[CGResidue],
) -> Tuple[List[Tuple[int, int]], int, int]:
    """
    构建可视化近似键：残基内 MST + 同链连续残基主链连接。

    不产生分子间键；超距边跳过并计数（避免周期镜像拉丝）。
    返回 (键对 model_idx, skipped_long, skipped_inter占位0)。
    """
    bonds: Set[Tuple[int, int]] = set()
    skipped_long = 0

    # 残基内
    for res in residues:
        edges, sk = _mst_edges(res.beads, INTRA_CAP)
        skipped_long += sk
        for e in edges:
            bonds.add(e)
        # 核酸：强制 BB1-BB2-BB3 路径（若在 cap 内）
        if res.kind == KIND_NA:
            path = _backbone_path_na(res)
            for i in range(len(path) - 1):
                d = _dist(path[i].coord, path[i + 1].coord)
                if d <= INTRA_CAP:
                    a, b = path[i].model_idx, path[i + 1].model_idx
                    bonds.add((a, b) if a < b else (b, a))
                else:
                    skipped_long += 1

    # 同链连续残基主链
    by_chain: Dict[Tuple[str, str, str], List[CGResidue]] = defaultdict(list)
    for res in residues:
        if res.kind not in (KIND_PROTEIN, KIND_NA):
            continue
        by_chain[(res.object, res.segi, res.chain)].append(res)

    for _, group in by_chain.items():
        # 按残基号排序；仅连接号差为 1 的邻居
        sortable: List[Tuple[int, CGResidue]] = []
        for res in group:
            ri = _resi_int(res.resi)
            if ri is not None:
                sortable.append((ri, res))
        sortable.sort(key=lambda x: x[0])
        for i in range(len(sortable) - 1):
            r0, res0 = sortable[i][0], sortable[i][1]
            r1, res1 = sortable[i + 1][0], sortable[i + 1][1]
            if r1 - r0 != 1:
                continue
            if res0.kind != res1.kind:
                continue
            if res0.kind == KIND_PROTEIN:
                b0s = _backbone_beads_protein(res0)
                b1s = _backbone_beads_protein(res1)
                if not b0s or not b1s:
                    continue
                # 取各残基第一个 BB
                d = _dist(b0s[0].coord, b1s[0].coord)
                if d <= PROTEIN_BB_CAP:
                    a, b = b0s[0].model_idx, b1s[0].model_idx
                    bonds.add((a, b) if a < b else (b, a))
                else:
                    skipped_long += 1
            elif res0.kind == KIND_NA:
                p0 = _backbone_path_na(res0)
                p1 = _backbone_path_na(res1)
                if not p0 or not p1:
                    continue
                # 连续残基：前残基末 BB → 后残基首 BB
                d = _dist(p0[-1].coord, p1[0].coord)
                if d <= NA_BB_CAP:
                    a, b = p0[-1].model_idx, p1[0].model_idx
                    bonds.add((a, b) if a < b else (b, a))
                else:
                    skipped_long += 1

    return sorted(bonds), skipped_long, 0


# ---------------------------------------------------------------------------
# 外部键对 / 简易拓扑
# ---------------------------------------------------------------------------

class TopologyError(ValueError):
    """拓扑或键对文件不可用（应明确失败，不得伪装成近似）。"""


def parse_bond_pairs_file(fp: str, n_atoms: int) -> List[Tuple[int, int]]:
    """
    解析全局零基键对文本。

    格式：每行两个整数 i j（可含注释 #）；索引针对整个所选对象原子顺序
    （与 cmd.get_model 顺序一致），不是 MartiniGlass 单分子局部列表。
    """
    if not fp or not os.path.isfile(fp):
        raise TopologyError("键对文件不存在：%s" % fp)
    pairs: List[Tuple[int, int]] = []
    with open(fp, "r", encoding="utf-8", errors="replace") as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#") or s.startswith(";"):
                continue
            # 允许逗号或空白分隔
            s = s.replace(",", " ")
            parts = s.split()
            if len(parts) < 2:
                raise TopologyError("第 %d 行无法解析为键对：%s" % (ln, line.rstrip()))
            try:
                i, j = int(parts[0]), int(parts[1])
            except ValueError as exc:
                raise TopologyError("第 %d 行索引非整数。" % ln) from exc
            if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
                raise TopologyError(
                    "第 %d 行索引越界（0..%d）：%d %d。索引须对整个所选对象全局编号。"
                    % (ln, n_atoms - 1, i, j)
                )
            if i == j:
                raise TopologyError("第 %d 行自环键被拒绝。" % ln)
            pairs.append((i, j) if i < j else (j, i))
    if not pairs:
        raise TopologyError("键对文件未包含任何有效键。")
    return pairs


def _reject_gromacs_preprocessor(text: str) -> None:
    """若含 #include / 条件编译则拒绝（未实现完整预处理器）。"""
    for pat, label in (
        (r"^\s*#\s*include\b", "#include"),
        (r"^\s*#\s*if(n?def)?\b", "条件编译"),
        (r"^\s*#\s*else\b", "条件编译"),
        (r"^\s*#\s*elif\b", "条件编译"),
    ):
        if re.search(pat, text, re.MULTILINE | re.IGNORECASE):
            raise TopologyError(
                "拓扑文件含 %s：本模块未实现 GROMACS 预处理器，"
                "请改用已展开的单文件，或提供全局零基键对文本。"
                % label
            )


def parse_itp_bonds(fp: str, n_atoms: int, atom_names=None) -> List[Tuple[int, int]]:
    """
    尝试从单个 .itp/.top 读取 [bonds]/[constraints]。

    要求 [atoms] 计数与对象原子数完全一致，否则明确失败。
    不支持多分子块自动匹配或 #include。
    """
    if not fp or not os.path.isfile(fp):
        raise TopologyError("拓扑文件不存在：%s" % fp)
    with open(fp, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    _reject_gromacs_preprocessor(text)

    section = ""
    atom_count = 0
    topology_names = []
    bonds: List[Tuple[int, int]] = []
    mol_blocks = 0
    for ln, raw in enumerate(text.splitlines(), 1):
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"^\[\s*(\w+)\s*\]$", line)
        if m:
            section = m.group(1).lower()
            if section == "moleculetype":
                mol_blocks += 1
                if mol_blocks > 1:
                    raise TopologyError(
                        "拓扑含多个 [moleculetype]，无法安全匹配到整个对象；"
                        "请提供与对象原子一一对应的单分子拓扑，或全局键对文件。"
                    )
            continue
        if section == "atoms":
            parts = line.split()
            if len(parts) < 5 or not parts[0].isdigit() or int(parts[0]) != atom_count + 1:
                raise TopologyError("拓扑原子记录必须按 1..N 连续排列，且含珠名称。")
            atom_count += 1
            topology_names.append(parts[4])
        elif section in ("bonds", "constraints"):
            parts = line.split()
            if len(parts) >= 2 and parts[0].lstrip("-").isdigit():
                try:
                    # GROMACS 为 1-based
                    i, j = int(parts[0]) - 1, int(parts[1]) - 1
                except ValueError as exc:
                    raise TopologyError("第 %d 行键索引无效。" % ln) from exc
                if i < 0 or j < 0:
                    raise TopologyError("第 %d 行键索引无效。" % ln)
                bonds.append((i, j) if i < j else (j, i))

    if atom_count == 0:
        raise TopologyError("拓扑中未解析到 [atoms]，无法校验。")
    if atom_count != n_atoms:
        raise TopologyError(
            "拓扑原子数 %d 与对象原子数 %d 不一致，拒绝当作精确拓扑。"
            % (atom_count, n_atoms)
        )
    if atom_names is not None and list(atom_names) != topology_names:
        raise TopologyError("拓扑珠名称或顺序与当前对象不匹配；请按对象原子顺序提供全局键对文件。")
    if not bonds:
        raise TopologyError("拓扑中未找到 [bonds]/[constraints]。")
    # 校验键下标
    for i, j in bonds:
        if i >= n_atoms or j >= n_atoms:
            raise TopologyError("拓扑键索引越界：%d %d" % (i, j))
    return bonds


def load_connectivity(
    fp: Optional[str], n_atoms: int, atom_names=None
) -> Tuple[Optional[List[Tuple[int, int]]], str]:
    """
    加载用户连通性；无文件则返回 (None, approximate)。

    按扩展名分流：.itp/.top 走拓扑解析；其余作全局键对。
    """
    if not fp:
        return None, PROVENANCE_APPROX
    ext = os.path.splitext(fp)[1].lower()
    if ext in (".itp", ".top"):
        return parse_itp_bonds(fp, n_atoms, atom_names), PROVENANCE_ITP
    return parse_bond_pairs_file(fp, n_atoms), PROVENANCE_BONDS_FILE


def filter_bonds_by_distance(
    pairs: Sequence[Tuple[int, int]],
    coords: Sequence[Sequence[float]],
    cap: float = PBC_REJECT,
) -> Tuple[List[Tuple[int, int]], int]:
    """丢弃超长键（周期拉丝），返回保留列表与跳过数。"""
    kept: List[Tuple[int, int]] = []
    skipped = 0
    cap2 = cap * cap
    for i, j in pairs:
        if _dist2(coords[i], coords[j]) > cap2:
            skipped += 1
            continue
        kept.append((i, j) if i < j else (j, i))
    return kept, skipped


# ---------------------------------------------------------------------------
# 模型构建与着色渲染
# ---------------------------------------------------------------------------

def _copy_model_with_bonds(cmd: Any, obj: str, pairs: Sequence[Tuple[int, int]]) -> Any:
    """复制源模型、清空自动键并写入可视化键。"""
    from chempy import Bond

    model = cmd.get_model(_quote_sel(obj))
    model.bond = []
    for i, j in pairs:
        b = Bond()
        b.index = [int(i), int(j)]
        b.order = 1
        model.add_bond(b)
    return model


def _owned_name(src: str, suffix: str) -> str:
    """生成自有对象名。"""
    safe = re.sub(r"[^\w\-]+", "_", src)[:40]
    return "%s%s_%s" % (OWNED_PREFIX, safe, suffix)


def _delete_if_exists(cmd: Any, name: str) -> None:
    """删除已存在对象（仅调用方保证为自有名）。"""
    try:
        if name in (cmd.get_names("all") or []):
            cmd.delete(name)
    except Exception:
        pass


def _copy_atom_sel(name: str, model_idx: int) -> str:
    """
    显示副本中按加载顺序选择原子。

    load_model 后原子顺序与源 model 一致，index 重编号为 1..N，
    故使用 model_idx+1，而非源对象上可能不连续的原始 index。
    """
    return "(%s and index %d)" % (_quote_sel(name), int(model_idx) + 1)


def _color_components(cmd, name, residues, mode):
    """批量按组分着色，避免对上万珠逐个调用选择器。"""
    groups = defaultdict(list)
    for res in residues:
        pal = COMPONENT_COLORS.get(res.component, COMPONENT_COLORS["other"])
        for bead in res.beads:
            nm = _norm_name(bead.name)
            if res.kind == KIND_NA:
                col = pal["base"] if nm.startswith("SC") else pal["backbone"]
            elif res.kind == KIND_PROTEIN:
                col = pal["sidechain"] if nm.startswith("SC") else pal["backbone"]
            elif res.kind in (KIND_LIPID, KIND_STEROL):
                is_tail = nm.startswith("C") or nm.startswith("EC") or bool(re.match(r"^\d", nm))
                col = pal.get("tail" if is_tail else "head", pal["default"])
            else:
                col = pal.get("default", "white")
            groups[col].append(bead.model_idx + 1)
    for color, indices in groups.items():
        cmd.color(color, "(%s and index %s)" % (_quote_sel(name), "+".join(map(str, indices))))


def _apply_reps(
    cmd: Any,
    name: str,
    residues: Sequence[CGResidue],
    mode: str,
    isolated_indices: Optional[Sequence[int]] = None,
) -> None:
    """设置表示：连续形态或球棍；孤立虚拟点以小球显示（不伪造键）。"""
    cmd.hide("everything", _quote_sel(name))
    indices = [b.model_idx + 1 for r in residues for b in r.beads]
    visible = "(%s and index %s)" % (_quote_sel(name), "+".join(map(str, indices)))
    if mode == MODE_BALLSTICK:
        cmd.show("sticks", visible)
        cmd.show("spheres", visible)
        try:
            cmd.set("stick_radius", 0.18, name)
            cmd.set("sphere_scale", 0.22, name)
        except Exception:
            pass
        return

    # 连续形态：全局细棍 + 头基小球 + 模板孤立虚拟点小球
    cmd.show("sticks", visible)
    try:
        cmd.set("stick_radius", 0.25, name)
        cmd.set("stick_ball", 0, name)
    except Exception:
        pass
    sphere_parts: List[str] = []
    for res in residues:
        if res.kind not in (KIND_LIPID, KIND_STEROL):
            continue
        for bead in res.beads:
            nm = _norm_name(bead.name)
            if nm in _HEAD_NAMES or nm == "ROH":
                sphere_parts.append(str(bead.model_idx + 1))
    if isolated_indices:
        for mi in isolated_indices:
            sphere_parts.append(str(int(mi) + 1))
    if sphere_parts:
        # 去重保序
        seen: Set[str] = set()
        uniq: List[str] = []
        for p in sphere_parts:
            if p not in seen:
                seen.add(p)
                uniq.append(p)
        hs = "(%s and index %s)" % (_quote_sel(name), "+".join(uniq))
        cmd.show("spheres", hs)
        try:
            cmd.set("sphere_scale", 0.28, name)
        except Exception:
            pass


def _cgo_tube_along(
    points: Sequence[Tuple[float, float, float]],
    rgb: Tuple[float, float, float],
    radius: float = 0.35,
    distance_cap: float = NA_BB_CAP,
) -> List[float]:
    """
    沿原始珠坐标生成 CGO 圆柱管（仅显示插值，不改坐标）。

    相邻点超距则断开，避免周期拉丝。
    """
    from pymol import cgo

    out: List[float] = []
    if len(points) < 2:
        return out
    r, g, b = rgb
    for i in range(len(points) - 1):
        p0, p1 = points[i], points[i + 1]
        if _dist(p0, p1) > distance_cap:
            continue
        out.extend([
            cgo.CYLINDER,
            p0[0], p0[1], p0[2],
            p1[0], p1[1], p1[2],
            radius,
            r, g, b,
            r, g, b,
        ])
    return out


def _cgo_base_plate(sc_coords: Sequence[Tuple[float, float, float]], rgb: Tuple[float, float, float]) -> List[float]:
    """由 SC 珠几何中心附近生成简易三角板（显示用）。"""
    from pymol import cgo

    if len(sc_coords) < 3:
        return []
    # 取前三个 SC 构成三角
    a, b, c = sc_coords[0], sc_coords[1], sc_coords[2]
    r, g, bl = rgb
    return [
        cgo.BEGIN, cgo.TRIANGLES,
        cgo.COLOR, r, g, bl,
        cgo.VERTEX, a[0], a[1], a[2],
        cgo.VERTEX, b[0], b[1], b[2],
        cgo.VERTEX, c[0], c[1], c[2],
        cgo.END,
    ]


def _build_guide_cgo(
    residues: Sequence[CGResidue],
) -> List[float]:
    """为蛋白/核酸主链与碱基板构建可选 CGO 引导。"""
    obj: List[float] = []
    # 蛋白 / 核酸按链串联 BB
    by_chain: Dict[Tuple[str, str, str, str], List[CGResidue]] = defaultdict(list)
    for res in residues:
        if res.kind not in (KIND_PROTEIN, KIND_NA):
            continue
        by_chain[(res.object, res.segi, res.chain, res.kind)].append(res)

    for key, group in by_chain.items():
        kind = key[3]
        sortable: List[Tuple[int, CGResidue]] = []
        for res in group:
            ri = _resi_int(res.resi)
            if ri is not None:
                sortable.append((ri, res))
        sortable.sort(key=lambda x: x[0])
        pts = []
        previous_resi = None
        rgb = (0.55, 0.55, 0.55) if kind == KIND_PROTEIN else (1.0, 0.45, 0.1)
        radius = 0.40 if kind == KIND_NA else 0.45
        cap = NA_BB_CAP if kind == KIND_NA else PROTEIN_BB_CAP
        for ri, res in sortable:
            if previous_resi is not None and ri != previous_resi + 1:
                obj.extend(_cgo_tube_along(pts, rgb, radius, cap))
                pts = []
            previous_resi = ri
            if kind == KIND_PROTEIN:
                bbs = _backbone_beads_protein(res)
                if bbs:
                    pts.append(bbs[0].coord)
            else:
                pts.extend(b.coord for b in _backbone_path_na(res))
                scs = [b.coord for b in res.beads if _norm_name(b.name).startswith("SC")]
                obj.extend(_cgo_base_plate(scs, (0.0, 0.55, 0.55)))
        obj.extend(_cgo_tube_along(pts, rgb, radius, cap))
    return obj


# ---------------------------------------------------------------------------
# 控制器
# ---------------------------------------------------------------------------

class CGDisplayController:
    """
    粗粒化显示控制器：仅跟踪自有生成对象，不删除用户对象。

    apply 成功后禁用源对象；restore 隐藏自有输出并恢复源先前可见性。
    """

    def __init__(self, cmd: Any = None) -> None:
        self.cmd = cmd or _global_cmd
        self.owned: List[str] = []
        self.source: Optional[str] = None
        self._source_was_enabled: bool = True
        self.last_detection: Optional[DetectionReport] = None
        self.last_apply: Optional[ApplyReport] = None

    def list_sources(self) -> List[str]:
        """可选源对象列表。"""
        return list_source_objects(self.cmd)

    def detect(self, obj: str) -> DetectionReport:
        """检测指定源对象。"""
        rep = detect_cg(obj, self.cmd)
        self.last_detection = rep
        return rep

    def _clean_owned(self, names: Optional[Sequence[str]] = None) -> None:
        """删除自有对象；默认清理控制器记录的全部。"""
        targets = list(names) if names is not None else list(self.owned)
        for n in targets:
            if is_owned_name(n):
                _delete_if_exists(self.cmd, n)
        if names is None:
            self.owned = []
        else:
            drop = set(targets)
            self.owned = [x for x in self.owned if x not in drop]

    def restore(self) -> str:
        """
        隐藏/删除自有显示对象，并恢复源对象先前启用状态。

        不改动源坐标与键。
        """
        self._clean_owned()
        msg = "已恢复。"
        if self.source:
            try:
                if self._source_was_enabled:
                    self.cmd.enable(self.source)
                else:
                    self.cmd.disable(self.source)
            except Exception:
                pass
            msg = "已隐藏粗粒化显示并恢复源对象“%s”可见性。" % self.source
        self.source = None
        self.last_apply = None
        return msg

    def apply(
        self,
        obj: str,
        mode: str = MODE_CONTINUOUS,
        bonds_file: Optional[str] = None,
        auto_templates: bool = True,
    ) -> ApplyReport:
        """
        生成独立显示副本并渲染；成功后禁用源对象。

        未提供 bonds_file 且 auto_templates=True 时，优先匹配已保存连接模板，
        未匹配残基再走几何近似；二者之间不猜测跨界连接。外部文件优先。
        失败时保持源可见，并仅清理本次新建的自有对象。
        """
        c = self.cmd
        report = ApplyReport(mode=validate_mode(mode))
        new_owned: List[str] = []

        if c is None or getattr(c, "_COb", None) is None:
            report.message = "PyMOL 引擎不可用。"
            return report

        # 关机/退出保护：基本探测
        try:
            c.get_names("public_objects", enabled_only=0)
        except Exception:
            report.message = "PyMOL 引擎已停止，取消操作。"
            return report

        det = detect_cg(obj, c)
        self.last_detection = det
        report.detection = det
        if det.empty:
            report.message = det.note or "无 CG 内容。"
            return report

        # 记录源可见性（在任何禁用之前）
        try:
            enabled_names = set(c.get_names("public_objects", enabled_only=1) or [])
            src_enabled = obj in enabled_names
        except Exception:
            src_enabled = True
        original_enabled = self._source_was_enabled if obj == self.source else src_enabled
        old_owned = list(self.owned)
        previous_source = self.source
        previous_source_enabled = self._source_was_enabled

        try:
            model = c.get_model(_quote_sel(obj))
        except Exception as exc:
            report.message = "读取源模型失败：%s" % exc
            return report

        n_atoms = len(model.atom)
        coords = [tuple(float(x) for x in at.coord[:3]) for at in model.atom]

        # 快照源状态用于不变性（调用方/测试可再验）
        try:
            src_bonds_before = len(getattr(model, "bond", []) or [])
        except Exception:
            src_bonds_before = -1

        provenance = PROVENANCE_APPROX
        if bonds_file:
            ext = os.path.splitext(bonds_file)[1].lower()
            provenance = PROVENANCE_ITP if ext in (".itp", ".top") else PROVENANCE_BONDS_FILE
        skipped_long = 0
        isolated_indices: List[int] = []
        try:
            conn, provenance = load_connectivity(
                bonds_file, n_atoms, [at.name for at in model.atom]
            )
            if conn is not None:
                # 外部键/拓扑优先，不做模板匹配
                pairs, skipped_long = filter_bonds_by_distance(conn, coords, PBC_REJECT)
                report.approximate_beads = 0
                report.matched_beads = 0
            elif auto_templates:
                from pymol import cg_templates as cgt

                match = cgt.match_templates(det.residues)
                note_parts: List[str] = []
                if match.warnings:
                    note_parts.extend(match.warnings[:6])
                matched_keys = match.matched_keys
                unmatched = [r for r in det.residues if r.key not in matched_keys]
                approx_pairs: List[Tuple[int, int]] = []
                sk_approx = 0
                if unmatched:
                    approx_pairs, sk_approx, _ = build_approximate_bonds(unmatched)
                # 仅合并，不在匹配/未匹配之间猜测连接
                merged: Set[Tuple[int, int]] = set(match.pairs)
                merged.update(approx_pairs)
                pairs, sk_dist = filter_bonds_by_distance(
                    sorted(merged), coords, PBC_REJECT
                )
                skipped_long = sk_approx + sk_dist
                isolated_indices = list(match.isolated_indices)
                report.matched_beads = int(match.matched_beads)
                report.approximate_beads = sum(len(r.beads) for r in unmatched)
                report.template_counts = dict(match.molecule_counts)
                report.isolated_indices = isolated_indices
                if match.matched_beads > 0 and report.approximate_beads == 0:
                    provenance = PROVENANCE_SAVED
                elif match.matched_beads > 0:
                    provenance = PROVENANCE_MIXED
                else:
                    provenance = PROVENANCE_APPROX
                    if not note_parts:
                        note_parts.append("未匹配到已保存连接模板，使用几何近似。")
                report.template_note = " ".join(note_parts)
            else:
                pairs, skipped_long, _ = build_approximate_bonds(det.residues)
                provenance = PROVENANCE_APPROX
                report.approximate_beads = det.bead_count
                report.matched_beads = 0
                report.template_note = "已关闭自动模板，纯几何近似。"
        except TopologyError as exc:
            report.message = str(exc)
            report.provenance = provenance
            report.provenance_label = PROVENANCE_LABELS.get(
                provenance, PROVENANCE_LABELS[PROVENANCE_APPROX]
            )
            return report
        except Exception as exc:
            report.message = "构建连通失败：%s" % exc
            return report

        report.provenance = provenance
        report.provenance_label = PROVENANCE_LABELS.get(
            provenance, PROVENANCE_LABELS[PROVENANCE_APPROX]
        )
        report.skipped_long = skipped_long
        report.bond_count = len(pairs)

        # 事务式创建：不覆盖同名前缀的用户对象，成功后再清理旧输出。
        mol_name = c.get_unused_name(_owned_name(obj, "mol"))
        cgo_name = c.get_unused_name(_owned_name(obj, "guide"))
        new_owned = [mol_name]

        try:
            bonded = _copy_model_with_bonds(c, obj, pairs)
            _delete_if_exists(c, mol_name)
            c.load_model(bonded, mol_name, zoom=0, finish=0)
            _apply_reps(
                c, mol_name, det.residues, report.mode, isolated_indices=isolated_indices
            )
            _color_components(c, mol_name, det.residues, report.mode)

            if report.mode == MODE_CONTINUOUS:
                guide = _build_guide_cgo(det.residues)
                _delete_if_exists(c, cgo_name)
                if guide:
                    # 必须在 load_cgo 前登记，否则 restore 不会清理 guide
                    new_owned.append(cgo_name)
                    c.load_cgo(guide, cgo_name, zoom=0)
            else:
                _delete_if_exists(c, cgo_name)

            # 成功后再替换旧输出，并恢复上一个源对象的可见性。
            self._clean_owned(old_owned)
            if previous_source and previous_source != obj:
                (c.enable if previous_source_enabled else c.disable)(previous_source)
            self.owned = [n for n in new_owned if n in (c.get_names("all") or [])]
            self.source = obj
            self._source_was_enabled = original_enabled
            try:
                c.disable(obj)
            except Exception:
                pass
            for n in self.owned:
                try:
                    c.enable(n)
                except Exception:
                    pass

            # 验证源未被改写
            try:
                model2 = c.get_model(_quote_sel(obj))
                if len(model2.atom) != n_atoms:
                    raise RuntimeError("源原子数被改变")
                if src_bonds_before >= 0 and len(model2.bond) != src_bonds_before:
                    # 某些 PyMOL 版本 get_model 键数不稳定时仅警告
                    pass
            except Exception:
                pass

            report.ok = True
            report.owned_objects = list(self.owned)
            lim = report.provenance_label
            extra = ""
            if skipped_long:
                extra = " 有 %d 处连接超过距离限制，已断开显示；可检查周期边界或提供拓扑。" % skipped_long
            tmpl = ""
            if report.template_note:
                tmpl = " %s" % report.template_note
            report.message = (
                "已生成粗粒化连通显示（%s）。溯源：%s。保留 CG 分辨率与原始坐标，"
                "未生成全原子坐标。%s%s"
                % (MODE_META[report.mode]["label"], lim, extra, tmpl)
            )
            self.last_apply = report
            return report

        except Exception as exc:
            # 失败：清理本次自有，保持源可见
            self._clean_owned(new_owned)
            try:
                if src_enabled:
                    c.enable(obj)
            except Exception:
                pass
            report.ok = False
            report.message = "生成失败，已保留源对象和先前显示：%s" % exc
            self.last_apply = report
            return report


def summary_text(det: DetectionReport, apply: Optional[ApplyReport] = None) -> str:
    """生成 GUI 用中文摘要。"""
    if det.empty:
        return det.note or "无检测结果。"
    lines = [
        "对象：%s" % det.object,
        "CG 珠数：%d　残基数：%d　排除溶剂/离子珠：%d"
        % (det.bead_count, det.residue_count, det.excluded_beads),
    ]
    if det.by_resn:
        top = list(det.by_resn.items())[:12]
        lines.append(
            "残基名：" + "，".join("%s×%d" % (k, v) for k, v in top)
        )
    if apply is not None and apply.detection is not None and apply.detection.object != det.object:
        apply = None
    if apply is not None:
        if not apply.ok:
            lines.append("状态：" + (apply.message or "失败"))
        lines.append("连接来源：" + apply.provenance_label)
        if apply.template_counts:
            counts = "，".join(
                "%s×%d" % (k, v) for k, v in sorted(apply.template_counts.items())
            )
            lines.append("模板分子：" + counts)
        lines.append("模板覆盖：%s / %s 珠；几何近似：%s 珠" % (
            format(apply.matched_beads, ","), format(det.bead_count, ","),
            format(apply.approximate_beads, ",")))
        if apply.template_note:
            lines.append("提示：" + apply.template_note)
        if apply.ok:
            lines.append(
                "显示连接：%s 对；超长跳过：%d"
                % (format(apply.bond_count, ","), apply.skipped_long)
            )
    lines.append("保留粗粒化分辨率与原始坐标；不生成全原子坐标。")
    return "\n".join(lines)
