# ==================================================
# 功能说明：科研出图后端——多配体启发式检测、三种智能主题样式与会话快照恢复。
# 使用方法：由 pmg_qt.publication_gui 或测试导入；通过 PublicationFigure 调用。
# 依赖环境：已安装的 PyMOL（pymol.cmd）；标准库 copy/dataclasses。
# 生成时间：2026-11-14 13:27
# ==================================================

"""一键科研出图纯后端逻辑（无 Qt 依赖）。"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from pymol import cmd as _global_cmd
from .publication_colors import apply_colors, normalize_color

# 布局样式标识（与视觉主题正交）
STYLE_OVERVIEW = "overview"
STYLE_POCKET = "pocket"

# 视觉主题标识（恰好三种）
THEME_JOURNAL = "journal"
THEME_ILLUSTRATION = "illustration"
THEME_SURFACE = "surface"

# 配体重原子下限（启发式，非绝对标准）
MIN_HEAVY_ATOMS = 5

# 口袋距离（Å）
POCKET_RADIUS = 4.5

# 临时选择名前缀，避免与用户命名冲突
_TMP_PREFIX = "_pubfig_"

# 常见水、离子、结晶添加剂（按残基名排除；启发式，需用户确认）
_EXCLUDED_RESN = frozenset({
    # 水
    "HOH", "WAT", "H2O", "TIP", "TIP3", "TIP4", "T3P", "T4P", "SPC", "SOL",
    # 常见离子
    "NA", "NA+", "NAT", "CL", "CL-", "CLA", "K", "K+", "CA", "CA2", "MG",
    "MG2", "ZN", "ZN2", "FE", "FE2", "FE3", "CU", "MN", "CO", "NI", "CD",
    "BR", "IOD", "I", "LI", "SR", "BA", "AL",
    # 常见缓冲/结晶添加剂
    "SO4", "PO4", "NO3", "NH4", "CO3",
    "GOL", "EDO", "BME", "DMS", "DMSO", "MPD", "ACT", "ACY", "ACE", "FMT",
    "IPA", "MOH", "EOH", "TRS", "EPE", "HEP", "MES", "CIT", "TAR", "PEG",
    "PGE", "PG4", "PE3", "PE4", "PE5", "1PE", "2PE", "P6G", "1PG", "2PG",
    "DOD", "D2O",
})

# 主题元数据（界面下拉/菜单与校验共用；恰好三项）
THEME_META: Dict[str, Dict[str, str]] = {
    THEME_JOURNAL: {
        "id": THEME_JOURNAL,
        "label": "期刊简约",
        "description": (
            "白底、低饱和度链色、细线卡通和棒状配体。"
        ),
    },
    THEME_ILLUSTRATION: {
        "id": THEME_ILLUSTRATION,
        "label": "卡通插画",
        "description": (
            "暖色背景、鲜明链色、加粗卡通和球棍配体。光线追踪时显示深色轮廓。"
        ),
    },
    THEME_SURFACE: {
        "id": THEME_SURFACE,
        "label": "通透表面",
        "description": (
            "蛋白半透明冷色表面叠加卡通上下文；配体不透明棒状/球棍突出。"
        ),
    },
}

# 主题顺序固定，供 UI 枚举
THEME_ORDER: Tuple[str, ...] = (
    THEME_JOURNAL,
    THEME_ILLUSTRATION,
    THEME_SURFACE,
)

# 主题切换时恢复表示的基准设置，再应用当前主题。
_OBJECT_THEME_DEFAULTS: Dict[str, float] = {
    "cartoon_transparency": 0.0,
    "transparency": 0.0,
    "sphere_scale": 1.0,  # 复位为默认球半径；表面主题再写小球棍
    "stick_radius": 0.25,
    "cartoon_oval_length": 1.15,
    "cartoon_oval_width": 0.4,
    "cartoon_loop_radius": 0.2,
    "cartoon_tube_radius": 0.5,
    "cartoon_rect_length": 1.4,
    "cartoon_rect_width": 0.4,
    "surface_quality": 0,
    "cartoon_color": -1,
    "cartoon_highlight_color": -1,
    "surface_color": -1,
    "stick_color": -1,
    "sphere_color": -1,
    "stick_ball": 0,
    "stick_transparency": 0.0,
    "sphere_transparency": 0.0,
}


def list_themes() -> List[Dict[str, str]]:
    """返回恰好三种主题的元数据列表（含 id/label/description）。"""
    return [dict(THEME_META[k]) for k in THEME_ORDER]


def validate_theme(theme: str) -> str:
    """
    校验主题标识；非法则抛出 ValueError，调用方应在突变场景前调用。

    返回规范化后的主题 id。
    """
    key = (theme or "").strip().lower()
    if key not in THEME_META:
        raise ValueError(
            "未知智能风格“%s”，可选：%s。"
            % (theme, "、".join(THEME_ORDER))
        )
    return key


@dataclass
class LigandCandidate:
    """单个配体候选（按 object/segi/chain/resi/resn 分组）。"""

    object: str
    segi: str
    chain: str
    resi: str
    resn: str
    indices: List[int]
    heavy_atom_count: int
    selection: str
    label: str


@dataclass(frozen=True)
class ProteinChain:
    object: str
    segi: str
    chain: str

    @property
    def key(self) -> Tuple[str, str, str]:
        return self.object, self.segi, self.chain

    @property
    def label(self) -> str:
        chain = f"链 {self.chain}" if self.chain else "无链名"
        segment = f"／{self.segi}" if self.segi else ""
        return f"{self.object}{segment} · {chain}"


@dataclass
class SceneInfo:
    """当前启用分子对象的检测结果。"""

    objects: List[str] = field(default_factory=list)
    has_protein: bool = False
    protein_atom_count: int = 0
    ligands: List[LigandCandidate] = field(default_factory=list)
    empty: bool = True
    note: str = ""
    protein_chains: List[ProteinChain] = field(default_factory=list)


def _quote_sel(s: str) -> str:
    """对选择表达式中的对象名做安全包裹。"""
    return 'model "%s"' % s.replace('\\', '\\\\').replace('"', '\\"')


def selection_from_indices(obj: str, indices: Sequence[int]) -> str:
    """
    用原子 index 构造稳定选择串，避免仅靠 resi/chain 歧义。

    设计思路：index 在对象内唯一；导出前不增删原子时可靠。
    """
    if not indices:
        return "none"
    uniq = sorted({int(i) for i in indices})
    # PyMOL 支持 index 1+2+3 语法
    joined = "+".join(str(i) for i in uniq)
    return "(%s and index %s)" % (_quote_sel(obj), joined)


def enabled_molecule_objects(cmd: Any = None) -> List[str]:
    """返回当前启用的分子对象名列表。"""
    c = cmd or _global_cmd
    names = c.get_names("public_objects", enabled_only=1) or []
    out: List[str] = []
    for n in names:
        try:
            if c.get_type(n) == "object:molecule":
                out.append(n)
        except Exception:
            continue
    return out


def detect_scene(cmd: Any = None) -> SceneInfo:
    """
    检测启用分子对象中的蛋白与配体候选。

    配体规则（启发式，需用户选择确认）：
    - organic 且非 polymer
    - 排除溶剂与常见添加剂残基名
    - 重原子数 >= MIN_HEAVY_ATOMS
    - 按 object/segi/chain/resi/resn 分组
    """
    c = cmd or _global_cmd
    objs = enabled_molecule_objects(c)
    info = SceneInfo(objects=list(objs))
    if not objs:
        info.empty = True
        info.note = "场景中没有启用的分子对象。"
        return info

    scope = " or ".join("(%s)" % _quote_sel(o) for o in objs)
    keys = set()
    info.protein_atom_count = int(c.iterate(
        f"({scope}) and polymer.protein", "keys.add((model, segi, chain))", space={"keys": keys}) or 0)
    info.protein_chains = [ProteinChain(*key) for key in sorted(keys)]
    info.has_protein = info.protein_atom_count > 0

    # 收集 organic 非聚合物原子，按残基键分组
    buckets: Dict[Tuple[str, str, str, str, str], List[Tuple[int, str]]] = {}
    space: Dict[str, Any] = {"buckets": buckets, "excluded": _EXCLUDED_RESN}

    def _collect(
        model: str,
        segi: str,
        chain: str,
        resi: str,
        resn: str,
        index: int,
        elem: str,
    ) -> None:
        """iterate 回调：按残基键累积原子。"""
        rn = (resn or "").strip().upper()
        if rn in space["excluded"]:
            return
        key = (model, segi or "", chain or "", str(resi), rn or resn or "")
        space["buckets"].setdefault(key, []).append((int(index), (elem or "").upper()))

    space["collect"] = _collect
    sele = "(%s) and organic and not polymer and not solvent" % scope
    c.iterate(
        sele,
        "collect(model, segi, chain, resi, resn, index, elem)",
        space=space,
    )

    ligands: List[LigandCandidate] = []
    for (model, segi, chain, resi, resn), atoms in sorted(
        buckets.items(), key=lambda x: (x[0][0], x[0][2], x[0][3], x[0][4])
    ):
        idxs = [a[0] for a in atoms]
        heavy = sum(1 for _, e in atoms if e and e not in {"H", "D"})
        if heavy < MIN_HEAVY_ATOMS:
            continue
        sel = selection_from_indices(model, idxs)
        chain_disp = chain if chain else "_"
        segi_disp = segi if segi else ""
        label = "%s/%s/%s/%s %s（重原子 %d）" % (
            model,
            segi_disp,
            chain_disp,
            resi,
            resn,
            heavy,
        )
        ligands.append(
            LigandCandidate(
                object=model,
                segi=segi,
                chain=chain,
                resi=resi,
                resn=resn,
                indices=idxs,
                heavy_atom_count=heavy,
                selection=sel,
                label=label,
            )
        )

    info.ligands = ligands
    info.empty = (not info.has_protein) and (not ligands)
    if info.empty:
        info.note = "未检测到蛋白或可用配体候选，请先加载并启用分子对象。"
    elif not info.has_protein and ligands:
        info.note = "未检测到蛋白；可将独立配体用于小分子出图（启发式候选，请核对）。"
    elif info.has_protein and not ligands:
        info.note = "检测到蛋白，但无配体候选；可使用无配体结构总览。"
    else:
        info.note = (
            "已识别蛋白及候选配体，请选择要展示的配体。"
        )
    return info


def candidate_key(lig: Optional[LigandCandidate]) -> Optional[Tuple[Any, ...]]:
    """配体候选稳定键，供确认状态比对。"""
    if lig is None:
        return None
    return (lig.object, lig.segi, lig.chain, lig.resi, lig.resn, tuple(lig.indices))


def normalize_ligands(
    ligand: Optional[LigandCandidate] = None,
    ligands: Optional[Sequence[LigandCandidate]] = None,
) -> List[LigandCandidate]:
    """
    合并单配体与多配体参数并去重。

    设计思路：ligands 非 None 时以其为准（含空列表=仅蛋白总览）；
    否则回退到 legacy 的单个 ligand。不做假合并候选对象。
    """
    if ligands is not None:
        items = list(ligands)
    elif ligand is not None:
        items = [ligand]
    else:
        items = []
    out: List[LigandCandidate] = []
    seen: Set[Tuple[Any, ...]] = set()
    for lig in items:
        if lig is None:
            continue
        key = candidate_key(lig)
        if key is None or key in seen:
            continue
        seen.add(key)
        out.append(lig)
    return out


def ligand_selection(ligands: Sequence[LigandCandidate]) -> Optional[str]:
    """
    将多个候选的 index 选择安全并集为一个选择串。

    跨对象时用 or 连接各候选既有 selection；空列表返回 None。
    """
    parts: List[str] = []
    for lig in ligands:
        sel = (lig.selection or "").strip()
        if not sel or sel == "none":
            continue
        parts.append("(%s)" % sel)
    if not parts:
        return None
    return " or ".join(parts)


def _candidate_in_scene(lig: LigandCandidate, scene_ligands: Sequence[LigandCandidate]) -> bool:
    """判断候选是否仍出现在当前 detect_scene 结果中。"""
    key = candidate_key(lig)
    return any(candidate_key(c) == key for c in scene_ligands)


class PublicationFigure:
    """
    科研出图控制器：会话快照、样式应用、PNG 导出。

    设计思路：首次改样式前内存快照 get_session；恢复时 deepcopy 再 set_session，
    保证多次恢复行为可预期。临时选择统一前缀并在 finally 中清理。
    布局（overview/pocket）与视觉主题（journal/illustration/surface）正交。
    """

    def __init__(self, cmd: Any = None, *, track_undo: bool = False) -> None:
        self.cmd = cmd or _global_cmd
        self._snapshot: Optional[dict] = None
        self._tmp_names: List[str] = []
        self._track_undo = track_undo
        self._undo = None

    def has_snapshot(self) -> bool:
        """是否已保存原始会话快照。"""
        return self._snapshot is not None

    def clear_snapshot(self) -> None:
        """丢弃快照（测试或关闭对话框时用）。"""
        self.clear_undo()
        self._snapshot = None

    def clear_undo(self):
        if self._undo is not None:
            undo, self._undo = self._undo, None
            undo.discard()

    def can_undo(self) -> bool:
        return self._undo is not None

    def undo_style(self) -> bool:
        if self._undo is None:
            return False
        undo, self._undo = self._undo, None
        return undo.restore()

    def ensure_snapshot(self) -> None:
        """在首次修改场景前保存内存会话快照。"""
        if self._snapshot is None:
            # deepcopy 避免后续 set_session 影响缓存副本
            self._snapshot = copy.deepcopy(self.cmd.get_session(quiet=1))

    def restore_original(self) -> bool:
        """
        恢复首次样式前的会话。

        返回是否成功执行恢复。
        """
        if self._snapshot is None:
            return False
        self.clear_undo()
        self.cmd.set_session(copy.deepcopy(self._snapshot), quiet=1)
        self.clear_snapshot()
        return True

    def _tmp(self, tag: str) -> str:
        """分配唯一临时选择名并登记以便清理。"""
        name = self.cmd.get_unused_name(_TMP_PREFIX + tag)
        self._tmp_names.append(name)
        return name

    def _cleanup_temps(self) -> None:
        """删除本轮创建的临时选择，不触碰用户对象。"""
        for n in self._tmp_names:
            try:
                if n in (self.cmd.get_names("all") or []):
                    self.cmd.delete(n)
            except Exception:
                pass
        self._tmp_names = []

    def _scope_sele(self, objs: Optional[Sequence[str]] = None) -> str:
        """构造启用分子对象的并集选择。"""
        names = list(objs) if objs is not None else enabled_molecule_objects(self.cmd)
        if not names:
            return "none"
        return " or ".join("(%s)" % _quote_sel(o) for o in names)

    def _reset_object_theme_settings(self, objs: Sequence[str]) -> None:
        """
        复位对象级透明度/半径等主题残留，防止表面或加粗卡通泄漏到期刊简约。

        设计思路：先统一清零/默认，再由当前主题写入；仅作用于启用分子对象名。
        """
        c = self.cmd
        for obj in objs:
            # 同时移除对象、原子与键级设置，避免配体小球半径等覆盖残留。
            c.unset_deep(" ".join(_OBJECT_THEME_DEFAULTS), obj, updates=0)
            for key, value in _OBJECT_THEME_DEFAULTS.items():
                c.set(key, value, obj)

    def _apply_theme_appearance(self, theme: str) -> None:
        """按主题写入全局背景/光照/轮廓（keep_appearance=False 时调用）。"""
        c = self.cmd
        c.bg_color("0xF3EDE3" if theme == THEME_ILLUSTRATION else "white")
        c.set("bg_gradient", 0)
        c.set("ray_opaque_background", 1)
        c.set("opaque_background", 1)
        c.set("ambient_occlusion_mode", 0)
        c.set("depth_cue", 0)
        c.set("ray_trace_fog", 0)
        c.set("ray_shadows", 0)
        if theme == THEME_ILLUSTRATION:
            c.set("ray_trace_mode", 1)
            c.set("ray_trace_color", "0x252336")
            c.set("ray_trace_gain", 0.13)
            c.set("direct", 0.55)
            c.set("reflect", 0.42)
            c.set("light_count", 3)
            c.set("shininess", 26)
            c.set("specular", 0.28)
            c.set("ambient", 0.32)
        elif theme == THEME_SURFACE:
            # 通透表面：柔和环境光，轮廓关闭以免表面发脏
            c.set("ray_trace_mode", 0)
            c.set("direct", 0.35)
            c.set("reflect", 0.5)
            c.set("light_count", 3)
            c.set("shininess", 20)
            c.set("specular", 0.1)
            c.set("ambient", 0.58)
        else:
            # 期刊简约：克制白底出版光
            c.set("ray_trace_mode", 0)
            c.set("direct", 0.36)
            c.set("reflect", 0.38)
            c.set("light_count", 2)
            c.set("shininess", 40)
            c.set("specular", 0.08)
            c.set("ambient", 0.56)

    def _apply_theme_reps(self, theme, objs, protein_sele, has_protein, lig_sele, pocket_sele):
        """Apply the preset geometry before assigning the independent color schemes."""
        c = self.cmd
        illustration = theme == THEME_ILLUSTRATION
        surface = theme == THEME_SURFACE
        shape = {
            THEME_JOURNAL: (0.16, 0.95, 0.30, 0.16, 0.42, 1.25, 0.31),
            THEME_ILLUSTRATION: (0.22, 1.55, 0.58, 0.30, 0.65, 1.58, 0.50),
            THEME_SURFACE: (0.20, 1.15, 0.40, 0.20, 0.50, 1.40, 0.40),
        }[theme]
        settings = ("stick_radius", "cartoon_oval_length", "cartoon_oval_width",
                    "cartoon_loop_radius", "cartoon_tube_radius", "cartoon_rect_length",
                    "cartoon_rect_width")
        pocket_transparency = {THEME_JOURNAL: 0.78, THEME_ILLUSTRATION: 0.66, THEME_SURFACE: 0.55}
        for obj in objs:
            for key, value in zip(settings, shape):
                c.set(key, value, obj)
            if surface:
                c.set("surface_quality", 1, obj)
                c.set("transparency", 0.55 if pocket_sele else 0.45, obj)
            if pocket_sele:
                c.set("cartoon_transparency", pocket_transparency[theme], obj)
            elif surface:
                c.set("cartoon_transparency", 0.35, obj)

        c.set("cartoon_fancy_helices", 1)
        c.set("cartoon_smooth_loops", 1)
        c.set("cartoon_flat_sheets", 1)
        c.set("cartoon_side_chain_helper", 1)
        if has_protein:
            c.show("cartoon", protein_sele)
            if surface:
                c.show("surface", protein_sele)
        if pocket_sele:
            c.show("sticks", pocket_sele)
            c.set_bond("stick_radius", 0.14, pocket_sele, pocket_sele)
        if lig_sele:
            c.show("sticks", lig_sele)
            c.set_bond("stick_radius", 0.24 if illustration else 0.20, lig_sele, lig_sele)
            if illustration or surface:
                heavy = f"({lig_sele}) and not elem H+D"
                c.show("spheres", heavy)
                c.set("sphere_scale", 0.26 if illustration else 0.22, heavy)

    def apply_style(
        self,
        style: str = STYLE_OVERVIEW,
        ligand: Optional[LigandCandidate] = None,
        keep_view: bool = False,
        objs: Optional[Sequence[str]] = None,
        keep_appearance: bool = False,
        theme: str = THEME_JOURNAL,
        ligands: Optional[Sequence[LigandCandidate]] = None,
        protein_color=None,
        ligand_color=None,
        protein_mode=None,
        ligand_mode=None,
        chain_colors=None,
        view=None,
        appearance=None,
    ) -> SceneInfo:
        """
        将出版启发样式应用到当前可编辑场景。

        style: overview | pocket（布局维度）
        theme: journal | illustration | surface（视觉主题，默认 journal 兼容旧调用）
        ligand: 单个配体（legacy）；与 ligands 二选一时 ligands 优先。
        ligands: 多配体列表；None/空表示无配体蛋白总览；独立多配体亦可。
        keep_view: True 时不改相机。
        keep_appearance: True 时保留当前背景、光照与轮廓/阴影，仍应用表示样式。
        protein_mode: chain（默认）| single；chain_colors 以 ProteinChain.key 为键。
        ligand_mode: element（默认）| carbon；ligand_color 只覆盖配体碳原子。
        未指定模式时，显式传入蛋白或配体颜色分别启用 single 或 carbon。
        """
        # 校验主题须在任何场景突变之前
        theme = validate_theme(theme)
        protein_color = normalize_color(protein_color, self.cmd)
        ligand_color = normalize_color(ligand_color, self.cmd)
        protein_mode = protein_mode or ("single" if protein_color is not None else "chain")
        ligand_mode = ligand_mode or ("carbon" if ligand_color is not None else "element")
        if protein_mode not in ("chain", "single") or ligand_mode not in ("element", "carbon"):
            raise ValueError("无法识别所选着色方式。")
        chain_colors = {key: normalize_color(color, self.cmd)
                        for key, color in (chain_colors or {}).items() if color is not None}
        # 归一化/去重须在快照与突变之前完成
        lig_list = normalize_ligands(ligand=ligand, ligands=ligands)

        info = detect_scene(self.cmd)
        if info.empty and not lig_list:
            raise RuntimeError(info.note or "场景为空，无法出图。")

        if style == STYLE_POCKET:
            if not info.has_protein or not lig_list:
                raise RuntimeError("结合口袋样式需要同时存在蛋白与所选配体。")

        # 任一过期/无效候选整单拒绝，不改场景
        for lig in lig_list:
            if not _candidate_in_scene(lig, info.ligands):
                raise RuntimeError("配体或场景已变化，请刷新候选后重试。")

        lig_sele = ligand_selection(lig_list)
        if style == STYLE_POCKET and lig_sele and not self.cmd.count_atoms(
            "((%s) and polymer.protein) within %s of (%s)"
            % (self._scope_sele(info.objects), POCKET_RADIUS, lig_sele)
        ):
            raise RuntimeError("配体附近 4.5 Å 内没有蛋白原子，无法生成结合口袋图。")
        self.ensure_snapshot()
        c = self.cmd
        obj_list = list(objs) if objs is not None else list(info.objects)
        scope = self._scope_sele(obj_list)

        undo = None
        if self._track_undo:
            from .publication_recipe import APPEARANCE_SETTINGS
            from .publication_undo import FigureUndo
            changed = (
                "orthoscopic", "antialias", "opaque_background", "ray_opaque_background",
                "cartoon_fancy_helices", "cartoon_smooth_loops", "cartoon_flat_sheets",
                "cartoon_side_chain_helper",
            )
            affected = list(dict.fromkeys([*obj_list, *(lig.object for lig in lig_list)]))
            undo = FigureUndo(c, affected, self._scope_sele(affected), changed + APPEARANCE_SETTINGS,
                              tuple(_OBJECT_THEME_DEFAULTS))
        try:
            if not keep_appearance:
                self._apply_theme_appearance(theme)
            if appearance is not None:
                for key, value in appearance.items():
                    c.set(key, value)
            c.set("orthoscopic", 1)
            c.set("antialias", 2)

            # 切换主题前复位对象级残留，避免表面/加粗卡通泄漏
            self._reset_object_theme_settings(obj_list)

            c.hide("everything", scope)
            c.hide("everything", "(%s) and solvent" % scope)

            protein_sele = "(%s) and polymer.protein" % scope
            has_protein = c.count_atoms(protein_sele) > 0

            if lig_sele and c.count_atoms(lig_sele) < 1:
                raise RuntimeError("所选配体选择为空，请刷新候选后重试。")

            pocket_sele = None
            if style == STYLE_POCKET and lig_sele and has_protein:
                # 口袋围绕全部所选配体并集
                pocket = self._tmp("pocket")
                c.select(
                    pocket,
                    "byres ((%s) within %s of (%s))"
                    % (protein_sele, POCKET_RADIUS, lig_sele),
                    0,
                )
                pocket_sele = "?" + pocket

            self._apply_theme_reps(
                theme=theme,
                objs=obj_list,
                protein_sele=protein_sele,
                has_protein=has_protein,
                lig_sele=lig_sele,
                pocket_sele=pocket_sele,
            )

            chains = [chain for chain in info.protein_chains if chain.object in obj_list]
            apply_colors(c, theme, protein_sele, chains, lig_sele, protein_mode=protein_mode,
                         ligand_mode=ligand_mode, protein_color=protein_color,
                         ligand_color=ligand_color, chain_colors=chain_colors)

            if view is not None:
                c.set_view(view, animate=0)
            elif not keep_view:
                if style == STYLE_POCKET and lig_sele and pocket_sele:
                    c.orient("(%s) or (%s)" % (lig_sele, pocket_sele), animate=-1)
                    c.zoom("(%s) or (%s)" % (lig_sele, pocket_sele), buffer=3.0, animate=-1)
                elif lig_sele and has_protein:
                    c.zoom("(%s) or (%s)" % (protein_sele, lig_sele), buffer=3.0 if theme == THEME_SURFACE else 2.0, animate=-1)
                elif lig_sele:
                    c.zoom(lig_sele, buffer=2.0, animate=-1)
                elif has_protein:
                    c.zoom(protein_sele, buffer=3.0 if theme == THEME_SURFACE else 2.0, animate=-1)
                else:
                    c.zoom(scope, animate=-1)
        except Exception:
            if undo is not None:
                undo.restore()
            raise
        finally:
            self._cleanup_temps()

        if undo is not None:
            self.clear_undo()
            self._undo = undo
        return info

    def export_png(
        self,
        fp: str,
        style: str = STYLE_OVERVIEW,
        ligand: Optional[LigandCandidate] = None,
        width: int = 1800,
        height: int = 1400,
        dpi: float = 300.0,
        transparent: bool = False,
        keep_view: bool = False,
        objs: Optional[Sequence[str]] = None,
        keep_appearance: bool = False,
        theme: str = THEME_JOURNAL,
        ligands: Optional[Sequence[LigandCandidate]] = None,
        protein_color=None,
        ligand_color=None,
        protein_mode=None,
        ligand_mode=None,
        chain_colors=None,
    ) -> str:
        """
        应用样式并以光线追踪导出 PNG。

        调用方应在此之前已选定保存路径（GUI 在突变前询问文件名）。
        keep_appearance: 保留背景/光照/轮廓阴影；透明导出仍会关闭不透明背景。
        theme: 须传入所选主题，不得静默回落到 journal。
        ligands: 多配体列表，完整传给 apply_style；与 ligand 兼容。
        """
        theme = validate_theme(theme)
        self.apply_style(
            style=style,
            ligand=ligand,
            keep_view=keep_view,
            objs=objs,
            keep_appearance=keep_appearance,
            theme=theme,
            ligands=ligands,
            protein_color=protein_color,
            ligand_color=ligand_color,
            protein_mode=protein_mode,
            ligand_mode=ligand_mode,
            chain_colors=chain_colors,
        )
        c = self.cmd
        # 透明与否由本次导出选项决定，不沿用上一次导出的透明状态。
        opaque = 0 if transparent else 1
        c.set("opaque_background", opaque)
        c.set("ray_opaque_background", opaque)
        c.png(str(fp), int(width), int(height), dpi=float(dpi), ray=1, quiet=0)
        return str(fp)


def find_candidate_by_label(
    ligands: Sequence[LigandCandidate], label: str
) -> Optional[LigandCandidate]:
    """按界面显示标签查找候选。"""
    for lig in ligands:
        if lig.label == label:
            return lig
    return None
