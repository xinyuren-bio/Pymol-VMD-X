# ==================================================
# 功能说明：已保存粗粒化拓扑连接模板的加载、校验与保守残基/链匹配。
# 使用方法：load_registry() 懒加载 JSON；match_templates(residues) 返回全局 model_idx 键对。
# 依赖环境：仅标准库；模板位于 $PYMOL_DATA/cg_templates 或源码 data/cg_templates。
# 生成时间：2026-11-14 16:30
# ==================================================

"""已保存 CG 连接模板库（无 Qt、不读原始力场 ITP）。

匹配结果仅表示「已保存连接模板」与源珠名/残基名一致，
不能鉴定力场版本、电荷或质子化状态；虚拟位点记录不作化学键。
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

# 残基键类型与 coarse_grained.CGResidue.key 一致
ResidueKey = Tuple[str, str, str, str, str]


@dataclass
class CGTemplate:
    """单个已保存连接模板。"""

    id: str
    molecule: str
    label: str
    scope: str
    component: str
    package_id: str
    resn_aliases: frozenset
    atoms: List[Dict[str, Any]]
    bonds: List[Tuple[int, int]]
    constraints: List[Tuple[int, int]]
    virtual_sites: List[Any]
    atom_aliases: Dict[str, List[str]]
    # 接受的珠名（含显式别名）→ 模板内 zero-based 下标（residue scope）
    name_to_local: Dict[str, int] = field(default_factory=dict)
    # (resnr, 珠名) → 下标（chain scope；珠名含别名展开）
    res_name_to_local: Dict[Tuple[int, str], int] = field(default_factory=dict)
    connected_local: frozenset = field(default_factory=frozenset)
    isolated_local: frozenset = field(default_factory=frozenset)


@dataclass
class TemplateLibrary:
    """模板注册库。"""

    templates: List[CGTemplate] = field(default_factory=list)
    package_ids: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    ok: bool = False
    paths: List[str] = field(default_factory=list)


@dataclass
class TemplateMatch:
    """一次匹配结果（键对为全局 model_idx）。"""

    pairs: List[Tuple[int, int]] = field(default_factory=list)
    matched_keys: Set[ResidueKey] = field(default_factory=set)
    matched_beads: int = 0
    molecule_counts: Dict[str, int] = field(default_factory=dict)
    template_atom_types: Dict[int, str] = field(default_factory=dict)
    isolated_indices: List[int] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


_REGISTRY: Optional[TemplateLibrary] = None


def _norm(s: str) -> str:
    """规范化名称。"""
    return (s or "").strip().upper()


def _resi_int(r: str) -> Optional[int]:
    """只接受完整整数残基号，插入码不能冒充连续残基。"""
    try:
        return int(str(r).strip())
    except (TypeError, ValueError):
        return None


def _pair(a: int, b: int) -> Tuple[int, int]:
    """有序无向边。"""
    return (a, b) if a < b else (b, a)


def _candidate_dirs() -> List[str]:
    """
    模板根目录候选：先 PYMOL_DATA/cg_templates，再源码 data/cg_templates。

    设计思路：安装后依赖环境变量；开发树用相对 modules/pymol 的回退路径。
    """
    out: List[str] = []
    seen: Set[str] = set()

    def add(p: str) -> None:
        n = os.path.normpath(p)
        if n not in seen:
            seen.add(n)
            out.append(n)

    env = os.environ.get("PYMOL_DATA") or ""
    if env:
        add(os.path.join(env, "cg_templates"))
    here = os.path.dirname(os.path.abspath(__file__))
    # modules/pymol -> <repo>/data/cg_templates
    add(os.path.join(here, "..", "..", "data", "cg_templates"))
    # 若包内嵌 pymol_path/data
    add(os.path.join(here, "pymol_path", "data", "cg_templates"))
    try:
        pymol_path = os.environ.get("PYMOL_PATH") or ""
        if pymol_path:
            add(os.path.join(pymol_path, "data", "cg_templates"))
    except Exception:
        pass
    return out


def _iter_json_paths() -> List[str]:
    """枚举可读的 templates.json。"""
    found: List[str] = []
    for root in _candidate_dirs():
        if not os.path.isdir(root):
            continue
        try:
            names = sorted(os.listdir(root))
        except OSError:
            continue
        for name in names:
            fp = os.path.join(root, name, "templates.json")
            if os.path.isfile(fp):
                found.append(fp)
        # 若根下直接放 json
        direct = os.path.join(root, "templates.json")
        if os.path.isfile(direct):
            found.append(direct)
    # 去重保序
    uniq: List[str] = []
    seen: Set[str] = set()
    for fp in found:
        n = os.path.normpath(fp)
        if n not in seen:
            seen.add(n)
            uniq.append(n)
    return uniq


def _edge_list(raw: Any, n: int, label: str, warnings: List[str]) -> Optional[List[Tuple[int, int]]]:
    """校验并规范化边列表；失败返回 None。"""
    if raw is None:
        return []
    if not isinstance(raw, list):
        warnings.append("模板边字段“%s”不是列表。" % label)
        return None
    out: List[Tuple[int, int]] = []
    for item in raw:
        if not (isinstance(item, (list, tuple)) and len(item) >= 2):
            warnings.append("模板“%s”含非法边记录。" % label)
            return None
        if type(item[0]) is not int or type(item[1]) is not int:
            warnings.append("模板“%s”边索引非整数。" % label)
            return None
        i, j = item[0], item[1]
        if i == j or i < 0 or j < 0 or i >= n or j >= n:
            warnings.append("模板“%s”边索引越界或自环：%s %s。" % (label, i, j))
            return None
        out.append(_pair(i, j))
    return out


def _build_name_maps(
    atoms: Sequence[Dict[str, Any]],
    aliases: Dict[str, List[str]],
    scope: str,
    tid: str,
    warnings: List[str],
) -> Optional[Tuple[Dict[str, int], Dict[Tuple[int, str], int]]]:
    """
    构建珠名映射；别名碰撞或缺珠名则失败。

    仅使用 JSON 中显式 atom_aliases，不做广泛猜测。
    """
    name_to_local: Dict[str, int] = {}
    res_name_to_local: Dict[Tuple[int, str], int] = {}
    canonical = {_norm(a.get("name", "")) for a in atoms}

    for key, vals in aliases.items():
        nk = _norm(key)
        if nk not in canonical:
            warnings.append("模板 %s 的别名键“%s”不是模板珠名。" % (tid, key))
            return None
        if not isinstance(vals, (list, tuple)):
            warnings.append("模板 %s 别名值非法。" % tid)
            return None

    for i, at in enumerate(atoms):
        nm = _norm(str(at.get("name", "")))
        if not nm:
            warnings.append("模板 %s 第 %d 珠缺少名称。" % (tid, i))
            return None
        accepted = [nm]
        for al in aliases.get(nm, []) or []:
            an = _norm(str(al))
            if an and an not in accepted:
                accepted.append(an)

        try:
            rnr = int(at.get("resnr", 1))
        except (TypeError, ValueError):
            warnings.append("模板 %s 残基号非法。" % tid)
            return None

        for an in accepted:
            if scope == "residue" and an in name_to_local and name_to_local[an] != i:
                warnings.append(
                    "模板 %s 珠名/别名“%s”映射冲突。" % (tid, an)
                )
                return None
            if scope == "residue":
                name_to_local[an] = i
            rk = (rnr, an)
            if rk in res_name_to_local and res_name_to_local[rk] != i:
                warnings.append(
                    "模板 %s 链内珠名“%s”在残基 %d 冲突。" % (tid, an, rnr)
                )
                return None
            res_name_to_local[rk] = i

    if scope == "residue":
        # 每个模板珠必须恰有一个规范名可反查
        for i, at in enumerate(atoms):
            nm = _norm(str(at.get("name", "")))
            if name_to_local.get(nm) != i:
                warnings.append("模板 %s 规范珠名映射失败。" % tid)
                return None
    return name_to_local, res_name_to_local


def _parse_template(raw: Dict[str, Any], package_id: str, warnings: List[str]) -> Optional[CGTemplate]:
    """解析并校验单个模板对象。"""
    tid = str(raw.get("id") or raw.get("molecule") or "").strip()
    if not tid:
        warnings.append("存在缺少 id 的模板，已跳过。")
        return None
    scope = str(raw.get("scope") or "residue").strip().lower()
    if scope not in ("residue", "chain"):
        warnings.append("模板 %s 的 scope 非法，已跳过。" % tid)
        return None
    atoms = raw.get("atoms")
    if not isinstance(atoms, list) or not atoms:
        warnings.append("模板 %s 缺少 atoms，已跳过。" % tid)
        return None
    n = len(atoms)
    for i, at in enumerate(atoms):
        if not isinstance(at, dict):
            warnings.append("模板 %s 原子记录损坏。" % tid)
            return None
        try:
            aid = int(at.get("id", i + 1))
        except (TypeError, ValueError):
            warnings.append("模板 %s 原子 id 非法。" % tid)
            return None
        if aid != i + 1:
            warnings.append("模板 %s 原子 id 必须为 1..N 连续，已跳过。" % tid)
            return None
        if not isinstance(at.get("type"), str) or not at["type"].strip():
            warnings.append("模板 %s 缺少珠子类型，已跳过。" % tid)
            return None

    bonds = _edge_list(raw.get("bonds"), n, tid + ".bonds", warnings)
    constraints = _edge_list(raw.get("constraints"), n, tid + ".constraints", warnings)
    if bonds is None or constraints is None:
        return None

    aliases_raw = raw.get("atom_aliases") or {}
    if not isinstance(aliases_raw, dict):
        warnings.append("模板 %s 的 atom_aliases 非法。" % tid)
        return None
    # 规范化别名表：键用规范大写珠名
    aliases: Dict[str, List[str]] = {}
    for k, vals in aliases_raw.items():
        if not isinstance(vals, list) or any(not isinstance(v, str) for v in vals):
            warnings.append("模板 %s 别名值必须为名称列表。" % tid)
            return None
        aliases[_norm(str(k))] = vals

    maps = _build_name_maps(atoms, aliases, scope, tid, warnings)
    if maps is None:
        return None
    name_to_local, res_name_to_local = maps

    connected: Set[int] = set()
    for i, j in bonds + constraints:
        connected.add(i)
        connected.add(j)
    isolated = frozenset(i for i in range(n) if i not in connected)

    resn_values = raw.get("resn_aliases")
    if not isinstance(resn_values, list) or any(not isinstance(x, str) for x in resn_values):
        warnings.append("模板 %s 的 resn_aliases 必须为名称列表。" % tid)
        return None
    resn_aliases = frozenset(_norm(x) for x in resn_values)
    if not resn_aliases:
        warnings.append("模板 %s 缺少 resn_aliases，已跳过。" % tid)
        return None

    return CGTemplate(
        id=tid,
        molecule=str(raw.get("molecule") or tid),
        label=str(raw.get("label") or tid),
        scope=scope,
        component=str(raw.get("component") or "other"),
        package_id=package_id,
        resn_aliases=resn_aliases,
        atoms=list(atoms),
        bonds=bonds,
        constraints=constraints,
        virtual_sites=list(raw.get("virtual_sites") or []),
        atom_aliases=aliases,
        name_to_local=name_to_local,
        res_name_to_local=res_name_to_local,
        connected_local=frozenset(connected),
        isolated_local=isolated,
    )


def _load_package(fp: str) -> Tuple[List[CGTemplate], List[str], str]:
    """从单个 templates.json 加载；损坏时返回空列表与警告。"""
    warnings: List[str] = []
    try:
        with open(fp, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as exc:
        return [], ["无法读取模板库 %s：%s" % (fp, exc)], ""

    if not isinstance(data, dict):
        return [], ["模板库根节点不是对象：%s" % fp], ""

    ver = data.get("schema_version", 1)
    try:
        ver_i = int(ver)
    except (TypeError, ValueError):
        ver_i = -1
    if ver_i != 1:
        return [], ["模板库 schema_version=%s 不受支持：%s" % (ver, fp)], ""

    package_id = str(data.get("id") or os.path.basename(os.path.dirname(fp)) or "unknown")
    raw_list = data.get("templates")
    if not isinstance(raw_list, list):
        return [], ["模板库缺少 templates 数组：%s" % fp], package_id

    # 包内 id 碰撞检查
    by_id: Dict[str, CGTemplate] = {}
    conflicted: Set[str] = set()
    for raw in raw_list:
        if not isinstance(raw, dict):
            warnings.append("包 %s 含非对象模板，已跳过。" % package_id)
            continue
        try:
            tpl = _parse_template(raw, package_id, warnings)
        except (TypeError, ValueError, AttributeError) as exc:
            warnings.append("包 %s 含损坏模板，已跳过：%s" % (package_id, exc))
            continue
        if tpl is None:
            continue
        if tpl.id in conflicted:
            continue
        if tpl.id in by_id:
            if by_id[tpl.id] == tpl:
                continue
            warnings.append("包 %s 中模板 id“%s”内容冲突，已拒绝该 id。" % (package_id, tpl.id))
            by_id.pop(tpl.id)
            conflicted.add(tpl.id)
            continue
        by_id[tpl.id] = tpl
    return list(by_id.values()), warnings, package_id


def load_registry(force_reload: bool = False) -> TemplateLibrary:
    """
    懒加载只读模板库。

    缺失或损坏时不抛异常：返回 ok=False 的空库并附带警告，供调用方回退几何近似。
    """
    global _REGISTRY
    if _REGISTRY is not None and not force_reload:
        return _REGISTRY

    lib = TemplateLibrary()
    paths = _iter_json_paths()
    if not paths:
        lib.warnings.append(
            "未找到已保存拓扑模板库（$PYMOL_DATA/cg_templates 或源码 data/cg_templates）。"
        )
        lib.ok = False
        _REGISTRY = lib
        return lib

    # 跨包模板 id 碰撞：同 id 不同内容则拒绝该 id
    by_id: Dict[str, CGTemplate] = {}
    conflicted: Set[str] = set()
    for fp in paths:
        templates, warns, pkg = _load_package(fp)
        lib.warnings.extend(warns)
        if pkg:
            lib.package_ids.append(pkg)
        lib.paths.append(fp)
        for tpl in templates:
            if tpl.id in conflicted:
                continue
            if tpl.id in by_id:
                if by_id[tpl.id] == tpl:
                    continue
                lib.warnings.append(
                    "模板 id“%s”在多个库中内容冲突，已拒绝该 id。" % tpl.id
                )
                by_id.pop(tpl.id)
                conflicted.add(tpl.id)
                continue
            by_id[tpl.id] = tpl

    lib.templates = list(by_id.values())
    lib.ok = len(lib.templates) > 0
    if not lib.ok and not lib.warnings:
        lib.warnings.append("模板库为空。")
    _REGISTRY = lib
    return lib


def reset_registry_cache() -> None:
    """清除懒加载缓存（测试用）。"""
    global _REGISTRY
    _REGISTRY = None


def _map_residue_beads(
    res: Any, tpl: CGTemplate
) -> Optional[Dict[int, Any]]:
    """
    将单个残基珠按名称（含显式别名）映射到模板下标。

    要求珠数完全一致、一一对应、无重复；不依赖原子顺序。
    """
    beads = list(getattr(res, "beads", []) or [])
    if len(beads) != len(tpl.atoms):
        return None
    used_local: Dict[int, Any] = {}
    for b in beads:
        nm = _norm(getattr(b, "name", ""))
        li = tpl.name_to_local.get(nm)
        if li is None:
            return None
        if li in used_local:
            return None
        used_local[li] = b
    if len(used_local) != len(tpl.atoms):
        return None
    return used_local


def _try_match_residue(
    res: Any, tpl: CGTemplate
) -> Optional[Tuple[List[Tuple[int, int]], Dict[int, str], List[int], ResidueKey]]:
    """尝试用 residue 模板匹配单个残基。"""
    mapping = _map_residue_beads(res, tpl)
    if mapping is None:
        return None
    pairs: List[Tuple[int, int]] = []
    for i, j in tpl.bonds + tpl.constraints:
        bi, bj = mapping[i], mapping[j]
        pairs.append(_pair(int(bi.model_idx), int(bj.model_idx)))
    types: Dict[int, str] = {}
    for li, bead in mapping.items():
        types[int(bead.model_idx)] = str(tpl.atoms[li].get("type") or "")
    isolated = [int(mapping[li].model_idx) for li in sorted(tpl.isolated_local)]
    key: ResidueKey = tuple(res.key)  # type: ignore[assignment]
    return pairs, types, isolated, key


def _chain_groups(
    residues: Sequence[Any],
) -> Dict[Tuple[str, str, str], List[Any]]:
    """按 object/segi/chain 分组。"""
    groups: Dict[Tuple[str, str, str], List[Any]] = defaultdict(list)
    for res in residues:
        groups[(res.object, res.segi, res.chain)].append(res)
    return groups


def _try_match_chain(
    group: Sequence[Any], tpl: CGTemplate
) -> Optional[Tuple[List[Tuple[int, int]], Dict[int, str], List[int], List[ResidueKey]]]:
    """
    尝试用 chain 模板匹配一整条链候选。

    仅完整精确匹配（如 poly-U10）：残基数、连续残基号、每残基珠名均须一致。
    不把超长链切成多个模板窗口。
    """
    # 仅保留别名内残基；同链若另有其他核酸残基则拒绝（其他序列）
    alias_hits = [r for r in group if _norm(r.resn) in tpl.resn_aliases]
    other_na = [
        r for r in group
        if _norm(r.resn) not in tpl.resn_aliases
        and getattr(r, "kind", "") == "nucleic"
    ]
    if other_na:
        return None
    if not alias_hits:
        return None

    # 模板残基号集合
    tpl_resnrs = sorted({int(a.get("resnr", 1)) for a in tpl.atoms})
    n_res = len(tpl_resnrs)
    if len(alias_hits) != n_res:
        return None
    if tpl_resnrs != list(range(tpl_resnrs[0], tpl_resnrs[0] + n_res)):
        return None

    sortable: List[Tuple[int, Any]] = []
    for res in alias_hits:
        ri = _resi_int(res.resi)
        if ri is None:
            return None
        sortable.append((ri, res))
    sortable.sort(key=lambda x: x[0])
    # 残基号不得重复
    nums = [x[0] for x in sortable]
    if len(set(nums)) != len(nums):
        return None
    # 必须整数连续（允许整体偏移）
    for i in range(len(nums) - 1):
        if nums[i + 1] - nums[i] != 1:
            return None

    # 按相对顺序映射到模板 resnr
    base = tpl_resnrs[0]
    mapping: Dict[int, Any] = {}
    keys: List[ResidueKey] = []
    for offset, (_ri, res) in enumerate(sortable):
        rnr = base + offset
        beads = list(res.beads)
        # 该残基在模板中的珠
        local_for_res = [
            i for i, a in enumerate(tpl.atoms) if int(a.get("resnr", 1)) == rnr
        ]
        if len(beads) != len(local_for_res):
            return None
        used: Set[int] = set()
        for b in beads:
            nm = _norm(b.name)
            li = tpl.res_name_to_local.get((rnr, nm))
            if li is None or li in used:
                return None
            used.add(li)
            mapping[li] = b
        if len(used) != len(local_for_res):
            return None
        keys.append(tuple(res.key))  # type: ignore[arg-type]

    if len(mapping) != len(tpl.atoms):
        return None

    pairs: List[Tuple[int, int]] = []
    for i, j in tpl.bonds + tpl.constraints:
        bi, bj = mapping[i], mapping[j]
        pairs.append(_pair(int(bi.model_idx), int(bj.model_idx)))
    types = {
        int(mapping[li].model_idx): str(tpl.atoms[li].get("type") or "")
        for li in mapping
    }
    isolated = [int(mapping[li].model_idx) for li in sorted(tpl.isolated_local)]
    return pairs, types, isolated, keys


def match_templates(
    residues: Sequence[Any],
    registry: Optional[TemplateLibrary] = None,
) -> TemplateMatch:
    """
    将源残基保守匹配到已保存连接模板。

    规则摘要：
    - residue scope：残基名别名 + 珠数 + 珠名（含显式别名）一一对应；
    - chain scope：整链精确长度与连续残基号，不切片超长链；
    - 同一候选命中多个不同模板时拒绝并警告；
    - 不修改源原子；template_atom_types 仅记录模板类型。
    """
    result = TemplateMatch()
    lib = registry if registry is not None else load_registry()
    if lib.warnings:
        result.warnings.extend(lib.warnings)
    if not lib.ok or not lib.templates:
        if not any("未找到" in w or "空" in w or "无法读取" in w for w in result.warnings):
            result.warnings.append("无可用已保存连接模板，回退几何近似。")
        return result

    residue_tpls = [t for t in lib.templates if t.scope == "residue"]
    chain_tpls = [t for t in lib.templates if t.scope == "chain"]

    claimed: Set[ResidueKey] = set()
    pairs: Set[Tuple[int, int]] = set()
    types: Dict[int, str] = {}
    isolated: List[int] = []
    counts: Dict[str, int] = defaultdict(int)

    # ---------- chain scope 优先（完整分子）----------
    groups = _chain_groups(residues)
    for _gk, group in groups.items():
        # 跳过已无剩余候选的组
        pending = [r for r in group if tuple(r.key) not in claimed]
        if not pending:
            continue
        hits: List[Tuple[CGTemplate, Any]] = []
        for tpl in chain_tpls:
            got = _try_match_chain(pending, tpl)
            if got is not None:
                hits.append((tpl, got))
        if len(hits) > 1:
            ids = "、".join(sorted({t.id for t, _ in hits}))
            result.warnings.append(
                "链 %s/%s/%s 同时匹配多个模板（%s），已跳过以免误选。"
                % (pending[0].object, pending[0].segi, pending[0].chain, ids)
            )
            continue
        if len(hits) == 1:
            tpl, (pp, tp, iso, keys) = hits[0]
            if any(k in claimed for k in keys):
                continue
            for k in keys:
                claimed.add(k)
            for e in pp:
                pairs.add(e)
            types.update(tp)
            isolated.extend(iso)
            counts[tpl.molecule] += 1

    # ---------- residue scope ----------
    for res in residues:
        key: ResidueKey = tuple(res.key)  # type: ignore[assignment]
        if key in claimed:
            continue
        hits_r: List[Tuple[CGTemplate, Any]] = []
        for tpl in residue_tpls:
            if _norm(res.resn) not in tpl.resn_aliases:
                continue
            got = _try_match_residue(res, tpl)
            if got is not None:
                hits_r.append((tpl, got))
        if len(hits_r) > 1:
            ids = "、".join(sorted({t.id for t, _ in hits_r}))
            result.warnings.append(
                "残基 %s:%s 同时匹配多个模板（%s），已跳过。"
                % (res.resn, res.resi, ids)
            )
            continue
        if len(hits_r) == 1:
            tpl, (pp, tp, iso, k) = hits_r[0]
            claimed.add(k)
            for e in pp:
                pairs.add(e)
            types.update(tp)
            isolated.extend(iso)
            counts[tpl.molecule] += 1

    result.pairs = sorted(pairs)
    result.matched_keys = claimed
    result.matched_beads = sum(
        len(r.beads) for r in residues if tuple(r.key) in claimed
    )
    result.molecule_counts = dict(sorted(counts.items()))
    result.template_atom_types = types
    # 去重保序
    seen_iso: Set[int] = set()
    iso_out: List[int] = []
    for i in isolated:
        if i not in seen_iso:
            seen_iso.add(i)
            iso_out.append(i)
    result.isolated_indices = iso_out
    return result


def library_summary(lib: Optional[TemplateLibrary] = None) -> str:
    """中文摘要，供 GUI/调试。"""
    reg = lib if lib is not None else load_registry()
    if not reg.ok:
        return "已保存模板：不可用。" + (
            (" " + "；".join(reg.warnings[:3])) if reg.warnings else ""
        )
    parts = [
        "%s(%d珠,%s)"
        % (t.molecule, len(t.atoms), t.scope)
        for t in reg.templates
    ]
    return "已保存连接模板 %d 个：%s" % (len(reg.templates), "，".join(parts))
