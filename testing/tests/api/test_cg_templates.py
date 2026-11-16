# ==================================================
# 功能说明：已保存粗粒化连接模板加载、保守匹配与 apply 集成测试。
# 使用方法：由仓库测试运行器或 python -m unittest 发现；禁止依赖 pytest。
# 依赖环境：已安装的 pymol；可选 PYMOLX_CG_SAMPLE 样本 PDB。
# 生成时间：2026-11-14 16:30
# ==================================================

"""cg_templates 与 coarse_grained 自动模板匹配单元测试。"""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from typing import List, Optional

from chempy import Atom
from chempy.models import Indexed

from pymol import cmd
from pymol import cg_templates as cgt
from pymol import coarse_grained as cg

_SAMPLE_PDB = os.environ.get("PYMOLX_CG_SAMPLE", "")


def _bead(
    name: str,
    resn: str,
    chain: str,
    resi: str,
    x: float,
    y: float,
    z: float,
    segi: str = "",
) -> Atom:
    """构造 CG 珠原子。"""
    at = Atom()
    at.name = name
    at.resn = resn
    at.chain = chain
    at.resi = str(resi)
    at.segi = segi
    at.symbol = "C"
    at.coord = [float(x), float(y), float(z)]
    return at


def _load_beads(name: str, atoms: list) -> None:
    """加载 chempy 模型对象。"""
    model = Indexed()
    for at in atoms:
        model.add_atom(at)
    if name in (cmd.get_names("all") or []):
        cmd.delete(name)
    cmd.load_model(model, name, zoom=0)


def _tpl_by_id(lib: cgt.TemplateLibrary, tid: str) -> cgt.CGTemplate:
    """按 id 取模板。"""
    for t in lib.templates:
        if t.id == tid:
            return t
    raise KeyError(tid)


def _atoms_from_template(
    tpl: cgt.CGTemplate,
    chain: str,
    resi: str,
    resn: Optional[str] = None,
    name_map: Optional[dict] = None,
    shuffle: bool = False,
    x0: float = 0.0,
    y0: float = 0.0,
) -> List[Atom]:
    """
    由模板 atoms 生成残基珠；可选打乱顺序与珠名别名替换。

    设计思路：坐标沿 x 均匀排布，保证模板键在距离过滤内。
    """
    atoms: List[Atom] = []
    rn = resn or str(tpl.atoms[0].get("resn") or tpl.molecule)
    order = list(range(len(tpl.atoms)))
    if shuffle:
        order = list(reversed(order))
    for k, li in enumerate(order):
        at = tpl.atoms[li]
        nm = str(at["name"])
        if name_map and nm in name_map:
            nm = name_map[nm]
        atoms.append(
            _bead(nm, rn, chain, resi, x0 + float(k) * 1.5, y0, 0.0)
        )
    return atoms


def _rna_chain_atoms(
    tpl: cgt.CGTemplate,
    chain: str,
    start_resi: int = 1,
    n_res: Optional[int] = None,
    skip_resi: Optional[int] = None,
    x0: float = 0.0,
) -> List[Atom]:
    """由 RNA 模板生成一条链；可截断/跳号以构造拒绝用例。"""
    by_rnr = {}
    for at in tpl.atoms:
        by_rnr.setdefault(int(at["resnr"]), []).append(at)
    resnrs = sorted(by_rnr)
    if n_res is not None:
        resnrs = resnrs[:n_res]
    out: List[Atom] = []
    emit_i = 0
    for rnr in resnrs:
        resi = start_resi + emit_i
        if skip_resi is not None and resi == skip_resi:
            emit_i += 1
            resi = start_resi + emit_i
        for j, at in enumerate(by_rnr[rnr]):
            out.append(
                _bead(
                    str(at["name"]),
                    str(at["resn"]),
                    chain,
                    str(resi),
                    x0 + float(emit_i) * 8.0 + float(j) * 1.2,
                    0.0,
                    0.0,
                )
            )
        emit_i += 1
    return out


class TestRegistrySchema(unittest.TestCase):
    """注册库加载与边/类型完整性。"""

    def setUp(self) -> None:
        cgt.reset_registry_cache()

    def test_load_registry_edges_and_types(self) -> None:
        """五个选用模板应可加载，边在范围内且含类型。"""
        lib = cgt.load_registry()
        self.assertTrue(lib.ok, lib.warnings)
        by_mol = {t.molecule: t for t in lib.templates}
        self.assertEqual(
            set(by_mol), {"AL15", "AL31", "DPPC", "CHOL", "RNA"}
        )
        self.assertEqual(len(by_mol["AL15"].atoms), 56)
        self.assertEqual(len(by_mol["AL15"].bonds), 55)
        self.assertEqual(len(by_mol["AL31"].atoms), 14)
        self.assertEqual(len(by_mol["AL31"].bonds), 13)
        self.assertEqual(len(by_mol["DPPC"].atoms), 12)
        self.assertEqual(len(by_mol["DPPC"].bonds), 11)
        self.assertEqual(len(by_mol["CHOL"].atoms), 8)
        self.assertEqual(len(by_mol["CHOL"].bonds) + len(by_mol["CHOL"].constraints), 6)
        self.assertEqual(len(by_mol["RNA"].atoms), 80)
        self.assertEqual(
            len(by_mol["RNA"].bonds) + len(by_mol["RNA"].constraints), 89
        )
        # CHOL 孤立虚拟点（无键/约束）
        self.assertEqual(len(by_mol["CHOL"].isolated_local), 3)
        for t in lib.templates:
            n = len(t.atoms)
            for i, j in t.bonds + t.constraints:
                self.assertGreaterEqual(i, 0)
                self.assertLess(j, n)
                self.assertNotEqual(i, j)
            for at in t.atoms:
                self.assertTrue(at.get("type"))
                self.assertTrue(at.get("name"))

    def test_corrupt_library_and_conflicting_ids_are_rejected(self):
        """坏 JSON 可恢复，同 id 的不同连接不能任取一个。"""
        template = _tpl_by_id(cgt.load_registry(), "DPPC")
        raw = {
            "id": template.id, "scope": "residue", "atoms": template.atoms,
            "resn_aliases": list(template.resn_aliases),
            "bonds": template.bonds, "constraints": template.constraints,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "templates.json"
            path.write_text('{invalid')
            templates, warnings, _ = cgt._load_package(str(path))
            self.assertFalse(templates)
            self.assertTrue(warnings)
            other = copy.deepcopy(raw)
            other["bonds"] = other["bonds"][:-1]
            path.write_text(json.dumps({"schema_version": 1, "templates": [raw, other]}))
            templates, warnings, _ = cgt._load_package(str(path))
            self.assertFalse(templates)
            self.assertTrue(any("冲突" in w for w in warnings))
        for field, bad in (("bonds", [[0.5, 1]]), ("bonds", [[0, 999]]),
                           ("atom_aliases", {"NC3": "PO4"})):
            broken = dict(raw, **{field: bad})
            warnings = []
            self.assertIsNone(cgt._parse_template(broken, "test", warnings))
            self.assertTrue(warnings)


class TestResidueMatching(unittest.TestCase):
    """残基 scope 匹配与拒绝。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cgt.reset_registry_cache()
        self.lib = cgt.load_registry()

    def test_shuffled_atoms_still_match_edges(self) -> None:
        """珠顺序打乱仍应按名映射得到正确边。"""
        tpl = _tpl_by_id(self.lib, "DPPC")
        atoms = _atoms_from_template(tpl, "L", "1", shuffle=True)
        _load_beads("dppc1", atoms)
        det = cg.detect_cg("dppc1", cmd)
        m = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m.matched_beads, 12)
        self.assertEqual(m.molecule_counts.get("DPPC"), 1)
        self.assertEqual(len(m.pairs), 11)
        # 边端点须为全局 model_idx
        name_to_idx = {}
        for res in det.residues:
            for b in res.beads:
                name_to_idx[cg._norm_name(b.name)] = b.model_idx
        # 模板边 0-1 为 NC3-PO4
        e = tuple(sorted((name_to_idx["NC3"], name_to_idx["PO4"])))
        self.assertIn(e, set(m.pairs))
        # 模板类型写入记录，不改源
        self.assertTrue(m.template_atom_types)
        model = cmd.get_model("dppc1")
        self.assertTrue(all(not hasattr(at, "cg_type") for at in model.atom))

    def test_ec_rotation_alias(self) -> None:
        """显式 EC10→0EC1 别名可匹配；重复珠名拒绝。"""
        tpl = _tpl_by_id(self.lib, "AL15")
        atoms = _atoms_from_template(
            tpl, "A", "1", name_map={"EC10": "0EC1"}, shuffle=True
        )
        _load_beads("al15a", atoms)
        det = cg.detect_cg("al15a", cmd)
        m = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m.matched_beads, 56)
        self.assertEqual(m.molecule_counts.get("AL15"), 1)
        self.assertEqual(len(m.pairs), 55)

        # 同时存在 EC10 与其别名 0EC1 → 同一模板珠重复映射，拒绝
        atoms2 = _atoms_from_template(tpl, "A", "2")
        for a in atoms2:
            if a.name == "C1A":
                a.name = "0EC1"
                break
        self.assertIn("EC10", [a.name for a in atoms2])
        _load_beads("al15bad", atoms2)
        det2 = cg.detect_cg("al15bad", cmd)
        m2 = cgt.match_templates(det2.residues, self.lib)
        self.assertEqual(m2.matched_beads, 0)

    def test_incomplete_duplicate_wrong_resn_rejected(self) -> None:
        """不完整、重复珠名、错误残基名均不匹配。"""
        tpl = _tpl_by_id(self.lib, "DPPC")
        # 不完整
        atoms = _atoms_from_template(tpl, "L", "1")[:8]
        _load_beads("inc", atoms)
        det = cg.detect_cg("inc", cmd)
        m = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m.matched_beads, 0)

        # 重复珠名
        atoms2 = _atoms_from_template(tpl, "L", "2")
        atoms2[1].name = atoms2[0].name
        _load_beads("dup", atoms2)
        det2 = cg.detect_cg("dup", cmd)
        m2 = cgt.match_templates(det2.residues, self.lib)
        self.assertEqual(m2.matched_beads, 0)

        # 错误残基名
        atoms3 = _atoms_from_template(tpl, "L", "3", resn="XXXX")
        _load_beads("badn", atoms3)
        det3 = cg.detect_cg("badn", cmd)
        # 可能被 classify 为 other 而丢弃；即使留下也不应匹配模板
        m3 = cgt.match_templates(det3.residues, self.lib)
        self.assertEqual(m3.molecule_counts.get("DPPC", 0), 0)


class TestRnaChainMatching(unittest.TestCase):
    """RNA 链边界与精确 10mer 规则。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cgt.reset_registry_cache()
        self.lib = cgt.load_registry()
        self.rna = _tpl_by_id(self.lib, "RNA")

    def test_exact_polyu10_with_offset(self) -> None:
        """恰好 10 个连续 U（可偏移）应匹配。"""
        atoms = _rna_chain_atoms(self.rna, "R", start_resi=100)
        _load_beads("u10", atoms)
        det = cg.detect_cg("u10", cmd)
        m = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m.molecule_counts.get("RNA"), 1)
        self.assertEqual(m.matched_beads, 80)
        self.assertEqual(len(m.pairs), 89)
        self.assertEqual(len(m.matched_keys), 10)

    def test_reject_gap_twenty_and_cross_chain(self) -> None:
        """缺号、20mer、跨链均拒绝或不跨连。"""
        # 10 个 U 但残基号有缺口（1–4,6–11）
        atoms_gap = _rna_chain_atoms(self.rna, "A", start_resi=1)
        for at in atoms_gap:
            ri = int(at.resi)
            if ri >= 5:
                at.resi = str(ri + 1)
        _load_beads("gap", atoms_gap)
        m_gap = cgt.match_templates(cg.detect_cg("gap", cmd).residues, self.lib)
        self.assertEqual(m_gap.molecule_counts.get("RNA", 0), 0)

        # 20mer：两份 10mer 同链连续 → 共 20 残基，不得切成两个
        a1 = _rna_chain_atoms(self.rna, "B", start_resi=1, x0=0.0)
        a2 = _rna_chain_atoms(self.rna, "B", start_resi=11, x0=200.0)
        _load_beads("u20", a1 + a2)
        m20 = cgt.match_templates(cg.detect_cg("u20", cmd).residues, self.lib)
        self.assertEqual(m20.molecule_counts.get("RNA", 0), 0)

        # 两链各 10mer：各自匹配，不跨链连边
        c1 = _rna_chain_atoms(self.rna, "C", start_resi=1, x0=0.0)
        c2 = _rna_chain_atoms(self.rna, "D", start_resi=1, x0=500.0)
        _load_beads("two", c1 + c2)
        det = cg.detect_cg("two", cmd)
        m2 = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m2.molecule_counts.get("RNA"), 2)
        chain_of = {}
        for res in det.residues:
            for b in res.beads:
                chain_of[b.model_idx] = res.chain
        for i, j in m2.pairs:
            self.assertEqual(chain_of[i], chain_of[j])

    def test_insertion_codes_do_not_match_integer_chain_template(self):
        atoms = _rna_chain_atoms(self.rna, "R")
        for atom in atoms:
            if atom.resi == "5":
                atom.resi = "5A"
        _load_beads("insertion", atoms)
        match = cgt.match_templates(cg.detect_cg("insertion", cmd).residues, self.lib)
        self.assertEqual(match.matched_beads, 0)


class TestAmbiguousAndMixed(unittest.TestCase):
    """模棱两可拒绝；部分模板 + 几何混合；关闭自动模板。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cgt.reset_registry_cache()
        self.lib = cgt.load_registry()

    def test_ambiguous_duplicate_templates_rejected(self) -> None:
        """同一残基命中两个不同模板时拒绝。"""
        base = _tpl_by_id(self.lib, "DPPC")
        twin = copy.deepcopy(base)
        twin.id = "DPPC_TWIN"
        twin.molecule = "DPPC_TWIN"
        lib = cgt.TemplateLibrary(
            templates=[base, twin],
            package_ids=["test"],
            ok=True,
        )
        atoms = _atoms_from_template(base, "L", "1")
        _load_beads("amb", atoms)
        det = cg.detect_cg("amb", cmd)
        m = cgt.match_templates(det.residues, lib)
        self.assertEqual(m.matched_beads, 0)
        self.assertTrue(any("多个模板" in w for w in m.warnings))

    def test_partial_template_plus_geometry_mixed(self) -> None:
        """完整 DPPC 走模板，不完整脂走几何；溯源为混合。"""
        dppc = _tpl_by_id(self.lib, "DPPC")
        atoms = _atoms_from_template(dppc, "L", "1", x0=0.0)
        # 迷你脂（不完整）放远处
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A", "C2A"]):
            atoms.append(
                _bead(nm, "DPPC", "L", "2", 100.0 + float(i) * 2.0, 0.0, 0.0)
            )
        _load_beads("mix", atoms)
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("mix", mode=cg.MODE_BALLSTICK, auto_templates=True)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(rep.provenance, cg.PROVENANCE_MIXED)
        self.assertEqual(rep.template_counts.get("DPPC"), 1)
        self.assertGreater(rep.matched_beads, 0)
        self.assertGreater(rep.approximate_beads, 0)
        # 匹配与未匹配之间无跨残基猜测键：检查无跨 resi 边
        det = ctrl.last_detection
        idx_resi = {}
        for res in det.residues:
            for b in res.beads:
                idx_resi[b.model_idx] = res.resi
        model = cmd.get_model(rep.owned_objects[0])
        for b in model.bond:
            i, j = int(b.index[0]), int(b.index[1])
            self.assertEqual(idx_resi[i], idx_resi[j])

    def test_auto_templates_off_pure_geometry(self) -> None:
        """关闭自动模板时完整 DPPC 也走纯近似。"""
        dppc = _tpl_by_id(self.lib, "DPPC")
        _load_beads("geo", _atoms_from_template(dppc, "L", "1"))
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("geo", auto_templates=False)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(rep.provenance, cg.PROVENANCE_APPROX)
        self.assertEqual(rep.matched_beads, 0)
        self.assertIn("非力场拓扑", rep.provenance_label)


class TestLifecycleWithTemplates(unittest.TestCase):
    """源不变与 restore 清理含 guide。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cgt.reset_registry_cache()
        self.lib = cgt.load_registry()

    def test_source_unchanged_and_restore_cleans_guide(self) -> None:
        """连续模式 RNA guide 须登记到 owned，restore 全部清理。"""
        rna = _tpl_by_id(self.lib, "RNA")
        _load_beads("src", _rna_chain_atoms(rna, "R", start_resi=1))
        n0 = cmd.count_atoms("src")
        c0 = [
            (a.name, tuple(a.coord))
            for a in cmd.get_model("src").atom
        ]
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("src", mode=cg.MODE_CONTINUOUS, auto_templates=True)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(rep.provenance, cg.PROVENANCE_SAVED)
        self.assertEqual(cmd.count_atoms("src"), n0)
        c1 = [
            (a.name, tuple(a.coord))
            for a in cmd.get_model("src").atom
        ]
        self.assertEqual(c0, c1)
        owned = list(rep.owned_objects)
        self.assertTrue(any(cg.is_owned_name(n) for n in owned))
        # RNA 连续模式应生成 guide，且必须在 owned 中以便 restore 清理
        guide_names = [
            n for n in (cmd.get_names("all") or [])
            if n.startswith(cg.OWNED_PREFIX) and "guide" in n
        ]
        self.assertTrue(guide_names, "应生成 guide CGO")
        for n in guide_names:
            self.assertIn(n, owned)
        ctrl.restore()
        for n in owned:
            self.assertNotIn(n, cmd.get_names("all") or [])
        enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertIn("src", enabled)

    def test_chol_isolated_indices_reported(self) -> None:
        """CHOL 虚拟孤立珠应出现在 isolated_indices。"""
        chol = _tpl_by_id(self.lib, "CHOL")
        _load_beads("chol", _atoms_from_template(chol, "", "1"))
        det = cg.detect_cg("chol", cmd)
        m = cgt.match_templates(det.residues, self.lib)
        self.assertEqual(m.matched_beads, 8)
        self.assertEqual(len(m.isolated_indices), 3)
        self.assertEqual(len(m.pairs), 6)


@unittest.skipUnless(os.path.isfile(_SAMPLE_PDB), "样本 PDB 不存在，跳过")
class TestSampleFullTemplateCoverage(unittest.TestCase):
    """用户样本应被模板完整覆盖。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cgt.reset_registry_cache()

    def test_sample_counts_and_edges(self) -> None:
        """13320 珠全覆盖；分子计数与边数符合预期。"""
        cmd.load(_SAMPLE_PDB, "solute", quiet=1)
        det = cg.detect_cg("solute", cmd)
        self.assertEqual(det.bead_count, 13320)
        m = cgt.match_templates(det.residues)
        self.assertEqual(m.matched_beads, 13320)
        self.assertEqual(
            m.molecule_counts,
            {
                "AL15": 20,
                "AL31": 500,
                "CHOL": 380,
                "DPPC": 100,
                "RNA": 12,
            },
        )
        # 距离过滤前：55*20+13*500+11*100+6*380+89*12 = 12048
        self.assertEqual(len(m.pairs), 12048)
        # 虚拟点未伪造为键：CHOL 每分子 3，RNA 每核苷酸 1。
        self.assertEqual(len(m.isolated_indices), 380 * 3 + 120)

        n0 = cmd.count_atoms("solute")
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("solute", mode=cg.MODE_CONTINUOUS, auto_templates=True)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(cmd.count_atoms("solute"), n0)
        self.assertEqual(rep.provenance, cg.PROVENANCE_SAVED)
        self.assertIn("已保存", rep.provenance_label)
        self.assertNotIn("Martini3", rep.provenance_label)
        self.assertNotIn("Martini 3", rep.message)
        self.assertEqual(rep.matched_beads, 13320)
        self.assertEqual(rep.approximate_beads, 0)
        ctrl.restore()
        self.assertIn("solute", cmd.get_names("objects", enabled_only=1) or [])


if __name__ == "__main__":
    unittest.main()
