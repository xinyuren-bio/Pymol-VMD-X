# ==================================================
# 功能说明：粗粒化连通显示后端检测、键图、Apply/Restore 与键文件校验测试。
# 使用方法：由仓库测试运行器或 python -m unittest 发现；禁止依赖 pytest。
# 依赖环境：已安装的 pymol；可选样本 PDB 路径存在时做冒烟测试。
# 生成时间：2026-11-14 13:36
# ==================================================

"""coarse_grained 模块单元测试（标准 unittest）。"""

from __future__ import annotations

import os
import tempfile
import unittest

from chempy import Atom, Bond
from chempy.models import Indexed

from pymol import cmd
from pymol import coarse_grained as cg

# 用户提供的只读样本（存在时才跑冒烟）
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
    """构造保留完整残基名的 CG 珠原子。"""
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
    """用 chempy 模型加载对象（避免 PDB 三字符残基名截断）。"""
    model = Indexed()
    for at in atoms:
        model.add_atom(at)
    if name in (cmd.get_names("all") or []):
        cmd.delete(name)
    cmd.load_model(model, name, zoom=0)


def _coords(obj: str):
    """返回 [(index, name, x,y,z), ...]。"""
    out = []
    space = {"out": out}

    def collect(index, name, x, y, z):
        out.append((int(index), str(name), float(x), float(y), float(z)))

    space["collect"] = collect
    cmd.iterate_state(1, obj, "collect(index, name, x, y, z)", space=space)
    return out


def _bond_set(obj: str):
    """源对象键的 frozenset（model 0-based 下标对）。"""
    model = cmd.get_model(obj)
    s = set()
    for b in model.bond:
        i, j = int(b.index[0]), int(b.index[1])
        s.add((i, j) if i < j else (j, i))
    return s


class TestCGSnapshotImport(unittest.TestCase):
    def test_marked_pdb_preserves_conect_and_displays_coarse_beads(self):
        cmd.reinitialize()
        cmd.set('connect_mode', 0)
        cmd.set('pdb_unbond_cations', 1)
        cmd.set('auto_show_classified', 1)

        rows = [
            ('BB1', 'U', 'A', 0.0, 'C'),
            ('BB2', 'U', 'A', 4.5, 'C'),
            ('SC3', 'U', 'A', 8.0, 'C'),
            ('SC4', 'U', 'A', 9.4, 'C'),
            ('GL1', 'AL31', 'B', 20.0, 'C'),
            ('K2A', 'AL31', 'B', 24.5, ''),
            ('C3A', 'AL31', 'B', 29.0, 'C'),
        ]
        records = [
            f'ATOM  {serial:5d} {name:>4} {resn:<4}{chain}   1    '
            f'{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {elem:>2}'
            for serial, (name, resn, chain, x, elem) in enumerate(rows, 1)
        ]
        contents = '\n'.join(records + ['CONECT    1    2', 'CONECT    5    6',
                                      'CONECT    6    7', 'END', ''])
        cmd.read_pdbstr('REMARK   EXPORTED CG SOLUTE BEADS\n' + contents, 'cg_snapshot')
        model = cmd.get_model('cg_snapshot')
        serial_pairs = {tuple(sorted(model.atom[i].id for i in b.index)) for b in model.bond}
        self.assertEqual(serial_pairs, {(1, 2), (5, 6), (6, 7)})
        self.assertEqual(cmd.count_atoms('cg_snapshot and resn U and rep sticks'), 2)
        self.assertEqual(cmd.count_atoms('cg_snapshot and resn U and rep nb_spheres'), 2)
        self.assertEqual(cmd.count_atoms('cg_snapshot and rep cartoon'), 0)

        cmd.read_pdbstr(contents, 'ordinary_pdb')
        model = cmd.get_model('ordinary_pdb')
        serial_pairs = {tuple(sorted(model.atom[i].id for i in b.index)) for b in model.bond}
        self.assertNotIn((5, 6), serial_pairs)
        self.assertIn((3, 4), serial_pairs)


class TestClassifyAndDetect(unittest.TestCase):
    """残基分类与检测。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_classify_lipid_rna_protein_sterol(self) -> None:
        """已知珠名应分到正确种类。"""
        k, c = cg.classify_residue("DPPC", ["NC3", "PO4", "GL1", "C1A", "C1B"])
        self.assertEqual(k, cg.KIND_LIPID)
        self.assertEqual(c, "dppc")
        k, c = cg.classify_residue("AL31", ["POH", "NC4", "GL1", "C1A", "K2A"])
        self.assertEqual(k, cg.KIND_LIPID)
        self.assertEqual(c, "al31")
        k, c = cg.classify_residue("AL15M", ["LNK", "C1A", "EC1", "0EC1"])
        self.assertEqual(k, cg.KIND_LIPID)
        self.assertEqual(c, "al15")
        k, c = cg.classify_residue("U", ["BB1", "BB2", "BB3", "SC1", "SC2"])
        self.assertEqual(k, cg.KIND_NA)
        self.assertEqual(c, "rna")
        k, c = cg.classify_residue("ALA", ["BB", "SC1"])
        self.assertEqual(k, cg.KIND_PROTEIN)
        k, c = cg.classify_residue("CHOL", ["ROH", "R1", "R2", "C1", "C2"])
        self.assertEqual(k, cg.KIND_STEROL)
        self.assertEqual(c, "chol")

    def test_detect_synthetic_mixture(self) -> None:
        """合成脂+RNA+蛋白应被检出；水被排除。"""
        atoms = []
        for ri, x0 in ((1, 0.0), (2, 5.0)):
            for nm, dx in (("BB1", 0), ("BB2", 1.5), ("BB3", 3.0), ("SC1", 1.0), ("SC2", 2.0)):
                atoms.append(_bead(nm, "U", "A", ri, x0 + dx, 0.0, 0.0))
        for nm, x in (
            ("NC3", 20), ("PO4", 22), ("GL1", 24), ("GL2", 24.5),
            ("C1A", 26), ("C2A", 28), ("C1B", 26.5), ("C2B", 28.5),
        ):
            atoms.append(_bead(nm, "DPPC", "", 10, x, 10.0, 0.0))
        for ri, x0 in ((1, 40.0), (2, 44.0)):
            atoms.append(_bead("BB", "ALA", "P", ri, x0, 0.0, 0.0))
            atoms.append(_bead("SC1", "ALA", "P", ri, x0 + 1.0, 1.0, 0.0))
        atoms.append(_bead("OW", "SOL", "W", 1, 100.0, 0.0, 0.0))
        _load_beads("mix", atoms)

        det = cg.detect_cg("mix", cmd)
        self.assertFalse(det.empty)
        self.assertGreaterEqual(det.by_kind.get(cg.KIND_NA, 0), 2)
        self.assertGreaterEqual(det.by_kind.get(cg.KIND_LIPID, 0), 1)
        self.assertGreaterEqual(det.by_kind.get(cg.KIND_PROTEIN, 0), 2)
        self.assertGreater(det.excluded_beads, 0)


class TestBondGraph(unittest.TestCase):
    """键图：连续、无分子间、无超长。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_lipid_mst_connected_no_intermolecular(self) -> None:
        """两脂分子各自连通，且无跨分子键。"""
        atoms = []
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A", "C2A"]):
            atoms.append(_bead(nm, "DPPC", "L", 1, float(i) * 2.0, 0.0, 0.0))
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A", "C2A"]):
            atoms.append(_bead(nm, "DPPC", "L", 2, 100.0 + float(i) * 2.0, 0.0, 0.0))
        _load_beads("lips", atoms)
        det = cg.detect_cg("lips", cmd)
        pairs, skipped, _ = cg.build_approximate_bonds(det.residues)
        self.assertGreaterEqual(len(pairs), 4)
        idx_resi = {}
        for res in det.residues:
            for b in res.beads:
                idx_resi[b.model_idx] = res.resi
        model = cmd.get_model("lips")
        for i, j in pairs:
            self.assertEqual(idx_resi[i], idx_resi[j], "不应出现分子间键")
            d = cg._dist(model.atom[i].coord, model.atom[j].coord)
            self.assertLessEqual(d, cg.INTRA_CAP + 1e-6)

    def test_rna_and_protein_backbone_continuity(self) -> None:
        """连续残基号应产生主链连接；间隔残基不连。"""
        atoms = []
        for ri, x0 in ((1, 0.0), (2, 4.0), (3, 8.0), (5, 20.0)):
            for nm, dx in (("BB1", 0.0), ("BB2", 1.2), ("BB3", 2.4), ("SC1", 1.0)):
                atoms.append(_bead(nm, "U", "A", ri, x0 + dx, 0.0, 0.0))
        for ri, x0 in ((10, 50.0), (11, 54.0)):
            atoms.append(_bead("BB", "GLY", "B", ri, x0, 5.0, 0.0))
            atoms.append(_bead("SC1", "GLY", "B", ri, x0, 6.0, 0.0))
        _load_beads("poly", atoms)
        det = cg.detect_cg("poly", cmd)
        pairs, _, _ = cg.build_approximate_bonds(det.residues)
        pairset = set(pairs)

        by_key = {}
        for res in det.residues:
            for b in res.beads:
                by_key[(res.chain, res.resi, cg._norm_name(b.name))] = b.model_idx

        def has_pair(a, b) -> bool:
            i, j = (a, b) if a < b else (b, a)
            return (i, j) in pairset

        self.assertTrue(has_pair(by_key[("A", "1", "BB3")], by_key[("A", "2", "BB1")]))
        self.assertFalse(has_pair(by_key[("A", "3", "BB3")], by_key[("A", "5", "BB1")]))
        self.assertTrue(has_pair(by_key[("B", "10", "BB")], by_key[("B", "11", "BB")]))

    def test_long_bonds_skipped(self) -> None:
        """超距珠不应被 MST 连上拉丝。"""
        atoms = [
            _bead("NC3", "DPPC", "X", 1, 0.0, 0.0, 0.0),
            _bead("PO4", "DPPC", "X", 1, 2.0, 0.0, 0.0),
            _bead("C4A", "DPPC", "X", 1, 50.0, 0.0, 0.0),
        ]
        _load_beads("far", atoms)
        det = cg.detect_cg("far", cmd)
        pairs, skipped, _ = cg.build_approximate_bonds(det.residues)
        model = cmd.get_model("far")
        for i, j in pairs:
            self.assertLessEqual(
                cg._dist(model.atom[i].coord, model.atom[j].coord),
                cg.INTRA_CAP + 1e-6,
            )
        self.assertGreater(skipped, 0)


class TestSourceUnchangedAndLifecycle(unittest.TestCase):
    """源不变性、重复 apply/restore、禁用对象安全。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def _mini_lipid(self, name: str = "src") -> None:
        """加载小型脂对象。"""
        atoms = []
        for i, nm in enumerate(
            ["NC3", "PO4", "GL1", "GL2", "C1A", "C2A", "C1B", "C2B"]
        ):
            atoms.append(_bead(nm, "DPPC", "A", 1, float(i) * 2.0, 0.0, 1.0))
        _load_beads(name, atoms)

    def test_source_counts_coords_bonds_unchanged(self) -> None:
        """apply 后源原子数、坐标、键不变。"""
        self._mini_lipid("src")
        # 人为加一条源键，确认不被改写
        model = cmd.get_model("src")
        b = Bond()
        b.index = [0, 1]
        b.order = 1
        model.bond = [b]
        cmd.delete("src")
        cmd.load_model(model, "src", zoom=0)

        n0 = cmd.count_atoms("src")
        c0 = _coords("src")
        b0 = _bond_set("src")
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("src", mode=cg.MODE_BALLSTICK)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(cmd.count_atoms("src"), n0)
        c1 = _coords("src")
        self.assertEqual(len(c0), len(c1))
        for a, b in zip(c0, c1):
            self.assertEqual(a[0], b[0])
            self.assertEqual(a[1], b[1])
            self.assertAlmostEqual(a[2], b[2], places=3)
            self.assertAlmostEqual(a[3], b[3], places=3)
            self.assertAlmostEqual(a[4], b[4], places=3)
        self.assertEqual(_bond_set("src"), b0)
        self.assertTrue(any(cg.is_owned_name(x) for x in rep.owned_objects))
        enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertNotIn("src", enabled)

    def test_repeat_apply_restore(self) -> None:
        """重复 apply 只清理自有输出；restore 恢复源可见。"""
        self._mini_lipid("src")
        ctrl = cg.CGDisplayController(cmd)
        r1 = ctrl.apply("src", mode=cg.MODE_CONTINUOUS)
        self.assertTrue(r1.ok, r1.message)
        owned1 = list(r1.owned_objects)
        r2 = ctrl.apply("src", mode=cg.MODE_BALLSTICK)
        self.assertTrue(r2.ok, r2.message)
        for n in cmd.get_names("public_objects", enabled_only=0) or []:
            if n != "src":
                self.assertTrue(cg.is_owned_name(n), n)
        msg = ctrl.restore()
        self.assertIn("恢复", msg)
        enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertIn("src", enabled)
        for n in owned1 + list(r2.owned_objects):
            self.assertNotIn(n, cmd.get_names("all") or [])

    def test_preexisting_output_name_is_not_overwritten(self):
        self._mini_lipid("src")
        cmd.pseudoatom("cgvis_src_mol", pos=[9, 9, 9])
        before = _coords("cgvis_src_mol")
        ctrl = cg.CGDisplayController(cmd)
        self.assertTrue(ctrl.apply("src").ok)
        ctrl.restore()
        self.assertEqual(_coords("cgvis_src_mol"), before)

    def test_switching_source_restores_previous_visibility(self):
        self._mini_lipid("source_a")
        self._mini_lipid("source_b")
        ctrl = cg.CGDisplayController(cmd)
        self.assertTrue(ctrl.apply("source_a").ok)
        self.assertTrue(ctrl.apply("source_b").ok)
        self.assertIn("source_a", cmd.get_names("objects", enabled_only=1))
        ctrl.restore()
        self.assertIn("source_b", cmd.get_names("objects", enabled_only=1))

    def test_failed_rebuild_keeps_previous_display(self):
        from unittest.mock import patch
        self._mini_lipid("src")
        ctrl = cg.CGDisplayController(cmd)
        self.assertTrue(ctrl.apply("src").ok)
        before = list(ctrl.owned)
        with patch.object(cmd, "load_model", side_effect=RuntimeError("test failure")):
            self.assertFalse(ctrl.apply("src").ok)
        self.assertEqual(ctrl.owned, before)
        self.assertTrue(all(name in cmd.get_names("objects", enabled_only=1) for name in before))

    def test_disabled_source_still_listed_and_safe(self) -> None:
        """禁用源仍可选；处理自有副本应被拒绝。"""
        self._mini_lipid("src")
        cmd.disable("src")
        srcs = cg.list_source_objects(cmd)
        self.assertIn("src", srcs)
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("src", mode=cg.MODE_BALLSTICK)
        self.assertTrue(rep.ok, rep.message)
        owned = rep.owned_objects[0]
        det = cg.detect_cg(owned, cmd)
        self.assertTrue(det.empty)


class TestBondFileValidation(unittest.TestCase):
    """畸形键文件拒绝；合法键对接受。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_malformed_and_oob_rejected(self) -> None:
        """越界/非整数/空文件应抛 TopologyError。"""
        with tempfile.TemporaryDirectory() as td:
            empty = os.path.join(td, "empty.txt")
            with open(empty, "w", encoding="utf-8") as f:
                f.write("# only comment\n")
            with self.assertRaises(cg.TopologyError):
                cg.parse_bond_pairs_file(empty, 10)

            bad = os.path.join(td, "bad.txt")
            with open(bad, "w", encoding="utf-8") as f:
                f.write("0 99\n")
            with self.assertRaises(cg.TopologyError):
                cg.parse_bond_pairs_file(bad, 10)

            junk = os.path.join(td, "junk.txt")
            with open(junk, "w", encoding="utf-8") as f:
                f.write("foo bar\n")
            with self.assertRaises(cg.TopologyError):
                cg.parse_bond_pairs_file(junk, 10)

    def test_valid_pairs_and_itp_mismatch(self) -> None:
        """合法键对可读；含 include 的 itp 与原子数不匹配应失败。"""
        with tempfile.TemporaryDirectory() as td:
            ok = os.path.join(td, "ok.txt")
            with open(ok, "w", encoding="utf-8") as f:
                f.write("0 1\n1 2\n")
            pairs = cg.parse_bond_pairs_file(ok, 5)
            self.assertEqual(len(pairs), 2)

            itp = os.path.join(td, "x.itp")
            with open(itp, "w", encoding="utf-8") as f:
                f.write('#include "other.itp"\n[ atoms ]\n')
            with self.assertRaises(cg.TopologyError):
                cg.parse_itp_bonds(itp, 5)

            itp2 = os.path.join(td, "y.itp")
            with open(itp2, "w", encoding="utf-8") as f:
                f.write(
                    "[ moleculetype ]\nLIP 3\n"
                    "[ atoms ]\n"
                    "1 C 1 DPPC NC3 1 0 0\n"
                    "2 C 1 DPPC PO4 1 0 0\n"
                    "[ bonds ]\n"
                    "1 2 1 0.4 0\n"
                )
            with self.assertRaises(cg.TopologyError):
                cg.parse_itp_bonds(itp2, 99)

    def test_apply_rejects_bad_bonds_keeps_source_visible(self) -> None:
        """apply 遇坏键文件失败时源仍可见且无残留自有对象。"""
        atoms = [
            _bead("NC3", "DPPC", "A", 1, 0.0, 0.0, 0.0),
            _bead("PO4", "DPPC", "A", 1, 2.0, 0.0, 0.0),
            _bead("GL1", "DPPC", "A", 1, 4.0, 0.0, 0.0),
        ]
        _load_beads("src", atoms)
        with tempfile.TemporaryDirectory() as td:
            bad = os.path.join(td, "bad.txt")
            with open(bad, "w", encoding="utf-8") as f:
                f.write("0 100\n")
            ctrl = cg.CGDisplayController(cmd)
            rep = ctrl.apply("src", bonds_file=bad)
            self.assertFalse(rep.ok)
            enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
            self.assertIn("src", enabled)
            for n in cmd.get_names("all") or []:
                self.assertFalse(cg.is_owned_name(n), n)


@unittest.skipUnless(os.path.isfile(_SAMPLE_PDB), "样本 PDB 不存在，跳过冒烟")
class TestSamplePdbSmoke(unittest.TestCase):
    """对用户样本做检测/渲染冒烟（只读加载，不改原文件）。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_detect_and_render_sample(self) -> None:
        """样本应检出 RNA/脂/固醇量级并成功 apply/restore。"""
        cmd.load(_SAMPLE_PDB, "solute", quiet=1)
        det = cg.detect_cg("solute", cmd)
        self.assertFalse(det.empty)
        self.assertGreater(det.bead_count, 10000)
        self.assertGreaterEqual(det.by_component.get("rna", 0), 100)
        self.assertGreaterEqual(det.by_component.get("al31", 0), 400)
        self.assertGreaterEqual(det.by_component.get("al15", 0), 10)
        self.assertGreaterEqual(det.by_component.get("dppc", 0), 50)
        self.assertGreaterEqual(det.by_component.get("chol", 0), 300)

        n0 = cmd.count_atoms("solute")
        ctrl = cg.CGDisplayController(cmd)
        rep = ctrl.apply("solute", mode=cg.MODE_CONTINUOUS, auto_templates=True)
        self.assertTrue(rep.ok, rep.message)
        self.assertEqual(cmd.count_atoms("solute"), n0)
        # 完整样本应由已保存连接模板覆盖，不再强制断言「非力场拓扑」
        self.assertEqual(rep.provenance, cg.PROVENANCE_SAVED)
        self.assertIn("已保存", rep.provenance_label)
        self.assertEqual(rep.matched_beads, det.bead_count)
        self.assertEqual(rep.approximate_beads, 0)
        ctrl.restore()
        enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertIn("solute", enabled)


class TestModeValidation(unittest.TestCase):
    """模式校验与摘要。"""

    def test_validate_mode(self) -> None:
        """中英文模式名可归一。"""
        self.assertEqual(cg.validate_mode("连续形态"), cg.MODE_CONTINUOUS)
        self.assertEqual(cg.validate_mode("球棍连接"), cg.MODE_BALLSTICK)
        with self.assertRaises(ValueError):
            cg.validate_mode("magic")


if __name__ == "__main__":
    unittest.main()
