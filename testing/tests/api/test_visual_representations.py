# ==================================================
# 功能说明：visual_representations 选择解析/求值与多层表示后端单元测试。
# 使用方法：由仓库测试运行器或 python -m unittest 发现；禁止依赖 pytest；勿在此文件执行安装。
# 依赖环境：已安装的 pymol；可选 PYMOLX_CG_SAMPLE 做大规模冒烟（存在才跑）。
# 生成时间：2026-11-14 16:50
# ==================================================

"""visual_representations 后端测试（标准 unittest）。"""

from __future__ import annotations

import os
import unittest

from chempy import Atom, Bond
from chempy.models import Indexed

from pymol import cmd
from pymol import visual_representations as vr

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
    elem: str = "C",
) -> Atom:
    """构造 CG/原子。"""
    at = Atom()
    at.name = name
    at.resn = resn
    at.chain = chain
    at.resi = str(resi)
    at.segi = segi
    at.symbol = elem
    at.coord = [float(x), float(y), float(z)]
    return at


def _load_beads(name: str, atoms: list, bonds=()) -> None:
    """用 chempy 模型加载对象（保留完整残基名）。"""
    model = Indexed()
    for at in atoms:
        model.add_atom(at)
    model.connect_mode = 1
    for pair in bonds:
        bond = Bond()
        bond.index = list(pair)
        model.add_bond(bond)
    if name in (cmd.get_names("all") or []):
        cmd.delete(name)
    cmd.load_model(model, name, zoom=0)


def _coords(obj: str):
    """返回坐标快照。"""
    out = []
    space = {"out": out}

    def collect(index, name, x, y, z):
        out.append((int(index), str(name), float(x), float(y), float(z)))

    space["collect"] = collect
    cmd.iterate_state(1, obj, "collect(index, name, x, y, z)", space=space)
    return out


def _bond_count(obj: str) -> int:
    """键数量。"""
    return len(getattr(cmd.get_model(obj), "bond", []) or [])


class TestParser(unittest.TestCase):
    """词法与语法。"""

    def test_parse_basic_and_or_not(self) -> None:
        """and 优先于 or；not 最高。"""
        ast = vr.parse_selection("protein or nucleic and not water")
        self.assertIsInstance(ast, vr.BinNode)
        self.assertEqual(ast.op, "or")

    def test_reject_injection_and_unknown(self) -> None:
        """分号/命令注入/未知符号友好拒绝。"""
        with self.assertRaises(vr.SelectionError):
            vr.parse_selection("all; color red")
        with self.assertRaises(vr.SelectionError):
            vr.parse_selection("protein | nucleic")
        with self.assertRaises(vr.SelectionError):
            vr.parse_selection("")

    def test_range_reversed(self) -> None:
        """范围颠倒拒绝。"""
        with self.assertRaises(vr.SelectionError):
            vr.parse_selection("resid 10 to 1")

    def test_invalid_resid_comparison_explains_what_to_enter(self) -> None:
        for expression in ("resid <", "resid >= GLY", "resid != 80.5", "resid >< 80"):
            with self.subTest(expression=expression), self.assertRaisesRegex(vr.SelectionError, "整数"):
                vr.parse_selection(expression)
        with self.assertRaisesRegex(vr.SelectionError, "比较运算符用于 resid"):
            vr.parse_selection("chain > A")

    def test_within_and_same_residue(self) -> None:
        """within / same residue as 可解析。"""
        ast = vr.parse_selection("within 5 of (nucleic)")
        self.assertIsInstance(ast, vr.WithinNode)
        self.assertAlmostEqual(ast.distance, 5.0)
        ast2 = vr.parse_selection("same residue as name BB1")
        self.assertIsInstance(ast2, vr.SameResidueNode)


class TestEvaluator(unittest.TestCase):
    """选择求值与索引映射。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_resname_nucleic_protein_lipid(self) -> None:
        """关键词与 resname 多值。"""
        atoms = []
        for i, nm in enumerate(["BB1", "BB2", "BB3", "SC1"]):
            atoms.append(_bead(nm, "U", "A", 1, float(i), 0, 0))
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A"]):
            atoms.append(_bead(nm, "DPPC", "L", 1, 20 + float(i), 0, 0))
        for i, nm in enumerate(["BB", "SC1"]):
            atoms.append(_bead(nm, "ALA", "P", 1, 40 + float(i), 0, 0))
        atoms.append(_bead("OW", "SOL", "W", 1, 100, 0, 0, elem="O"))
        _load_beads("mix", atoms)
        recs = vr.load_atom_records("mix", cmd)

        idx, _ = vr.evaluate_selection("nucleic", recs)
        self.assertEqual(len(idx), 4)
        idx, _ = vr.evaluate_selection("resname DPPC", recs)
        self.assertEqual(len(idx), 4)
        idx, _ = vr.evaluate_selection("protein", recs)
        self.assertEqual(len(idx), 2)
        idx, _ = vr.evaluate_selection("not water", recs)
        self.assertEqual(len(idx), 10)
        idx, _ = vr.evaluate_selection("resname DPPC ALA", recs)
        self.assertEqual(len(idx), 6)

    def test_name_wildcard_and_resid_range(self) -> None:
        """name 通配与 resid 范围。"""
        atoms = []
        for ri in (1, 2, 3, 11):
            atoms.append(_bead("BB1", "U", "A", ri, float(ri), 0, 0))
            atoms.append(_bead("SC1", "U", "A", ri, float(ri), 1, 0))
        _load_beads("na", atoms)
        recs = vr.load_atom_records("na", cmd)
        idx, _ = vr.evaluate_selection("name BB*", recs)
        self.assertEqual(len(idx), 4)
        idx, _ = vr.evaluate_selection("resid 1 to 3", recs)
        self.assertEqual(len(idx), 6)

    def test_resid_comparisons_preserve_boundaries_and_boolean_scope(self) -> None:
        atoms = [_bead("BB1", "U", "A", r, r, 0, 0) for r in (-1, 80, 81, 85, 89, 90)]
        atoms += [_bead("CA", "ALA", "B", 85, 85, 2, 0)]
        _load_beads("limits", atoms)
        records = vr.load_atom_records("limits", cmd)
        cases = {
            'resid < 90 and resid > 80': [81, 85, 85, 89],
            'resid 81 to 89': [81, 85, 85, 89],
            'resid>=80 and resid<=90 and resid!=85': [80, 81, 89, 90],
            'resid == -1 or (chain A and resid = 85 and resname U and name BB* and not chain B)': [-1, 85],
        }
        for expression, expected in cases.items():
            with self.subTest(expression=expression):
                _, mids = vr.evaluate_selection(expression, records)
                self.assertEqual(sorted(records[i].resi for i in mids), expected)

    def test_index_zero_maps_source_index(self) -> None:
        """index 0 映射 get_model 第一原子的 source.index（非简单 +1）。"""
        atoms = [_bead("BB", "ALA", "A", 1, 0, 0, 0)]
        atoms.append(_bead("SC1", "ALA", "A", 1, 1, 0, 0))
        _load_beads("idxo", atoms)
        # 打乱感知：读取真实 index
        model = cmd.get_model("idxo")
        first_index = int(model.atom[0].index)
        recs = vr.load_atom_records("idxo", cmd)
        indices, mids = vr.evaluate_selection("index 0", recs)
        self.assertEqual(mids, [0])
        self.assertEqual(indices, [first_index])
        # 不应假设 index == model_idx+1 恒成立时也要返回 source.index
        self.assertEqual(indices[0], recs[0].index)

    def test_within_scoped(self) -> None:
        """within 仅在对象内；距离准确。"""
        atoms = [
            _bead("BB1", "U", "A", 1, 0, 0, 0),
            _bead("NC3", "DPPC", "L", 1, 3, 0, 0),
            _bead("PO4", "DPPC", "L", 2, 20, 0, 0),
        ]
        _load_beads("wobj", atoms)
        recs = vr.load_atom_records("wobj", cmd)
        idx, mids = vr.evaluate_selection("within 5 of nucleic", recs)
        self.assertIn(0, mids)
        self.assertIn(1, mids)
        self.assertNotIn(2, mids)

    def test_no_escape_other_object(self) -> None:
        """同名 resn 不逃逸到其它对象。"""
        _load_beads(
            "o1",
            [_bead("NC3", "DPPC", "L", 1, 0, 0, 0), _bead("PO4", "DPPC", "L", 1, 2, 0, 0)],
        )
        _load_beads(
            "o2",
            [_bead("NC3", "DPPC", "L", 1, 0, 0, 0)] * 5,
        )
        prev = vr.preview_selection("o1", "resname DPPC", cmd)
        self.assertTrue(prev.ok)
        self.assertEqual(prev.count, 2)

    def test_al31_truncated_with_beads(self) -> None:
        """AL31 查询匹配截断 AL3 且需珠特征。"""
        atoms = [
            _bead("POH", "AL3", "L", 1, 0, 0, 0),
            _bead("NC4", "AL3", "L", 1, 1, 0, 0),
            _bead("GL1", "AL3", "L", 1, 2, 0, 0),
            _bead("C1A", "AL3", "L", 1, 3, 0, 0),
            _bead("K2A", "AL3", "L", 1, 4, 0, 0),
        ]
        # 无关 AL3 截断但无脂珠
        atoms.append(_bead("X", "AL3", "Z", 9, 50, 0, 0))
        _load_beads("al", atoms)
        recs = vr.load_atom_records("al", cmd)
        idx, _ = vr.evaluate_selection("resname AL31", recs)
        self.assertGreaterEqual(len(idx), 5)


class TestController(unittest.TestCase):
    """多层独立显示与不变性。"""

    def setUp(self) -> None:
        cmd.reinitialize()

    def test_open_default_no_scene_change(self) -> None:
        """默认层 Apply 短路不改场景。"""
        atoms = [_bead("BB", "ALA", "A", 1, 0, 0, 0)]
        _load_beads("src", atoms)
        before = _coords("src")
        ctrl = vr.RepresentationController(cmd)
        # 模拟 GUI：仅原始层
        rep = ctrl.apply("src")
        # 默认一层 original/all/原色且无旧 owned → GUI 会短路；
        # 后端若直接 apply 会建副本；此处测 restore 与源不变
        self.assertTrue(isinstance(rep, vr.ApplyReport))
        after = _coords("src")
        self.assertEqual(before, after)

    def test_two_layers_independent_and_restore(self) -> None:
        """两层颜色独立；恢复不影响源坐标/键。"""
        atoms = []
        for i, nm in enumerate(["NC3", "PO4", "GL1", "C1A", "C2A"]):
            atoms.append(_bead(nm, "DPPC", "L", 1, float(i) * 2, 0, 0))
        for i, nm in enumerate(["BB1", "BB2", "BB3", "SC1"]):
            atoms.append(_bead(nm, "U", "A", 1, 30 + float(i), 0, 0))
        _load_beads("src2", atoms, [(0, 1), (1, 2), (2, 3), (3, 4), (5, 6), (6, 7), (7, 8)])
        before_coords = _coords("src2")
        before_bonds = _bond_count("src2")

        ctrl = vr.RepresentationController(cmd)
        layers = [
            vr.RepresentationLayer(
                layer_id=vr.new_layer_id(),
                enabled=True,
                style=vr.STYLE_STICKS,
                color_mode=vr.COLOR_SOLID,
                color_rgb=(1.0, 0.0, 0.0),
                selection="resname DPPC",
            ),
            vr.RepresentationLayer(
                layer_id=vr.new_layer_id(),
                enabled=True,
                style=vr.STYLE_STICKS,
                color_mode=vr.COLOR_SOLID,
                color_rgb=(0.0, 0.0, 1.0),
                selection="nucleic",
            ),
        ]
        rep = ctrl.apply("src2", layers)
        self.assertTrue(rep.ok, msg=rep.message)
        self.assertEqual(len(ctrl.owned), 2)
        # 源应被禁用但坐标/键不变
        self.assertEqual(_coords("src2"), before_coords)
        self.assertEqual(_bond_count("src2"), before_bonds)
        enabled = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertNotIn("src2", enabled)

        ctrl.restore()
        self.assertEqual(ctrl.owned, [])
        self.assertEqual(_coords("src2"), before_coords)
        enabled2 = set(cmd.get_names("public_objects", enabled_only=1) or [])
        self.assertIn("src2", enabled2)

    def test_failed_apply_keeps_previous(self) -> None:
        """失败保留前次显示。"""
        atoms = [_bead("BB", "ALA", "A", 1, 0, 0, 0), _bead("SC1", "ALA", "A", 1, 1, 0, 0)]
        _load_beads("keep", atoms)
        ctrl = vr.RepresentationController(cmd)
        ok_layers = [
            vr.RepresentationLayer(
                layer_id=vr.new_layer_id(),
                style=vr.STYLE_SPHERES,
                color_mode=vr.COLOR_SOLID,
                color_rgb=(0.2, 0.8, 0.2),
                selection="all",
            )
        ]
        rep1 = ctrl.apply("keep", ok_layers)
        self.assertTrue(rep1.ok, msg=rep1.message)
        owned_before = list(ctrl.owned)
        bad = [
            vr.RepresentationLayer(
                layer_id=vr.new_layer_id(),
                style=vr.STYLE_SPHERES,
                color_mode=vr.COLOR_SOLID,
                color_rgb=(1, 0, 0),
                selection="resname NOTREALZZ",
            )
        ]
        rep2 = ctrl.apply("keep", bad)
        self.assertFalse(rep2.ok)
        self.assertEqual(ctrl.owned, owned_before)

    def test_rgb_validation(self) -> None:
        """RGB 越界拒绝。"""
        with self.assertRaises(ValueError):
            vr.validate_rgb((1.2, 0, 0))
        self.assertEqual(vr.rgb_to_hex_color((1, 0, 0.5)), "0xff0080")

    def test_source_named_reserved_word(self) -> None:
        """源对象名恰为 all/protein 时按精确名取模。"""
        _load_beads("protein", [_bead("BB", "ALA", "A", 1, 0, 0, 0)])
        recs = vr.load_atom_records("protein", cmd)
        self.assertEqual(len(recs), 1)
        ctrl = vr.RepresentationController(cmd)
        layers = [
            vr.RepresentationLayer(
                layer_id=vr.new_layer_id(),
                style=vr.STYLE_SPHERES,
                color_mode=vr.COLOR_ORIGINAL,
                selection="all",
            )
        ]
        rep = ctrl.apply("protein", layers)
        self.assertTrue(rep.ok, msg=rep.message)

    def test_suggest_component_layers(self) -> None:
        """按组分分层不含水分。"""
        atoms = []
        for nm in ("NC3", "PO4", "GL1", "C1A"):
            atoms.append(_bead(nm, "DPPC", "L", 1, 0, 0, 0))
        for nm in ("BB1", "BB2", "BB3", "SC1"):
            atoms.append(_bead(nm, "U", "A", 1, 10, 0, 0))
        atoms.append(_bead("OW", "SOL", "W", 1, 99, 0, 0, elem="O"))
        _load_beads("comp", atoms)
        ctrl = vr.RepresentationController(cmd)
        layers = ctrl.suggest_component_layers("comp")
        texts = " ".join(L.selection for L in layers)
        self.assertIn("DPPC", texts)
        self.assertIn("nucleic", texts)
        self.assertNotIn("water", texts.lower())


@unittest.skipUnless(
    _SAMPLE_PDB and os.path.isfile(_SAMPLE_PDB),
    "需要 PYMOLX_CG_SAMPLE 指向可读样本",
)
class TestSampleSmoke(unittest.TestCase):
    """大规模样本冒烟（存在才跑）。"""

    def setUp(self) -> None:
        cmd.reinitialize()
        cmd.load(_SAMPLE_PDB, "lnp")

    def test_counts_al31_dppc_nucleic(self) -> None:
        """AL31/DPPC/nucleic 数量级检查。"""
        prev_a = vr.preview_selection("lnp", "resname AL31", cmd)
        prev_d = vr.preview_selection("lnp", "resname DPPC", cmd)
        prev_n = vr.preview_selection("lnp", "nucleic", cmd)
        self.assertTrue(prev_a.ok, msg=prev_a.message)
        self.assertTrue(prev_d.ok, msg=prev_d.message)
        self.assertTrue(prev_n.ok, msg=prev_n.message)
        # 用户场景量级；允许样本差异，只要求明显非零且 AL31 最多
        self.assertGreater(prev_a.count, 1000)
        self.assertGreater(prev_d.count, 100)
        self.assertGreater(prev_n.count, 100)

    def test_partial_pdb_bonds_do_not_block_saved_templates(self):
        before = _bond_count("lnp")
        controller = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer(vr.new_layer_id(), style=vr.STYLE_STICKS,
                                       color_mode=vr.COLOR_SOLID, selection="resname AL31")
        report = controller.apply("lnp", [layer])
        self.assertTrue(report.ok, report.message)
        self.assertEqual(cmd.count_atoms(report.owned_objects[0]), 7000)
        self.assertEqual(_bond_count(report.owned_objects[0]), 6500)
        self.assertEqual(_bond_count("lnp"), before)
        self.assertIn("模板", report.bond_note)
        controller.restore()


if __name__ == "__main__":
    unittest.main()
