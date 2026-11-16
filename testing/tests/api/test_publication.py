# ==================================================
# 功能说明：科研出图后端检测、多配体样式、三主题与会话恢复测试。
# 使用方法：在仓库根目录用本机 PyMOL/测试框架运行本文件。
# 依赖环境：本地已编译/安装的 pymol；标准 unittest（禁止依赖 pytest）。
# 生成时间：2026-11-14 13:30
# ==================================================

"""publication 模块检测与恢复测试。"""

from __future__ import annotations

import unittest
from pymol import cmd, publication as pub


def _pdb_het(resn: str, resi: int, chain: str, x: float, elem: str = "C") -> str:
    """生成单原子 HETATM 行，用于水/离子/添加剂占位。"""
    return (
        "HETATM%5d  %s   %s %s%4d    %8.3f%8.3f%8.3f  1.00 20.00           %s\n"
        % (resi, elem.ljust(3)[:3], resn.ljust(3)[:3], chain, resi, x, 0.0, 0.0, elem)
    )


def _load_additive_soup(name: str = "additives") -> None:
    """加载常见水、离子、添加剂（应被排除）。"""
    lines = [
        _pdb_het("HOH", 1, "A", 0.0, "O"),
        _pdb_het("HOH", 2, "A", 1.0, "O"),
        _pdb_het("NA", 3, "A", 2.0, "NA"),
        _pdb_het("CL", 4, "A", 3.0, "CL"),
        _pdb_het("GOL", 5, "A", 4.0, "C"),
        _pdb_het("EDO", 6, "A", 5.0, "C"),
        _pdb_het("SO4", 7, "A", 6.0, "S"),
        "END\n",
    ]
    cmd.read_pdbstr("".join(lines), name)


def _merge(dst: str, src: str) -> None:
    """将 src 原子并入 dst 对象后删除 src。"""
    tmp = cmd.get_unused_name("_merge_tmp")
    cmd.create(tmp, "(%s) or (%s)" % (dst, src))
    cmd.delete(dst)
    cmd.delete(src)
    cmd.set_name(tmp, dst)


class TestPublicationDetect(unittest.TestCase):
    """配体/蛋白启发式检测。"""

    def test_empty_scene(self):
        """空场景应标记 empty 并给出中文说明。"""
        cmd.reinitialize()
        info = pub.detect_scene(cmd)
        self.assertTrue(info.empty)
        self.assertFalse(info.has_protein)
        self.assertEqual(info.ligands, [])
        self.assertTrue(info.note)

    def test_protein_only(self):
        """仅蛋白：无配体候选，可用无配体总览。"""
        cmd.reinitialize()
        cmd.fab("ACDEFGHIK", "prot")
        info = pub.detect_scene(cmd)
        self.assertFalse(info.empty)
        self.assertTrue(info.has_protein)
        self.assertEqual(len(info.ligands), 0)

    def test_standalone_ligand(self):
        """独立小分子（无蛋白）应作为配体候选。"""
        cmd.reinitialize()
        cmd.fragment("benzene", "lig")
        info = pub.detect_scene(cmd)
        self.assertFalse(info.empty)
        self.assertFalse(info.has_protein)
        self.assertGreaterEqual(len(info.ligands), 1)
        self.assertGreaterEqual(info.ligands[0].heavy_atom_count, 5)
        # 选择串应基于 index
        self.assertIn("index", info.ligands[0].selection)
        self.assertGreater(cmd.count_atoms(info.ligands[0].selection), 0)

    def test_multiple_ligands_exclude_additives_water(self):
        """多配体保留；水/离子/常见添加剂排除。"""
        cmd.reinitialize()
        cmd.fab("ACDEFGHIKLMNPQRSTVWY", "complex")
        cmd.fragment("benzene", "lig1")
        cmd.alter("lig1", 'resn="BNZ"; resi="901"; chain="L"')
        cmd.translate([20.0, 0.0, 0.0], "lig1")
        cmd.fragment("indane", "lig2")
        cmd.alter("lig2", 'resn="IND"; resi="902"; chain="L"')
        cmd.translate([0.0, 20.0, 0.0], "lig2")
        _merge("complex", "lig1")
        _merge("complex", "lig2")
        _load_additive_soup("junk")
        _merge("complex", "junk")

        info = pub.detect_scene(cmd)
        self.assertTrue(info.has_protein)
        resns = {lig.resn.upper() for lig in info.ligands}
        self.assertIn("BNZ", resns)
        self.assertIn("IND", resns)
        for bad in ("HOH", "NA", "CL", "GOL", "EDO", "SO4"):
            self.assertNotIn(bad, resns)
        self.assertGreaterEqual(len(info.ligands), 2)

    def test_disabled_object_ignored(self):
        """未启用对象不参与检测。"""
        cmd.reinitialize()
        cmd.fab("AAAA", "on")
        cmd.fragment("benzene", "offlig")
        cmd.disable("offlig")
        info = pub.detect_scene(cmd)
        self.assertTrue(info.has_protein)
        self.assertEqual(len(info.ligands), 0)


class TestPublicationStyleRestore(unittest.TestCase):
    """样式应用与会话恢复。"""

    def test_overview_and_restore_repeatable(self):
        """预览后恢复应可重复，且不删除用户对象。"""
        cmd.reinitialize()
        cmd.fab("ACDEFGHIK", "prot")
        cmd.fragment("benzene", "bnz")
        cmd.alter("bnz", 'resn="BNZ"; resi="1"; chain="A"')
        n_atoms_before = cmd.count_atoms()
        names_before = set(cmd.get_names("public_objects"))

        fig = pub.PublicationFigure(cmd)
        info = pub.detect_scene(cmd)
        self.assertTrue(info.has_protein)
        self.assertGreaterEqual(len(info.ligands), 1)
        lig = info.ligands[0]

        fig.apply_style(pub.STYLE_OVERVIEW, ligand=lig, keep_view=True)
        self.assertTrue(fig.has_snapshot())
        self.assertIn("prot", cmd.get_names("public_objects"))
        self.assertIn("bnz", cmd.get_names("public_objects"))

        # 第一次恢复
        self.assertTrue(fig.restore_original())
        self.assertEqual(cmd.count_atoms(), n_atoms_before)
        self.assertEqual(set(cmd.get_names("public_objects")), names_before)

        # 再次应用并恢复，行为应可预期
        fig.apply_style(pub.STYLE_OVERVIEW, ligand=None, keep_view=True)
        self.assertTrue(fig.restore_original())
        self.assertEqual(cmd.count_atoms(), n_atoms_before)
        self.assertEqual(set(cmd.get_names("public_objects")), names_before)

        # 临时选择不应残留
        leftover = [n for n in cmd.get_names("all") if n.startswith("_pubfig_")]
        self.assertEqual(leftover, [])

    def test_pocket_requires_protein_and_ligand(self):
        """口袋样式在缺少蛋白或配体时应报错。"""
        cmd.reinitialize()
        cmd.fab("ACDEFGHIK", "prot")
        fig = pub.PublicationFigure(cmd)
        with self.assertRaises(RuntimeError):
            fig.apply_style(pub.STYLE_POCKET, ligand=None, keep_view=True)

        cmd.reinitialize()
        cmd.fragment("benzene", "lig")
        info = pub.detect_scene(cmd)
        fig = pub.PublicationFigure(cmd)
        with self.assertRaises(RuntimeError):
            fig.apply_style(pub.STYLE_POCKET, ligand=info.ligands[0], keep_view=True)

    def test_pocket_style_hides_solvent_keeps_atoms(self):
        """口袋样式隐藏溶剂但不删原子。"""
        cmd.reinitialize()
        cmd.fab("ACDEFGHIKLMN", "prot")
        cmd.fragment("benzene", "lig")
        cmd.alter("lig", 'resn="BNZ"; resi="900"; chain="A"')
        # 把配体挪近蛋白，便于口袋选择
        cmd.origin("prot")
        cmd.translate([5.0, 0.0, 0.0], "lig")
        _load_additive_soup("wat")
        n_before = cmd.count_atoms()

        info = pub.detect_scene(cmd)
        lig = next(x for x in info.ligands if x.resn.upper() == "BNZ")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(pub.STYLE_POCKET, ligand=lig, keep_view=True)
        self.assertEqual(cmd.count_atoms(), n_before)
        self.assertAlmostEqual(cmd.get_setting_float("cartoon_transparency", "prot"), 0.78, places=2)
        # 溶剂应无可见表示（隐藏），原子仍在
        self.assertGreater(cmd.count_atoms("wat"), 0)
        self.assertEqual(cmd.count_atoms("wat and (rep sticks or rep lines or rep nonbonded or rep spheres or rep cartoon)"), 0)


class TestPublicationValidation(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()

    def test_stale_candidate_does_not_mutate_scene(self):
        cmd.fragment("benzene", "lig")
        candidate = pub.detect_scene().ligands[0]
        cmd.alter("lig", 'resn="NEW"')
        fig = pub.PublicationFigure()
        with self.assertRaises(RuntimeError):
            fig.apply_style(ligand=candidate)
        self.assertFalse(fig.has_snapshot())

    def test_pocket_without_nearby_protein_rejected(self):
        cmd.fab("AAAA", "prot")
        cmd.fragment("benzene", "lig")
        cmd.translate([1000, 1000, 1000], "lig")
        fig = pub.PublicationFigure()
        with self.assertRaises(RuntimeError):
            fig.apply_style(pub.STYLE_POCKET, pub.detect_scene().ligands[0])
        self.assertFalse(fig.has_snapshot())

    def test_user_selection_and_disabled_representation_preserved(self):
        cmd.fragment("benzene", "lig")
        cmd.fragment("benzene", "hidden")
        cmd.disable("hidden")
        cmd.select("_pubfig_pocket", "lig")
        before = cmd.get_names("all")
        fig = pub.PublicationFigure()
        fig.apply_style(ligand=pub.detect_scene().ligands[0])
        self.assertEqual(cmd.get_names("all"), before)
        self.assertNotIn("hidden", cmd.get_names("objects", enabled_only=1))
        self.assertEqual(cmd.count_atoms("hidden and rep lines"), 12)


class TestPublicationThemes(unittest.TestCase):
    """三种智能风格：差异、切换复位、作用域与导出传参。"""

    def setUp(self):
        cmd.reinitialize()

    def _complex(self):
        """蛋白 + 近距配体。"""
        cmd.fab("ACDEFGHIKLMN", "prot")
        cmd.fragment("benzene", "lig")
        cmd.alter("lig", 'resn="BNZ"; resi="900"; chain="A"')
        cmd.origin("prot")
        cmd.translate([5.0, 0.0, 0.0], "lig")
        return pub.detect_scene(cmd).ligands[0]

    def test_exactly_three_themes_metadata(self):
        """后端恰好暴露三种主题元数据。"""
        themes = pub.list_themes()
        self.assertEqual(len(themes), 3)
        ids = [t["id"] for t in themes]
        self.assertEqual(ids, list(pub.THEME_ORDER))
        for t in themes:
            self.assertIn("label", t)
            self.assertIn("description", t)
            self.assertTrue(t["label"])
            self.assertTrue(t["description"])

    def test_invalid_theme_rejected_before_mutation(self):
        """非法主题在突变前拒绝，且不写快照。"""
        cmd.fab("AAAA", "prot")
        fig = pub.PublicationFigure(cmd)
        with self.assertRaises(ValueError):
            fig.apply_style(theme="neon_disco")
        self.assertFalse(fig.has_snapshot())

    def test_three_themes_visually_distinct(self):
        """三主题在表示/透明度/轮廓上可区分。"""
        lig = self._complex()
        fig = pub.PublicationFigure(cmd)

        fig.apply_style(theme=pub.THEME_JOURNAL, ligand=lig, keep_view=True)
        self.assertGreater(cmd.count_atoms("prot and polymer.protein and rep cartoon"), 0)
        self.assertEqual(cmd.count_atoms("prot and polymer.protein and rep surface"), 0)
        self.assertAlmostEqual(cmd.get_setting_float("cartoon_transparency", "prot"), 0.0, places=2)
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 0)
        j_loop = float(cmd.get_setting_float("cartoon_loop_radius", "prot"))

        fig.apply_style(theme=pub.THEME_ILLUSTRATION, ligand=lig, keep_view=True)
        self.assertGreater(cmd.count_atoms("prot and polymer.protein and rep cartoon"), 0)
        self.assertEqual(cmd.count_atoms("prot and polymer.protein and rep surface"), 0)
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 1)
        i_loop = float(cmd.get_setting_float("cartoon_loop_radius", "prot"))
        self.assertGreater(i_loop, j_loop)

        fig.apply_style(theme=pub.THEME_SURFACE, ligand=lig, keep_view=True)
        self.assertGreater(cmd.count_atoms("prot and polymer.protein and rep surface"), 0)
        self.assertGreater(cmd.count_atoms("prot and polymer.protein and rep cartoon"), 0)
        self.assertGreater(cmd.get_setting_float("transparency", "prot"), 0.2)
        self.assertGreater(cmd.count_atoms("lig and rep sticks"), 0)
        self.assertEqual(cmd.count_atoms("lig and rep surface"), 0)

    def test_theme_switch_resets_surface_and_radii(self):
        """从通透表面切回期刊简约时，表面与加粗卡通不得泄漏。"""
        lig = self._complex()
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(theme=pub.THEME_SURFACE, ligand=lig, keep_view=True)
        fig.apply_style(theme=pub.THEME_ILLUSTRATION, ligand=lig, keep_view=True)
        fig.apply_style(theme=pub.THEME_JOURNAL, ligand=lig, keep_view=True)
        self.assertEqual(cmd.count_atoms("prot and polymer.protein and rep surface"), 0)
        self.assertAlmostEqual(cmd.get_setting_float("transparency", "prot"), 0.0, places=2)
        self.assertAlmostEqual(cmd.get_setting_float("cartoon_transparency", "prot"), 0.0, places=2)
        self.assertAlmostEqual(
            float(cmd.get_setting_float("cartoon_loop_radius", "prot")),
            float(pub._OBJECT_THEME_DEFAULTS["cartoon_loop_radius"]),
            places=2,
        )
        self.assertEqual(int(cmd.get_setting_int("ray_trace_mode")), 0)

    def test_disabled_objects_out_of_scope(self):
        """禁用对象不被主题改写表示。"""
        cmd.fab("AAAA", "on")
        cmd.fragment("benzene", "offlig")
        cmd.disable("offlig")
        cmd.show("lines", "offlig")
        off_lines = cmd.count_atoms("offlig and rep lines")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(theme=pub.THEME_SURFACE, ligand=None, keep_view=True)
        self.assertEqual(cmd.count_atoms("offlig and rep lines"), off_lines)
        self.assertNotIn("offlig", cmd.get_names("objects", enabled_only=1))

    def test_restore_after_theme(self):
        """主题预览后 get_session 恢复可重复。"""
        lig = self._complex()
        names_before = set(cmd.get_names("public_objects"))
        n_before = cmd.count_atoms()
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(theme=pub.THEME_ILLUSTRATION, ligand=lig, keep_view=True)
        self.assertTrue(fig.restore_original())
        self.assertEqual(set(cmd.get_names("public_objects")), names_before)
        self.assertEqual(cmd.count_atoms(), n_before)

    def test_export_png_passes_theme(self):
        """export_png 须把所选主题传给 apply_style，不得静默回落 journal。"""
        from unittest.mock import patch

        lig = self._complex()
        fig = pub.PublicationFigure(cmd)
        seen = {}

        real_apply = fig.apply_style

        def _wrap(*a, **k):
            seen["theme"] = k.get("theme", a[5] if len(a) > 5 else None)
            return real_apply(*a, **k)

        with patch.object(fig, "apply_style", side_effect=_wrap):
            with patch.object(cmd, "png"):
                fig.export_png(
                    "ignored.png",
                    ligand=lig,
                    theme=pub.THEME_SURFACE,
                    keep_view=True,
                )
        self.assertEqual(seen.get("theme"), pub.THEME_SURFACE)

    def test_protein_only_and_standalone_ligand_themes(self):
        """无配体蛋白与独立配体均可应用主题。"""
        cmd.fab("ACDEFG", "prot")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(theme=pub.THEME_JOURNAL, ligand=None, keep_view=True)
        self.assertGreater(cmd.count_atoms("prot and rep cartoon"), 0)

        cmd.reinitialize()
        cmd.fragment("benzene", "lig")
        lig = pub.detect_scene(cmd).ligands[0]
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(theme=pub.THEME_ILLUSTRATION, ligand=lig, keep_view=True)
        self.assertGreater(cmd.count_atoms("lig and rep sticks"), 0)


class TestPublicationMultiLigand(unittest.TestCase):
    """多配体归一化、并集着色/口袋与校验。"""

    def setUp(self):
        cmd.reinitialize()

    def _two_in_one_object(self):
        """同一对象内两个配体候选。"""
        cmd.fab("ACDEFGHIKLMNPQRSTVWY", "complex")
        cmd.fragment("benzene", "lig1")
        cmd.alter("lig1", 'resn="BNZ"; resi="901"; chain="L"')
        cmd.translate([20.0, 0.0, 0.0], "lig1")
        cmd.fragment("indane", "lig2")
        cmd.alter("lig2", 'resn="IND"; resi="902"; chain="L"')
        cmd.translate([0.0, 20.0, 0.0], "lig2")
        _merge("complex", "lig1")
        _merge("complex", "lig2")
        info = pub.detect_scene(cmd)
        self.assertGreaterEqual(len(info.ligands), 2)
        return info.ligands

    def _two_across_objects(self):
        """跨对象两个配体。"""
        cmd.fab("ACDEFGHIKLMN", "prot")
        cmd.fragment("benzene", "lig_a")
        cmd.alter("lig_a", 'resn="BNZ"; resi="1"; chain="A"')
        cmd.translate([5.0, 0.0, 0.0], "lig_a")
        cmd.fragment("indane", "lig_b")
        cmd.alter("lig_b", 'resn="IND"; resi="2"; chain="B"')
        cmd.translate([0.0, 5.0, 0.0], "lig_b")
        info = pub.detect_scene(cmd)
        self.assertGreaterEqual(len(info.ligands), 2)
        return info.ligands

    def test_normalize_dedupe_and_legacy_single(self):
        """ligands 优先；去重；legacy 单配体仍可用。"""
        ligs = self._two_across_objects()
        a, b = ligs[0], ligs[1]
        out = pub.normalize_ligands(ligand=a, ligands=None)
        self.assertEqual(len(out), 1)
        self.assertEqual(pub.candidate_key(out[0]), pub.candidate_key(a))
        out2 = pub.normalize_ligands(ligands=[a, a, b])
        self.assertEqual(len(out2), 2)
        empty = pub.normalize_ligands(ligands=[])
        self.assertEqual(empty, [])
        none = pub.normalize_ligands()
        self.assertEqual(none, [])

    def test_two_ligands_same_object_shown(self):
        """同一对象两配体均显示 sticks。"""
        ligs = self._two_in_one_object()
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(ligands=ligs[:2], keep_view=True, theme=pub.THEME_JOURNAL)
        sel = pub.ligand_selection(ligs[:2])
        self.assertIsNotNone(sel)
        self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % sel), 0)
        for lig in ligs[:2]:
            self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % lig.selection), 0)

    def test_two_ligands_across_objects_shown_and_export(self):
        """跨对象两配体着色，且 export_png 传入完整列表。"""
        from unittest.mock import patch

        ligs = self._two_across_objects()
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(ligands=ligs[:2], keep_view=True, theme=pub.THEME_SURFACE)
        for lig in ligs[:2]:
            self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % lig.selection), 0)

        seen = {}
        real_apply = fig.apply_style

        def _wrap(*a, **k):
            seen["ligands"] = k.get("ligands")
            seen["theme"] = k.get("theme")
            return real_apply(*a, **k)

        with patch.object(fig, "apply_style", side_effect=_wrap):
            with patch.object(cmd, "png"):
                fig.export_png(
                    "ignored.png",
                    ligands=ligs[:2],
                    theme=pub.THEME_ILLUSTRATION,
                    keep_view=True,
                )
        self.assertEqual(seen.get("theme"), pub.THEME_ILLUSTRATION)
        self.assertEqual(len(seen.get("ligands") or []), 2)

    def test_pocket_union_includes_both_sites(self):
        """口袋 byres 并集覆盖两个近距配体位点。"""
        cmd.fab("ACDEFGHIKLMNPQRSTVWYACDEFG", "prot")
        cmd.fragment("benzene", "lig1")
        cmd.alter("lig1", 'resn="BNZ"; resi="901"; chain="L"')
        cmd.origin("prot")
        cmd.translate([4.0, 0.0, 0.0], "lig1")
        cmd.fragment("indane", "lig2")
        cmd.alter("lig2", 'resn="IND"; resi="902"; chain="L"')
        cmd.origin("prot")
        cmd.translate([-4.0, 0.0, 0.0], "lig2")
        ligs = pub.detect_scene(cmd).ligands
        self.assertGreaterEqual(len(ligs), 2)
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(
            style=pub.STYLE_POCKET,
            ligands=ligs[:2],
            keep_view=True,
            theme=pub.THEME_JOURNAL,
        )
        # 两配体均应可见；口袋残基 sticks 非空
        for lig in ligs[:2]:
            self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % lig.selection), 0)
        self.assertGreater(
            cmd.count_atoms("prot and polymer.protein and rep sticks"), 0
        )

    def test_stale_secondary_candidate_rejected_before_mutation(self):
        """任一过期候选整单拒绝，且不写快照。"""
        ligs = self._two_across_objects()
        a, b = ligs[0], ligs[1]
        # 使第二个候选过期
        cmd.alter(b.object, 'resn="ZZZ"')
        stale_b = b
        fig = pub.PublicationFigure(cmd)
        with self.assertRaises(RuntimeError):
            fig.apply_style(ligands=[a, stale_b], keep_view=True)
        self.assertFalse(fig.has_snapshot())

    def test_empty_ligands_protein_overview(self):
        """空 ligands 列表支持仅蛋白总览。"""
        cmd.fab("ACDEFG", "prot")
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(ligands=[], keep_view=True, theme=pub.THEME_JOURNAL)
        self.assertGreater(cmd.count_atoms("prot and rep cartoon"), 0)

    def test_standalone_multiple_ligands(self):
        """无蛋白时多个独立配体可出图。"""
        cmd.fragment("benzene", "lig_a")
        cmd.alter("lig_a", 'resn="BNZ"; resi="1"; chain="A"')
        cmd.fragment("indane", "lig_b")
        cmd.alter("lig_b", 'resn="IND"; resi="2"; chain="B"')
        ligs = pub.detect_scene(cmd).ligands
        self.assertGreaterEqual(len(ligs), 2)
        fig = pub.PublicationFigure(cmd)
        fig.apply_style(ligands=ligs, keep_view=True, theme=pub.THEME_ILLUSTRATION)
        sel = pub.ligand_selection(ligs)
        self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % sel), 0)


def _qt_available() -> bool:
    """检测 pymol.Qt 是否可用。"""
    try:
        from pymol.Qt import QtWidgets  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_qt_available(), "pymol.Qt 不可用")
class TestSmartThemeLigandConfirm(unittest.TestCase):
    """智能主题多配体确认/取消（Qt mock choose_ligands）。"""

    def setUp(self):
        cmd.reinitialize()
        from pymol.Qt import QtWidgets
        app = QtWidgets.QApplication.instance()
        if app is None:
            self.app = QtWidgets.QApplication([])
        else:
            self.app = app
        from pmg_qt import publication_gui as pg
        pg.clear_publication_dialog_ref()

    def tearDown(self):
        from pmg_qt import publication_gui as pg
        pg.clear_publication_dialog_ref()

    def _two_ligands(self):
        cmd.fab("ACDEFGHIKLMNPQRSTVWY", "complex")
        cmd.fragment("benzene", "lig1")
        cmd.alter("lig1", 'resn="BNZ"; resi="901"; chain="L"')
        cmd.translate([20.0, 0.0, 0.0], "lig1")
        cmd.fragment("indane", "lig2")
        cmd.alter("lig2", 'resn="IND"; resi="902"; chain="L"')
        cmd.translate([0.0, 20.0, 0.0], "lig2")
        _merge("complex", "lig1")
        _merge("complex", "lig2")

    def test_single_candidate_auto_select(self):
        """单候选自动选用，无需多选对话框。"""
        from unittest.mock import patch
        from pymol.Qt import QtWidgets
        from pmg_qt import publication_gui as pg

        cmd.fab("AAAA", "prot")
        cmd.fragment("benzene", "lig")
        cmd.alter("lig", 'resn="BNZ"; resi="1"; chain="A"')
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        with patch.object(pg, "choose_ligands") as chooser:
            ok = pg.apply_smart_theme(parent, pub.THEME_JOURNAL)
            chooser.assert_not_called()
        self.assertTrue(ok)
        self.assertGreater(cmd.count_atoms("prot and rep cartoon"), 0)

    def test_multi_candidate_cancel_leaves_scene(self):
        """多选对话框取消时场景/快照不变。"""
        from unittest.mock import patch
        from pymol.Qt import QtWidgets
        from pmg_qt import publication_gui as pg

        self._two_ligands()
        before_cartoon = cmd.count_atoms("complex and rep cartoon")
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        with patch.object(pg, "choose_ligands", return_value=([], False)):
            ok = pg.apply_smart_theme(parent, pub.THEME_ILLUSTRATION)
        self.assertFalse(ok)
        self.assertEqual(cmd.count_atoms("complex and rep cartoon"), before_cartoon)
        dlg = pg.get_or_create_publication_dialog(parent, show=False)
        self.assertFalse(dlg.fig.has_snapshot())

    def test_multi_candidate_confirm_subset_reused(self):
        """确认子集后后续主题复用，不再弹窗；两配体均着色。"""
        from unittest.mock import patch
        from pymol.Qt import QtWidgets
        from pmg_qt import publication_gui as pg

        self._two_ligands()
        info = pub.detect_scene(cmd)
        self.assertGreaterEqual(len(info.ligands), 2)
        subset = info.ligands[:2]
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        with patch.object(pg, "choose_ligands", return_value=(subset, True)) as chooser:
            self.assertTrue(pg.apply_smart_theme(parent, pub.THEME_JOURNAL))
            self.assertEqual(chooser.call_count, 1)
            self.assertTrue(pg.apply_smart_theme(parent, pub.THEME_SURFACE))
            self.assertEqual(chooser.call_count, 1)
        for lig in subset:
            self.assertGreater(cmd.count_atoms("(%s) and rep sticks" % lig.selection), 0)

    def test_checkbox_selection_preserved_across_refresh_and_export_path(self):
        """对话框勾选在刷新与主题/导出路径间保持。"""
        from unittest.mock import patch
        from pymol.Qt import QtWidgets, QtCore
        from pmg_qt import publication_gui as pg

        self._two_ligands()
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        dlg = pg.get_or_create_publication_dialog(parent, show=False)
        self.assertGreaterEqual(dlg.ligand_list.count(), 2)
        # 显式勾选前两个
        dlg.ligand_list.blockSignals(True)
        for i in range(dlg.ligand_list.count()):
            state = (
                QtCore.Qt.CheckState.Checked
                if i < 2
                else QtCore.Qt.CheckState.Unchecked
            )
            dlg.ligand_list.item(i).setCheckState(state)
        dlg.ligand_list.blockSignals(False)
        dlg._mark_ligands_explicit()
        keys_before = [pub.candidate_key(x) for x in dlg._current_ligands()]
        self.assertEqual(len(keys_before), 2)
        theme_before = dlg._current_theme()

        dlg.refresh_candidates()
        keys_after = [pub.candidate_key(x) for x in dlg._current_ligands()]
        self.assertEqual(keys_after, keys_before)

        # 智能主题应复用勾选，不弹 chooser；主题与配体保留
        with patch.object(pg, "choose_ligands") as chooser:
            self.assertTrue(pg.apply_smart_theme(parent, pub.THEME_SURFACE))
            chooser.assert_not_called()
        self.assertEqual(
            [pub.candidate_key(x) for x in dlg._current_ligands()], keys_before
        )
        self.assertEqual(dlg._current_theme(), pub.THEME_SURFACE)

        # 打开导出对话框路径（再次 get）不得把勾选重置为“仅第一个”
        dlg2 = pg.open_publication_dialog(parent)
        self.assertIs(dlg2, dlg)
        self.assertEqual(
            [pub.candidate_key(x) for x in dlg2._current_ligands()], keys_before
        )
        self.assertEqual(dlg2._current_theme(), pub.THEME_SURFACE)
        self.assertNotEqual(theme_before, pub.THEME_SURFACE)

        # 预览/导出内部必须传完整列表
        seen = {}
        real_apply = dlg.fig.apply_style

        def _wrap(*a, **k):
            seen["ligands"] = k.get("ligands")
            return real_apply(*a, **k)

        with patch.object(dlg.fig, "apply_style", side_effect=_wrap):
            dlg.on_preview()
        self.assertEqual(len(seen.get("ligands") or []), 2)
    def test_empty_scene_message(self):
        """空场景给出中文提示且不改场景。"""
        from unittest.mock import patch
        from pymol.Qt import QtWidgets
        from pmg_qt import publication_gui as pg

        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        with patch.object(QtWidgets.QMessageBox, "information") as info:
            ok = pg.apply_smart_theme(parent, pub.THEME_JOURNAL)
        self.assertFalse(ok)
        info.assert_called()
        args = info.call_args[0]
        self.assertTrue(any(isinstance(a, str) and a for a in args))


@unittest.skipUnless(_qt_available(), "Qt unavailable")
class TestExplicitProteinOnlyChoice(unittest.TestCase):
    def test_clear_single_ligand_is_respected(self):
        from pymol.Qt import QtWidgets
        from pmg_qt import publication_gui as pg
        cmd.reinitialize()
        cmd.fab("AAAA", "prot")
        cmd.fragment("benzene", "lig")
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        parent = QtWidgets.QMainWindow()
        parent.cmd = cmd
        dialog = pg.get_or_create_publication_dialog(parent)
        dialog._on_clear_ligands()
        selected, accepted = dialog.resolve_ligands_for_smart_theme(parent)
        self.assertTrue(accepted)
        self.assertEqual(selected, [])
        pg.clear_publication_dialog_ref()
        parent.deleteLater()
