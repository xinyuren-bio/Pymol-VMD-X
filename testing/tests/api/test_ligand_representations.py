"""Connectivity checks and scoped sticks/ball-and-stick rendering regressions."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import unittest

from chempy import Atom, Bond
from chempy.models import Indexed
from pymol import cmd, visual_representations as vr, cg_templates


def load_ligand(name, bonds=()):
    model = Indexed()
    model.connect_mode = 1
    for i in range(3):
        atom = Atom()
        atom.name = 'C%d' % i
        atom.resn = 'LIG'
        atom.resi = '1'
        atom.symbol = 'C'
        atom.coord = [i * 1.4, 0, 0]
        atom.hetatm = 1
        model.add_atom(atom)
    for i, j, order in bonds:
        bond = Bond()
        bond.index = [i, j]
        bond.order = order
        model.add_bond(bond)
    cmd.load_model(model, name, zoom=0)


def snapshot(name):
    model = cmd.get_model(name)
    colors = []
    cmd.iterate(name, 'colors.append(color)', space={'colors': colors})
    return ([tuple(a.coord) for a in model.atom],
            [(tuple(b.index), b.order) for b in model.bond], colors,
            cmd.count_atoms(name + ' and rep sticks'),
            cmd.count_atoms(name + ' and rep spheres'))


class LigandRepresentations(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()

    def test_missing_bonds_preserves_previous_display_and_other_object(self):
        load_ligand('source')
        load_ligand('other', [(0, 1, 2)])
        cmd.show('spheres', 'other')
        before = snapshot('source'), snapshot('other')
        ctrl = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer('ligand', selection='resname LIG', style='spheres')
        self.assertTrue(ctrl.apply('source', [layer]).ok)
        owned = list(ctrl.owned)
        visible = set(cmd.get_names('objects', enabled_only=1))
        for style in ('sticks', 'ballstick', '棍状', '球棍'):
            layer.style = style
            result = ctrl.apply('source', [layer])
            self.assertFalse(result.ok, result.message)
            self.assertIn('0 条内部键', result.message)
            self.assertIn('未自动添加化学键', result.message)
            self.assertEqual(ctrl.owned, owned)
            self.assertEqual(set(cmd.get_names('objects', enabled_only=1)), visible)
        self.assertEqual((snapshot('source'), snapshot('other')), before)

    def test_sticks_and_ballstick_keep_source_bonds_and_orders(self):
        load_ligand('source', [(0, 1, 2), (1, 2, 1)])
        cmd.set('stick_radius', 0, 'source')
        cmd.set_bond('stick_radius', 0, 'source', 'source')
        cmd.set('stick_transparency', 1, 'source')
        cmd.set_bond('stick_transparency', 1, 'source', 'source')
        cmd.set('sphere_transparency', 1, 'source')
        before = snapshot('source')
        ctrl = vr.RepresentationController(cmd)
        preview = ctrl.preview('source', 'resname LIG')
        self.assertEqual((preview.count, preview.bond_count), (3, 2))
        layer = vr.RepresentationLayer('ligand', selection='resname LIG')
        for style, spheres in [('sticks', 0), ('ballstick', 3), ('sticks', 0)]:
            layer.style = style
            result = ctrl.apply('source', [layer])
            self.assertTrue(result.ok, result.message)
            output = ctrl.owned[0]
            self.assertEqual(cmd.count_atoms(output + ' and rep sticks'), 3)
            self.assertEqual(cmd.count_atoms(output + ' and rep spheres'), spheres)
            self.assertEqual([b.order for b in cmd.get_model(output).bond], [2, 1])
            self.assertAlmostEqual(cmd.get_setting_float('stick_radius', output), .18, places=5)
            self.assertEqual(cmd.get_setting_float('stick_transparency', output), 0)
            self.assertEqual(result.layer_bond_counts[output], 2)
            overrides = cmd.get_bond('stick_radius', output, output)
            self.assertTrue(all(value is None for _, values in overrides
                                for _, _, value in values))
        self.assertEqual(snapshot('source'), before)
        self.assertEqual(cmd.get_setting_float('stick_radius', 'source'), 0)
        self.assertEqual(cmd.get_setting_float('stick_transparency', 'source'), 1)

    def test_boundary_bonds_do_not_expand_rule(self):
        load_ligand('source', [(0, 1, 1), (1, 2, 1)])
        ctrl = vr.RepresentationController(cmd)
        preview = ctrl.preview('source', 'index 1')
        self.assertEqual((preview.count, preview.bond_count, preview.external_bond_count), (1, 0, 2))
        layer = vr.RepresentationLayer('ligand', selection='index 1', style='sticks')
        result = ctrl.apply('source', [layer])
        self.assertFalse(result.ok)
        self.assertIn('未扩大选择范围', result.message)
        self.assertFalse(ctrl.owned)
        layer.selection = 'index 0 to 1'
        self.assertTrue(ctrl.apply('source', [layer]).ok)
        self.assertEqual(cmd.count_atoms(ctrl.owned[0]), 2)
        self.assertEqual(len(cmd.get_model(ctrl.owned[0]).bond), 1)

    def test_mixed_cg_template_retains_ligand_bond_order(self):
        load_ligand('source', [(0, 1, 2), (1, 2, 1)])
        model = cmd.get_model('source')
        template = next(t for t in cg_templates.load_registry().templates if t.id == 'DPPC')
        for i, spec in enumerate(template.atoms):
            atom = Atom()
            atom.name = spec['name']; atom.resn = 'DPPC'; atom.resi = '2'
            atom.symbol = 'C'; atom.coord = [i * .3, 20, 0]
            model.add_atom(atom)
        model.connect_mode = 1
        cmd.delete('source')
        cmd.load_model(model, 'source', finish=0, zoom=0)
        before = snapshot('source')
        ctrl = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer('ligand', selection='resname LIG', style='sticks')
        result = ctrl.apply('source', [layer])
        self.assertTrue(result.ok, result.message)
        self.assertEqual(cmd.count_atoms(ctrl.owned[0]), 3)
        self.assertEqual([b.order for b in cmd.get_model(ctrl.owned[0]).bond], [2, 1])
        self.assertEqual(snapshot('source'), before)

    def test_unknown_cg_does_not_guess_connections(self):
        load_ligand('source')
        cmd.alter('source', "resn='UNKNOWN'; name='BB'")
        ctrl = vr.RepresentationController(cmd)
        result = ctrl.apply('source', [vr.RepresentationLayer('cg', style='sticks')])
        self.assertFalse(result.ok)
        self.assertEqual(len(cmd.get_model('source').bond), 0)


SAMPLE = os.environ.get('PYMOLX_LIGAND_SAMPLE', '')


@unittest.skipUnless(os.path.isfile(SAMPLE), '需要 PYMOLX_LIGAND_SAMPLE 指向 dpol-24.pdb')
class ActualDpolLigand(unittest.TestCase):
    def test_actual_object_selection_switch_and_scope(self):
        cmd.reinitialize()
        cmd.load(SAMPLE, 'dpol-24', zoom=0)
        load_ligand('other', [(0, 1, 1)])
        before = snapshot("'dpol-24'"), snapshot('other')
        ctrl = vr.RepresentationController(cmd)
        preview = ctrl.preview('dpol-24', 'resname LIG')
        self.assertTrue(preview.ok, preview.message)
        self.assertEqual((preview.count, preview.bond_count), (22, 0))
        layer = vr.RepresentationLayer('ligand', selection='resname LIG', style='spheres')
        self.assertTrue(ctrl.apply('dpol-24', [layer]).ok)
        self.assertEqual(cmd.count_atoms(ctrl.owned[0]), 22)
        previous = list(ctrl.owned)
        for style in ('sticks', 'ballstick'):
            layer.style = style
            result = ctrl.apply('dpol-24', [layer])
            self.assertFalse(result.ok, result.message)
            self.assertIn('22 个原子/珠', result.message)
            self.assertEqual(ctrl.owned, previous)
        # Existing protein connectivity in the same real object renders both styles.
        layer.selection = 'protein and resid 1'
        for style in ('sticks', 'ballstick'):
            layer.style = style
            result = ctrl.apply('dpol-24', [layer])
            self.assertTrue(result.ok, result.message)
            self.assertGreater(len(cmd.get_model(ctrl.owned[0]).bond), 0)
        self.assertEqual((snapshot("'dpol-24'"), snapshot('other')), before)
        ctrl.restore()


class LigandDialog(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pymol.Qt import QtWidgets
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def test_manual_and_auto_apply_show_missing_bonds_without_scene_change(self):
        from pymol.Qt import QtWidgets
        from pmg_qt.visual_representations_gui import VisualRepresentationsDialog
        cmd.reinitialize()
        load_ligand('source')
        parent = QtWidgets.QWidget(); parent.cmd = cmd
        dialog = VisualRepresentationsDialog(parent, cmd)
        dialog.select_target('source')
        dialog.ed_sele.setText('resname LIG')
        dialog._preview()
        self.assertIn('源内部键：0', dialog.lbl_preview.text())
        dialog.cmb_style.setCurrentIndex(dialog.cmb_style.findData('spheres'))
        dialog.btn_apply.click()
        previous = list(dialog.ctrl.owned)
        self.assertTrue(previous)
        dialog.cmb_style.setCurrentIndex(dialog.cmb_style.findData('sticks'))
        self.assertIn('只显示', dialog.lbl_style_help.text())
        dialog.btn_apply.click()
        self.assertIn('0 条内部键', dialog.status.text())
        self.assertEqual(dialog.ctrl.owned, previous)
        dialog.chk_auto.setChecked(True)
        dialog.cmb_style.setCurrentIndex(dialog.cmb_style.findData('ballstick'))
        self.assertIn('小原子球', dialog.lbl_style_help.text())
        self.assertIn('0 条内部键', dialog.status.text())
        self.assertEqual(dialog.ctrl.owned, previous)
        dialog.close(); parent.close()

    @unittest.skipUnless(os.path.isfile(SAMPLE), '需要真实 dpol-24 样本')
    def test_actual_dpol_dialog_select_and_switch_styles(self):
        from pymol.Qt import QtCore, QtWidgets
        from pmg_qt.visual_representations_gui import VisualRepresentationsDialog
        cmd.reinitialize()
        cmd.load(SAMPLE, 'dpol-24', zoom=0)
        before = snapshot("'dpol-24'")
        parent = QtWidgets.QWidget(); parent.cmd = cmd
        dialog = VisualRepresentationsDialog(parent, cmd)
        dialog.select_target('dpol-24')
        item = next(dialog.list_resnames.item(i) for i in range(dialog.list_resnames.count())
                    if dialog.list_resnames.item(i).data(QtCore.Qt.ItemDataRole.UserRole) == 'LIG')
        dialog._choose_resname(item)
        self.assertEqual(dialog.ed_sele.text(), 'resname LIG')
        self.assertIn('22 珠/原子；源内部键：0', dialog.lbl_preview.text())
        dialog.cmb_style.setCurrentIndex(dialog.cmb_style.findData('spheres'))
        dialog.btn_apply.click()
        previous = list(dialog.ctrl.owned)
        self.assertEqual(cmd.count_atoms(previous[0]), 22)
        for style in ('sticks', 'ballstick'):
            dialog.cmb_style.setCurrentIndex(dialog.cmb_style.findData(style))
            dialog.btn_apply.click()
            self.assertIn('22 个原子/珠，但有 0 条内部键', dialog.status.text())
            self.assertEqual(dialog.ctrl.owned, previous)
            self.assertIn('22 珠/原子；源内部键：0', dialog.lbl_preview.text())
        self.assertEqual(snapshot("'dpol-24'"), before)
        dialog.ctrl.restore()
        dialog.close(); parent.close()


if __name__ == '__main__':
    unittest.main()
