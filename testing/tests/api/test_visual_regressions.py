"""Regression checks for selection boundaries, transactional copies and row editing."""
import copy
import os
import unittest
from unittest.mock import patch

from chempy import Atom
from chempy.models import Indexed
from pymol import cmd, visual_representations as vr, cg_templates


def _make_source(name):
    model=Indexed()
    for i in range(12):
        atom=Atom();atom.name='C'+str(i);atom.resn='LIG';atom.resi='1'
        atom.symbol='C';atom.coord=[i,0,0];model.add_atom(atom)
    cmd.load_model(model,name,zoom=0)


class VisualRegressions(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()

    def test_all_scope_reapply_and_restore(self):
        _make_source('source_a')
        _make_source('source_b')
        _make_source('hidden_source')
        cmd.disable('hidden_source')
        ctrl = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer(vr.new_layer_id(), style='spheres',
                                       color_mode='solid', color_rgb=(.2,.4,.6))
        before = cmd.count_atoms('source_a or source_b')
        for _ in range(2):
            result = ctrl.apply('all', [layer])
            self.assertTrue(result.ok, result.message)
            self.assertEqual(sum(result.layer_counts.values()), before)
            self.assertEqual(len(ctrl.owned), 2)
        ctrl.restore()
        self.assertEqual(set(cmd.get_names('objects', enabled_only=1)), {'source_a','source_b'})

    def test_failed_style_does_not_leak_new_object(self):
        _make_source('source')
        ctrl = vr.RepresentationController(cmd)
        result = ctrl.apply('source')
        self.assertTrue(result.ok, result.message)
        before = set(cmd.get_names('all'))
        with patch.object(vr, '_apply_color', side_effect=RuntimeError('test failure')):
            result = ctrl.apply('source')
        self.assertFalse(result.ok)
        self.assertEqual(set(cmd.get_names('all')), before)

    def test_cg_coordinates_states_colors_and_bonds(self):
        template = next(t for t in cg_templates.load_registry().templates if t.id=='DPPC')
        model = Indexed()
        for i, spec in enumerate(template.atoms):
            atom = Atom()
            atom.name=spec['name']; atom.resn='DPPC'; atom.resi='1'; atom.chain='A'
            atom.symbol='C'; atom.coord=[i*.3,0,0]
            model.add_atom(atom)
        cmd.load_model(model, 'cg', finish=0)
        cmd.load_coordset(cmd.get_coords('cg')+[0,2,0], 'cg', state=2)
        cmd.color('0x224466','cg')
        original = cmd.get_coords('cg',state=2).copy()
        ctrl = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer(vr.new_layer_id(),style='sticks')
        result = ctrl.apply('cg',[layer])
        self.assertTrue(result.ok,result.message)
        output = ctrl.owned[0]
        self.assertEqual(cmd.count_states(output),2)
        self.assertTrue((original==cmd.get_coords(output,state=2)).all())
        self.assertEqual(len(cmd.get_model(output).bond),11)
        colors=[]; cmd.iterate(output,'colors.append(color)',space={'colors':colors})
        self.assertEqual(len(set(colors)),1)
        self.assertAlmostEqual(cmd.get_color_tuple(colors[0])[0],34/255,places=5)
        self.assertEqual(len(cmd.get_model('cg').bond),0)

    def test_selection_target_and_sparse_indices(self):
        _make_source('source')
        cmd.remove('source and index 2')
        cmd.select('subset','source and index 5+7+9')
        ctrl=vr.RepresentationController(cmd)
        self.assertEqual(ctrl.preview('subset','all').count,3)
        self.assertEqual(ctrl.preview('subset','same residue as all').count,3)
        self.assertEqual(ctrl.preview('subset','within 0 of all').count,3)
        result=ctrl.apply('subset')
        self.assertTrue(result.ok,result.message)
        self.assertEqual(cmd.count_atoms(ctrl.owned[0]),3)

    def test_integer_fields_and_keyword_like_names(self):
        for text in ('index 1.5','resid 1.5 to 2.5','serial -1','index -1 to 2','name CA; delete all'):
            with self.assertRaises(vr.SelectionError,msg=text): vr.parse_selection(text)
        cmd.pseudoatom('scope',name='123HG',resn='NAME1')
        for text in ('name 123HG','resname NAME1','name "123HG"'):
            self.assertEqual(vr.preview_selection('scope',text,cmd).count,1,text)


if __name__=='__main__':
    unittest.main()
