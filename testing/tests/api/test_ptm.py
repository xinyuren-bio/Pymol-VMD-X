import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from pymol import cmd, ptm


def atoms(model):
    return {(*tuple_for(a), a.name): a for a in model.atom}


def tuple_for(atom):
    return atom.segi, atom.chain, atom.resi, atom.resn


def bonds(model):
    return {frozenset((*tuple_for(model.atom[i]), model.atom[i].name)
                      for i in bond.index): bond.order for bond in model.bond}


class PTMModels(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()
        cmd.fab('ASTYKK', 'source', chain='A')
        self.source = ptm.object_model(cmd, 'source')
        self.assertEqual(len(ptm.list_sites(self.source)), 5)

    def build_all(self):
        requests = zip(ptm.list_sites(self.source), ('SEP', 'TPO', 'PTR', 'XRW', 'D8R'))
        return ptm.build_model(self.source, requests)

    def test_existing_coordinates_stereochemistry_and_pdb_round_trip(self):
        result, sites, _ = self.build_all()
        self.assertEqual([site.resn for site in sites], ['SEP', 'TPO', 'PTR', 'XRW', 'D8R'])
        original_coords = {(a.chain, a.resi, a.name): a.coord for a in self.source.atom
                           if a.symbol != 'H'}
        changed_coords = {(a.chain, a.resi, a.name): a.coord for a in result.atom}
        self.assertTrue(all(changed_coords[key] == pos for key, pos in original_coords.items()))
        self.assertEqual(atoms(cmd.get_model('source')).keys(), atoms(self.source).keys())
        ptm.load_result(cmd, 'modified', result, sites, preview=False)

        written = cmd.get_model('modified')
        self.assertEqual(bonds(result), bonds(written))
        original_setting = cmd.get_setting_boolean('pdb_conect_nodup')
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'modified.pdb'
            ptm.export_pdb(cmd, 'modified', path)
            self.assertIn('MD TOPOLOGY AND FORCE FIELD PARAMETERS: NONE', path.read_text())
            cmd.load(str(path), 'round_trip', zoom=0)
        self.assertEqual(cmd.get_setting_boolean('pdb_conect_nodup'), original_setting)
        read = cmd.get_model('round_trip')
        self.assertEqual(atoms(written).keys(), atoms(read).keys())
        self.assertEqual(bonds(written), bonds(read))
        for site in sites:
            if site.resn not in ('XRW', 'D8R'):
                continue
            names = {a.name: i for i, a in enumerate(read.atom) if site == ptm.Site.from_atom(a)}
            stereo = ptm._lactyl_configuration(np.asarray(read.get_coord_list()), names, site.resn)
            self.assertEqual(stereo, '2S' if site.resn == 'XRW' else '2R')

    def test_repair_only_selected_imported_residue(self):
        result, sites, _ = self.build_all()
        first, second = sites[:2]
        old = bonds(result)
        missing = {frozenset(((site.segment, site.chain, site.resi, site.resn, 'P'),
                              (site.segment, site.chain, site.resi, site.resn, atom)))
                   for site, atom in ((first, 'OG'), (second, 'OG1'))}
        result.bond = [b for b in result.bond if
                       frozenset((*tuple_for(result.atom[i]), result.atom[i].name)
                                 for i in b.index) not in missing]
        repaired, _, _ = ptm.build_model(result, [(first, first.resn)])
        first_edge = frozenset(((first.segment, first.chain, first.resi, first.resn, 'P'),
                                (first.segment, first.chain, first.resi, first.resn, 'OG')))
        self.assertEqual(bonds(repaired), {edge: order for edge, order in old.items()
                                         if edge not in missing or edge == first_edge})
        self.assertEqual(repaired.get_coord_list(), result.get_coord_list())

    def test_missing_anchor_leaves_source_intact(self):
        cmd.remove('source & resn SER & name OG')
        incomplete = cmd.get_model('source')
        site = next(s for s in ptm.list_sites(incomplete) if s.resn == 'SER')
        with self.assertRaisesRegex(ValueError, 'OG'):
            ptm.build_model(incomplete, [(site, 'SEP')])
        self.assertEqual(atoms(cmd.get_model('source')).keys(), atoms(incomplete).keys())


class PTMDialogWorkflow(unittest.TestCase):
    def test_preview_apply_and_close(self):
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        try:
            from pymol.Qt import QtWidgets
        except ImportError:
            self.skipTest('Qt unavailable')
        from pmg_qt.ptm_gui import PTMDialog

        cmd.reinitialize()
        cmd.fab('ASTYKK', 'source', chain='A')
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        parent = QtWidgets.QWidget()
        parent.cmd = cmd
        dialog = PTMDialog(parent)
        dialog.show()
        app.processEvents()
        try:
            dialog._add_row_selection(0)
            last = dialog.table.rowCount() - 1
            dialog.table.cellWidget(last, 2).setCurrentIndex(1)
            dialog._add_row_selection(last)
            dialog._build(preview=True)
            preview = dialog.preview
            self.assertIn(preview, cmd.get_names('objects'))
            self.assertEqual(cmd.count_atoms('%source & resn SEP+D8R'), 0)
            dialog._build(preview=False)
            self.assertNotIn(preview, cmd.get_names('objects'))
            self.assertEqual(cmd.count_atoms('%' + dialog.result + ' & name P'), 1)
            self.assertGreater(cmd.count_atoms('%' + dialog.result + ' & resn D8R & rep sticks'), 0)
            dialog._build(preview=True)
            preview = dialog.preview
            dialog.reject()
            self.assertNotIn(preview, cmd.get_names('objects'))
            self.assertIn(dialog.result, cmd.get_names('objects'))
        finally:
            dialog.close()
            parent.close()
