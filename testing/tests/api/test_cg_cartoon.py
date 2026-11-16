"""Cartoon paths follow the selected CG topology and update across states."""

import unittest
from unittest.mock import patch

from chempy import Atom, Bond
from chempy.models import Indexed
from pymol import cgo, cmd, visual_representations as vr
from pymol.cg_cartoon import BASE_RADIUS, TUBE_RADIUS, geometry


def _rna_model():
    model = Indexed()
    model.connect_mode = 1
    for chain, resi, x, y in [('A', '1', 0, 0), ('A', '2', 5, 1),
                             ('A', '3', 31, 1), ('B', '1', 5, -1)]:
        for name, dx, dy in [('BB1', 0, 0), ('BB2', 2, 0), ('BB3', 4, 0),
                             ('SC1', 4, 2), ('SC2', 4, 3)]:
            atom = Atom()
            atom.index = len(model.atom) + 1
            atom.name, atom.chain, atom.resi = name, chain, resi
            atom.resn, atom.symbol = 'U', 'C'
            atom.coord = [x + dx, y + dy, 0]
            model.add_atom(atom)
    # Missing A1→A2, overlong A2→A3, and an explicit cross-chain contact.
    pairs = [(5*i+j, 5*i+j+1) for i in range(4) for j in (0, 1)]
    pairs += [(2, 3), (3, 4), (7, 8), (7, 10), (2, 15)]
    for i, j in pairs:
        bond = Bond()
        bond.index = [i, j]
        model.add_bond(bond)
    return model


class CGCartoonTests(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()

    def test_breaks_and_base_directions_follow_existing_bonds(self):
        model = _rna_model()
        indices = {a.index for a in model.atom}
        data, connected = geometry(model, indices, dict.fromkeys(indices, (0.2, 0.4, 0.6)))
        self.assertTrue(connected)
        self.assertEqual(len(data) % 14, 0)
        tubes = []
        bases = []
        for i in range(0, len(data), 14):
            self.assertEqual(data[i], cgo.SAUSAGE)
            a, b, radius = data[i+1:i+4], data[i+4:i+7], data[i+7]
            if radius == TUBE_RADIUS:
                tubes.append((a, b))
            elif radius == BASE_RADIUS:
                bases.append((a, b))
        self.assertGreater(len(tubes), 8)  # more than the raw backbone edges
        runs = [(0, 4, 0), (5, 9, 1), (31, 35, 1), (5, 9, -1)]
        for a, b in tubes:
            self.assertTrue(any(lo <= a[0] <= hi and lo <= b[0] <= hi
                                and abs(a[1]-y) < 1e-6 and abs(b[1]-y) < 1e-6
                                for lo, hi, y in runs), (a, b))
        self.assertEqual(bases, [([4, 0, 0], [4, 2.5, 0])])

    def test_layer_builds_and_replaces_colored_geometry_for_each_state(self):
        cmd.load_model(_rna_model(), 'source', finish=0, zoom=0)
        cmd.load_coordset(cmd.get_coords('source') + [0, 0, 10], 'source', state=2)
        ctrl = vr.RepresentationController(cmd)
        layer = vr.RepresentationLayer('rna', selection='nucleic and chain B',
                                       style='cartoon', color_mode='solid', color_rgb=(1, 0, 0))
        with patch.object(cmd, 'load_cgo', wraps=cmd.load_cgo) as calls:
            report = ctrl.apply('source', [layer])
        self.assertTrue(report.ok, report.message)
        self.assertEqual(cmd.count_atoms(layer.owned_object), 5)
        self.assertEqual(cmd.count_atoms(layer.owned_object + ' and visible'), 0)
        self.assertEqual(calls.call_count, 2)
        loads = [call.args for call in calls.call_args_list]
        self.assertEqual(cmd.count_states(loads[0][1]), 2)
        for z, (data, _) in zip((0, 10), loads):
            self.assertEqual(data[3], z)
            self.assertEqual(data[8:11], [1, 0, 0])

        old_names = set(ctrl.owned)
        layer.color_rgb = (0, 1, 0)
        with patch.object(cmd, 'load_cgo', wraps=cmd.load_cgo) as calls:
            report = ctrl.apply('source', [layer])
        self.assertTrue(report.ok, report.message)
        self.assertEqual(calls.call_args.args[0][8:11], [0, 1, 0])
        self.assertFalse(old_names.intersection(cmd.get_names('objects')))
        ctrl.restore()
        self.assertEqual(cmd.get_names('objects', enabled_only=1), ['source'])
