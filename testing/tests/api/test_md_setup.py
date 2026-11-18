from io import BytesIO
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest

from pymol import cmd, md_charmm_gui, md_parameters as parameters, md_setup, ptm


class CharmmGuiArchiveTests(unittest.TestCase):
    def write_archive(self, path, files):
        with tarfile.open(path, 'w:gz') as archive:
            for name, data in files.items():
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                archive.addfile(entry, BytesIO(data))

    def test_preserves_topology_and_numeric_stage_order(self):
        with TemporaryDirectory() as temp:
            work = Path(temp)
            job = work / 'job.tgz'
            inputs = {
                'topol.top': b'#include "toppar/forcefield.itp"\n',
                'toppar/forcefield.itp': b'[ defaults ]\n; columns\n1 2 yes 1.0 1.0\n',
                'step3_input.gro': b'coordinates\n',
                'index.ndx': b'[ System ]\n1\n',
                'step4.0_minimization.mdp': b'integrator = steep\n',
                'step4.10_equilibration.mdp': b'nsteps = 10\n',
                'step4.2_equilibration.mdp': b'nsteps = 2\r\n',
                'step5_production.mdp': b'nsteps = 100\n',
            }
            self.write_archive(job, {'job/gromacs/' + k: v for k, v in inputs.items()})
            md_charmm_gui._read_inputs(job, work)
            files = md_charmm_gui._topology_files(work)
            self.assertEqual(files, ['topol.top', 'toppar/forcefield.itp'])
            stages = md_charmm_gui._write_stages(work)
            self.assertEqual(list(stages), ['em', 'eq', 'eq2', 'md'])
            self.assertEqual(stages['eq'], 'step4.2_equilibration.mdp')
            for stage, source in stages.items():
                self.assertEqual((work / 'mdp' / f'{stage}.mdp').read_bytes(), inputs[source])
            (work / 'toppar' / 'forcefield.itp').write_text('[ defaults ]\n1 2 yes 0.5 0.8333\n')
            with self.assertRaisesRegex(ValueError, '不符'):
                md_charmm_gui._topology_files(work)

    def test_rejects_unsafe_members_and_missing_dependencies(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / 'job.tgz'
            self.write_archive(archive, {'job/gromacs/../../unexpected.itp': b'bad'})
            with self.assertRaisesRegex(ValueError, '无效路径'):
                md_charmm_gui._read_inputs(archive, root)
            (root / 'topol.top').write_text('#include "toppar/missing.itp"\n')
            with self.assertRaisesRegex(ValueError, '缺少.*toppar/missing.itp'):
                md_charmm_gui._topology_files(root)


class MDSetupTests(unittest.TestCase):
    def setUp(self):
        cmd.reinitialize()
        cmd.fab('ASTYKKG', 'peptide', chain='A')

    def modified_pdb(self, codes):
        source = cmd.get_model('peptide')
        requests = [(site, code) for site, code in zip(ptm.list_sites(source), codes)]
        result, sites, _ = ptm.build_model(source, requests)
        ptm.load_result(cmd, 'ptm_fixture', result, sites, False)
        return cmd.get_pdbstr('%ptm_fixture', state=1)

    def test_detects_phosphorylated_and_lactyl_residues(self):
        pdb = self.modified_pdb(('SEP', 'TPO', 'PTR', 'XRW', 'D8R'))
        residues, waters = md_setup.inspect_pdb(pdb)
        self.assertEqual(waters, 0)
        self.assertEqual([key[2] for key in residues],
                         ['ALA', 'SEP', 'TPO', 'PTR', 'XRW', 'D8R', 'GLY'])
        _, mapping = md_setup._normalized_pdb(residues, {})
        self.assertEqual([entry['charge'] for entry in mapping], [None, -2, -2, -2, 0, 0, None])
        with self.assertRaisesRegex(ValueError, '构型不符'):
            # Reflecting the whole input changes lactyl chirality while
            # leaving bond lengths, atom names and residue order unchanged.
            flipped = []
            for row in pdb.splitlines():
                if row.startswith(('ATOM  ', 'HETATM')):
                    row = row[:30] + f'{-float(row[30:38]):8.3f}' + row[38:]
                flipped.append(row)
            md_setup.inspect_pdb('\n'.join(flipped))

    def test_missing_lactyl_template_stops_before_running_tools(self):
        pdb = self.modified_pdb(('SEP', 'TPO', 'PTR', 'XRW', 'D8R'))
        with TemporaryDirectory() as temp:
            root = Path(temp)
            output = root / 'blocked.zip'
            with self.assertRaisesRegex(ValueError, '缺少 ff14SB 蛋白内残基模板'):
                md_setup.build_package(pdb, output, md_setup.Runtime(root), library=root / 'library')
            self.assertFalse(output.exists())

    def test_monoanion_tyrosine_accepts_ccd_and_amber_oxygen_names(self):
        pdb = self.modified_pdb(('SEP', 'TPO', 'PTR')).replace(' PTR A', ' Y1P A')
        renamed = pdb.replace(' OH  Y1P A', ' OG  Y1P A')
        self.assertNotEqual(renamed, pdb)
        outputs = []
        for source in (pdb, renamed):
            residues, _ = md_setup.inspect_pdb(source)
            result, mapping = md_setup._normalized_pdb(residues, parameters.PHOSPHATE_ATOM_NAMES)
            entry = next(row for row in mapping if row['amber_residue'] == 'Y1P')
            self.assertEqual(entry['charge'], -1)
            self.assertIn('OG', entry['atom_names'])
            outputs.append(result)
        self.assertEqual(*outputs)

    def test_archived_parameters_reused_and_tampering_rejected(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            stock = root / 'installed'
            for name in (*parameters.PHOS_FILES, *parameters.BASE_FILES):
                file = stock / name
                file.parent.mkdir(exist_ok=True, parents=True)
                file.write_text(name)
            archive, record = parameters.archived_phosphate(stock, root / 'archive')
            manifest = archive / 'template.json'
            legacy = dict(record)
            del legacy['base_sha256']
            manifest.write_text(json.dumps(legacy))
            self.assertEqual(parameters.archived_phosphate(stock, root / 'archive'), (archive, record))
            for name in parameters.PHOS_FILES:
                (stock / name).unlink()
            self.assertEqual(parameters.archived_phosphate(stock, root / 'archive'), (archive, record))
            base = stock / parameters.BASE_FILES[0]
            base.write_text('changed base')
            with self.assertRaisesRegex(ValueError, '基础文件'):
                parameters.archived_phosphate(stock, root / 'archive')
            base.write_text(parameters.BASE_FILES[0])
            (archive / 'frcmod.phosaa14SB').write_text('changed parameters')
            with self.assertRaisesRegex(ValueError, 'SHA'):
                parameters.archived_phosphate(stock, root / 'archive')

    def test_unknown_residue_and_incomplete_chain_rejected(self):
        pdb = cmd.get_pdbstr('%peptide', state=1)
        unknown = pdb.replace('ALA A', 'LIG A')
        with self.assertRaisesRegex(ValueError, 'LIG'):
            md_setup.inspect_pdb(unknown)
        multichain = pdb.replace('GLY A', 'GLY B')
        with self.assertRaisesRegex(ValueError, '每条蛋白链'):
            md_setup.inspect_pdb(multichain)

    def test_combine_keeps_scene_coordinates_and_maps_duplicate_chain_names(self):
        cmd.fab('LKT', 'other_peptide', chain='A')
        cmd.translate([26, 0, 0], object='other_peptide', camera=0)
        objects = ('peptide', 'other_peptide')
        pdb, origins = md_setup.combine_protein_sources(
            (name, 1, cmd.get_pdbstr('%' + name, state=1)) for name in objects)
        residues, _ = md_setup.inspect_pdb(pdb)
        self.assertEqual(list(origins), ['A', 'B'])
        self.assertEqual(origins['B'], {'object': 'other_peptide', 'state': 1, 'source_chain': 'A'})
        self.assertEqual(len(residues), 10)
        self.assertEqual(residues[('B', '1', 'LEU')]['N'][2],
                         [round(float(c), 3) for c in cmd.get_coords('%other_peptide & index 1')[0]])
        normalized, mapping = md_setup._normalized_pdb(residues, {}, origins)
        self.assertEqual(normalized.splitlines().count('TER'), 2)
        self.assertEqual(mapping[-1]['object'], 'other_peptide')

    def test_cross_chain_conect_reports_distance_and_deduplicates(self):
        cmd.fab('LKT', 'other_peptide', chain='B')
        pdb, _ = md_setup.combine_protein_sources(
            (name, 1, cmd.get_pdbstr('%' + name, state=1))
            for name in ('peptide', 'other_peptide'))
        atoms = [line for line in pdb.splitlines() if line.startswith('ATOM  ')]
        a = next(line for line in atoms if line[12:27] == ' N   ALA A   1 ')
        b = next(line for line in atoms if line[12:27] == ' N   LEU B   1 ')
        first, second = (line[6:11] for line in (a, b))
        pdb += f'CONECT{first}{second}{second}\nCONECT{second}{first}\n'
        with self.assertRaises(ValueError) as raised:
            md_setup.inspect_pdb(pdb, conect_source=md_setup.OBJECT_CONECT)
        message = str(raised.exception)
        self.assertIn('1 个跨链原子对', message)
        self.assertEqual(message.count('A/1 ALA N ↔ B/1 LEU N'), 1)
        self.assertIn('Å · CONECT（由当前对象键表导出', message)

    def test_ssbond_annotation_retained_through_chain_remapping(self):
        cmd.fab('ACCA', 'sulfur_peptide', chain='A')
        cmd.translate([30, 0, 0], object='sulfur_peptide', camera=0)
        source = cmd.get_pdbstr('%sulfur_peptide', state=1)
        source += 'SSBOND   1 CYS A    2    CYS A    3                          1555   1555  2.06\n'
        pdb, _ = md_setup.combine_protein_sources(
            (('peptide', 1, cmd.get_pdbstr('%peptide', state=1)), ('sulfur_peptide', 1, source)))
        with self.assertRaisesRegex(ValueError, r'B/2 CYS SG ↔ B/3 CYS SG .* PDB SSBOND'):
            md_setup.inspect_pdb(pdb)

    def test_lactyl_fragment_keeps_selected_atom_order_and_connectivity(self):
        previous = cmd.get_setting_boolean('pdb_conect_all')
        try:
            cmd.set('pdb_conect_all', 1)
            model = cmd.get_model('peptide')
            lysine = next(site for site in ptm.list_sites(model) if site.resn == 'LYS')
            result, sites, _ = ptm.build_model(model, [(lysine, 'XRW')])
            ptm.load_result(cmd, 'lactyl_fixture', result, sites, False)
            pdb = cmd.get_pdbstr('%lactyl_fixture', state=1)
        finally:
            cmd.set('pdb_conect_all', previous)
        residues, _ = md_setup.inspect_pdb(pdb)
        source = list(residues)
        site = next(key for key in source if key[2] == 'XRW')
        i = source.index(site)
        fragment = md_setup.lactyl_validation_peptide(pdb, site)
        selected, _ = md_setup.inspect_pdb(fragment)
        self.assertEqual(list(selected), source[i-1:i+2])
        self.assertEqual(list(selected[site]), list(residues[site]))
        self.assertIn('FORCE FIELD PARAMETERS: NONE', fragment)
        cmd.read_pdbstr(fragment, 'validation_fragment')
        for left, right in zip(source[i-1:i+1], source[i:i+2]):
            model = cmd.get_model(f'validation_fragment & ((resi {left[1]} & name C) | '
                                  f'(resi {right[1]} & name N))')
            self.assertEqual(len(model.bond), 1)

    def test_register_template_checks_base_files_and_never_replaces(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            base, library, candidate = (root / name for name in ('base', 'library', 'candidate'))
            for name in parameters.BASE_FILES:
                file = base / name
                file.parent.mkdir(exist_ok=True, parents=True)
                file.write_text(name)
            candidate.mkdir()
            atoms = [(name, {'C': 6, 'N': 7, 'O': 8}[element])
                     for name, (element, _) in parameters.TEMPLATES['XRW'].atoms.items()
                     if name != 'OXT']
            index = {name: i for i, (name, _) in enumerate(atoms, 1)}
            # A format-only fixture. No force-field calculation uses these files.
            off = ['!entry.XRW.unit.atoms table']
            off.extend(f'"{name}" "type" 0 1 0 {i} {number} 0'
                       for i, (name, number) in enumerate(atoms, 1))
            off.extend(('!entry.XRW.unit.connect array int', str(index['N']), str(index['C'])))
            (candidate / 'residue.lib').write_text('\n'.join(off) + '\n')
            (candidate / 'residue.frcmod').write_text('fixture')
            record = {
                'residue': 'XRW', 'force_field': 'ff14SB', 'configuration': '2S', 'charge': 0,
                'source': 'format-test fixture',
                'sha256': parameters.checksums(candidate, ('residue.lib', 'residue.frcmod')),
                'base_sha256': parameters.checksums(base, parameters.BASE_FILES),
            }
            manifest = candidate / 'template.json'
            manifest.write_text(json.dumps(record))
            target = parameters.register_lactyl_template(candidate, library, 'XRW', base)
            self.assertEqual(target, library / 'XRW')
            self.assertEqual(parameters.register_lactyl_template(candidate, library, 'XRW', base), target)
            record['source'] = 'different record'
            manifest.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, '已有 XRW'):
                parameters.register_lactyl_template(candidate, library, 'XRW', base)
            (base / parameters.BASE_FILES[0]).write_text('different base')
            with self.assertRaisesRegex(ValueError, '基础文件'):
                parameters.register_lactyl_template(candidate, root / 'other', 'XRW', base)
            self.assertFalse((root / 'other').exists())
