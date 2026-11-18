"""Prepare a portable GROMACS package for a solvated protein assembly."""

from collections import defaultdict
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import platform
import re
from string import ascii_lowercase, ascii_uppercase, digits
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from . import md_parameters as parameters
from .md_protocol import RUN_SCRIPT, STAGE_LABELS, stage_inputs
from .ptm import LACTYL_STEREO_ATOMS, _lactyl_configuration
from .ptm_templates import TEMPLATES

AMINO_ACIDS = frozenset('ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP HIS ILE LEU '
                        'LYS MET PHE PRO SER THR TRP TYR VAL ASH GLH LYN'.split())
WATER = {'HOH', 'WAT', 'SOL'}
PDB_CHAINS = ascii_uppercase + ascii_lowercase + digits
OBJECT_CONECT = 'CONECT（由当前对象键表导出，可能来自载入时的距离推断）'


@dataclass(frozen=True)
class Runtime:
    prefix: Path

    @property
    def leap_data(self):
        return self.prefix / 'dat' / 'leap'

    def executable(self, name):
        path = self.prefix / 'bin' / name
        # The Conda macOS SIMD wrapper loses DYLD_LIBRARY_PATH when its
        # system shell starts. Invoke its actual executable in that case.
        if name == 'gmx' and platform.system() == 'Darwin' and platform.machine() == 'arm64':
            native = self.prefix / 'bin.ARM_NEON_ASIMD' / name
            if native.is_file():
                path = native
        if not os.access(path, os.X_OK):
            raise ValueError(f'运行环境缺少 {name}，请选包含 AmberTools、ParmEd 和 GROMACS 的环境。')
        return path

    def environment(self):
        env = dict(os.environ)
        env['LC_ALL'] = 'C'
        env['PATH'] = str(self.prefix / 'bin') + os.pathsep + env.get('PATH', '')
        if platform.system() == 'Darwin':
            dirs = [str(self.prefix / 'lib'), str(self.prefix / 'lib.ARM_NEON_ASIMD')]
            env['DYLD_LIBRARY_PATH'] = os.pathsep.join(dirs + ([env['DYLD_LIBRARY_PATH']]
                                                                if env.get('DYLD_LIBRARY_PATH') else []))
        return env


def default_runtime():
    for path in (os.environ.get('PYMOLX_MD_RUNTIME'), os.environ.get('CONDA_PREFIX'),
                 Path.home() / 'miniconda3' / 'envs' / 'md_web'):
        if path and (Path(path) / 'bin' / 'tleap').is_file() and (Path(path) / 'bin' / 'gmx').is_file():
            return Path(path)
    return None


def combine_protein_sources(sources):
    """Combine (object, state, PDB) records without changing coordinates."""
    entries = []
    reserved = set()
    for name, state, pdb in sources:
        lines = pdb.splitlines()
        if sum(line.startswith('MODEL ') for line in lines) > 1:
            raise ValueError(f'对象 {name} 包含多个坐标状态。')
        chains = dict.fromkeys(line[21] for line in lines
                               if line.startswith(('ATOM  ', 'HETATM')) and len(line) > 21)
        if not chains:
            raise ValueError(f'对象 {name} 在所选状态没有原子。')
        reserved.update(chains)
        entries.append((name, state, lines, chains))
    available = iter(chain for chain in PDB_CHAINS if chain not in reserved)
    origins = {}
    coordinates, connections = [], []
    for name, state, lines, chains in entries:
        for chain in chains:
            label = chain if chain not in origins else next(available, None)
            if label is None:
                raise ValueError('所选链数超出 PDB 链标识容量。')
            chains[chain] = label
            origins[label] = {'object': name, 'state': state, 'source_chain': chain.strip()}
        serials = {}
        for line in lines:
            if line.startswith(('ATOM  ', 'HETATM')):
                serial = str(len(coordinates) + 1)
                if len(serial) > 5 or line[6:11] in serials or len(line) < 54:
                    raise ValueError(f'对象 {name} 的 PDB 原子编号或记录无效。')
                serials[line[6:11]] = serial.rjust(5)
                coordinates.append(line[:6] + serials[line[6:11]] + line[11:21]
                                   + chains[line[21]] + line[22:])
        for line in lines:
            if line.startswith('CONECT') and line[6:11] in serials:
                neighbors = [serials[line[i:i+5]] for i in range(11, len(line), 5)
                             if line[i:i+5] in serials]
                if neighbors:
                    connections.append('CONECT' + serials[line[6:11]] + ''.join(neighbors))
            elif line.startswith(('LINK  ', 'SSBOND')):
                first, second = (21, 51) if line.startswith('LINK  ') else (15, 29)
                if len(line) > second and line[first] in chains and line[second] in chains:
                    connections.append(line[:first] + chains[line[first]] + line[first+1:second]
                                       + chains[line[second]] + line[second+1:])
    return '\n'.join([*coordinates, *connections, 'END', '']), origins


def _check_protein_connections(lines, residues, conect_source):
    atoms_by_id = {}
    atoms_by_site = {}
    for (chain, number, code), atoms in residues.items():
        for row, _, xyz in atoms.values():
            site = (chain, number, code, row[12:16].strip())
            atoms_by_id[row[6:11]] = site
            atoms_by_site[site] = xyz
    pairs = defaultdict(set)
    for line in lines:
        if line.startswith('CONECT'):
            left = atoms_by_id.get(line[6:11])
            for i in range(11, len(line), 5):
                right = atoms_by_id.get(line[i:i+5])
                if left and right and (left[0] != right[0] or (
                        left[2:] == right[2:] == ('CYS', 'SG') and left != right)):
                    pairs[tuple(sorted((left, right)))].add(conect_source)
        elif line.startswith('LINK  ') and len(line) >= 57:
            left = (line[21], line[22:27].strip(), line[17:20].strip(), line[12:16].strip())
            right = (line[51], line[52:57].strip(), line[47:50].strip(), line[42:46].strip())
            if left[0] != right[0]:
                pairs[tuple(sorted((left, right)))].add('PDB LINK')
        elif line.startswith('SSBOND') and len(line) >= 36:
            left = (line[15], line[17:22].strip(), 'CYS', 'SG')
            right = (line[29], line[31:36].strip(), 'CYS', 'SG')
            pairs[tuple(sorted((left, right)))].add('PDB SSBOND')
    cysteines = [site for site in atoms_by_site if site[2:] == ('CYS', 'SG')]
    for i, left in enumerate(cysteines):
        for right in cysteines[i+1:]:
            pair = tuple(sorted((left, right)))
            if pair not in pairs and math.dist(atoms_by_site[left], atoms_by_site[right]) < 2.4:
                pairs[pair].add('SG–SG 近接触，小于 2.40 Å，需核对二硫键')
    if pairs:
        cross_count = sum(left[0] != right[0] for left, right in pairs)
        headline = f'检测到 {cross_count} 个跨链原子对，需要核查连接。' if cross_count else '检测到可能的二硫键。'
        details = [headline]
        for (left, right), evidence in sorted(pairs.items()):
            distance = (f'{math.dist(atoms_by_site[left], atoms_by_site[right]):.3f} Å'
                        if left in atoms_by_site and right in atoms_by_site else '坐标缺失')
            first, second = (f'{chain.strip() or "无链"}/{number} {code} {atom}'
                             for chain, number, code, atom in (left, right))
            details.append(f'{first} ↔ {second} · {distance} · {"、".join(sorted(evidence))}')
        details.append('请对照原始结构核查。真实共价连接需要单独准备交联拓扑。')
        raise ValueError('\n'.join(details))


def inspect_pdb(text, *, conect_source='PDB CONECT'):
    """Check the input and retain the heavy atoms of each protein chain."""
    lines = text.splitlines()
    if sum(line.startswith('MODEL ') for line in lines) > 1:
        raise ValueError('请选择单个坐标状态。')
    residues = defaultdict(dict)
    discarded_water = set()
    for line in lines:
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        name, code = line[12:16].strip(), line[17:20].strip()
        if code == 'Y1P' and name == 'OG':
            name = 'OH'
        key = (line[21:22], line[22:27].strip(), code)
        if code in WATER:
            discarded_water.add(key)
            continue
        if code not in AMINO_ACIDS | parameters.PTM_CODES:
            raise ValueError(f'{key[0].strip()}/{key[1]} 的 {code} 不在蛋白打包范围。请先确认组分。')
        element = line[76:78].strip().upper() or name.lstrip('0123456789')[:1].upper()
        if element in ('H', 'D'):
            if parameters.PHOSPHATE.get(code) == -2 and name in (
                'H1P', 'HOP2', 'HOP3', 'HO2P', 'HO3P'
            ):
                raise ValueError(f'{key[1]} {code} 含磷酸氢，请先核对单负或双负状态。')
            continue
        if len(line) < 54 or line[16].strip() or not name or name in residues[key]:
            raise ValueError(f'{key[0].strip()}/{key[1]} 有替代构象、重名或不完整的原子记录。')
        try:
            xyz = [float(line[offset:offset+8]) for offset in (30, 38, 46)]
        except ValueError as exc:
            raise ValueError(f'{key[1]} {name} 的 PDB 坐标字段无效。') from exc
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError(f'{key[1]} {name} 的坐标不是有限数值。')
        residues[key][name] = (line, element, xyz)
    chains = defaultdict(list)
    for key in residues:
        chains[key[0]].append(key)
    keys = [key for chain in chains.values() for key in chain]
    if not keys or any(len(chain) < 2 for chain in chains.values()):
        raise ValueError('每条蛋白链至少需要两个残基。')
    if len(keys) > 9999 or len({key[:2] for key in keys}) != len(keys):
        raise ValueError('残基编号有歧义或超出 PDB 容量。')

    _check_protein_connections(lines, residues, conect_source)
    termini = {end for chain in chains.values() for end in (chain[0], chain[-1])}
    for i, key in enumerate(keys):
        chain, number, code = key
        atoms = residues[key]
        if not {'N', 'CA', 'C', 'O'} <= atoms.keys():
            raise ValueError(f'{number} {code} 缺少蛋白主链重原子。')
        if i and keys[i-1][0] == chain and np.linalg.norm(
                np.subtract(residues[keys[i-1]]['C'][2], atoms['N'][2])) > 1.9:
            raise ValueError(f'{keys[i-1][1]}–{number} 主链断开，请先处理缺失残基或链边界。')
        if code not in parameters.PTM_CODES:
            continue
        if key in termini:
            raise ValueError(f'{number} {code} 位于链端，需要独立核验的末端模板。')
        ccd = parameters.PHOSPHATE_PARENT.get(code, code)
        template = TEMPLATES[ccd]
        expected = set(template.atoms) - {'OXT'}
        if (atoms.keys() != expected or
                any(atoms[name][1] != template.atoms[name][0] for name in expected)):
            raise ValueError(f'{number} {code} 的重原子名称或元素不符合 PTM 模板。')
        if code in parameters.LACTYL:
            names = {name: j for j, name in enumerate(atoms)}
            geometry = np.asarray([entry[2] for entry in atoms.values()])
            if _lactyl_configuration(geometry, names, code) != parameters.LACTYL[code]:
                raise ValueError(f'{number} {code} 的乳酸侧链坐标与预期构型不符。')
    return {key: residues[key] for key in keys}, len(discarded_water)


def lactyl_validation_peptide(pdb, site):
    """Extract the modified site and its two original neighbors for validation."""
    residues, _ = inspect_pdb(pdb)
    sites = [key for key in residues if key[0] == site[0]]
    if site not in residues or site[2] not in parameters.LACTYL:
        raise ValueError('所选乳酸化位点已变化，请刷新。')
    index = sites.index(site)
    if index == 0 or index == len(sites) - 1:
        raise ValueError('位点需要同时具备 N 端和 C 端相邻残基。')
    atoms = [entry[0] for key in sites[index-1:index+2] for entry in residues[key].values()]
    serials = {int(line[6:11]) for line in atoms}
    if len(serials) != len(atoms):
        raise ValueError('短肽范围内的 PDB 原子序号有重复。')
    connections = []
    for row in pdb.splitlines():
        if not row.startswith('CONECT') or int(row[6:11]) not in serials:
            continue
        numbers = [int(row[i:i+5]) for i in range(11, len(row), 5) if row[i:i+5].strip()]
        retained = ''.join(f'{number:5d}' for number in numbers if number in serials)
        if retained:
            connections.append(row[:11] + retained)
    chain, number, code = site
    header = [
        'REMARK 900 PYMOL-VMD-X INTERNAL PTM STRUCTURE VALIDATION FRAGMENT',
        f'REMARK 900 SOURCE SITE {chain.strip() or "_"} {number} {code} ({parameters.LACTYL[code]})',
        'REMARK 900 MD TOPOLOGY AND FORCE FIELD PARAMETERS: NONE',
    ]
    return '\n'.join([*header, *atoms, 'TER', *connections, 'END', ''])


def _normalized_pdb(residues, aliases, origins=None):
    """Renumber in input order and keep a mapping to the original chain and IDs."""
    lines = []
    mapping = []
    index = 0
    origins = origins or {}
    for i, ((chain, resi, code), atoms) in enumerate(residues.items(), 1):
        if mapping and mapping[-1]['chain'] != chain.strip():
            lines.append('TER')
        amber_code = 'HIE' if code == 'HIS' else code
        names = {name: aliases.get(code, {}).get(name, name) for name in atoms}
        for name, (line, _, _) in atoms.items():
            atom = names[name]
            field = atom if len(atom) == 4 else f' {atom:<3}'
            index += 1
            if index > 99999:
                raise ValueError('超出 PDB 的原子序号容量。')
            lines.append(f'ATOM  {index:5d} {field} {amber_code:>3} {chain}{i:4d} {line[27:76]}{line[76:78]:>2}')
        record = {
            'chain': chain.strip(), 'resi': resi, 'source_residue': code,
            'amber_residue': amber_code, 'atom_names': sorted(names.values()),
            'charge': parameters.PHOSPHATE.get(code, 0 if code in parameters.LACTYL else None),
        }
        record.update(origins.get(chain, {'source_chain': chain.strip()}))
        if code in parameters.LACTYL:
            record['lactyl_stereo_atoms'] = [names[name] for name in LACTYL_STEREO_ATOMS[code]]
            record['lactyl_configuration'] = parameters.LACTYL[code]
            record['lactyl_bonds'] = [[names[a], names[b]]
                                     for a, b, _ in TEMPLATES[code].bonds if a in names and b in names]
        mapping.append(record)
    return '\n'.join([*lines, 'TER', 'END', '']), mapping


def _run(runtime, work, tool, args, log, timeout=300):
    command = [str(runtime.executable(tool)), *args]
    result = subprocess.run(command, cwd=work, env=runtime.environment(),
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, errors='replace',
                            timeout=timeout)
    (work / log).write_text(result.stdout, encoding='utf-8')
    if result.returncode:
        lines = result.stdout.splitlines()
        relevant = {j for i, line in enumerate(lines) if re.match(r'^\s*(?:ERROR|WARNING|Line)', line)
                    for j in range(i, min(i + 5, len(lines)))}
        relevant.update(range(max(0, len(lines) - 17), len(lines)))
        raise ValueError(f'{tool} 执行失败（{result.returncode}）。\n' +
                         '\n'.join(lines[i] for i in sorted(relevant)))
    return result.stdout


_CONVERT = r"""import json
import re
from pathlib import Path
import numpy as np
import parmed

mapping = json.loads(Path('residue-map.json').read_text())
source_xyz = {(int(line[22:26]), line[12:16].strip()): [float(line[i:i+8]) for i in (30, 38, 46)]
              for line in Path('protein.pdb').read_text().splitlines() if line.startswith('ATOM  ')}
structure = parmed.load_file('system.prmtop', xyz='system.rst7')
solute = structure.residues[:len(mapping)]
if len(solute) != len(mapping) or abs(sum(atom.charge for atom in structure.atoms)) > 1e-3:
    raise ValueError('Protein residue count or neutralized system charge is inconsistent.')
bonds = {frozenset((bond.atom1.idx, bond.atom2.idx)) for bond in structure.bonds}
previous_c = None
chain_indices = {}
displacements = []
for i, (residue, original) in enumerate(zip(solute, mapping)):
    atoms = {a.name: a for a in residue.atoms if a.atomic_number > 1}
    expected = set(original['atom_names'])
    same_chain = i > 0 and original['chain'] == mapping[i-1]['chain']
    terminal = i == len(solute) - 1 or original['chain'] != mapping[i+1]['chain']
    if (residue.name != original['amber_residue']
            or expected - atoms.keys()
            or atoms.keys() - expected - ({'OXT'} if terminal else set())):
        raise ValueError(f"Unresolved heavy atoms at {original['chain']}/{original['resi']} "
                         f"{original['source_residue']}: missing {sorted(expected - atoms.keys())}, "
                         f"added {sorted(atoms.keys() - expected)}; repair the source before export.")
    if same_chain and frozenset((previous_c, atoms['N'].idx)) not in bonds:
        raise ValueError(f"Missing peptide bond at position {i+1}.")
    charge = original['charge']
    if charge is not None and abs(sum(a.charge for a in residue.atoms) - charge) > 1e-4:
        raise ValueError(f"Parameter charge mismatch at position {i+1}.")
    if 'lactyl_configuration' in original:
        for a, b in original['lactyl_bonds']:
            if frozenset((atoms[a].idx, atoms[b].idx)) not in bonds:
                raise ValueError(f'Missing XRW/D8R bond {a}–{b} at position {i+1}.')
        stereo_atoms = [atoms[name] for name in original['lactyl_stereo_atoms']]
        vec = np.array([[a.xx, a.xy, a.xz] for a in stereo_atoms])
        signed_volume = np.linalg.det(vec[1:] - vec[0])
        stereochemistry = '2S' if signed_volume > 0 else '2R'
        if not abs(signed_volume) >= 0.1 or stereochemistry != original['lactyl_configuration']:
            raise ValueError(f'XRW/D8R stereochemistry mismatch at position {i+1}.')
    previous_c = atoms['C'].idx
    chain_indices.setdefault(original['chain'], []).extend(atom.idx + 1 for atom in residue.atoms)
    displacements.extend(np.subtract([atoms[name].xx, atoms[name].xy, atoms[name].xz],
                                     source_xyz[(i+1, name)]) for name in expected)
if not np.allclose(displacements, displacements[0], atol=0.003, rtol=0):
    raise ValueError('Solvation changed the input protein geometry or relative chain positions.')

split = structure.split()
actual_molecules = {}
for atom in structure.atoms:
    actual_molecules.setdefault(atom.marked - 1, []).append(atom.idx + 1)
if list(actual_molecules.values())[:len(chain_indices)] != list(chain_indices.values()):
    raise ValueError('Covalent topology does not match the input chain boundaries.')
if [i for indices in actual_molecules.values() for i in indices] != list(range(1, len(structure.atoms) + 1)):
    raise ValueError('Coordinate order does not match the topology molecule order.')
parmed.gromacs.GromacsGroFile.write(structure, 'system.gro', precision=5, combine='all')
structure.save('system.top', overwrite=True)
topology = Path('system.top')
contents = topology.read_text()
sections = list(re.finditer(r'^\[\s*(?:moleculetype|system)\s*\]', contents, re.M))
if len(sections) != len(split) + 1:
    raise ValueError('Unexpected number of GROMACS molecule definitions.')
protein_types = [(index, molecule) for index, (molecule, copies) in enumerate(split)
                 if min(copies) < len(chain_indices)]
for index, molecule in reversed(protein_types):
    name = 'posre.itp' if len(protein_types) == 1 else f'posre-{index+1}.itp'
    Path(name).write_text('[ position_restraints ]\n' +
        ''.join(f'{atom.idx+1:7d} 1 1000 1000 1000\n' for atom in molecule.atoms if atom.atomic_number > 1))
    insert = sections[index+1].start()
    contents = contents[:insert] + f'#ifdef POSRES\n#include "{name}"\n#endif\n\n' + contents[insert:]
topology.write_text(contents)
natoms = sum(map(len, chain_indices.values()))
groups = {
    'System': list(range(1, len(structure.atoms) + 1)),
    'Protein': list(range(1, natoms + 1)),
    'NonProtein': list(range(natoms + 1, len(structure.atoms) + 1)),
    'PTM': [a.idx + 1 for r, m in zip(solute, mapping) if m['charge'] is not None for a in r.atoms],
    **{f'Chain_{chain or "blank"}': indices for chain, indices in chain_indices.items()},
}
with Path('index.ndx').open('w') as stream:
    for group, indices in groups.items():
        if indices:
            stream.write(f'[ {group} ]\n')
            for i in range(0, len(indices), 15):
                stream.write(' '.join(map(str, indices[i:i+15])) + '\n')
print(f'ParmEd {parmed.__version__}; {len(chain_indices)} chains; '
      f'{len(solute)} protein residues; {len(structure.atoms)} system atoms')
"""


def _preflight_stage(runtime, work, stage, *, label=None, topology='system.top', coordinates='system.gro'):
    name = f'{stage}-preflight'
    try:
        log = _run(runtime, work, 'gmx', [
            'grompp', '-f', f'mdp/{stage}.mdp', '-c', coordinates, '-r', coordinates,
            '-n', 'index.ndx', '-p', topology, '-o', f'{name}.tpr', '-po', f'{name}.mdp',
        ], f'{name}.log')
        ignored = [line for line in log.splitlines() if 'Ignoring obsolete mdp entry' in line]
        if ignored:
            raise ValueError('\n'.join(ignored))
        processed = {}
        for line in (work / f'{name}.mdp').read_text().splitlines():
            key, sep, value = line.partition(';')[0].partition('=')
            if sep:
                processed[key.strip()] = value.strip()
        minimizers = {'steep', 'cg', 'l-bfgs'}
        dynamics = {'md', 'md-vv', 'md-vv-avek', 'sd', 'bd'}
        integrators = minimizers if stage.startswith('em') else dynamics
        if processed.get('integrator') not in integrators:
            raise ValueError(f'integrator = {processed.get("integrator")} 无法作为此阶段运行。')
        if int(processed['nsteps']) <= 0:
            raise ValueError('nsteps 必须为正整数，启动脚本据此判断阶段是否完成。')
        (work / f'{name}.tpr').unlink()
        return processed
    except (ValueError, subprocess.TimeoutExpired) as exc:
        stage_name = label or f'{STAGE_LABELS[stage]}（{stage}.mdp）'
        raise ValueError(f'{stage_name} 预检失败。\n{exc}') from exc


def _check_zip_path(output):
    path = Path(output).expanduser().absolute()
    if path.suffix.lower() != '.zip' or not path.parent.is_dir() or path.exists():
        raise ValueError('请选择已存在的输出目录及一个尚未使用的 ZIP 文件名。')
    return path


def _write_zip(work, output):
    archive_path = work.parent / 'package.zip'
    with ZipFile(archive_path, 'w', compression=ZIP_DEFLATED) as archive:
        for item in sorted(work.rglob('*')):
            if item.is_file():
                archive.write(item, item.relative_to(work))
    os.link(archive_path, output)


def build_package(pdb, output, runtime, *, library=parameters.LIBRARY, production_ns=10,
                  temperature=300, padding=12, mdp_overrides=None, chain_origins=None):
    """Prepare and validate a ZIP with preflight logs and a target-side launcher."""
    output = _check_zip_path(output)
    library = Path(library).expanduser()
    if not (math.isfinite(production_ns) and 0.001 <= production_ns <= 10000
            and math.isfinite(temperature) and 273 <= temperature <= 400
            and math.isfinite(padding) and 11 <= padding <= 30):
        raise ValueError('请核对模拟时长、温度和盒子边距。')
    residues, waters = inspect_pdb(pdb, conect_source=OBJECT_CONECT if chain_origins else 'PDB CONECT')
    chain_origins = chain_origins or {}
    stages = stage_inputs(temperature, round(production_ns * 500000))
    overrides = mdp_overrides or {}
    for stage, content in overrides.items():
        if stage not in stages:
            raise ValueError(f'未知 MDP 阶段 {stage}。')
        if not isinstance(content, str) or not content.strip() or '\0' in content:
            raise ValueError(f'{STAGE_LABELS[stage]}（{stage}.mdp）为空或含无效字符。')
    stages.update(overrides)
    codes = {key[2] for key in residues} & parameters.PTM_CODES
    with TemporaryDirectory(prefix='.pymolx-md-', dir=output.parent) as temporary:
        work = Path(temporary) / 'bundle'
        work.mkdir()
        commands, records, aliases = parameters.copy_parameters(codes, runtime.leap_data, library, work)
        base_hashes = parameters.checksums(runtime.leap_data, parameters.BASE_FILES)
        for code in codes & parameters.LACTYL.keys():
            parameters.check_lactyl_force_field(records[code], base_hashes)
        normalized, mapping = _normalized_pdb(residues, aliases, chain_origins)
        (work / 'source.pdb').write_text(pdb, encoding='utf-8')
        (work / 'protein.pdb').write_text(normalized, encoding='utf-8')
        (work / 'residue-map.json').write_text(json.dumps(mapping, indent=2) + '\n')
        chains = {}
        for record in mapping:
            chain = record['chain']
            chains.setdefault(chain, {
                key: record[key] for key in ('object', 'state', 'source_chain') if key in record
            }).setdefault('index_group', 'Chain_' + (chain or 'blank'))
        (work / 'chain-map.json').write_text(json.dumps(chains, indent=2, ensure_ascii=False) + '\n',
                                            encoding='utf-8')
        leap = [
            'source leaprc.protein.ff14SB', 'source leaprc.water.tip3p', *commands,
            'protein = loadpdb protein.pdb',
            'check protein', f'solvatebox protein TIP3PBOX {padding:g}',
            'addionsrand protein Na+ 0', 'addionsrand protein Cl- 0',
            'saveamberparm protein system.prmtop system.rst7', 'quit',
        ]
        (work / 'prepare.leap').write_text('\n'.join(leap) + '\n')
        (work / 'convert.py').write_text(_CONVERT)
        gmx = _run(runtime, work, 'gmx', ['--version'], 'gromacs-preparation.txt', timeout=20)
        leap_log = _run(runtime, work, 'tleap', ['-f', 'prepare.leap'], 'prepare.log')
        if not re.search(r'Exiting LEaP: Errors = 0;', leap_log) or not (work / 'system.prmtop').is_file():
            raise ValueError('AmberTools 未生成有效拓扑。\n' + '\n'.join(leap_log.splitlines()[-16:]))
        conversion = _run(runtime, work, 'python', ['convert.py'], 'conversion.log')
        mdps = work / 'mdp'
        mdps.mkdir()
        processed = {}
        for stage, settings in stages.items():
            (mdps / f'{stage}.mdp').write_bytes(settings.encode('utf-8'))
            processed[stage] = _preflight_stage(runtime, work, stage)
        production = processed['md']
        steps = int(production['nsteps'])
        timestep = float(production['dt'])
        launcher = work / 'run_md.sh'
        launcher.write_text(RUN_SCRIPT)
        launcher.chmod(0o755)
        conda_record = next((runtime.prefix / 'conda-meta').glob('ambertools-*.json'), None)
        amber_version = 'see prepare.log'
        if conda_record:
            installed = json.loads(conda_record.read_text(encoding='utf-8'))
            amber_version = f"{installed['version']} ({installed['build']})"
        manifest = {
            'base_force_field': 'ff14SB', 'water': 'TIP3P', 'ambertools': amber_version,
            'base_sha256': base_hashes, 'ptm_parameters': records,
            'residue_map': 'residue-map.json', 'chain_map': 'chain-map.json',
            'removed_input_waters': waters,
            'ions': 'Na+/Cl- neutralization; no background salt',
            'termini': 'charged N/C termini; HIS defaults to HIE',
            'default_temperature_K': temperature,
            'production_steps': steps, 'production_ns': steps * timestep / 1000,
            'box_padding_A': padding, 'timestep_ps': timestep, 'custom_mdp_stages': list(overrides),
            'preflight_gromacs': next((line.split(':', 1)[-1].strip() for line in gmx.splitlines()
                                       if line.startswith('GROMACS version:')), 'gromacs-preparation.txt'),
            'validation': 'Topology/coordinate assembly and grompp for all stages on the prepared system. '
                          'The launcher regenerates each TPR with the actual prior-stage files.',
            'conversion': conversion.strip(),
        }
        for option, coupling, key in (('ref-t', 'tcoupl', 'temperature_K'),
                                       ('ref-p', 'pcoupl', 'pressure_bar')):
            targets = production.get(option, '').split()
            if production.get(coupling, 'no').lower() != 'no' and len(targets) == 1:
                manifest[key] = float(targets[0])
        (work / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (work / 'README.txt').write_text(
            'Unzip on the target and run: bash run_md.sh\n'
            'Requires GROMACS >=2021 with GPU support. Set GMX_BIN to select an executable.\n'
            'For a CPU run: bash run_md.sh --cpu\n'
            'Stages: CPU minimization, NVT/NPT equilibration, then production. See mdp/ for actual settings.\n'
            'To resume, run the same command in this directory. Keep all stage files.\n'
            'Review equilibration and model suitability before interpreting production data.\n'
            'Force-field identity, residue mapping and preprocessing logs are in this package.\n'
        )
        _write_zip(work, output)
    return manifest


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdb', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, default=default_runtime())
    parser.add_argument('--library', type=Path, default=parameters.LIBRARY)
    parser.add_argument('--ns', type=float, default=10)
    args = parser.parse_args()
    if args.runtime is None:
        parser.error('provide --runtime for AmberTools, ParmEd and GROMACS')
    try:
        build_package(args.pdb.read_text(), args.output, Runtime(args.runtime.expanduser().absolute()),
                      library=args.library, production_ns=args.ns)
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        parser.exit(1, f'{exc}\n')
    print(args.output)


if __name__ == '__main__':
    main()
