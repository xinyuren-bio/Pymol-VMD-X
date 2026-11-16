"""Build and inspect single-state PTM structure models with explicit bonds."""

from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from chempy import Bond

from . import CmdException
from .ptm_templates import TEMPLATES


PARENT_CODES = {
    'SER': ('SEP',),
    'THR': ('TPO',),
    'TYR': ('PTR',),
    'LYS': ('XRW', 'D8R'),
}
BACKBONE_EXCEPT_CA = {'N', 'C', 'O', 'OXT'}


@dataclass(frozen=True)
class Site:
    segment: str
    chain: str
    resi: str
    resn: str

    @classmethod
    def from_atom(cls, atom):
        return cls(atom.segi, atom.chain, atom.resi, atom.resn)

    def label(self):
        chain = self.chain or '无链'
        return f'{chain}/{self.resi} {self.resn}' + (f' [{self.segment}]' if self.segment else '')

    def choices(self):
        return (self.resn,) if self.resn in TEMPLATES else PARENT_CODES.get(self.resn, ())


def _residues(model):
    records = defaultdict(list)
    for index, atom in enumerate(model.atom):
        records[Site.from_atom(atom)].append(index)
    return records


def list_sites(model):
    """List addressable standard and already modified residues."""
    return [site for site in _residues(model) if site.choices()]


def _index_pair(first, second):
    return tuple(sorted((first, second)))


def _unit(vector):
    length = np.linalg.norm(vector)
    if not np.isfinite(length) or length < 1e-6:
        raise ValueError('侧链定位原子的坐标退化。')
    return vector / length


def _frame(anchor, parent, reference):
    axis = _unit(parent - anchor)
    normal = _unit(np.cross(axis, reference - parent))
    return np.stack((axis, np.cross(normal, axis), normal))


def _place(template, names, coords, environment):
    """Align to the target side chain and sample rotations about its last bond."""
    reference = np.asarray([template.atoms[name][1] for name in template.frame])
    existing = coords[[names[name] for name in template.frame]]
    fragment = np.asarray([template.atoms[name][1] for name in template.fragment])
    fragment = (fragment - reference[0]) @ _frame(*reference).T
    angles = np.deg2rad(np.arange(0, 360, 15))
    x, y, z = fragment.T
    candidates = np.stack((
        np.broadcast_to(x, (len(angles), len(x))),
        np.cos(angles)[:, None] * y - np.sin(angles)[:, None] * z,
        np.sin(angles)[:, None] * y + np.cos(angles)[:, None] * z,
    ), axis=-1) @ _frame(*existing) + existing[0]

    # Score only nearby heavy atoms. The attachment atom and its immediate
    # parent are excluded from the score, which is a placement heuristic.
    neighbors = environment[np.linalg.norm(environment - existing[0], axis=1) < 8.0]
    distances = np.linalg.norm(candidates[:, :, None] - neighbors, axis=-1)
    scores = np.maximum(2.6 - distances, 0).reshape(len(angles), -1)
    best = np.argmin(np.sum(scores ** 2, axis=1))
    return candidates[best], int(np.count_nonzero(distances[best] < 2.2))


def _named_atoms(model, indices, site):
    names = {}
    alternates = set()
    for index in indices:
        atom = model.atom[index]
        if atom.symbol in ('H', 'D'):
            continue
        if atom.alt:
            alternates.add(atom.alt)
        if len(alternates) > 1 or atom.name in names:
            raise ValueError(f'{site.label()} 有替代构象或重名原子，请先选择单一构象。')
        names[atom.name] = index
    return names, next(iter(alternates), '')


def _add_bond(bonds, first, second, order):
    key = _index_pair(first, second)
    bond = Bond()
    bond.index = list(key)
    bond.order = order
    bonds[key] = bond


def _lactyl_configuration(coords, names, code):
    # CIP priority at the lactyl carbon is OH > carbonyl C > methyl C > H.
    center, oxygen, carbonyl, methyl = {
        'XRW': ('C24', 'O26', 'C22', 'C25'),
        'D8R': ('C2', 'O2', 'C1', 'C3'),
    }[code]
    stereo = (center, oxygen, carbonyl, methyl)
    if not all(name in names for name in stereo):
        return None
    center, *ligands = coords[[names[name] for name in stereo]]
    volume = np.linalg.det(ligands - center)
    if abs(volume) < 0.1:
        return '未能判定'
    return '2S' if volume > 0 else '2R'


def build_model(source_model, requests):
    """Return a ChemPy model, the output sites, and placement diagnostics.

    requests is an iterable of (Site, CCD code). Existing modified sites use
    their own CCD code and keep their coordinates while repairing local bonds.
    """
    requests = dict(requests)
    if not requests:
        raise ValueError('请至少选择一个位点。')
    model = deepcopy(source_model)
    records = _residues(model)
    bonds = {_index_pair(*bond.index): bond for bond in model.bond}
    omitted = set()
    modified = []
    warnings = []

    for site, code in requests.items():
        if site not in records or code not in site.choices():
            raise ValueError(f'{site.label()} 当前不存在或不支持所选修饰。')
        template = TEMPLATES[code]
        indices = records[site]
        names, alt = _named_atoms(model, indices, site)
        coords = np.asarray(model.get_coord_list(), dtype=float)
        existing = site.resn == code
        if existing:
            # Match only exact CCD atom names and elements. Absent atoms are
            # never fabricated when displaying a structure imported from disk.
            heavy_names = {name for name, index in names.items() if name in template.atoms
                           and model.atom[index].symbol == template.atoms[name][0]}
            expected = [(a, b, order) for a, b, order in template.bonds
                        if a in heavy_names and b in heavy_names]
            valid = {_index_pair(names[a], names[b]) for a, b, _ in expected}
            local = {names[name] for name in heavy_names}
            for key in tuple(bonds):
                if key[0] in local and key[1] in local and key not in valid:
                    del bonds[key]
            for first, second, order in expected:
                i, j = names[first], names[second]
                distance = np.linalg.norm(coords[i] - coords[j])
                if not 0.7 <= distance <= 2.3:
                    raise ValueError(f'{site.label()} 的 {first}–{second} 距离为 {distance:.2f} Å，请核对结构。')
                _add_bond(bonds, i, j, order)
            missing = set(template.fragment + template.frame[:1]) - heavy_names
            if missing:
                warnings.append(f'{site.label()} 缺少 {", ".join(sorted(missing))}，已连接现有原子。')
            if template.parent == 'LYS':
                found = _lactyl_configuration(coords, {n: names[n] for n in heavy_names}, code)
                expected_stereo = '2S' if code == 'XRW' else '2R'
                if found is not None and found != expected_stereo:
                    warnings.append(f'{site.label()} 坐标构型为 {found}，与 CCD 的 {expected_stereo} 不一致。')
        else:
            required = set(template.atoms) - set(template.fragment) - {'OXT'}
            missing = {name for name in required if name not in names
                       or model.atom[names[name]].symbol != template.atoms[name][0]}
            if missing:
                raise ValueError(f'{site.label()} 缺少建模所需原子 {", ".join(sorted(missing))}。')
            if set(template.fragment) & names.keys():
                raise ValueError(f'{site.label()} 已有修饰原子，请核对残基名称。')
            anchor, parent = map(names.__getitem__, template.frame[:2])
            attached = {key[1] if key[0] == anchor else key[0] for key in bonds if anchor in key}
            hydrogens = {index for index in attached if model.atom[index].symbol in ('H', 'D')}
            if attached - hydrogens != {parent}:
                raise ValueError(f'{site.label()} 的 {template.frame[0]} 已有附加键或缺少侧链连接。')
            # Also recognize named protons in inputs without explicit H bonds.
            proton_names = {
                'SER': {'HG', 'HG1'}, 'THR': {'HG1'}, 'TYR': {'HH'},
                'LYS': {'HZ1', 'HZ2', 'HZ3', '1HZ', '2HZ', '3HZ'},
            }[template.parent]
            hydrogens.update(i for i in indices if model.atom[i].symbol in ('H', 'D')
                             and model.atom[i].name in proton_names)
            omitted.update(hydrogens)
            environment = coords[[i for i, atom in enumerate(model.atom)
                                  if i not in (anchor, parent) and atom.symbol not in ('H', 'D')]]
            fragment, close_contacts = _place(template, names, coords, environment)
            if close_contacts:
                warnings.append(f'{site.label()} 有 {close_contacts} 处新增重原子近接触（<2.2 Å）。')
            for i in indices:
                model.atom[i].resn = code
                model.atom[i].partial_charge = 0.0
                model.atom[i].text_type = '??'
                model.atom[i].numeric_type = -9999
            original = model.atom[anchor]
            original.formal_charge = 0
            for name, coord in zip(template.fragment, fragment):
                atom = original.new_in_residue()
                atom.name = name
                atom.symbol = template.atoms[name][0]
                atom.coord = coord.tolist()
                atom.q, atom.b, atom.ss, atom.alt = original.q, original.b, original.ss, alt
                atom.formal_charge = -1 if name in ('O2P', 'O3P') else 0
                names[name] = len(model.atom)
                model.atom.append(atom)
            for first, second, order in template.bonds:
                if first in template.fragment or second in template.fragment:
                    _add_bond(bonds, names[first], names[second], order)
        modified.append(Site(site.segment, site.chain, site.resi, code))

    keep = [i for i in range(len(model.atom)) if i not in omitted]
    remap = {old: new for new, old in enumerate(keep)}
    model.atom = [model.atom[i] for i in keep]
    model.bond = []
    for (i, j), bond in bonds.items():
        if i in remap and j in remap:
            bond.index = [remap[i], remap[j]]
            model.bond.append(bond)
    model.index = None
    return model, tuple(modified), tuple(warnings)


def object_model(cmd, name, state=1):
    if name not in cmd.get_names('objects') or cmd.get_type(name) != 'object:molecule':
        raise ValueError('源分子对象已删除或更名，请刷新。')
    model = cmd.get_model('%' + name, state)
    if not model.atom:
        raise ValueError('当前状态没有坐标，请检查对象的状态。')
    return model


def select_sites(cmd, name, sites, selection, state=1, sidechain=False):
    keys = set(sites)
    atoms = object_model(cmd, name, state).atom
    indices = [a.index for a in atoms if Site.from_atom(a) in keys
               and (not sidechain or a.name not in BACKBONE_EXCEPT_CA and a.symbol not in ('H', 'D'))]
    if not indices:
        raise ValueError('所选位点当前没有可显示的原子。')
    cmd.select_list(selection, name, indices, state=state, mode='index')


def load_result(cmd, name, model, sites, preview):
    marker = cmd.get_unused_name('_ptm_atoms')
    try:
        cmd.load_model(model, name, state=1, zoom=0)
        if not cmd.get_setting_boolean('static_singletons'):
            cmd.set('state', 1, name)
        select_sites(cmd, name, sites, marker, sidechain=True)
        cmd.hide('everything', '%' + name)
        if not preview:
            cmd.show('cartoon', '%' + name)
        cmd.show('sticks', '%' + marker)
        cmd.color('cyan', '%' + marker + ' & elem C')
    except (CmdException, ValueError):
        cmd.delete(name)
        raise
    finally:
        cmd.delete(marker)


def export_pdb(cmd, name, filename):
    """Export coordinates and ordered CONECT records from the applied object."""
    model = object_model(cmd, name)
    if len(model.atom) > 99999 or any(
        len(a.chain) > 1 or len(a.segi) > 4 or len(a.resn) > 3 or len(a.name) > 4
        or not -999 <= a.resi_number <= 9999 for a in model.atom
    ):
        raise ValueError('对象超过 PDB 的链名、编号或原子名限制，无法无损导出。')
    # The native PDB writer can preserve bond orders through repeated CONECT.
    # Restore the writer settings after capturing the structure.
    with cmd.lockcm:
        conect_all = cmd.get_setting_boolean('pdb_conect_all')
        nodup = cmd.get_setting_boolean('pdb_conect_nodup')
        retain_ids = cmd.get_setting_boolean('pdb_retain_ids')
        try:
            cmd.set('pdb_conect_all', 1, updates=0)
            cmd.set('pdb_conect_nodup', 0, updates=0)
            cmd.set('pdb_retain_ids', 0, updates=0)
            pdb = cmd.get_pdbstr('%' + name, state=1)
        finally:
            cmd.set('pdb_conect_all', conect_all, updates=0)
            cmd.set('pdb_conect_nodup', nodup, updates=0)
            cmd.set('pdb_retain_ids', retain_ids, updates=0)
    if 'CONECT' not in pdb:
        raise ValueError('PDB 导出未写入连接记录。')
    codes = sorted({site.resn for site in list_sites(model) if site.resn in TEMPLATES})
    remarks = (
        'REMARK 900 PYMOL-VMD-X PTM STRUCTURE MODEL\n'
        'REMARK 900 MD TOPOLOGY AND FORCE FIELD PARAMETERS: NONE\n'
        f'REMARK 900 WWPDB CCD: {" ".join(codes)}\n'
        'REMARK 900 CCD CONFIGURATION XRW: L-LACTYL (2S); D8R: D-LACTYL (2R)\n'
    )
    Path(filename).write_text(remarks + pdb, encoding='utf-8')
