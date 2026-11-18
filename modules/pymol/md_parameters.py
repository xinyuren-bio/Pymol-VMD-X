"""Reusable Amber PTM parameter records for the protein MD builder."""

from hashlib import sha256
import json
from pathlib import Path
import shlex
import shutil
from tempfile import TemporaryDirectory

from .ptm_templates import TEMPLATES

FF = 'ff14SB'
LIBRARY = Path.home() / '.pymolx' / 'md-parameters' / FF
PHOSPHATE = {'SEP': -2, 'TPO': -2, 'PTR': -2, 'S1P': -1, 'T1P': -1, 'Y1P': -1}
PHOSPHATE_PARENT = {'S1P': 'SEP', 'T1P': 'TPO', 'Y1P': 'PTR'}
# phosaa14SB calls the Tyr phenolic oxygen OG in Y1P and OH in PTR.
PHOSPHATE_ATOM_NAMES = {'Y1P': {'OH': 'OG'}}
LACTYL = {'XRW': '2S', 'D8R': '2R'}
PTM_CODES = PHOSPHATE.keys() | LACTYL.keys()
PHOS_SOURCE = 'https://doi.org/10.1021/acs.jctc.4c00732'
LACTYL_SOURCE = 'https://doi.org/10.1038/s41589-025-02107-8'
PHOS_FILES = {
    'lib/phosaa14SB.lib': 'phosaa14SB.lib',
    'parm/frcmod.phosaa14SB': 'frcmod.phosaa14SB',
    'cmd/leaprc.phosaa14SB': 'leaprc.phosaa14SB',
}
BASE_FILES = (
    'cmd/leaprc.protein.ff14SB', 'cmd/leaprc.water.tip3p', 'parm/parm10.dat',
    'parm/frcmod.ff14SB', 'lib/amino12.lib', 'lib/aminont12.lib',
    'lib/aminoct12.lib', 'lib/atomic_ions.lib', 'lib/solvents.lib',
    'parm/frcmod.tip3p', 'parm/frcmod.ionsjc_tip3p',
    'parm/frcmod.ions1lm_126_tip3p', 'parm/frcmod.ions234lm_126_tip3p',
)
LACTYL_FILES = ('residue.lib', 'residue.frcmod', 'template.json')


def digest(file):
    return sha256(file.read_bytes()).hexdigest()


def checksums(directory, names):
    return {name: digest(directory / name) for name in names}


def read_record(directory, files):
    record = json.loads((directory / 'template.json').read_text(encoding='utf-8'))
    if (not isinstance(record, dict) or record.get('force_field') != FF
            or record.get('sha256') != checksums(directory, files)):
        raise ValueError(f'参数库「{directory.name}」的力场或 SHA‑256 校验不匹配。')
    return record


def archived_phosphate(leap_data, library):
    """Seed the archive once from the installed AmberTools distribution."""
    base_hashes = checksums(leap_data, BASE_FILES)
    target = library / 'phosaa14SB'
    if not target.exists():
        sources = [leap_data / filename for filename in PHOS_FILES]
        if any(not path.is_file() for path in sources):
            raise ValueError('找不到 AmberTools 的 phosaa14SB 参数文件，请检查运行环境。')
        library.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=library) as temporary:
            pending = Path(temporary) / 'phosaa14SB'
            pending.mkdir()
            for path in sources:
                shutil.copyfile(path, pending / path.name)
            record = {
                'force_field': FF, 'source': PHOS_SOURCE,
                'residue_charges': PHOSPHATE,
                'base_sha256': base_hashes,
                'sha256': checksums(pending, PHOS_FILES.values()),
            }
            (pending / 'template.json').write_text(json.dumps(record, indent=2) + '\n')
            try:
                pending.rename(target)
            except OSError:
                if not (target / 'template.json').is_file():
                    raise
    record = read_record(target, PHOS_FILES.values())
    if record.get('residue_charges') != PHOSPHATE or record.get('source') != PHOS_SOURCE:
        raise ValueError('phosaa14SB 参数记录的来源或残基电荷与程序支持的配置不一致。')
    if 'base_sha256' not in record:
        # Bind archives created before the base-force-field fingerprints were recorded.
        for installed, stored in PHOS_FILES.items():
            path = leap_data / installed
            if not path.is_file() or digest(path) != record['sha256'][stored]:
                raise ValueError('旧版 phosaa14SB 归档与所选工具环境无法核对。请另选参数库。')
        record['base_sha256'] = base_hashes
        with TemporaryDirectory(dir=target) as temporary:
            pending = Path(temporary) / 'template.json'
            pending.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
            pending.replace(target / 'template.json')
    if record['base_sha256'] != base_hashes:
        raise ValueError('phosaa14SB 归档与所选 ff14SB/TIP3P 基础文件的 SHA‑256 不匹配。')
    return target, record


def _off_entries(file, code, subsection):
    lines = file.read_text(encoding='utf-8').splitlines()
    marker = f'!entry.{code}.unit.{subsection} '
    for i, line in enumerate(lines):
        if line.startswith(marker):
            for item in lines[i + 1:]:
                if item.startswith('!'):
                    break
                yield shlex.split(item)
            return


def _validate_lactyl(path, code):
    missing = [name for name in LACTYL_FILES if not (path / name).is_file()]
    if missing:
        raise ValueError(f'{code}（{LACTYL[code]}）缺少 ff14SB 蛋白内残基模板。'
                         f'待补文件 {", ".join(missing)}。')
    record = read_record(path, LACTYL_FILES[:2])
    aliases = record.get('atom_names', {})
    nonterminal = set(TEMPLATES[code].atoms) - {'OXT'}
    if (record.get('residue') != code or record.get('configuration') != LACTYL[code]
            or record.get('charge') != 0 or not isinstance(record.get('source'), str)
            or not record['source'].strip()
            or not isinstance(aliases, dict) or not aliases.keys() <= nonterminal
            or not all(isinstance(name, str) and name.isascii() and name.isalnum()
                       and len(name) <= 4 for name in aliases.values())):
        raise ValueError(f'{code} 参数记录的构型、电荷、来源或原子名映射无效。')
    library = path / 'residue.lib'
    names = [line.split()[0] for line in library.read_text(encoding='utf-8').splitlines()
             if line.startswith('!entry.') and '.unit.atoms ' in line]
    if names != [f'!entry.{code}.unit.atoms']:
        raise ValueError(f'OFF 库必须只定义一个名为 {code} 的单残基模板。')
    atoms = list(_off_entries(library, code, 'atoms'))
    connection = list(_off_entries(library, code, 'connect'))
    try:
        heavy = {row[0]: int(row[6]) for row in atoms if int(row[6]) != 1}
        mapped = {aliases.get(name, name): {'C': 6, 'N': 7, 'O': 8}[TEMPLATES[code].atoms[name][0]]
                  for name in nonterminal}
        total = sum(float(row[7]) for row in atoms)
        ends = [int(row[0]) for row in connection]
    except (IndexError, ValueError) as exc:
        raise ValueError(f'{code} OFF 库的原子或连接位点字段无效。') from exc
    if (len(mapped) != len(nonterminal) or mapped != heavy or not abs(total) <= 1e-4
            or len({row[0] for row in atoms}) != len(atoms)
            or len(ends) != 2 or not all(1 <= i <= len(atoms) for i in ends)
            or [atoms[i-1][0] for i in ends] != ['N', 'C']):
        raise ValueError(f'{code} OFF 库与 CCD 原子名、元素、净电荷或 N/C 肽键连接不匹配。')
    return record


def check_lactyl_force_field(record, base_hashes):
    if record.get('base_sha256') != base_hashes:
        raise ValueError(f'{record["residue"]} 参数记录与所选 ff14SB/TIP3P 基础文件的 SHA‑256 不匹配。')


def lactyl_template(library, code):
    """Read an independently validated, internal-residue Amber OFF template."""
    path = library / code
    record = _validate_lactyl(path, code)
    return path, record


def register_lactyl_template(path, library, code, leap_data):
    """Register one compatible external template with an atomic directory move."""
    record = _validate_lactyl(path, code)
    base_hashes = checksums(leap_data, BASE_FILES)
    check_lactyl_force_field(record, base_hashes)
    target = library / code
    if target.exists():
        if all((target / name).is_file() and digest(target / name) == digest(path / name)
               for name in LACTYL_FILES):
            return target
        raise ValueError(f'参数库已有 {code} 目录。请选择另一个参数库，或先核对并备份已有文件。')
    library.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=library) as temporary:
        pending = Path(temporary) / code
        pending.mkdir()
        for name in LACTYL_FILES:
            shutil.copyfile(path / name, pending / name)
        check_lactyl_force_field(_validate_lactyl(pending, code), base_hashes)
        pending.rename(target)
    return target


def copy_parameters(codes, leap_data, library, work):
    """Stage exact parameter files and return LEaP commands and provenance."""
    custom = {code: lactyl_template(library, code) for code in sorted(codes & LACTYL.keys())}
    records = {}
    commands = []
    if codes & PHOSPHATE.keys():
        source, record = archived_phosphate(leap_data, library)
        dest = work / 'parameters' / source.name
        dest.mkdir(parents=True)
        for name in (*PHOS_FILES.values(), 'template.json'):
            shutil.copyfile(source / name, dest / name)
        prefix = 'parameters/phosaa14SB/'
        script = (source / 'leaprc.phosaa14SB').read_text(encoding='utf-8')
        for name in ('phosaa14SB.lib', 'frcmod.phosaa14SB'):
            script = script.replace(name, prefix + name)
        (dest / 'load.leap').write_text(script, encoding='utf-8')
        commands.append('source ' + prefix + 'load.leap')
        records['phosaa14SB'] = record
    aliases = {code: PHOSPHATE_ATOM_NAMES[code] for code in codes & PHOSPHATE_ATOM_NAMES.keys()}
    for code, (source, record) in custom.items():
        dest = work / 'parameters' / code
        dest.mkdir(parents=True)
        for name in ('residue.lib', 'residue.frcmod', 'template.json'):
            shutil.copyfile(source / name, dest / name)
        commands.extend((f'loadamberparams parameters/{code}/residue.frcmod',
                         f'loadoff parameters/{code}/residue.lib'))
        records[code] = record
        aliases[code] = record.get('atom_names', {})
    return commands, records, aliases
