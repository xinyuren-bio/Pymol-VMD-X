"""Import a CHARMM-GUI Solution Builder project without changing its force field."""

from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
from tempfile import TemporaryDirectory

from .md_protocol import charmm_gui_run_script
from .md_setup import _check_zip_path, _preflight_stage, _run, _write_zip


SOLUTION_BUILDER_URL = 'https://charmm-gui.org/?doc=input/solution'
_INPUTS = ('topol.top', 'step3_input.gro', 'index.ndx')
_STAGE = re.compile(r'^step(\d+(?:\.\d+)*)_(minimization(?:_\w+)?|equilibration(?:_\w+)?|production)\.mdp$')
_INCLUDE = re.compile(r'^\s*#\s*include\s*["<]([^">]+)[">]', re.M)


def _read_inputs(source, work):
    """Copy only data from one GROMACS output directory in a TAR or TGZ."""
    try:
        with tarfile.open(source, mode='r:*') as archive:
            selected = {}
            for entry in archive:
                path = PurePosixPath(entry.name)
                if 'gromacs' not in path.parts or path.name.startswith('._'):
                    continue
                root = path.parts.index('gromacs') + 1
                if path.is_absolute() or '..' in path.parts or '\\' in entry.name:
                    raise ValueError(f'CHARMM-GUI 归档含无效路径 {entry.name}。')
                if entry.isdir() or path.suffix.lower() not in {'.top', '.itp', '.gro', '.ndx', '.mdp'}:
                    continue
                if not entry.isfile():
                    raise ValueError(f'拓扑或输入必须是普通文件，发现 {entry.name}。')
                relative = PurePosixPath(*path.parts[root:])
                files = selected.setdefault(path.parts[:root], {})
                if relative in files:
                    raise ValueError(f'CHARMM-GUI 归档含重复的输入 {entry.name}。')
                files[relative] = entry
            required = {PurePosixPath(name) for name in _INPUTS}
            projects = [files for files in selected.values() if required <= files.keys()]
            if len(projects) != 1:
                raise ValueError('需要一个完整的 Solution Builder GROMACS 作业，其中含 '
                                 'topol.top、step3_input.gro 和 index.ndx。')
            for relative, entry in projects[0].items():
                path = work.joinpath(*relative.parts)
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(entry) as data, path.open('xb') as output:
                    shutil.copyfileobj(data, output)
    except tarfile.TarError as exc:
        raise ValueError('请选择 CHARMM-GUI 下载的 .tgz、.tar.gz 或 .tar 归档。') from exc


def _topology_files(work):
    visited = {}
    pending = [work / 'topol.top']
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        contents = path.read_text(encoding='utf-8')
        visited[path] = contents
        for name in _INCLUDE.findall(contents):
            relative = PurePosixPath(name)
            if relative.is_absolute() or '..' in relative.parts or '\\' in name:
                raise ValueError(f'{path.name} 引用了归档目录外的拓扑 {name}。')
            match = next((folder / name for folder in (path.parent, work)
                          if (folder / name).is_file()), None)
            if match is None:
                raise ValueError(f'CHARMM-GUI 归档缺少 {path.name} 所引用的 {name}。')
            pending.append(match)
    for contents in visited.values():
        defaults = re.search(r'^\s*\[\s*defaults\s*\]([^\[]+)', contents, re.M)
        if defaults:
            row = next((s for line in defaults[1].splitlines()
                        if (s := line.partition(';')[0].split())), [])
            if len(row) >= 5 and row[:3] == ['1', '2', 'yes'] and tuple(map(float, row[3:5])) == (1., 1.):
                return sorted(path.relative_to(work).as_posix() for path in visited)
    raise ValueError('拓扑的 [ defaults ] 与 CHARMM 的非键和 1–4 约定不符。'
                     '请在 CHARMM-GUI 中选择 CHARMM36m 并重新导出。')


def _write_stages(work):
    ordered = []
    for path in work.glob('step*.mdp'):
        match = _STAGE.fullmatch(path.name)
        if match:
            order, kind = match.groups()
            ordered.append((tuple(map(int, order.split('.'))), kind.split('_')[0], path))
    ordered.sort()
    kinds = [kind for _, kind, _ in ordered]
    if ('minimization' not in kinds or 'equilibration' not in kinds or kinds.count('production') != 1
            or kinds != sorted(kinds, key=('minimization', 'equilibration', 'production').index)):
        raise ValueError('找不到完整的 CHARMM-GUI 最小化、平衡和生产阶段 MDP。')
    mdp_dir = work / 'mdp'
    mdp_dir.mkdir()
    prefixes = {'minimization': 'em', 'equilibration': 'eq', 'production': 'md'}
    stages = {}
    for _, kind, path in ordered:
        prefix = prefixes[kind]
        count = sum(stage.startswith(prefix) for stage in stages)
        stage = prefix + (str(count + 1) if count else '')
        settings = path.read_bytes()
        text = settings.decode('utf-8')
        if re.search(r'^\s*include\s*=\s*[^\s;]', text, re.M | re.I):
            raise ValueError(f'{path.name} 指定了额外的拓扑搜索路径，请提供自包含的 GROMACS 下载包。')
        (mdp_dir / f'{stage}.mdp').write_bytes(settings)
        stages[stage] = path.name
    return stages


def build_package(source, output, runtime):
    """Validate an imported GROMACS system and attach a one-command launcher."""
    output = _check_zip_path(output)
    source = Path(source).expanduser()
    with TemporaryDirectory(prefix='.pymolx-charmm-', dir=output.parent) as temporary:
        work = Path(temporary) / 'bundle'
        work.mkdir()
        _read_inputs(source, work)
        files = _topology_files(work)
        stages = _write_stages(work)
        version = _run(runtime, work, 'gmx', ['--version'], 'gromacs-preparation.txt', timeout=20)
        for stage, original in stages.items():
            _preflight_stage(runtime, work, stage, label=original,
                             topology='topol.top', coordinates='step3_input.gro')
        launcher = work / 'run_md.sh'
        launcher.write_text(charmm_gui_run_script(stages), encoding='utf-8')
        launcher.chmod(0o755)
        digest = sha256()
        with source.open('rb') as archive:
            for block in iter(lambda: archive.read(1048576), b''):
                digest.update(block)
        manifest = {
            'origin': 'CHARMM-GUI Solution Builder GROMACS archive',
            'archive': source.name, 'archive_sha256': digest.hexdigest(),
            'force_field': 'CHARMM (original imported topology)',
            'topology_files': files, 'stages': stages,
            'validation': 'GROMACS grompp for all stages on the starting coordinates',
            'preflight_gromacs': next((line.split(':', 1)[1].strip() for line in version.splitlines()
                                       if line.startswith('GROMACS version:')), 'gromacs-preparation.txt'),
        }
        (work / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        (work / 'PYMOLX-README.txt').write_text(
            'Run on the target machine: bash run_md.sh\n'
            'Use bash run_md.sh --cpu for testing without GPU; select GROMACS with GMX_BIN.\n'
            'The original GROMACS input files from CHARMM-GUI are included. See manifest.json for '
            'the stage name mapping. The runner uses the stage copies under mdp/.\n',
            encoding='utf-8',
        )
        _write_zip(work, output)
    return manifest
