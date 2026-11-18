"""Qt export dialog for portable GROMACS protein simulations."""

from pathlib import Path

from pymol import CmdException, md_charmm_gui, md_parameters, md_setup
from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.Qt.utils import AsyncFunc
from pymol.md_protocol import STAGE_LABELS, stage_inputs


class StageMDPDialog(QtWidgets.QDialog):
    def __init__(self, parent, defaults, overrides):
        super().__init__(parent)
        self.setWindowTitle('GROMACS 阶段参数')
        self.resize(700, 560)
        self.defaults = defaults
        self.overrides = dict(overrides)
        self.editors = {}
        root = QtWidgets.QVBoxLayout(self)
        self.tabs = QtWidgets.QTabWidget()
        root.addWidget(self.tabs, 1)
        for stage, label in STAGE_LABELS.items():
            editor = QtWidgets.QPlainTextEdit()
            editor.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont))
            editor.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
            editor.setPlainText(overrides.get(stage, defaults[stage]))
            editor.textChanged.connect(lambda stage=stage: self._edited(stage))
            self.editors[stage] = editor
            self.tabs.addTab(editor, label)
        self.status = QtWidgets.QLabel()
        self.status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        controls = QtWidgets.QHBoxLayout()
        for title, handler in (('导入本阶段…', self._import), ('导出本阶段…', self._export),
                               ('恢复本阶段默认', self._reset)):
            button = QtWidgets.QPushButton(title)
            button.clicked.connect(handler)
            controls.addWidget(button)
        controls.addStretch(1)
        buttons = QtWidgets.QDialogButtonBox()
        buttons.addButton('应用参数', QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.addButton('取消', QtWidgets.QDialogButtonBox.ButtonRole.RejectRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        controls.addWidget(buttons)
        root.addLayout(controls)
        self.tabs.currentChanged.connect(self._update_status)
        self._update_status()

    def _stage(self):
        return tuple(STAGE_LABELS)[self.tabs.currentIndex()]

    def _update_status(self):
        stage = self._stage()
        origin = '自定义文件内容' if stage in self.overrides else '由表单生成默认值'
        self.status.setText(f'{stage}.mdp · {origin}。生成模拟包时逐阶段运行 grompp。')

    def _edited(self, stage):
        self.overrides[stage] = self.editors[stage].toPlainText()
        self._update_status()

    def _set_text(self, stage, content):
        editor = self.editors[stage]
        with QtCore.QSignalBlocker(editor):
            editor.setPlainText(content)
        self._update_status()

    def _import(self):
        stage = self._stage()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, '导入本阶段 MDP', '', 'MDP (*.mdp);;所有文件 (*)')
        if filename:
            try:
                content = Path(filename).read_bytes().decode('utf-8')
                if not content.strip() or '\0' in content:
                    raise ValueError('文件为空或含无效字符。')
            except (OSError, ValueError) as exc:
                self.status.setText(f'{stage}.mdp 导入失败。{exc}')
                return
            self.overrides[stage] = content
            self._set_text(stage, content)

    def _export(self):
        stage = self._stage()
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, '导出本阶段 MDP', f'{stage}.mdp', 'MDP (*.mdp)')
        if filename:
            try:
                content = self.overrides.get(stage, self.defaults[stage])
                path = Path(filename)
                if not path.suffix:
                    path = path.with_suffix('.mdp')
                path.write_bytes(content.encode('utf-8'))
            except OSError as exc:
                self.status.setText(f'{stage}.mdp 导出失败。{exc}')
            else:
                self.status.setText(f'已导出 {path}。')

    def _reset(self):
        stage = self._stage()
        self.overrides.pop(stage, None)
        self._set_text(stage, self.defaults[stage])


class MDSetupDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.cmd = parent.cmd
        self.worker = None
        self.lactyl_sites = []
        self.mdp_overrides = {}
        self.setWindowTitle('生成 GROMACS 模拟包')
        self.resize(640, 445)
        root = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        root.addLayout(form)

        self.objects = QtWidgets.QListWidget()
        self.objects.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.objects.setMaximumHeight(110)
        self.state = QtWidgets.QSpinBox()
        self.state.setMinimum(1)
        refresh = QtWidgets.QPushButton('刷新对象')
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.objects, 1)
        row.addWidget(QtWidgets.QLabel('状态'))
        row.addWidget(self.state)
        row.addWidget(refresh)
        form.addRow('蛋白对象', row)
        form.addRow('', QtWidgets.QLabel('按住 ⌘／Ctrl 单独选中多个对象，或用 Shift 连续选择。'))

        prefix = md_setup.default_runtime()
        self.runtime = self._folder_row(form, '工具环境', prefix or '',
                                        'GROMACS 环境，AMBER 打包另需 AmberTools 和 ParmEd')
        self.library = self._folder_row(form, 'PTM 参数库', md_parameters.LIBRARY,
                                        '保存各修饰的 ff14SB 参数文件和来源')
        form.addRow('模拟体系', QtWidgets.QLabel('AMBER ff14SB · TIP3P · Na⁺／Cl⁻ 中和'))
        charmm_row = QtWidgets.QHBoxLayout()
        self.charmm = QtWidgets.QPushButton('导入 CHARMM-GUI 作业…')
        site = QtWidgets.QPushButton('打开 Solution Builder')
        charmm_row.addWidget(self.charmm)
        charmm_row.addWidget(site)
        charmm_row.addStretch()
        form.addRow('CHARMM36m', charmm_row)

        self.production = QtWidgets.QDoubleSpinBox()
        self.production.setDecimals(3)
        self.production.setRange(0.001, 10000)
        self.production.setValue(10)
        self.production.setSuffix(' ns')
        form.addRow('生产时长', self.production)
        self.temperature = QtWidgets.QSpinBox()
        self.temperature.setRange(273, 400)
        self.temperature.setValue(300)
        self.temperature.setSuffix(' K')
        form.addRow('温度', self.temperature)
        self.stage_button = QtWidgets.QPushButton('编辑四阶段 MDP…')
        self.stage_status = QtWidgets.QLabel()
        self.stage_status.setWordWrap(True)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.stage_button)
        row.addWidget(self.stage_status, 1)
        form.addRow('阶段参数', row)
        self._update_stage_status()

        self.summary = QtWidgets.QPlainTextEdit()
        self.summary.setReadOnly(True)
        root.addWidget(self.summary, 1)
        self.status = QtWidgets.QLabel()
        self.status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.status.setWordWrap(True)
        root.addWidget(self.status)

        row = QtWidgets.QHBoxLayout()
        self.import_parameters = QtWidgets.QPushButton('管理 PTM 参数…')
        self.peptide = QtWidgets.QPushButton('导出位点短肽…')
        row.addWidget(self.import_parameters)
        row.addWidget(self.peptide)
        row.addStretch(1)
        self.export_button = QtWidgets.QPushButton('生成模拟包…')
        close = QtWidgets.QPushButton('关闭')
        row.addWidget(self.export_button)
        row.addWidget(close)
        root.addLayout(row)

        self.objects.itemSelectionChanged.connect(self._change_object)
        self.state.valueChanged.connect(self._inspect)
        self.library.editingFinished.connect(self._inspect)
        self.runtime.editingFinished.connect(self._inspect)
        refresh.clicked.connect(self.refresh)
        self.stage_button.clicked.connect(self._edit_stages)
        self.import_parameters.clicked.connect(self._register_parameters)
        self.charmm.clicked.connect(self._import_charmm)
        site.clicked.connect(lambda: QtGui.QDesktopServices.openUrl(
            QtCore.QUrl(md_charmm_gui.SOLUTION_BUILDER_URL)))
        self.peptide.clicked.connect(self._export_validation_peptide)
        self.export_button.clicked.connect(self._export)
        close.clicked.connect(self.close)
        self.refresh()

    def _folder_row(self, form, title, value, hint):
        field = QtWidgets.QLineEdit(str(value))
        field.setPlaceholderText(hint)
        browse = QtWidgets.QPushButton('浏览…')
        row = QtWidgets.QHBoxLayout()
        row.addWidget(field, 1)
        row.addWidget(browse)

        def choose():
            path = QtWidgets.QFileDialog.getExistingDirectory(self, title, field.text())
            if path:
                field.setText(path)
                self._inspect()

        browse.clicked.connect(choose)
        form.addRow(title, row)
        return field

    def refresh(self):
        selected = set(self._selected_objects())
        with QtCore.QSignalBlocker(self.objects):
            self.objects.clear()
            self.objects.addItems([name for name in self.cmd.get_names('public_objects')
                                   if self.cmd.get_type(name) == 'object:molecule'])
            for row in range(self.objects.count()):
                item = self.objects.item(row)
                item.setSelected(item.text() in selected)
            if self.objects.count() and not self.objects.selectedItems():
                self.objects.setCurrentRow(0)
        self._change_object()

    def _selected_objects(self):
        return [self.objects.item(row).text() for row in range(self.objects.count())
                if self.objects.item(row).isSelected()]

    def _change_object(self):
        names = self._selected_objects()
        dynamic = [(name, count) for name in names if (count := self.cmd.count_states('%' + name)) > 1]
        try:
            state = self.cmd.get_object_state(dynamic[0][0]) if dynamic else 1
        except CmdException:
            state = 1
        with QtCore.QSignalBlocker(self.state):
            self.state.setMaximum(min((count for _, count in dynamic), default=1))
            self.state.setValue(state)
        self._inspect()

    def _source(self, explicit_bonds=False):
        names = self._selected_objects()
        if not names or not set(names) <= set(self.cmd.get_names('public_objects')):
            raise ValueError('请选择至少一个蛋白分子对象，或刷新已变化的对象列表。')

        def files():
            for name in names:
                selection = '%' + name
                state = 1 if self.cmd.count_states(selection) == 1 else self.state.value()
                if any(len(chain) > 1 for chain in self.cmd.get_chains(selection, state=state)):
                    raise ValueError(f'{name} 包含多字符链标识，请先将链转换为单字符 PDB 标识。')
                yield name, state, self.cmd.get_pdbstr(selection, state=state)

        if not explicit_bonds:
            return md_setup.combine_protein_sources(files())
        with self.cmd.lockcm:
            previous = self.cmd.get_setting_boolean('pdb_conect_all')
            try:
                self.cmd.set('pdb_conect_all', 1, updates=0)
                return md_setup.combine_protein_sources(files())
            finally:
                self.cmd.set('pdb_conect_all', previous, updates=0)

    def _update_stage_status(self):
        labels = [label for stage, label in STAGE_LABELS.items() if stage in self.mdp_overrides]
        status = '已自定义 ' + '、'.join(labels) + '。' if labels else '所有阶段使用默认参数。'
        self.stage_status.setText(status + ' 表单值只用于默认阶段。')

    def _edit_stages(self):
        if self.is_running():
            return
        defaults = stage_inputs(self.temperature.value(), round(self.production.value() * 500000))
        dialog = StageMDPDialog(self, defaults, self.mdp_overrides)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.mdp_overrides = dialog.overrides
            self._update_stage_status()

    def _runtime(self):
        value = self.runtime.text().strip()
        if not value:
            raise ValueError('请选择 AmberTools／ParmEd／GROMACS 工具环境。')
        return md_setup.Runtime(Path(value).expanduser().absolute())

    def _inspect(self):
        if self.is_running():
            return
        self.export_button.setEnabled(False)
        self.peptide.setEnabled(False)
        self.lactyl_sites = []
        try:
            pdb, origins = self._source(explicit_bonds=True)
            residues, waters = md_setup.inspect_pdb(pdb, conect_source=md_setup.OBJECT_CONECT)
            sites = [key for key in residues if key[2] in md_parameters.PTM_CODES]
            self.lactyl_sites = [key for key in sites if key[2] in md_parameters.LACTYL]
            self.peptide.setEnabled(bool(self.lactyl_sites))
            chains = dict.fromkeys(key[0] for key in residues)
            summary = [f'{len(self._selected_objects())} 个对象 · {len(chains)} 条蛋白链 · {len(residues)} 个残基']
            for chain in chains:
                origin = origins[chain]
                summary.append(f'{chain.strip() or "无链"} ← {origin["object"]}／'
                               f'{origin["source_chain"] or "无链"} · 状态 {origin["state"]}')
            if waters:
                summary.append(f'重新溶剂化时会替换 {waters} 个原有水分子。')
            library = Path(self.library.text()).expanduser()
            ready = True
            base_hashes = None
            archive_error = None
            if any(code in md_parameters.PHOSPHATE for _, _, code in sites):
                try:
                    runtime = self._runtime()
                    path = library / 'phosaa14SB'
                    if path.exists():
                        md_parameters.archived_phosphate(runtime.leap_data, library)
                    elif not all((runtime.leap_data / name).is_file() for name in md_parameters.PHOS_FILES):
                        raise ValueError('工具环境缺少 phosaa14SB 参数文件。')
                except (OSError, ValueError) as exc:
                    archive_error = str(exc)
                    ready = False
            for chain, number, code in sites:
                address = f'{chain.strip() or "无链"}/{number}  {code}'
                if code in md_parameters.PHOSPHATE:
                    state = md_parameters.PHOSPHATE[code]
                    summary.append(f'{address}  ·  phosaa14SB（{state}）')
                else:
                    try:
                        _, record = md_parameters.lactyl_template(library, code)
                        if base_hashes is None:
                            base_hashes = md_parameters.checksums(
                                self._runtime().leap_data, md_parameters.BASE_FILES)
                        md_parameters.check_lactyl_force_field(record, base_hashes)
                    except (ValueError, OSError) as exc:
                        summary.append(f'{address}  ·  {exc}')
                        summary.append(f'参数目录 {library / code}')
                        ready = False
                    else:
                        summary.append(f'{address}  ·  已找到 {md_parameters.LACTYL[code]} 模板')
            if archive_error:
                summary.append(archive_error)
            if not sites:
                summary.append('未发现 PTM 残基。')
            self.summary.setPlainText('\n'.join(summary))
            self.export_button.setEnabled(ready)
            self.status.setText('生成后解压，在 GPU 机器运行 bash run_md.sh。' if ready
                                else '可导入已核验的 ff14SB 链内参数，或导出位点短肽供参数化验证。')
        except (CmdException, ValueError) as exc:
            self.summary.setPlainText(str(exc))
            self.status.clear()

    def _choose_lactyl_site(self):
        if self.is_running() or not self.lactyl_sites:
            return None
        if len(self.lactyl_sites) == 1:
            return self.lactyl_sites[0]
        labels = [f'{chain.strip() or "无链"}/{number} {code}'
                  for chain, number, code in self.lactyl_sites]
        label, ok = QtWidgets.QInputDialog.getItem(self, '选择乳酸化位点', '位点', labels, 0, False)
        return self.lactyl_sites[labels.index(label)] if ok else None

    def _register_parameters(self):
        if self.is_running():
            return
        library = Path(self.library.text()).expanduser()
        try:
            runtime = self._runtime()
            base_hashes = md_parameters.checksums(runtime.leap_data, md_parameters.BASE_FILES)
        except (OSError, ValueError) as exc:
            self.status.setText(str(exc))
            return
        phosphate = 'phosaa14SB'
        types = [phosphate, *md_parameters.LACTYL]
        labels = []
        for code in types:
            try:
                if code == phosphate:
                    if not (library / code).exists():
                        raise FileNotFoundError
                    md_parameters.archived_phosphate(runtime.leap_data, library)
                else:
                    _, record = md_parameters.lactyl_template(library, code)
                    md_parameters.check_lactyl_force_field(record, base_hashes)
                status = '已归档'
            except (OSError, ValueError):
                status = '待归档或需核对'
            detail = 'SEP / TPO / PTR / S1P / T1P / Y1P' if code == phosphate else md_parameters.LACTYL[code]
            labels.append(f'{code} · {detail} · {status}')
        label, ok = QtWidgets.QInputDialog.getItem(
            self, '管理 PTM 参数', '选择要核对或导入的参数', labels, 0, False)
        if not ok:
            return
        code = types[labels.index(label)]
        try:
            if code == phosphate:
                target, _ = md_parameters.archived_phosphate(runtime.leap_data, library)
            else:
                source = QtWidgets.QFileDialog.getExistingDirectory(
                    self, f'选择含 {code} 的链内参数目录', str(library))
                if not source:
                    return
                target = md_parameters.register_lactyl_template(
                    Path(source), library, code, runtime.leap_data)
        except (ValueError, OSError) as exc:
            self.status.setText(str(exc))
            return
        self._inspect()
        self.status.setText(f'已登记 {target}。生成模拟包时会核对肽键、构型和电荷。')

    def _export_validation_peptide(self):
        site = self._choose_lactyl_site()
        if site is None:
            return
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, '导出位点及两侧相邻残基', f'{self._selected_objects()[0]}_{site[0].strip()}{site[1]}_peptide.pdb',
            'PDB (*.pdb)')
        if not filename:
            return
        try:
            source, _ = self._source(explicit_bonds=True)
            pdb = md_setup.lactyl_validation_peptide(source, site)
            path = Path(filename).with_suffix('.pdb')
            with path.open('x', encoding='utf-8') as output:
                output.write(pdb)
            self.status.setText(f'已导出短肽 {path}。这是待验证的结构输入。')
        except (CmdException, ValueError, OSError) as exc:
            self.status.setText(str(exc))

    def _export(self):
        if self.is_running():
            return
        try:
            pdb, origins = self._source(explicit_bonds=True)
            runtime = self._runtime()
        except (CmdException, ValueError) as exc:
            self.status.setText(str(exc))
            return
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, '保存 GROMACS 模拟包', self._selected_objects()[0] + '_md.zip', 'ZIP (*.zip)')
        if not filename:
            return
        path = Path(filename).with_suffix('.zip')
        library = Path(self.library.text()).expanduser()
        self._start_build(md_setup.build_package, path, pdb, path, runtime, library=library,
                          production_ns=self.production.value(), temperature=self.temperature.value(),
                          padding=12, mdp_overrides=dict(self.mdp_overrides), chain_origins=origins)

    def _import_charmm(self):
        if self.is_running():
            return
        try:
            runtime = self._runtime()
        except ValueError as exc:
            self.status.setText(str(exc))
            return
        source, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, '选择 CHARMM-GUI Solution Builder 的 GROMACS 作业归档',
            '', 'CHARMM-GUI (*.tgz *.tar.gz *.tar)')
        if not source:
            return
        source = Path(source)
        label = source.name.removesuffix('.gz').removesuffix('.tar').removesuffix('.tgz')
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, '保存 CHARMM-GUI 模拟包', f'{label}_md.zip', 'ZIP (*.zip)')
        if filename:
            path = Path(filename).with_suffix('.zip')
            self._start_build(md_charmm_gui.build_package, path, source, path, runtime)

    def _start_build(self, builder, output, *args, **kwargs):
        self.export_button.setEnabled(False)
        self.stage_button.setEnabled(False)
        self.import_parameters.setEnabled(False)
        self.charmm.setEnabled(False)
        status = ('正在读取 CHARMM-GUI 原始拓扑并预检各阶段…' if builder is md_charmm_gui.build_package
                  else '正在生成溶剂化体系、转换拓扑并运行 GROMACS 预检…')
        self.status.setText(status)
        self.output = output
        self.worker = AsyncFunc(builder, finishslot=self._finished)
        self.worker(*args, **kwargs)

    def is_running(self):
        return self.worker is not None and self.worker.isRunning()

    def _finished(self, payload):
        manifest, error = payload
        self.worker.wait()
        self.worker = None
        self.stage_button.setEnabled(True)
        self.import_parameters.setEnabled(True)
        self.charmm.setEnabled(True)
        if getattr(self.parent(), '_quit_requested', False):
            return
        self._inspect()
        if error is None:
            self.status.setText(f'模拟包已生成并通过预检\n{self.output}')
            if 'origin' in manifest:
                stages = '\n'.join(manifest['stages'].values())
                self.summary.setPlainText(f'CHARMM-GUI · {manifest["archive"]}\n'
                                          f'使用归档内的坐标、拓扑、力场和阶段 MDP。\n{stages}')
        else:
            self.status.setText('生成失败。')
            self.summary.setPlainText(str(error))

    def done(self, result):
        if not self.is_running():
            super().done(result)
