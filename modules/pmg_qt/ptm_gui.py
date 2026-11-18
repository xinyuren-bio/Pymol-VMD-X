"""Qt dialog for selecting, previewing and exporting PTM sites."""

from pathlib import Path

from pymol import CmdException, ptm
from pymol.Qt import QtCore, QtWidgets

Qt = QtCore.Qt


class PTMDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.cmd = parent.cmd
        self.preview = None
        self.result = None
        self.sites = []
        self.setWindowTitle('蛋白质翻译后修饰（PTM）')
        self.resize(690, 500)

        root = QtWidgets.QVBoxLayout(self)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('对象'))
        self.objects = QtWidgets.QComboBox()
        row.addWidget(self.objects, 1)
        row.addWidget(QtWidgets.QLabel('状态'))
        self.state = QtWidgets.QSpinBox()
        self.state.setMinimum(1)
        row.addWidget(self.state)
        refresh = QtWidgets.QPushButton('刷新')
        row.addWidget(refresh)
        root.addLayout(row)

        filters = QtWidgets.QHBoxLayout()
        self.chains = QtWidgets.QComboBox()
        filters.addWidget(self.chains)
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText('搜索残基编号、名称或链')
        self.search.setClearButtonEnabled(True)
        filters.addWidget(self.search, 1)
        self.imported_only = QtWidgets.QCheckBox('仅已有修饰')
        filters.addWidget(self.imported_only)
        root.addLayout(filters)

        self.table = QtWidgets.QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(('链 / 残基编号', '残基', '修饰'))
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        root.addWidget(self.table, 1)

        settings = QtWidgets.QLabel('结构模型 · 无 MD 拓扑和力场参数')
        settings.setToolTip('以 wwPDB CCD 理想重原子坐标建模。磷酸基按 −2 形式电荷显示。')
        root.addWidget(settings)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        root.addWidget(self.status)

        buttons = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton('选中位点')
        preview = QtWidgets.QPushButton('预览侧链')
        apply = QtWidgets.QPushButton('应用为新对象')
        self.restore_button = QtWidgets.QPushButton('恢复原氨基酸')
        self.restore_button.setEnabled(False)
        self.restore_button.setToolTip('将选中的已有修饰位点恢复成标准残基，生成独立对象。')
        self.export_button = QtWidgets.QPushButton('导出 PDB…')
        close = QtWidgets.QPushButton('关闭')
        for button in (choose, preview, apply, self.restore_button, self.export_button):
            buttons.addWidget(button)
        buttons.addStretch(1)
        buttons.addWidget(close)
        root.addLayout(buttons)
        self.export_button.setEnabled(False)

        refresh.clicked.connect(self.refresh)
        self.objects.currentIndexChanged.connect(self._change_object)
        self.state.valueChanged.connect(self._load_sites)
        self.chains.currentIndexChanged.connect(self._filter)
        self.search.textChanged.connect(self._filter)
        self.imported_only.toggled.connect(self._filter)
        choose.clicked.connect(self._select_source)
        preview.clicked.connect(lambda: self._build(preview=True))
        apply.clicked.connect(lambda: self._build(preview=False))
        self.restore_button.clicked.connect(lambda: self._build(preview=False, restore=True))
        self.table.itemSelectionChanged.connect(self._update_restore)
        self.export_button.clicked.connect(self._export)
        close.clicked.connect(self.close)
        self.refresh()

    def refresh(self):
        current = self.objects.currentText()
        self.objects.blockSignals(True)
        self.objects.clear()
        self.objects.addItems([name for name in self.cmd.get_names('public_objects')
                               if self.cmd.get_type(name) == 'object:molecule'])
        if self.objects.findText(current) >= 0:
            self.objects.setCurrentText(current)
        self.objects.blockSignals(False)
        self._change_object()

    def _change_object(self):
        name = self.objects.currentText()
        self.state.blockSignals(True)
        count = max(1, self.cmd.count_states('%' + name)) if name else 1
        try:
            current = self.cmd.get_object_state(name) if name else 1
        except CmdException:
            current = 1
        self.state.setMaximum(count)
        self.state.setValue(current)
        self.state.blockSignals(False)
        self._load_sites()

    def _load_sites(self):
        self._clear_preview()
        self.sites = []
        name = self.objects.currentText()
        if name:
            try:
                model = ptm.object_model(self.cmd, name, self.state.value())
                self.sites = ptm.list_sites(model)
                self.status.setText(f'{len(self.sites)} 个位点。可用 Shift 或 Command／Ctrl 多选。'
                                    if self.sites else '当前状态没有可用的 PTM 位点。')
            except (CmdException, ValueError) as exc:
                self.status.setText(str(exc))
        else:
            self.status.setText('请先打开蛋白结构。')

        self.chains.blockSignals(True)
        self.chains.clear()
        self.chains.addItem('全部链', None)
        for chain in sorted({site.chain for site in self.sites}):
            self.chains.addItem('链 ' + chain if chain else '无链', chain)
        self.chains.blockSignals(False)

        self.table.setRowCount(0)
        self.table.setRowCount(len(self.sites))
        for row, site in enumerate(self.sites):
            address = f'{site.chain or "无链"} / {site.resi}'
            if site.segment:
                address += f' [{site.segment}]'
            self.table.setItem(row, 0, QtWidgets.QTableWidgetItem(address))
            self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(site.resn))
            choices = site.choices()
            if len(choices) == 1:
                text = ptm.TEMPLATES[choices[0]].label
                if site.resn == choices[0]:
                    text += ' · 已有'
                self.table.setItem(row, 2, QtWidgets.QTableWidgetItem(text))
            else:
                combo = QtWidgets.QComboBox()
                for code in choices:
                    combo.addItem(ptm.TEMPLATES[code].label, code)
                combo.activated.connect(lambda _, r=row: self._add_row_selection(r))
                self.table.setCellWidget(row, 2, combo)
        self._filter()

    def _add_row_selection(self, row):
        selection = self.table.selectionModel()
        flags = QtCore.QItemSelectionModel.SelectionFlag
        selection.select(self.table.model().index(row, 0), flags.Select | flags.Rows)

    def _filter(self):
        chain = self.chains.currentData()
        query = self.search.text().strip().casefold()
        self.table.clearSelection()
        for row, site in enumerate(self.sites):
            self.table.setRowHidden(row,
                (chain is not None and site.chain != chain)
                or (self.imported_only.isChecked() and site.resn not in ptm.TEMPLATES)
                or (bool(query) and query not in site.label().casefold()))

    def _requests(self):
        for index in sorted(self.table.selectionModel().selectedRows(), key=lambda i: i.row()):
            row = index.row()
            site = self.sites[row]
            combo = self.table.cellWidget(row, 2)
            yield site, combo.currentData() if combo else site.choices()[0]

    def _update_restore(self):
        self.restore_button.setEnabled(any(site.resn in ptm.TEMPLATES for site, _ in self._requests()))

    def _select_source(self):
        try:
            selection = self.cmd.get_unused_name('ptm_sites', 0)
            ptm.select_sites(self.cmd, self.objects.currentText(),
                             [site for site, _ in self._requests()], selection,
                             state=self.state.value())
            self.status.setText(f'已创建源对象选择「{selection}」。')
        except (CmdException, ValueError) as exc:
            self.status.setText(str(exc))

    def _build(self, *, preview, restore=False):
        from .publication_gui import publication_export_running
        if publication_export_running():
            self.status.setText('正在导出图片，请等导出完成。')
            return
        source = self.objects.currentText()
        try:
            model = ptm.object_model(self.cmd, source, self.state.value())
            if restore:
                built, sites = ptm.restore_model(
                    model, (site for site, _ in self._requests() if site.resn in ptm.TEMPLATES))
                warnings = ()
                prefix = source + '_restored'
            else:
                built, sites, warnings = ptm.build_model(model, self._requests())
                prefix = '_ptm_preview' if preview else source + '_PTM'
            name = self.cmd.get_unused_name(prefix, 0)
            ptm.load_result(self.cmd, name, built, sites, preview)
            self._clear_preview()
            if preview:
                self.preview = name
                message = f'预览 {len(sites)} 个位点。'
            else:
                self.result = name
                self.export_button.setEnabled(True)
                if restore:
                    message = f'已在独立对象「{name}」中恢复 {len(sites)} 个位点。'
                else:
                    message = f'已生成独立对象「{name}」，包含 {len(sites)} 个位点。'
            self.status.setText(' '.join((message, *warnings)))
        except (CmdException, ValueError) as exc:
            self.status.setText(str(exc))

    def _export(self):
        if not self.result:
            return
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, '导出结构模型', self.result + '.pdb', 'PDB 结构 (*.pdb)')
        if not filename:
            return
        path = Path(filename)
        if path.suffix.lower() != '.pdb':
            path = path.with_suffix('.pdb')
        try:
            ptm.export_pdb(self.cmd, self.result, path)
            self.status.setText(f'已导出 {path.name}，包含显式连接记录。')
        except (CmdException, ValueError, OSError) as exc:
            self.status.setText(str(exc))

    def _clear_preview(self):
        if self.preview and getattr(self.cmd, '_COb', None) is not None:
            self.cmd.delete(self.preview)
        self.preview = None

    def done(self, result):
        self._clear_preview()
        super().done(result)
