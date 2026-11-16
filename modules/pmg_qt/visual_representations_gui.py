# ==================================================
# 功能说明：中文无模态 VMD 风格“选择与表示”窗口（多层独立显示）。
# 使用方法：open_visual_representations(parent, target='all')；退出清理 clear_visual_representations_ref()。
# 依赖环境：PyMOL + PySide6（pymol.Qt）；后端 pymol.visual_representations。
# 生成时间：2026-11-14 16:50
# ==================================================

"""中文图形表示对话框（类似 VMD Graphical Representations）。"""

from __future__ import annotations

from typing import Any, List, Optional, Tuple
from collections import Counter

from pymol.Qt import QtCore, QtGui, QtWidgets

from pymol import visual_representations as vr
from pmg_qt.color_palette import make_swatch_icon, open_color_palette_popup, rgb_to_hex

Qt = QtCore.Qt

_dialog_ref: Optional["VisualRepresentationsDialog"] = None

# 简明例子 chips
_EXAMPLE_CHIPS: Tuple[str, ...] = (
    "resname AL31",
    "resname DPPC",
    "nucleic",
    "protein",
    "name BB*",
    "not water",
    "resid 1 to 10",
)


def _parent_alive(parent: Optional[QtWidgets.QWidget]) -> bool:
    """判断父窗口是否仍可用于引擎访问。"""
    if parent is None:
        return False
    try:
        from pymol.Qt import sip
    except Exception:
        try:
            import shiboken6 as sip  # type: ignore
        except Exception:
            sip = None  # type: ignore
    if sip is not None:
        try:
            if sip.isdeleted(parent):
                return False
        except Exception:
            pass
    try:
        parent.objectName()
        return True
    except RuntimeError:
        return False


class LayerTable(QtWidgets.QTableWidget):
    """三列规则表：显示方式 / 颜色 / 选择表达式；行可开关。"""

    def __init__(self, parent: Optional[QtWidgets.QWidget] = None) -> None:
        super().__init__(0, 3, parent)
        self.setHorizontalHeaderLabels(["显示方式", "颜色", "选择表达式"])
        self.horizontalHeader().setStretchLastSection(True)
        self.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.verticalHeader().setVisible(True)


class VisualRepresentationsDialog(QtWidgets.QDialog):
    """
    无模态选择与表示窗口。

    对象下拉含 all / 可见分子 / 禁用源 / CG 连接副本，不含本窗口自有 visrep_* 副本。
    默认不自动应用；可选自动应用仅在 editingFinished / 控件更改时触发。
    """

    def __init__(self, parent: QtWidgets.QWidget, cmd: Any) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("选择与表示（V）· Graphical Representations")
        self.setMinimumSize(640, 480)
        self.cmd = cmd
        self.ctrl = vr.RepresentationController(cmd)
        self._busy = False
        self._closed = False
        self._loading = False
        self._editing_layer_id = None
        self._color_popups: List[Any] = []
        self._build_ui()
        self.refresh_objects()
        self._load_layers_to_table()

    def _build_ui(self) -> None:
        """构建中文界面。"""
        root = QtWidgets.QVBoxLayout(self)
        root.setSpacing(8)

        tip = QtWidgets.QLabel(
            "用 VMD 常用选择语法，分别设置磷脂、核酸和蛋白的显示与颜色。"
            "每行可独立显示；点击“应用”查看，“恢复源对象”返回原来的画面。"
        )
        tip.setWordWrap(True)
        tip.setStyleSheet("color:#444;font-size:11px;")
        root.addWidget(tip)

        # 对象
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("对象："))
        self.cmb_obj = QtWidgets.QComboBox()
        self.cmb_obj.setMinimumWidth(200)
        row.addWidget(self.cmb_obj, 1)
        self.btn_refresh = QtWidgets.QPushButton("刷新")
        row.addWidget(self.btn_refresh)
        root.addLayout(row)

        # 规则表
        self.table = LayerTable(self)
        root.addWidget(self.table, 1)

        # 层操作
        brow = QtWidgets.QHBoxLayout()
        self.btn_new = QtWidgets.QPushButton("新建表示")
        self.btn_del = QtWidgets.QPushButton("删除表示")
        self.btn_components = QtWidgets.QPushButton("按组分分层")
        self.btn_components.setToolTip(
            "按脂质残基类型、核酸整体、蛋白整体创建不同颜色规则（不含水分离子）。"
        )
        brow.addWidget(self.btn_new)
        brow.addWidget(self.btn_del)
        brow.addWidget(self.btn_components)
        brow.addStretch(1)
        root.addLayout(brow)

        # 编辑区
        edit = QtWidgets.QGroupBox("当前表示编辑")
        form = QtWidgets.QFormLayout(edit)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

        self.chk_enabled = QtWidgets.QCheckBox("启用此行")
        self.chk_enabled.setChecked(True)
        form.addRow("", self.chk_enabled)

        self.ed_sele = QtWidgets.QLineEdit("all")
        self.ed_sele.setPlaceholderText("例如 resid < 90 and resid > 80")
        self.ed_sele.setToolTip("支持 resid < 90、resid 81 to 89、resname、name、chain，及 and / or / not 和括号。")
        form.addRow("选择表达式：", self.ed_sele)

        resnames = QtWidgets.QWidget()
        res_layout = QtWidgets.QVBoxLayout(resnames)
        res_layout.setContentsMargins(0, 0, 0, 0)
        res_layout.setSpacing(4)
        search_row = QtWidgets.QHBoxLayout()
        self.ed_resname_search = QtWidgets.QLineEdit()
        self.ed_resname_search.setPlaceholderText("搜索当前对象的 resname…")
        self.lbl_resname_count = QtWidgets.QLabel("0 种")
        search_row.addWidget(self.ed_resname_search, 1)
        search_row.addWidget(self.lbl_resname_count)
        res_layout.addLayout(search_row)
        self.list_resnames = QtWidgets.QListWidget()
        self.list_resnames.setObjectName("resname_shortcuts")
        self.list_resnames.setViewMode(QtWidgets.QListView.ViewMode.IconMode)
        self.list_resnames.setMovement(QtWidgets.QListView.Movement.Static)
        self.list_resnames.setResizeMode(QtWidgets.QListView.ResizeMode.Adjust)
        self.list_resnames.setWrapping(True)
        self.list_resnames.setSpacing(4)
        self.list_resnames.setMinimumHeight(70)
        self.list_resnames.setMaximumHeight(110)
        self.list_resnames.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.list_resnames.setStyleSheet(
            "QListWidget { border: 1px solid #c6cbd1; border-radius: 5px; }"
            "QListWidget::item { padding: 5px 9px; background: #edf2f7; border-radius: 4px; }"
            "QListWidget::item:selected { background: #c6e3ff; color: #153d63; }")
        res_layout.addWidget(self.list_resnames)
        form.addRow("resname 快选：", resnames)

        # 例子 chips
        chip_row = QtWidgets.QHBoxLayout()
        chip_row.addWidget(QtWidgets.QLabel("例子："))
        for ex in _EXAMPLE_CHIPS:
            b = QtWidgets.QToolButton()
            b.setText(ex)
            b.setToolTip("填入：" + ex)
            b.clicked.connect(lambda _=False, t=ex: self._insert_example(t))
            chip_row.addWidget(b)
        chip_row.addStretch(1)
        form.addRow(chip_row)

        prev_row = QtWidgets.QHBoxLayout()
        self.btn_preview = QtWidgets.QPushButton("预览选择")
        self.lbl_preview = QtWidgets.QLabel("匹配：—")
        prev_row.addWidget(self.btn_preview)
        prev_row.addWidget(self.lbl_preview, 1)
        form.addRow(prev_row)

        self.cmb_style = QtWidgets.QComboBox()
        for meta in vr.list_styles():
            self.cmb_style.addItem(meta["label"], meta["id"])
        form.addRow("显示方式：", self.cmb_style)
        self.lbl_style_help = QtWidgets.QLabel()
        self.lbl_style_help.setWordWrap(True)
        form.addRow("", self.lbl_style_help)
        self._update_style_help()

        self.cmb_color_mode = QtWidgets.QComboBox()
        for meta in vr.list_color_modes():
            self.cmb_color_mode.addItem(meta["label"], meta["id"])
        form.addRow("颜色模式：", self.cmb_color_mode)

        crow = QtWidgets.QHBoxLayout()
        self.btn_color = QtWidgets.QToolButton()
        self.btn_color.setToolTip("打开色板")
        self._set_color_button((0.8, 0.2, 0.2))
        crow.addWidget(self.btn_color)
        crow.addWidget(QtWidgets.QLabel("单色色块（颜色模式=单色时生效）"))
        crow.addStretch(1)
        form.addRow("颜色：", crow)

        root.addWidget(edit)

        self.chk_auto = QtWidgets.QCheckBox("自动应用")
        self.chk_auto.setChecked(False)
        root.addWidget(self.chk_auto)

        self.status = QtWidgets.QLabel("就绪。打开窗口不会自动改场景。")
        self.status.setWordWrap(True)
        root.addWidget(self.status)

        # 底部按钮
        foot = QtWidgets.QHBoxLayout()
        self.btn_apply = QtWidgets.QPushButton("应用")
        self.btn_restore = QtWidgets.QPushButton("恢复源对象")
        self.btn_close = QtWidgets.QPushButton("关闭")
        foot.addWidget(self.btn_apply)
        foot.addWidget(self.btn_restore)
        foot.addStretch(1)
        foot.addWidget(self.btn_close)
        root.addLayout(foot)

        # 信号
        self.btn_refresh.clicked.connect(self.refresh_objects)
        self.btn_new.clicked.connect(self._new_layer)
        self.btn_del.clicked.connect(self._delete_layer)
        self.btn_components.clicked.connect(self._layers_by_component)
        self.btn_preview.clicked.connect(self._preview)
        self.btn_color.clicked.connect(self._pick_color)
        self.btn_apply.clicked.connect(self.apply_layers)
        self.btn_restore.clicked.connect(self.restore_source)
        self.btn_close.clicked.connect(self.close)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
        self.cmb_obj.currentIndexChanged.connect(self._on_source_changed)
        self.chk_enabled.toggled.connect(self._on_edit_changed)
        self.cmb_style.currentIndexChanged.connect(self._on_edit_changed)
        self.cmb_color_mode.currentIndexChanged.connect(self._on_edit_changed)
        self.ed_sele.editingFinished.connect(self._on_sele_finished)
        self.ed_resname_search.textChanged.connect(self._filter_resnames)
        self.list_resnames.itemClicked.connect(self._choose_resname)

    def _engine_ok(self) -> bool:
        """退出/导出 busy 保护。"""
        if self._closed or self._busy:
            return False
        parent = self.parent()
        if getattr(parent, "_quit_requested", False) or getattr(parent, "_shutting_down", False):
            return False
        if getattr(self.cmd, "_COb", None) is None:
            return False
        try:
            from pmg_qt.publication_gui import publication_export_running

            if publication_export_running():
                self.status.setText("正在导出图片，请完成后再调整表示。")
                return False
        except Exception:
            pass
        if not _parent_alive(parent):
            return False
        try:
            self.cmd.get_names("public_objects", enabled_only=0)
            return True
        except Exception:
            self.status.setText("PyMOL 引擎不可用或已退出。")
            return False

    def _set_color_button(self, rgb: Tuple[float, float, float]) -> None:
        """更新色块按钮外观与 tooltip。"""
        self._current_rgb = (
            float(rgb[0]),
            float(rgb[1]),
            float(rgb[2]),
        )
        self.btn_color.setIcon(make_swatch_icon(self._current_rgb))
        self.btn_color.setIconSize(QtCore.QSize(28, 18))
        self.btn_color.setText("")
        self.btn_color.setToolTip(
            "RGB (%.2f, %.2f, %.2f) / %s"
            % (self._current_rgb[0], self._current_rgb[1], self._current_rgb[2],
               rgb_to_hex(self._current_rgb))
        )

    def _related_for_source(self, obj: str) -> List[str]:
        """
        若源是粗粒化窗口拥有的 mol，收集对应 guide 以便暂隐。

        不导入 coarse_grained_gui；仅按命名约定推断，避免影响其它用户 CGO。
        """
        related: List[str] = []
        if not obj:
            return related
        import sys
        cg_gui = sys.modules.get("pmg_qt.coarse_grained_gui")
        dialog = getattr(cg_gui, "_dialog_ref", None)
        ctrl = getattr(dialog, "ctrl", None)
        if ctrl is not None and obj in ctrl.owned:
            for name in ctrl.owned:
                if name != obj and name in self.cmd.get_names("all"):
                    if self.cmd.get_type(name) == "object:cgo":
                        related.append(name)
        return related

    def refresh_objects(self) -> None:
        """刷新对象下拉。"""
        if not self._engine_ok():
            return
        cur = self.cmb_obj.currentData()
        self.cmb_obj.blockSignals(True)
        self.cmb_obj.clear()
        self.cmb_obj.addItem("all（当前可见分子）", "all")
        for n in vr.list_candidate_objects(self.cmd, exclude_owned=True):
            label = n
            try:
                enabled = n in (self.cmd.get_names("public_objects", enabled_only=1) or [])
                if not enabled:
                    label = "%s（已禁用）" % n
                if n.startswith("cgvis_"):
                    label = "%s（CG 连接副本）" % n
            except Exception:
                pass
            self.cmb_obj.addItem(label, n)
        self.cmb_obj.blockSignals(False)
        if cur:
            idx = self.cmb_obj.findData(cur)
            if idx >= 0:
                self.cmb_obj.setCurrentIndex(idx)
        self._refresh_resnames()

    def _refresh_resnames(self) -> None:
        """Read actual names within the current source; never retain stale entries."""
        self.list_resnames.clear()
        self.ed_resname_search.clear()
        self.lbl_resname_count.setText("0 种")
        self.list_resnames.setToolTip("")
        if not self._engine_ok():
            return
        try:
            source = self._current_source()
            counts = Counter()
            if source in (self.cmd.get_names("selections") or []):
                scopes = [vr._quote_sel(source)]
            else:
                scopes = [vr._quote_sel(obj) for obj in self.ctrl._resolve_sources(source)]
            for scope in scopes:
                self.cmd.iterate(scope, "counts[resn] += 1", space={"counts": counts})
            for name, count in sorted(counts.items(), key=lambda pair: pair[0].casefold()):
                label = name if name else "（空残基名）"
                item = QtWidgets.QListWidgetItem("%s · %s" % (label, format(count, ",")))
                item.setData(Qt.ItemDataRole.UserRole, name)
                item.setToolTip("resname：%s\n匹配 %s 个原子/珠" % (label, format(count, ",")))
                self.list_resnames.addItem(item)
            self.lbl_resname_count.setText("%d 种" % len(counts))
            if not counts:
                self.list_resnames.setToolTip("当前对象没有可用的原子或珠。")
        except Exception as exc:
            self.list_resnames.clear()
            self.list_resnames.setToolTip("读取 resname 失败：%s" % exc)

    def _filter_resnames(self, text) -> None:
        query = text.casefold().strip()
        visible = 0
        for index in range(self.list_resnames.count()):
            item = self.list_resnames.item(index)
            name = str(item.data(Qt.ItemDataRole.UserRole) or "")
            hidden = query not in name.casefold()
            item.setHidden(hidden)
            visible += not hidden
        total = self.list_resnames.count()
        self.lbl_resname_count.setText("%d / %d 种" % (visible, total) if query else "%d 种" % total)

    def _choose_resname(self, item) -> None:
        if not self._engine_ok():
            return
        name = str(item.data(Qt.ItemDataRole.UserRole) or "")
        # Quote names which could be parsed as keywords or contain whitespace.
        value = name if __import__('re').fullmatch(r"[A-Za-z0-9_+.*?\-]+", name) else '"' + name + '"'
        self.ed_sele.setText("resname " + value)
        self._preview()
        if self.chk_auto.isChecked():
            self.apply_layers()

    def select_target(self, target: str) -> None:
        """从对象行 V / 菜单传入 target 并选中。"""
        t = (target or "all").strip() or "all"
        owned_row = next((i for i, layer in enumerate(self.ctrl.layers)
                          if layer.owned_object == t or
                          layer.layer_id == self.ctrl.object_layers.get(t)), None)
        if owned_row is not None and self.ctrl.source:
            t = self.ctrl.source
        self.refresh_objects()
        # The quick V entry chooses the only available molecule automatically.
        if t.lower() == "all":
            candidates = self.ctrl._resolve_sources("all")
            if len(candidates) == 1:
                t = candidates[0]
            elif not candidates and self.ctrl.source:
                t = self.ctrl.source
        idx = self.cmb_obj.findData(t)
        if idx < 0 and t.lower() == "all":
            idx = 0
        if idx >= 0:
            self.cmb_obj.setCurrentIndex(idx)
        else:
            # 选择名或不在列表中的对象：临时加入
            self.cmb_obj.addItem(t, t)
            self.cmb_obj.setCurrentIndex(self.cmb_obj.count() - 1)
        self._sync_related()
        if owned_row is not None:
            self.table.selectRow(owned_row)

    def _current_source(self) -> str:
        """当前对象 data。"""
        data = self.cmb_obj.currentData()
        return str(data) if data is not None else "all"

    def _sync_related(self) -> None:
        """根据源更新 related_hide 名单。"""
        src = self._current_source()
        if src.lower() == "all":
            self.ctrl.set_related_objects([name for obj in self.ctrl._resolve_sources("all")
                                           for name in self._related_for_source(obj)])
        else:
            self.ctrl.set_related_objects(self._related_for_source(src))

    def _on_source_changed(self, _idx: int = 0) -> None:
        """切换源对象。"""
        if self._loading:
            return
        self._sync_related()
        self._refresh_resnames()
        self.status.setText("已选择对象“%s”。更改后需点击应用。" % self._current_source())

    def _load_layers_to_table(self) -> None:
        """把 controller 层刷到表格。"""
        selected_id = self._editing_layer_id
        self._loading = True
        self.table.setRowCount(0)
        for i, L in enumerate(self.ctrl.layers):
            self.table.insertRow(i)
            style_item = QtWidgets.QTableWidgetItem(
                vr.STYLE_LABELS.get(L.style, L.style)
            )
            color_item = QtWidgets.QTableWidgetItem(
                vr.COLOR_MODE_LABELS.get(L.color_mode, L.color_mode)
            )
            if L.color_mode == vr.COLOR_SOLID:
                color_item.setIcon(make_swatch_icon(L.color_rgb))
            sele_item = QtWidgets.QTableWidgetItem(L.selection)
            for it in (style_item, color_item, sele_item):
                it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if not L.enabled:
                    it.setForeground(QtGui.QBrush(QtGui.QColor(150, 150, 150)))
            self.table.setItem(i, 0, style_item)
            self.table.setItem(i, 1, color_item)
            self.table.setItem(i, 2, sele_item)
            # 行头显示开关状态
            hdr = "✓" if L.enabled else "·"
            self.table.setVerticalHeaderItem(
                i, QtWidgets.QTableWidgetItem(hdr)
            )
        if self.ctrl.layers:
            row = next((i for i, layer in enumerate(self.ctrl.layers)
                        if layer.layer_id == selected_id), 0)
            self.table.selectRow(row)
            self._fill_editor(self.ctrl.layers[row])
        else:
            self._editing_layer_id = None
        self._loading = False

    def _selected_layer_index(self) -> int:
        """当前选中行。"""
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return 0
        return int(rows[0].row())

    def _fill_editor(self, L: vr.RepresentationLayer) -> None:
        """编辑区填入层。"""
        self._loading = True
        self._editing_layer_id = L.layer_id
        self.chk_enabled.setChecked(L.enabled)
        self.ed_sele.setText(L.selection)
        si = self.cmb_style.findData(L.style)
        if si >= 0:
            self.cmb_style.setCurrentIndex(si)
        self._update_style_help()
        ci = self.cmb_color_mode.findData(L.color_mode)
        if ci >= 0:
            self.cmb_color_mode.setCurrentIndex(ci)
        self._set_color_button(L.color_rgb)
        self.lbl_preview.setText("匹配：—")
        self._loading = False

    def _on_row_selected(self) -> None:
        """表格行切换 → 同步编辑区。"""
        if self._loading:
            return
        # 先提交上一行
        self._commit_editor_to_layer()
        idx = self._selected_layer_index()
        if 0 <= idx < len(self.ctrl.layers):
            self._fill_editor(self.ctrl.layers[idx])
            self._preview()

    def _commit_editor_to_layer(self) -> None:
        """把编辑区写回当前层（不改场景）。"""
        if self._loading or not self.ctrl.layers:
            return
        idx = next((i for i, layer in enumerate(self.ctrl.layers)
                    if layer.layer_id == self._editing_layer_id), -1)
        if idx < 0 or idx >= len(self.ctrl.layers):
            return
        L = self.ctrl.layers[idx]
        L.enabled = bool(self.chk_enabled.isChecked())
        L.selection = (self.ed_sele.text() or "").strip()
        L.style = str(self.cmb_style.currentData() or vr.STYLE_ORIGINAL)
        L.color_mode = str(
            self.cmb_color_mode.currentData() or vr.COLOR_ORIGINAL
        )
        L.color_rgb = self._current_rgb
        # 刷新表行显示
        self._loading = True
        if self.table.item(idx, 0):
            self.table.item(idx, 0).setText(
                vr.STYLE_LABELS.get(L.style, L.style)
            )
            self.table.item(idx, 1).setText(
                vr.COLOR_MODE_LABELS.get(L.color_mode, L.color_mode)
            )
            self.table.item(idx, 1).setIcon(
                make_swatch_icon(L.color_rgb) if L.color_mode == vr.COLOR_SOLID else QtGui.QIcon())
            self.table.item(idx, 2).setText(L.selection)
            self.table.setVerticalHeaderItem(
                idx,
                QtWidgets.QTableWidgetItem("✓" if L.enabled else "·"),
            )
        self._loading = False

    def _insert_example(self, text: str) -> None:
        """填入例子表达式。"""
        self.ed_sele.setText(text)
        self._commit_editor_to_layer()
        if self.chk_auto.isChecked():
            self.apply_layers()

    def _new_layer(self) -> None:
        """新建表示行。"""
        self._commit_editor_to_layer()
        self.ctrl.add_layer(
            selection="all",
            style=vr.STYLE_STICKS,
            color_mode=vr.COLOR_BYRESN,
        )
        self._load_layers_to_table()
        self.table.selectRow(len(self.ctrl.layers) - 1)
        self.status.setText("已新建表示行，编辑后点击应用。")

    def _delete_layer(self) -> None:
        """删除当前表示。"""
        if not self.ctrl.layers:
            return
        idx = self._selected_layer_index()
        lid = self.ctrl.layers[idx].layer_id
        self.ctrl.remove_layer(lid)
        self._load_layers_to_table()
        self.status.setText("已删除表示行（需应用后更新显示；删尽可恢复源）。")
        if self.chk_auto.isChecked():
            self.apply_layers()

    def _layers_by_component(self) -> None:
        """一键按组分分层。"""
        if not self._engine_ok():
            return
        src = self._current_source()
        try:
            layers = self.ctrl.suggest_component_layers(src)
            self.ctrl.set_layers(layers)
            self._load_layers_to_table()
            self.status.setText(
                "已按组分生成 %d 层规则，点击应用生效。" % len(layers)
            )
        except Exception as exc:
            self.status.setText("按组分分层失败：%s" % exc)

    def _preview(self) -> None:
        """预览当前选择表达式。"""
        if not self._engine_ok():
            return
        self._commit_editor_to_layer()
        src = self._current_source()
        expr = (self.ed_sele.text() or "").strip()
        rep = self.ctrl.preview(src, expr)
        if rep.ok:
            text = "匹配：%s 珠/原子；源内部键：%s" % (
                format(rep.count, ","), format(rep.bond_count, ","))
            if rep.external_bond_count:
                text += "；跨选择边界键：%s（不显示）" % format(rep.external_bond_count, ",")
            self.lbl_preview.setText(text)
            self.status.setText(rep.message)
        else:
            self.lbl_preview.setText("匹配：—")
            self.status.setText(rep.message)

    def _pick_color(self) -> None:
        """打开色板编辑单色。"""
        if self._closed:
            return

        def on_color(r: float, g: float, b: float) -> None:
            self._set_color_button((r, g, b))
            self._commit_editor_to_layer()
            if self.cmb_color_mode.currentData() != vr.COLOR_SOLID:
                idx = self.cmb_color_mode.findData(vr.COLOR_SOLID)
                if idx >= 0:
                    self.cmb_color_mode.setCurrentIndex(idx)
            if self.chk_auto.isChecked():
                self.apply_layers()

        def on_orig() -> None:
            idx = self.cmb_color_mode.findData(vr.COLOR_ORIGINAL)
            if idx >= 0:
                self.cmb_color_mode.setCurrentIndex(idx)
            self._commit_editor_to_layer()
            if self.chk_auto.isChecked():
                self.apply_layers()

        def on_elem() -> None:
            idx = self.cmb_color_mode.findData(vr.COLOR_BYELEMENT)
            if idx >= 0:
                self.cmb_color_mode.setCurrentIndex(idx)
            self._commit_editor_to_layer()
            if self.chk_auto.isChecked():
                self.apply_layers()

        pop = open_color_palette_popup(
            self,
            current_rgb=self._current_rgb,
            on_color=on_color,
            on_original=on_orig,
            on_element=on_elem,
            show_legacy_actions=True,
        )
        self._color_popups.append(pop)

    def _update_style_help(self) -> None:
        style = self.cmb_style.currentData()
        self.lbl_style_help.setText({
            vr.STYLE_STICKS: "棍状：只显示选择范围内的键，不显示原子球；需要内部连接信息。",
            vr.STYLE_BALLSTICK: "球棍：小原子球＋棍状键；需要内部连接信息。仅有坐标时请选球体。",
            vr.STYLE_CARTOON: "粗粒化核酸：沿已有骨架键绘制平滑管，保留链间与缺键处的断点。碱基方向明确时显示短棒。",
        }.get(style, ""))
        self.lbl_style_help.setVisible(bool(self.lbl_style_help.text()))

    def _on_edit_changed(self, *_args: Any) -> None:
        """控件更改。"""
        self._update_style_help()
        if self._loading:
            return
        self._commit_editor_to_layer()
        if self.chk_auto.isChecked():
            self.apply_layers()

    def _on_sele_finished(self) -> None:
        """选择表达式编辑结束（非逐字）。"""
        if self._loading:
            return
        self._commit_editor_to_layer()
        if self.chk_auto.isChecked():
            self.apply_layers()

    def apply_layers(self) -> None:
        """应用所有启用层。"""
        if not self._engine_ok():
            return
        self._commit_editor_to_layer()
        src = self._current_source()
        self._sync_related()
        self._busy = True
        self.btn_apply.setEnabled(False)
        try:
            # 仅一层且原始样式/原色且未接管时：不改场景
            only = self.ctrl.layers
            if (
                len(only) == 1
                and only[0].style == vr.STYLE_ORIGINAL
                and only[0].color_mode == vr.COLOR_ORIGINAL
                and only[0].selection.strip().lower() == "all"
                and not self.ctrl.owned
            ):
                self.status.setText("当前为原始表示，未修改场景。")
                return
            rep = self.ctrl.apply(src)
            msg = rep.message or ""
            if rep.warnings:
                msg += " " + " ".join(rep.warnings[:3])
            self.status.setText(msg)
            self.refresh_objects()
            # 保持选中源
            if self.ctrl.source:
                idx = self.cmb_obj.findData(self.ctrl.source)
                if idx >= 0:
                    self.cmb_obj.setCurrentIndex(idx)
            self._load_layers_to_table()
        except Exception as exc:
            self.status.setText("应用失败：%s" % exc)
        finally:
            self._busy = False
            self.btn_apply.setEnabled(True)
            # Recompute after table refresh, including automatic resname selection.
            if self.ctrl.layers and not self._closed:
                status = self.status.text()
                self._preview()
                self.status.setText(status)

    def restore_source(self) -> None:
        """恢复源对象。"""
        if not self._engine_ok():
            return
        try:
            msg = self.ctrl.restore()
            self.status.setText(msg)
            self.refresh_objects()
            self._load_layers_to_table()
        except Exception as exc:
            self.status.setText("恢复失败：%s" % exc)

    def closeEvent(self, event: Any) -> None:
        """关闭保留当前显示；标记避免退出后访引擎。"""
        self._closed = True
        super().closeEvent(event)


def clear_visual_representations_ref() -> None:
    """清除对话框单例引用（应用退出清理用）。"""
    global _dialog_ref
    _dialog_ref = None


def open_visual_representations(
    parent: QtWidgets.QWidget, target: str = "all"
) -> VisualRepresentationsDialog:
    """打开（或激活）选择与表示窗口；固定入口供主 GUI / 测试调用。"""
    global _dialog_ref
    cmd = getattr(parent, "cmd", None)
    if cmd is None:
        from pymol import cmd as cmd  # noqa: N813
    if _dialog_ref is None or _dialog_ref.parent() is not parent:
        _dialog_ref = VisualRepresentationsDialog(parent, cmd)
        _dialog_ref.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
    else:
        _dialog_ref._closed = False
        if not _dialog_ref._busy:
            _dialog_ref.refresh_objects()
    _dialog_ref.select_target(target or "all")
    _dialog_ref.show()
    _dialog_ref.raise_()
    _dialog_ref.activateWindow()
    return _dialog_ref
