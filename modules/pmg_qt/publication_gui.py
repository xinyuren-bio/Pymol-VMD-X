# ==================================================
# 功能说明：科研出图无模态 Qt 对话框、多配体勾选与智能主题预览共享控制器。
# 使用方法：由 pymol_qt_gui / figure_toolbar 调用 open_publication_dialog 或 apply_smart_theme。
# 依赖环境：PyMOL + PySide6（pymol.Qt）；后端 pymol.publication。
# 生成时间：2026-11-14 13:28
# ==================================================

"""中文科研出图对话框与智能风格入口。"""

from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple

from pymol.Qt import QtCore, QtWidgets
from pymol.Qt.utils import AsyncFunc, PopupOnException, getSaveFileNameWithExt
from pymol import publication as pub
from .publication_color_controls import PublicationColorControls

from pmg_qt.figure_toolbar import is_appearance_customized

Qt = QtCore.Qt


def choose_ligands(
    parent: QtWidgets.QWidget,
    candidates: Sequence[pub.LigandCandidate],
    selected: Optional[Sequence[pub.LigandCandidate]] = None,
) -> Tuple[List[pub.LigandCandidate], bool]:
    """
    多选配体确认对话框；返回 (所选列表, 是否点了确定)。

    selected 为 None 时默认全选；否则按先前键勾选。取消不改写调用方状态。
    空列表在确定后表示仅蛋白总览（若场景有蛋白）。
    """
    dlg = QtWidgets.QDialog(parent)
    dlg.setWindowTitle("选择配体")
    dlg.setModal(True)
    root = QtWidgets.QVBoxLayout(dlg)
    tip = QtWidgets.QLabel(
        "检测到多个配体候选，请勾选要用于智能风格的配体"
        "（可全不选以做蛋白总览）："
    )
    tip.setWordWrap(True)
    root.addWidget(tip)

    lst = QtWidgets.QListWidget()
    lst.setMaximumHeight(160)
    prev_keys = None
    if selected is not None:
        prev_keys = {pub.candidate_key(x) for x in selected if x is not None}
    for lig in candidates:
        item = QtWidgets.QListWidgetItem(lig.label)
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsEnabled
        )
        item.setData(Qt.ItemDataRole.UserRole, lig)
        if prev_keys is None:
            item.setCheckState(Qt.CheckState.Checked)
        else:
            item.setCheckState(
                Qt.CheckState.Checked
                if pub.candidate_key(lig) in prev_keys
                else Qt.CheckState.Unchecked
            )
        lst.addItem(item)
    root.addWidget(lst)

    row = QtWidgets.QHBoxLayout()
    btn_all = QtWidgets.QPushButton("全选")
    btn_clear = QtWidgets.QPushButton("清空")
    row.addWidget(btn_all)
    row.addWidget(btn_clear)
    row.addStretch(1)
    root.addLayout(row)

    def _set_all(checked: bool) -> None:
        """全选或清空列表勾选。"""
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for i in range(lst.count()):
            lst.item(i).setCheckState(state)

    btn_all.clicked.connect(lambda: _set_all(True))
    btn_clear.clicked.connect(lambda: _set_all(False))

    buttons = QtWidgets.QDialogButtonBox(
        QtWidgets.QDialogButtonBox.StandardButton.Ok
        | QtWidgets.QDialogButtonBox.StandardButton.Cancel
    )
    buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("确定")
    buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Cancel).setText("取消")
    buttons.accepted.connect(dlg.accept)
    buttons.rejected.connect(dlg.reject)
    root.addWidget(buttons)

    if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
        return [], False
    chosen: List[pub.LigandCandidate] = []
    for i in range(lst.count()):
        item = lst.item(i)
        if item.checkState() == Qt.CheckState.Checked:
            lig = item.data(Qt.ItemDataRole.UserRole)
            if lig is not None:
                chosen.append(lig)
    return chosen, True


class PublicationDialog(QtWidgets.QDialog):
    """
    紧凑的科研出图对话框（无模态）。

    说明文字会标明样式为出版启发、候选为启发式识别，不声称期刊合规。
    与工具栏/菜单共享同一 fig 控制器，保证 get_session 快照可恢复。
    """

    def __init__(self, parent: QtWidgets.QWidget, cmd: Any) -> None:
        # 无模态，便于边调视角边出图
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("科研出图 · 预览与导出")
        self.setMinimumWidth(440)
        self.cmd = cmd
        self.fig = pub.PublicationFigure(cmd)
        self._scene = pub.SceneInfo()
        self._busy = False
        self._worker: Optional[AsyncFunc] = None
        # 用户已明确确认的配体键集合（智能主题复用；构造自动勾选不算确认）
        self._ligands_explicit: bool = False
        self._confirmed_keys: List[Tuple[Any, ...]] = []

        self._build_ui()
        self.refresh_candidates()

    def _build_ui(self) -> None:
        """构建紧凑中文界面（多配体可勾选列表）。"""
        root = QtWidgets.QVBoxLayout(self)
        root.setSpacing(8)

        tip = QtWidgets.QLabel(
            "布局可选结构总览或结合口袋；智能风格为期刊简约 / 卡通插画 / 通透表面。"
            "可勾选多个配体；若已用顶部工具栏调整背景/光照，可勾选“使用当前背景与光照”保留外观；"
            "预览后可再调视角并勾选“保持当前视角”导出。"
        )
        tip.setWordWrap(True)
        tip.setStyleSheet("color: #444; font-size: 11px;")
        root.addWidget(tip)

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        root.addWidget(self.status)

        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        lig_col = QtWidgets.QVBoxLayout()
        lig_col.setSpacing(4)
        self.ligand_list = QtWidgets.QListWidget()
        self.ligand_list.setMaximumHeight(120)
        self.ligand_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection
        )
        lig_col.addWidget(self.ligand_list)
        self.ligand_count_label = QtWidgets.QLabel("已选 0 个配体")
        self.ligand_count_label.setStyleSheet("color: #555; font-size: 11px;")
        lig_col.addWidget(self.ligand_count_label)
        row_lig = QtWidgets.QHBoxLayout()
        self.btn_select_all = QtWidgets.QPushButton("全选")
        self.btn_clear = QtWidgets.QPushButton("清空")
        self.btn_refresh = QtWidgets.QPushButton("刷新")
        self.btn_select_all.clicked.connect(self._on_select_all)
        self.btn_clear.clicked.connect(self._on_clear_ligands)
        self.btn_refresh.clicked.connect(self.refresh_candidates)
        row_lig.addWidget(self.btn_select_all)
        row_lig.addWidget(self.btn_clear)
        row_lig.addWidget(self.btn_refresh)
        row_lig.addStretch(1)
        lig_col.addLayout(row_lig)
        form.addRow("配体候选：", lig_col)

        # 布局样式（overview / pocket）
        self.style_overview = QtWidgets.QRadioButton("结构总览")
        self.style_pocket = QtWidgets.QRadioButton("结合口袋（4.5 Å）")
        self.style_overview.setChecked(True)
        style_box = QtWidgets.QHBoxLayout()
        style_box.addWidget(self.style_overview)
        style_box.addWidget(self.style_pocket)
        style_box.addStretch(1)
        form.addRow("布局：", style_box)

        # 智能风格（三种主题）
        self.theme_combo = QtWidgets.QComboBox()
        for meta in pub.list_themes():
            self.theme_combo.addItem(meta["label"], meta["id"])
        self.theme_combo.setCurrentIndex(0)
        self.theme_desc = QtWidgets.QLabel("")
        self.theme_desc.setWordWrap(True)
        self.theme_desc.setStyleSheet("color: #555; font-size: 11px;")
        theme_col = QtWidgets.QVBoxLayout()
        theme_col.setSpacing(2)
        theme_col.addWidget(self.theme_combo)
        theme_col.addWidget(self.theme_desc)
        form.addRow("智能风格：", theme_col)
        self.theme_combo.currentIndexChanged.connect(self._update_theme_desc)
        self._update_theme_desc()

        self.color_controls = PublicationColorControls(self.cmd, self._current_theme, self)
        self.color_controls.changed.connect(self._on_colors_changed)
        form.addRow("配色：", self.color_controls)

        # 导出参数
        self.spin_w = QtWidgets.QSpinBox()
        self.spin_w.setRange(100, 8000)
        self.spin_w.setValue(1800)
        self.spin_h = QtWidgets.QSpinBox()
        self.spin_h.setRange(100, 8000)
        self.spin_h.setValue(1400)
        self.spin_dpi = QtWidgets.QSpinBox()
        self.spin_dpi.setRange(72, 600)
        self.spin_dpi.setValue(300)
        wh = QtWidgets.QHBoxLayout()
        wh.addWidget(QtWidgets.QLabel("宽"))
        wh.addWidget(self.spin_w)
        wh.addWidget(QtWidgets.QLabel("高"))
        wh.addWidget(self.spin_h)
        wh.addWidget(QtWidgets.QLabel("DPI"))
        wh.addWidget(self.spin_dpi)
        form.addRow("尺寸：", wh)

        self.chk_transparent = QtWidgets.QCheckBox("透明背景")
        self.chk_keep_view = QtWidgets.QCheckBox("保持当前视角")
        self.chk_keep_appearance = QtWidgets.QCheckBox("使用当前背景与光照")
        self.chk_keep_appearance.setToolTip(
            "保留当前背景、光照以及轮廓/阴影设置；仍会应用蛋白/配体表示样式。"
            "透明背景选项始终可单独生效。"
        )
        self.chk_keep_appearance.setChecked(is_appearance_customized(self.cmd))
        self._appearance_user_set = False
        self.chk_keep_appearance.clicked.connect(self._appearance_choice_changed)
        opt = QtWidgets.QHBoxLayout()
        opt.addWidget(self.chk_transparent)
        opt.addWidget(self.chk_keep_view)
        opt.addWidget(self.chk_keep_appearance)
        opt.addStretch(1)
        form.addRow("选项：", opt)
        root.addLayout(form)

        btns = QtWidgets.QHBoxLayout()
        self.btn_preview = QtWidgets.QPushButton("预览样式")
        self.btn_restore = QtWidgets.QPushButton("恢复原始")
        self.btn_restore.setToolTip("恢复首次预览或出图前的整个会话，包括视角和对象。")
        self.btn_export = QtWidgets.QPushButton("导出 PNG")
        self.btn_export.setDefault(True)
        self.btn_export.setStyleSheet(
            "QPushButton { font-weight: bold; padding: 6px 12px;"
            " background: #2F6FED; color: white; border-radius: 4px; }"
            "QPushButton:disabled { background: #99A; color: #EEE; }"
        )
        self.btn_preview.clicked.connect(self.on_preview)
        self.btn_restore.clicked.connect(self.on_restore)
        self.btn_export.clicked.connect(self.on_export)
        btns.addWidget(self.btn_preview)
        btns.addWidget(self.btn_restore)
        btns.addStretch(1)
        btns.addWidget(self.btn_export)
        root.addLayout(btns)

        self.busy_bar = QtWidgets.QProgressBar()
        self.busy_bar.setRange(0, 0)
        self.busy_bar.setVisible(False)
        root.addWidget(self.busy_bar)

        self.ligand_list.itemChanged.connect(self._on_ligand_item_changed)
        self.style_overview.toggled.connect(self._update_pocket_enabled)

    def _update_theme_desc(self) -> None:
        """根据主题下拉刷新说明文字。"""
        tid = self.theme_combo.currentData()
        meta = pub.THEME_META.get(tid or "", {})
        self.theme_desc.setText(meta.get("description", ""))
        if hasattr(self, "color_controls"):
            self.color_controls.refresh()

    def _on_colors_changed(self):
        self.status.setText("配色已设置，点击“预览样式”应用；导出时会使用当前配色。")

    def _set_busy(self, on: bool, msg: str = "") -> None:
        """切换忙碌态，避免重复点击导致并发出图；含列表与主题下拉。"""
        self._busy = on
        self.busy_bar.setVisible(on)
        for w in (
            self.btn_preview,
            self.btn_restore,
            self.btn_export,
            self.btn_refresh,
            self.btn_select_all,
            self.btn_clear,
            self.ligand_list,
            self.style_overview,
            self.style_pocket,
            self.theme_combo,
            self.color_controls,
            self.spin_w,
            self.spin_h,
            self.spin_dpi,
            self.chk_transparent,
            self.chk_keep_view,
            self.chk_keep_appearance,
        ):
            w.setEnabled(not on)
        if on and msg:
            self.status.setText(msg)
        if not on:
            self._update_pocket_enabled()

    def _update_pocket_enabled(self) -> None:
        """无蛋白+至少一配体时禁用口袋样式。"""
        if self._busy:
            return
        ligs = self._current_ligands()
        ok = bool(self._scene.has_protein and ligs)
        self.style_pocket.setEnabled(ok)
        if not ok and self.style_pocket.isChecked():
            self.style_overview.setChecked(True)

    def _mark_ligands_explicit(self) -> None:
        """将当前勾选记为用户显式确认（供智能主题复用）。"""
        ligs = self._current_ligands()
        self._ligands_explicit = True
        self._confirmed_keys = [
            k for k in (pub.candidate_key(x) for x in ligs) if k is not None
        ]
        self._update_ligand_count_label()
        self._update_pocket_enabled()

    def _on_ligand_item_changed(self, _item: QtWidgets.QListWidgetItem = None) -> None:
        """列表勾选变化视为显式选择确认。"""
        if self.ligand_list.signalsBlocked():
            return
        self._mark_ligands_explicit()

    def _on_select_all(self) -> None:
        """全选所有配体候选（显式确认）。"""
        self.ligand_list.blockSignals(True)
        for i in range(self.ligand_list.count()):
            self.ligand_list.item(i).setCheckState(Qt.CheckState.Checked)
        self.ligand_list.blockSignals(False)
        self._mark_ligands_explicit()

    def _on_clear_ligands(self) -> None:
        """清空勾选（显式确认；有蛋白时可做仅蛋白总览）。"""
        self.ligand_list.blockSignals(True)
        for i in range(self.ligand_list.count()):
            self.ligand_list.item(i).setCheckState(Qt.CheckState.Unchecked)
        self.ligand_list.blockSignals(False)
        self._mark_ligands_explicit()

    def _update_ligand_count_label(self) -> None:
        """刷新已选数量标签与状态栏提示。"""
        ligs = self._current_ligands()
        n = len(ligs)
        self.ligand_count_label.setText("已选 %d 个配体" % n)
        names = "；".join(
            "%s/%s/%s" % (x.object, x.chain or "_", x.resi) for x in ligs[:6]
        )
        if n > 6:
            names += "…"
        self.ligand_count_label.setToolTip(names if names else "未选配体（蛋白总览）")

    def refresh_candidates(self) -> None:
        """重新检测启用对象中的配体候选；尽量保留仍有效的勾选。"""
        if self._busy:
            return
        previous_keys = {
            pub.candidate_key(x)
            for x in self._current_ligands()
            if pub.candidate_key(x) is not None
        }
        # 若尚无列表项，回退到已确认键
        if not previous_keys and self._ligands_explicit:
            previous_keys = set(self._confirmed_keys)

        had_list_items = self.ligand_list.count() > 0
        had_nonempty_explicit = bool(self._ligands_explicit and self._confirmed_keys)

        self._scene = pub.detect_scene(self.cmd)
        self.ligand_list.blockSignals(True)
        self.ligand_list.clear()
        # 首次填充可自动全选，但不把构造勾选记为显式确认
        auto_check_all = (
            not self._ligands_explicit and not previous_keys and not had_list_items
        )
        for lig in self._scene.ligands:
            item = QtWidgets.QListWidgetItem(lig.label)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
            )
            item.setData(Qt.ItemDataRole.UserRole, lig)
            key = pub.candidate_key(lig)
            if previous_keys and key in previous_keys:
                item.setCheckState(Qt.CheckState.Checked)
            elif self._ligands_explicit and key in set(self._confirmed_keys):
                item.setCheckState(Qt.CheckState.Checked)
            elif auto_check_all:
                item.setCheckState(Qt.CheckState.Checked)
            else:
                item.setCheckState(Qt.CheckState.Unchecked)
            self.ligand_list.addItem(item)
        self.ligand_list.blockSignals(False)

        # 同步确认键：去掉过期项；非空确认若全部消失则复位确认
        if self._ligands_explicit:
            alive = [
                pub.candidate_key(x)
                for x in self._current_ligands()
                if pub.candidate_key(x) is not None
            ]
            if had_nonempty_explicit and not alive:
                self._ligands_explicit = False
                self._confirmed_keys = []
            else:
                self._confirmed_keys = list(alive)

        n = len(self._scene.ligands)
        n_sel = len(self._current_ligands())
        if self._scene.empty:
            self.status.setText(self._scene.note or "场景为空。")
        else:
            self.status.setText(
                "%s  对象 %d 个；配体候选 %d 个；已选 %d 个。"
                % (self._scene.note, len(self._scene.objects), n, n_sel)
            )
        self._update_ligand_count_label()
        self._update_pocket_enabled()
        self.btn_restore.setEnabled(self.fig.has_snapshot())

    def _current_ligands(self) -> List[pub.LigandCandidate]:
        """当前勾选的全部配体候选。"""
        out: List[pub.LigandCandidate] = []
        for i in range(self.ligand_list.count()):
            item = self.ligand_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                lig = item.data(Qt.ItemDataRole.UserRole)
                if lig is not None:
                    out.append(lig)
        return out

    def _current_ligand(self) -> Optional[pub.LigandCandidate]:
        """兼容：返回勾选列表中的第一个配体。"""
        ligs = self._current_ligands()
        return ligs[0] if ligs else None

    def _current_style(self) -> str:
        """当前布局样式标识。"""
        if self.style_pocket.isChecked():
            return pub.STYLE_POCKET
        return pub.STYLE_OVERVIEW

    def _current_theme(self) -> str:
        """当前智能风格标识。"""
        tid = self.theme_combo.currentData()
        return tid if tid else pub.THEME_JOURNAL

    def set_theme_ui(self, theme: str) -> None:
        """同步主题下拉（不触发预览）。"""
        theme = pub.validate_theme(theme)
        for i in range(self.theme_combo.count()):
            if self.theme_combo.itemData(i) == theme:
                self.theme_combo.blockSignals(True)
                self.theme_combo.setCurrentIndex(i)
                self.theme_combo.blockSignals(False)
                self._update_theme_desc()
                break

    def _keep_view(self) -> bool:
        return bool(self.chk_keep_view.isChecked())

    def _appearance_choice_changed(self) -> None:
        """标记用户已手动改过“使用当前背景与光照”。"""
        self._appearance_user_set = True

    def force_theme_appearance_choice(self) -> None:
        """
        智能主题点击后：取消“使用当前背景与光照”并记为故意选择，
        避免工具栏定制标志在后续导出时覆盖主题外观。
        """
        from pmg_qt.figure_toolbar import clear_appearance_customized
        clear_appearance_customized(self.cmd)
        self.chk_keep_appearance.setChecked(False)
        self._appearance_user_set = False

    def _keep_appearance(self) -> bool:
        """是否保留当前背景与光照（及轮廓/阴影）。"""
        return bool(self.chk_keep_appearance.isChecked())

    def resolve_ligands_for_smart_theme(
        self, parent: QtWidgets.QWidget
    ) -> Tuple[List[pub.LigandCandidate], bool]:
        """
        为智能主题解析配体列表：单候选自动选；多候选弹多选框；取消不改场景。

        返回 (ligands, ok)。ok=False 表示用户取消。
        蛋白无配体时 ligands=[] 且 ok=True；空场景 ok=False。
        已确认子集（含显式空选）直接复用。
        """
        self.refresh_candidates()
        info = self._scene
        if info.empty:
            return [], False

        ligands = list(info.ligands)
        if not ligands:
            # 仅蛋白：允许无配体应用主题
            self._ligands_explicit = True
            self._confirmed_keys = []
            return [], True

        # 显式空选也必须保留；不能因恰好只有一个候选而擅自重新勾选。
        if self._ligands_explicit:
            keyset = set(self._confirmed_keys)
            chosen = [lig for lig in ligands if pub.candidate_key(lig) in keyset]
            self._sync_list_to_confirmed()
            return chosen, True

        if len(ligands) == 1:
            lig = ligands[0]
            key = pub.candidate_key(lig)
            self._ligands_explicit = True
            self._confirmed_keys = [key] if key is not None else []
            self._sync_list_to_confirmed()
            return [lig], True

        # 无确认：弹出可复用多选框（默认全选；selected=None）
        chosen, ok = choose_ligands(parent, ligands, selected=None)
        if not ok:
            return [], False
        self._ligands_explicit = True
        self._confirmed_keys = [
            k for k in (pub.candidate_key(x) for x in chosen) if k is not None
        ]
        self._sync_list_to_confirmed()
        return chosen, True

    def resolve_ligand_for_smart_theme(
        self, parent: QtWidgets.QWidget
    ) -> Tuple[Optional[pub.LigandCandidate], bool]:
        """兼容旧接口：返回所选列表的第一个配体。"""
        ligs, ok = self.resolve_ligands_for_smart_theme(parent)
        if not ok:
            return None, False
        return (ligs[0] if ligs else None), True

    def _sync_list_to_confirmed(self) -> None:
        """将对话框勾选同步到已确认键（不触发二次确认逻辑）。"""
        keyset = set(self._confirmed_keys)
        self.ligand_list.blockSignals(True)
        for i in range(self.ligand_list.count()):
            item = self.ligand_list.item(i)
            lig = item.data(Qt.ItemDataRole.UserRole)
            key = pub.candidate_key(lig) if lig is not None else None
            item.setCheckState(
                Qt.CheckState.Checked
                if key is not None and key in keyset
                else Qt.CheckState.Unchecked
            )
        self.ligand_list.blockSignals(False)
        self._update_ligand_count_label()
        self._update_pocket_enabled()

    def apply_theme_preview(
        self, theme: str, ligands: Optional[Sequence[pub.LigandCandidate]] = None
    ) -> str:
        """
        立即应用主题预览（总览布局），使用主题自身外观。

        返回简短中文状态文案；ligands 为完整所选列表。
        """
        theme = pub.validate_theme(theme)
        self.set_theme_ui(theme)
        self.force_theme_appearance_choice()
        self.style_overview.setChecked(True)
        lig_list = list(ligands) if ligands is not None else self._current_ligands()
        self.fig.apply_style(
            style=pub.STYLE_OVERVIEW,
            ligands=lig_list,
            keep_view=False,
            keep_appearance=False,
            theme=theme,
            **self.color_controls.overrides(theme),
        )
        self.btn_restore.setEnabled(True)
        label = pub.THEME_META[theme]["label"]
        n = len(lig_list)
        if n:
            target = "%d 个配体" % n
            tip_names = "；".join(x.resn for x in lig_list[:4])
            if n > 4:
                tip_names += "…"
        elif self._scene.has_protein:
            target = "蛋白总览（无配体）"
            tip_names = ""
        else:
            target = "当前分子"
            tip_names = ""
        msg = "已应用「%s」→ %s" % (label, target)
        self.status.setText(msg)
        self.status.setToolTip(tip_names)
        self._update_ligand_count_label()
        return msg

    @PopupOnException.decorator
    def on_preview(self) -> None:
        """预览：直接改当前可编辑场景。"""
        ligs = self._current_ligands()
        if self._scene.empty and not ligs:
            QtWidgets.QMessageBox.information(self, "科研出图", self._scene.note)
            return
        style = self._current_style()
        theme = self._current_theme()
        if style == pub.STYLE_POCKET and (not ligs or not self._scene.has_protein):
            QtWidgets.QMessageBox.warning(
                self, "科研出图", "结合口袋需要蛋白与至少一个配体，请勾选配体候选。"
            )
            return
        self.fig.apply_style(
            style=style,
            ligands=ligs,
            keep_view=self._keep_view(),
            keep_appearance=self._keep_appearance(),
            theme=theme,
            **self.color_controls.overrides(theme),
        )
        n = len(ligs)
        self.status.setText(
            "已应用预览样式（已选 %d 个配体）。可继续调整视角，或点“恢复原始”。" % n
        )
        self.btn_restore.setEnabled(True)

    @PopupOnException.decorator
    def on_restore(self) -> None:
        """恢复首次样式前的会话快照。"""
        if not self.fig.restore_original():
            QtWidgets.QMessageBox.information(
                self, "科研出图", "尚无快照（需先预览或出图一次）。"
            )
            return
        self.status.setText("已恢复原始会话。")
        self.refresh_candidates()

    def on_export(self) -> None:
        """导出 PNG：先选文件名，再异步应用样式并光线追踪。"""
        if self._busy:
            return
        ligs = self._current_ligands()
        if self._scene.empty and not ligs:
            QtWidgets.QMessageBox.information(self, "科研出图", self._scene.note)
            return
        style = self._current_style()
        theme = self._current_theme()
        if style == pub.STYLE_POCKET and (not ligs or not self._scene.has_protein):
            QtWidgets.QMessageBox.warning(
                self, "科研出图", "结合口袋需要蛋白与至少一个配体，请勾选配体候选。"
            )
            return

        parent = self.parent() or self
        initial = getattr(parent, "initialdir", "") or ""
        fp = getSaveFileNameWithExt(
            self,
            "保存出版图 PNG",
            initial,
            filter="PNG 图像 (*.png)",
        )
        if not fp:
            return
        if hasattr(parent, "initialdir"):
            import os

            parent.initialdir = os.path.dirname(fp)

        w = int(self.spin_w.value())
        h = int(self.spin_h.value())
        dpi = float(self.spin_dpi.value())
        transparent = bool(self.chk_transparent.isChecked())
        keep_view = self._keep_view()
        keep_appearance = self._keep_appearance()
        fig = self.fig
        # 捕获当前勾选列表，避免异步期间被刷新冲掉
        ligs_job = list(ligs)
        color_overrides = self.color_controls.overrides(theme)

        def _job() -> str:
            """后台线程：样式 + ray PNG（传入全部所选配体与主题）。"""
            return fig.export_png(
                fp,
                style=style,
                ligands=ligs_job,
                width=w,
                height=h,
                dpi=dpi,
                transparent=transparent,
                keep_view=keep_view,
                keep_appearance=keep_appearance,
                theme=theme,
                **color_overrides,
            )

        def _done(payload: tuple) -> None:
            """导出结束回调（主线程）。"""
            result, exc = payload
            self._set_busy(False)
            self.btn_restore.setEnabled(self.fig.has_snapshot())
            if getattr(self.parent(), "_quit_requested", False):
                return
            if exc is not None:
                QtWidgets.QMessageBox.warning(
                    self, "科研出图", "导出失败：%s" % (exc,)
                )
                self.status.setText("导出失败。")
                return
            self.refresh_candidates()
            self.status.setText(
                "已保存：%s（导出时已选 %d 个配体）" % (result, len(ligs_job))
            )

        try:
            self.fig.ensure_snapshot()
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "科研出图", "无法保存原始场景：%s" % exc)
            return
        self._set_busy(True, "正在光线追踪导出，请稍候…")
        self._worker = AsyncFunc(_job, finishslot=_done)
        self._worker()

    def is_export_running(self) -> bool:
        """导出工作线程是否仍在运行（勿在运行中销毁 QThread）。"""
        if self._busy:
            return True
        w = self._worker
        return w is not None and w.isRunning()

    def closeEvent(self, event):
        if self._busy:
            event.ignore()
            return
        super().closeEvent(event)


_dialog_ref: Optional[PublicationDialog] = None


def publication_export_running() -> bool:
    """全局：当前科研出图导出是否仍在进行。"""
    dlg = _dialog_ref
    return dlg is not None and dlg.is_export_running()


def clear_publication_dialog_ref() -> None:
    """清除出图对话框单例引用（退出清理用，不强制终止线程）。"""
    global _dialog_ref
    _dialog_ref = None


def get_or_create_publication_dialog(
    parent: QtWidgets.QWidget, show: bool = False
) -> PublicationDialog:
    """
    获取或创建共享 PublicationDialog（同一 fig 控制器），默认不显示。

    设计思路：工具栏智能风格与菜单/导出对话框必须共用快照，禁止另建控制器。
    打开时刷新候选但保留已选主题与已勾选配体，不重置为默认。
    """
    global _dialog_ref
    cmd = getattr(parent, "cmd", None)
    if cmd is None:
        from pymol import cmd as cmd  # noqa: N813
    if _dialog_ref is None or _dialog_ref.parent() is not parent:
        _dialog_ref = PublicationDialog(parent, cmd)
        _dialog_ref.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
    if not _dialog_ref.is_export_running():
        _dialog_ref.refresh_candidates()
        if not _dialog_ref._appearance_user_set:
            _dialog_ref.chk_keep_appearance.setChecked(is_appearance_customized(cmd))
    if show:
        _dialog_ref.show()
        _dialog_ref.raise_()
        _dialog_ref.activateWindow()
    return _dialog_ref


def open_publication_dialog(parent: QtWidgets.QWidget) -> PublicationDialog:
    """打开（或激活）无模态科研出图对话框。"""
    return get_or_create_publication_dialog(parent, show=True)


def _feedback(parent: QtWidgets.QWidget, msg: str) -> None:
    """优先用主窗口 statusBar，否则静默可用时弹信息。"""
    bar = getattr(parent, "statusBar", None)
    if callable(bar):
        try:
            bar().showMessage(msg, 8000)
            return
        except Exception:
            pass
    # 无 statusBar 时不打断；对话框自身 status 已写


def apply_smart_theme(parent: QtWidgets.QWidget, theme: str) -> bool:
    """
    工具栏/菜单共用的智能风格入口：立即预览，不打开导出对话框。

    多配体且无确认时弹出多选；取消则场景不变。导出进行中拒绝切换。
    返回是否成功应用。
    """
    if publication_export_running():
        _feedback(parent, "正在导出，请稍候再切换智能风格。")
        return False
    try:
        theme = pub.validate_theme(theme)
    except ValueError as exc:
        QtWidgets.QMessageBox.warning(parent, "科研出图", str(exc))
        return False

    dlg = get_or_create_publication_dialog(parent, show=False)
    if dlg.is_export_running():
        _feedback(parent, "正在导出，请稍候再切换智能风格。")
        return False

    info = pub.detect_scene(dlg.cmd)
    if info.empty:
        note = info.note or "场景中没有启用的分子对象。"
        QtWidgets.QMessageBox.information(parent, "科研出图", note)
        _feedback(parent, note)
        return False

    ligs, ok = dlg.resolve_ligands_for_smart_theme(parent)
    if not ok:
        # 用户取消：场景保持不变
        _feedback(parent, "已取消，场景未改动。")
        return False

    try:
        msg = dlg.apply_theme_preview(theme, ligs)
    except Exception as exc:
        QtWidgets.QMessageBox.warning(parent, "科研出图", "应用风格失败：%s" % (exc,))
        return False
    _feedback(parent, msg)
    # 同步工具栏勾选（主题会改背景/光照）
    bar = getattr(parent, "figure_toolbar", None)
    if bar is not None and hasattr(bar, "refresh_checked_states"):
        try:
            bar.refresh_checked_states()
        except Exception:
            pass
    return True
