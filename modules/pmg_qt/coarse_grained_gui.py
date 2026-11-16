# ==================================================
# 功能说明：粗粒化 Martini 风格连通显示无模态 Qt 对话框（检测/应用/恢复）。
# 使用方法：由工具栏或菜单调用 open_coarse_grained_dialog(parent)；退出可用 clear_coarse_grained_dialog_ref。
# 依赖环境：PyMOL + PySide6（pymol.Qt）；后端 pymol.coarse_grained / cg_templates。
# 生成时间：2026-11-14 16:30
# ==================================================

"""中文粗粒化连通显示对话框。"""

from __future__ import annotations

from typing import Any, Optional

from pymol.Qt import QtCore, QtWidgets

from pymol import coarse_grained as cg

Qt = QtCore.Qt

_dialog_ref: Optional["CoarseGrainedDialog"] = None


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


class CoarseGrainedDialog(QtWidgets.QDialog):
    """
    无模态粗粒化显示对话框。

    可选源对象（含已禁用）、刷新检测、显示模式、自动匹配已保存拓扑、
    可选键/拓扑文件、Apply/Restore。默认先模板后几何近似。
    """

    def __init__(self, parent: QtWidgets.QWidget, cmd: Any) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("粗粒化连通显示 · Martini 风格")
        self.setMinimumWidth(460)
        self.cmd = cmd
        self.ctrl = cg.CGDisplayController(cmd)
        self._busy = False
        self._closed = False
        self._build_ui()
        self.refresh_sources()

    def _build_ui(self) -> None:
        """构建中文界面。"""
        root = QtWidgets.QVBoxLayout(self)
        root.setSpacing(8)

        tip = QtWidgets.QLabel(
            "自动匹配已保存的脂质与 RNA 连接模板，生成连续或球棍显示。"
            "未匹配的组分使用几何近似，原始粗粒化坐标保持不变。"
        )
        tip.setWordWrap(True)
        tip.setStyleSheet("color: #444; font-size: 11px;")
        root.addWidget(tip)

        # 源对象
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("源对象："))
        self.cmb_obj = QtWidgets.QComboBox()
        self.cmb_obj.setMinimumWidth(180)
        row.addWidget(self.cmb_obj, 1)
        self.btn_refresh = QtWidgets.QPushButton("刷新检测")
        row.addWidget(self.btn_refresh)
        root.addLayout(row)

        # 模式
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(QtWidgets.QLabel("显示模式："))
        self.cmb_mode = QtWidgets.QComboBox()
        for m in cg.list_modes():
            self.cmb_mode.addItem(m["label"], m["id"])
        row2.addWidget(self.cmb_mode, 1)
        root.addLayout(row2)

        # 自动模板匹配（默认开启）
        self.chk_auto_templates = QtWidgets.QCheckBox("自动匹配已保存拓扑")
        self.chk_auto_templates.setChecked(True)
        self.chk_auto_templates.setToolTip(
            "按完整残基与珠名组合匹配。PDB 通常没有力场珠子类型，"
            "名称匹配不能鉴定力场版本或电荷状态。关闭则使用纯几何近似。"
            "显示保留粗粒化分辨率，不生成全原子坐标。"
        )
        root.addWidget(self.chk_auto_templates)

        # 可选键文件
        row3 = QtWidgets.QHBoxLayout()
        row3.addWidget(QtWidgets.QLabel("键/拓扑文件："))
        self.ed_bonds = QtWidgets.QLineEdit()
        self.ed_bonds.setPlaceholderText(
            "可选覆盖：全局零基键对 .txt，或原子数匹配的 .itp/.top"
        )
        row3.addWidget(self.ed_bonds, 1)
        self.btn_browse = QtWidgets.QPushButton("浏览…")
        row3.addWidget(self.btn_browse)
        self.btn_clear_file = QtWidgets.QPushButton("清除")
        row3.addWidget(self.btn_clear_file)
        root.addLayout(row3)

        # 摘要
        self.summary = QtWidgets.QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setMaximumHeight(180)
        root.addWidget(self.summary)

        self.status = QtWidgets.QLabel("就绪。")
        self.status.setWordWrap(True)
        root.addWidget(self.status)

        # 按钮
        brow = QtWidgets.QHBoxLayout()
        self.btn_apply = QtWidgets.QPushButton("应用显示")
        self.btn_restore = QtWidgets.QPushButton("恢复源对象")
        brow.addWidget(self.btn_apply)
        brow.addWidget(self.btn_restore)
        brow.addStretch(1)
        self.btn_close = QtWidgets.QPushButton("关闭")
        brow.addWidget(self.btn_close)
        root.addLayout(brow)

        self.btn_refresh.clicked.connect(self.refresh_detection)
        self.btn_browse.clicked.connect(self._browse_bonds)
        self.btn_clear_file.clicked.connect(lambda: self.ed_bonds.clear())
        self.btn_apply.clicked.connect(self.apply_display)
        self.btn_restore.clicked.connect(self.restore_display)
        self.btn_close.clicked.connect(self.close)
        self.cmb_obj.currentIndexChanged.connect(self.refresh_detection)

    def _engine_ok(self) -> bool:
        """退出/关机状态下禁止访问引擎。"""
        if self._closed or self._busy:
            return False
        parent = self.parent()
        if getattr(parent, "_quit_requested", False) or getattr(parent, "_shutting_down", False):
            return False
        if getattr(self.cmd, "_COb", None) is None:
            return False
        from pmg_qt.publication_gui import publication_export_running
        if publication_export_running():
            self.status.setText("正在导出图片，请完成后再调整粗粒化显示。")
            return False
        if not _parent_alive(self.parent()):
            return False
        try:
            self.cmd.get_names("public_objects", enabled_only=0)
            return True
        except Exception:
            self.status.setText("PyMOL 引擎不可用或已退出。")
            return False

    def refresh_sources(self) -> None:
        """刷新源对象下拉（含禁用对象，排除自有副本）。"""
        if not self._engine_ok():
            return
        cur = self.cmb_obj.currentData()
        self.cmb_obj.blockSignals(True)
        self.cmb_obj.clear()
        for n in self.ctrl.list_sources():
            self.cmb_obj.addItem(n, n)
        self.cmb_obj.blockSignals(False)
        if cur:
            idx = self.cmb_obj.findData(cur)
            if idx >= 0:
                self.cmb_obj.setCurrentIndex(idx)
        if self.cmb_obj.count() == 0:
            self.summary.setPlainText("场景中没有可用的分子对象。")
            self.status.setText("无源对象。")
        else:
            self.refresh_detection()

    def refresh_detection(self) -> None:
        """对当前源对象运行检测并更新摘要。"""
        if not self._engine_ok():
            return
        obj = self.cmb_obj.currentData()
        if not obj:
            self.summary.setPlainText("请选择源对象。")
            return
        try:
            det = self.ctrl.detect(str(obj))
            self.summary.setPlainText(cg.summary_text(det, self.ctrl.last_apply))
            self.status.setText(det.note or "检测完成。")
        except Exception as exc:
            self.status.setText("检测失败：%s" % exc)

    def _browse_bonds(self) -> None:
        """浏览键对或拓扑文件。"""
        fp, _flt = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "选择键对或拓扑文件",
            "",
            "键对与单文件拓扑 (*.txt *.bonds *.dat *.itp *.top);;所有文件 (*)",
        )
        if fp:
            self.ed_bonds.setText(fp)

    def apply_display(self) -> None:
        """应用粗粒化显示。"""
        if not self._engine_ok():
            return
        obj = self.cmb_obj.currentData()
        if not obj:
            self.status.setText("请选择源对象。")
            return
        mode = self.cmb_mode.currentData() or cg.MODE_CONTINUOUS
        bonds = (self.ed_bonds.text() or "").strip() or None
        auto_tpl = bool(self.chk_auto_templates.isChecked())
        self._busy = True
        self.btn_apply.setEnabled(False)
        try:
            rep = self.ctrl.apply(
                str(obj),
                mode=str(mode),
                bonds_file=bonds,
                auto_templates=auto_tpl,
            )
            det = self.ctrl.last_detection or cg.DetectionReport(empty=True)
            self.summary.setPlainText(cg.summary_text(det, rep))
            status = rep.message or ""
            if rep.ok:
                status = "已应用：模板覆盖 %s / %s 珠；几何近似 %s 珠。" % (
                    format(rep.matched_beads, ","), format(det.bead_count, ","),
                    format(rep.approximate_beads, ","),
                )
            self.status.setText(status)
            # 应用后源可能被禁用，仍保留在列表中
            self.refresh_sources()
        except Exception as exc:
            self.status.setText("应用失败：%s" % exc)
        finally:
            self._busy = False
            self.btn_apply.setEnabled(True)

    def restore_display(self) -> None:
        """恢复源对象可见性并清理自有输出。"""
        if not self._engine_ok():
            return
        try:
            msg = self.ctrl.restore()
            self.status.setText(msg)
            if self.ctrl.last_detection:
                self.summary.setPlainText(
                    cg.summary_text(self.ctrl.last_detection, None)
                )
            self.refresh_sources()
        except Exception as exc:
            self.status.setText("恢复失败：%s" % exc)

    def closeEvent(self, event: Any) -> None:
        """关闭时标记，避免退出后访问引擎。"""
        self._closed = True
        super().closeEvent(event)


def clear_coarse_grained_dialog_ref() -> None:
    """清除对话框单例引用（应用退出清理用）。"""
    global _dialog_ref
    _dialog_ref = None


def get_or_create_coarse_grained_dialog(
    parent: QtWidgets.QWidget, show: bool = False
) -> CoarseGrainedDialog:
    """获取或创建无模态粗粒化对话框。"""
    global _dialog_ref
    cmd = getattr(parent, "cmd", None)
    if cmd is None:
        from pymol import cmd as cmd  # noqa: N813
    if _dialog_ref is None or _dialog_ref.parent() is not parent:
        _dialog_ref = CoarseGrainedDialog(parent, cmd)
        _dialog_ref.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
    else:
        _dialog_ref._closed = False
        if not _dialog_ref._busy:
            _dialog_ref.refresh_sources()
    if show:
        _dialog_ref.show()
        _dialog_ref.raise_()
        _dialog_ref.activateWindow()
    return _dialog_ref


def open_coarse_grained_dialog(parent: QtWidgets.QWidget) -> CoarseGrainedDialog:
    """打开（或激活）粗粒化连通显示对话框，供工具栏/菜单集成。"""
    return get_or_create_coarse_grained_dialog(parent, show=True)
