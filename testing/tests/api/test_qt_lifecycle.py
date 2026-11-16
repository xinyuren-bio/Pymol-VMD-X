# ==================================================
# 功能说明：Qt 优雅退出与 _qt_quit 钩子的回归测试（含 offscreen 子进程）。
# 使用方法：由仓库测试运行器或 python -m unittest 发现本文件；勿在无显示环境强依赖窗口。
# 依赖环境：pymol + PySide6；子进程用例使用 QT_QPA_PLATFORM=offscreen。
# 生成时间：2026-11-14 12:55
# ==================================================

"""Qt 生命周期与 cmd._qt_quit 兼容性测试（标准 unittest）。"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
import unittest


def _qt_available() -> bool:
    """检测当前解释器是否可导入 pymol.Qt。"""
    try:
        from pymol.Qt import QtWidgets  # noqa: F401
        return True
    except Exception:
        return False


_LIFECYCLE_SCRIPT = textwrap.dedent(
    r"""
    import os
    import sys
    import threading

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    mode = os.environ.get("PYMOLX_TEST_QUIT_MODE", "request")
    code = int(os.environ.get("PYMOLX_TEST_QUIT_CODE", "0"))
    delay = int(os.environ.get("PYMOLX_TEST_QUIT_DELAY_MS", "300"))

    # 安静启动，缩短路径
    sys.argv = ["PymolX", "-q", "-Q"]

    from pymol import invocation
    invocation.parse_args(list(sys.argv))

    from pymol.Qt import QtCore
    from pmg_qt.pymol_qt_gui import PyMOLQtGUI, execapp

    _orig_show = PyMOLQtGUI.show

    def _show(self):
        _orig_show(self)

        def _fire():
            if mode == "close":
                self.close()
            elif mode == "confirm":
                self.confirm_quit()
            elif mode == "cmd":
                self.cmd.quit(code)
            elif mode == "pending_feedback":
                # 引擎停止后已经排队的刷新回调不得再访问释放的数据。
                QtCore.QTimer.singleShot(0, self._finish_shutdown)
                QtCore.QTimer.singleShot(0, self.update_feedback)
                QtCore.QTimer.singleShot(0, self.update_progress)
            elif mode == "app_quit":
                from pymol.Qt import QtWidgets
                QtWidgets.QApplication.instance().quit()
            elif mode == "export_close":
                import tempfile
                from pmg_qt import publication_gui as pg
                self.cmd.load(os.environ["PYMOLX_TEST_STRUCTURE"])
                self.open_publication_dialog()
                dialog = pg._dialog_ref
                dialog.spin_w.setValue(900)
                dialog.spin_h.setValue(700)
                dialog.style_pocket.setChecked(True)
                output = os.path.join(tempfile.mkdtemp(prefix="pymolx_quit_"), "figure.png")
                pg.getSaveFileNameWithExt = lambda *a, **kw: output
                dialog.on_export()
                assert dialog.is_export_running()
                self.close()
            elif mode == "cmd_bg":
                threading.Thread(
                    target=lambda: self.cmd.quit(code), daemon=True
                ).start()
            else:
                self.request_quit(code)

        QtCore.QTimer.singleShot(delay, _fire)

    PyMOLQtGUI.show = _show
    rc = execapp()
    raise SystemExit(int(rc) if rc is not None else 0)
    """
)


class TestQtQuitHandlerAPI(unittest.TestCase):
    """无 GUI：验证 optional _qt_quit 钩子与 headless 路径兼容性。"""

    def test_qt_quit_handler_short_circuits_internal_quit(self):
        """注册 _qt_quit 时 internal._quit 不得走到 _cmd.quit/PExit。"""
        from pymol import cmd, internal

        cmd.fragment("benzene", "quit_hook_test")

        seen = []

        def _handler(code=0, _self=None):
            seen.append(int(code))
            return "qt"

        prev = getattr(cmd, "_qt_quit", None)
        cmd._qt_quit = _handler
        try:
            r = internal._quit(7, cmd)
            self.assertEqual(seen, [7])
            self.assertEqual(r, "qt")
        finally:
            if prev is None and hasattr(cmd, "_qt_quit"):
                delattr(cmd, "_qt_quit")
            else:
                cmd._qt_quit = prev

    def test_quit_uses_qt_handler_when_registered(self):
        """cmd.quit 在存在 _qt_quit 时应直接委托并保留退出码。"""
        from pymol import cmd

        seen = []

        def _handler(code=0, _self=None):
            seen.append(int(code))
            return None

        prev = getattr(cmd, "_qt_quit", None)
        cmd._qt_quit = _handler
        try:
            cmd.quit(3)
            self.assertEqual(seen, [3])
        finally:
            if prev is None and hasattr(cmd, "_qt_quit"):
                delattr(cmd, "_qt_quit")
            else:
                cmd._qt_quit = prev

    def test_headless_native_quit_preserves_code(self):
        """没有 Qt 窗口时仍使用原有退出路径。"""
        child = subprocess.run(
            [sys.executable, "-c",
             "from pymol import cmd; cmd.fragment('benzene'); cmd.quit(3); print('UNREACHABLE')"],
            capture_output=True, text=True, timeout=15,
        )
        self.assertEqual(child.returncode, 3, child.stderr)
        self.assertNotIn("UNREACHABLE", child.stdout)


@unittest.skipUnless(_qt_available(), "pymol.Qt unavailable")
class TestQtLifecycleSubprocess(unittest.TestCase):
    """真实 Qt 子进程：offscreen 下 close/菜单/cmd.quit 应正常返回退出码。"""

    def _run(self, mode: str, code: int = 0, timeout: int = 90) -> subprocess.CompletedProcess:
        env = os.environ.copy()
        env["QT_QPA_PLATFORM"] = "offscreen"
        from pathlib import Path
        env["PYMOLX_TEST_STRUCTURE"] = str(Path(__file__).resolve().parents[2] / "data" / "1oky.pdb.gz")
        env["PYMOLX_TEST_QUIT_MODE"] = mode
        env["PYMOLX_TEST_QUIT_CODE"] = str(code)
        env["PYMOLX_TEST_QUIT_DELAY_MS"] = "400"
        # 避免 macOS 上 cocoa 插件隐藏旗标等问题干扰 headless
        env.setdefault("QT_LOGGING_RULES", "*.debug=false")
        return subprocess.run(
            [sys.executable, "-c", _LIFECYCLE_SCRIPT],
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    def test_request_quit_returns_zero(self):
        """request_quit(0) 子进程应正常结束且返回 0。"""
        p = self._run("request", 0)
        self.assertEqual(
            p.returncode,
            0,
            msg="stdout=%r stderr=%r" % (p.stdout[-2000:], p.stderr[-2000:]),
        )

    def test_request_quit_preserves_code(self):
        """request_quit / cmd.quit 应保留非零退出码。"""
        p = self._run("cmd", 5)
        self.assertEqual(
            p.returncode,
            5,
            msg="stdout=%r stderr=%r" % (p.stdout[-2000:], p.stderr[-2000:]),
        )

    def test_window_close_event(self):
        """窗口 close() 应走优雅退出而非崩溃（非 0/-11/-6 信号）。"""
        p = self._run("close", 0)
        self.assertEqual(
            p.returncode,
            0,
            msg="stdout=%r stderr=%r" % (p.stdout[-2000:], p.stderr[-2000:]),
        )
        self.assertNotIn("Segmentation fault", p.stderr)

    def test_confirm_quit_menu_path(self):
        """File/Quit（confirm_quit）应正常退出。"""
        p = self._run("confirm", 0)
        self.assertEqual(
            p.returncode,
            0,
            msg="stdout=%r stderr=%r" % (p.stdout[-2000:], p.stderr[-2000:]),
        )

    def test_pending_feedback_after_engine_stop(self):
        p = self._run("pending_feedback")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])

    def test_application_quit_event(self):
        p = self._run("app_quit")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])

    def test_close_during_real_export(self):
        p = self._run("export_close", timeout=20)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])

    def test_background_thread_quit(self):
        """后台线程 cmd.quit 应投递到主线程并保留退出码。"""
        p = self._run("cmd_bg", 4)
        self.assertEqual(
            p.returncode,
            4,
            msg="stdout=%r stderr=%r" % (p.stdout[-2000:], p.stderr[-2000:]),
        )


if __name__ == "__main__":
    unittest.main()
