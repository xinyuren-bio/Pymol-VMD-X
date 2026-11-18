"""Refresh an existing Pymol-VMD-X macOS installation from this checkout."""

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import tempfile
from zoneinfo import ZoneInfo


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_APP = Path("/Applications/Pymol-VMD-X.app")
METADATA = Path("pymol/_local_build.json")


def replace(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        staged = Path(stream.name)
        try:
            stream.write(data)
            stream.close()
            staged.chmod(0o644)
            staged.replace(path)
        finally:
            staged.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--app", type=Path, default=DEFAULT_APP)
    args = parser.parse_args()
    launcher = args.app / "Contents/MacOS/PymolX"
    site = Path(sysconfig.get_path("purelib")).resolve()
    if sys.platform != "darwin" or sys.prefix == sys.base_prefix:
        parser.error("请使用 macOS 安装版所用虚拟环境的 Python 运行。")
    if not launcher.is_file() or os.fsencode(str(site)) + b"\0" not in launcher.read_bytes():
        parser.error("所选应用的启动器没有指向这个虚拟环境。")
    if not list((site / "pymol").glob("_cmd*.so")):
        parser.error("虚拟环境里缺少已编译的 PyMOL 内核。")

    def git(*args):
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()

    build = {
        "commit": git("rev-parse", "--short", "HEAD"),
        "dirty": bool(git("status", "--porcelain")),
        "built_at": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds"),
    }
    patches = {METADATA: json.dumps(build, ensure_ascii=False).encode("utf-8")}
    for package in ("pymol", "pmg_qt"):
        for source in (ROOT / "modules" / package).rglob("*"):
            if not source.is_file() or source.suffix not in {".py", ".ui"}:
                continue
            relative = source.relative_to(ROOT / "modules")
            data = source.read_bytes()
            installed = site / relative
            if installed.is_file() and installed.read_bytes() == data:
                continue
            if source.suffix == ".py":
                compile(data, str(source), "exec")
            patches[relative] = data

    originals = {}
    try:
        for relative, data in patches.items():
            installed = site / relative
            installed.parent.mkdir(parents=True, exist_ok=True)
            originals[installed] = installed.read_bytes() if installed.is_file() else None
            replace(installed, data)
        with tempfile.TemporaryDirectory() as folder:
            marker = Path(folder) / "passed"
            script = Path(folder) / "check.py"
            checks = (
                "test_figure_toolbar.TestFigureToolbarQt.test_clicking_smart_styles_changes_the_viewport_settings",
                "test_visual_representations_gui.TestVisualRepGui.test_in_place_applies_colored_sticks_and_cg_cartoon_together",
            )
            script.write_text(
                "import sys, unittest\n"
                "from pathlib import Path\n"
                "from pymol import branding\n"
                f"assert Path(branding.__file__).resolve().is_relative_to(Path({str(site)!r}))\n"
                f"assert branding.local_build_info() == {build!r}\n"
                f"sys.path.insert(0, {str(ROOT / 'testing/tests/api')!r})\n"
                f"tests = {checks!r}\n"
                "result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromNames(tests))\n"
                f"if result.wasSuccessful() and result.testsRun == {len(checks)} and not result.skipped:\n"
                f"    Path({str(marker)!r}).touch()\n",
                encoding="utf-8",
            )
            env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
            subprocess.run([str(launcher), "-ckq", "-r", str(script)], cwd=ROOT, env=env,
                           timeout=45, check=True)
            if not marker.exists():
                raise RuntimeError("安装版的风格按钮或原位表示验证未通过，请核对所用的内核是否为新构建。")
    except BaseException:
        for installed, old in reversed(list(originals.items())):
            if old is None:
                installed.unlink(missing_ok=True)
            else:
                replace(installed, old)
        raise
    print(f"已更新 {len(patches) - 1} 个界面和 Python 文件。请完全退出并重新打开应用。")
    print(f"本机构建 {build['commit']}  {build['built_at']}")


if __name__ == "__main__":
    main()
