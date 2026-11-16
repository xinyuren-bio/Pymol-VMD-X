/* ==================================================
 * 功能说明：macOS 应用入口。可执行文件放在 PymolX.app 内，菜单栏和程序坞显示 PymolX。
 * 使用方法：由安装步骤编译进 PymolX.app/Contents/MacOS/PymolX，不要单独当脚本运行。
 * 依赖环境：已安装并编译好的 pymolx 虚拟环境，以及对应版本的 Python 开发头文件。
 * 生成时间：2026-11-14 11:20
 * ================================================== */

#include <Python.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef PYMOLX_SITE
#define PYMOLX_SITE ""
#endif

/* 从 SystemExit 里取出退出码。图形界面正常退出时会走到这里。 */
static int pymolx_exit_code(void)
{
    PyObject *typ = NULL;
    PyObject *val = NULL;
    PyObject *tb = NULL;
    int code = 1;

    PyErr_Fetch(&typ, &val, &tb);
    if (val && PyLong_Check(val)) {
        code = (int)PyLong_AsLong(val);
    } else if (val && PyExceptionInstance_Check(val)) {
        PyObject *attr = PyObject_GetAttrString(val, "code");
        if (attr == Py_None || attr == NULL) {
            code = 0;
        } else if (PyLong_Check(attr)) {
            code = (int)PyLong_AsLong(attr);
        }
        Py_XDECREF(attr);
    }
    Py_XDECREF(typ);
    Py_XDECREF(val);
    Py_XDECREF(tb);
    return code;
}

int main(int argc, char **argv)
{
    wchar_t **wargv = NULL;
    char snippet[2048];
    int i = 0;
    int rc = 0;
    int code = 0;

    setlocale(LC_ALL, "");
    if (PYMOLX_SITE[0] == '\0') {
        fprintf(stderr, "PymolX: 未指定软件安装目录\n");
        return 1;
    }

    wargv = (wchar_t **)calloc((size_t)argc + 1, sizeof(wchar_t *));
    if (wargv == NULL) {
        fprintf(stderr, "PymolX: 内存不足\n");
        return 1;
    }
    for (i = 0; i < argc; i++) {
        wargv[i] = Py_DecodeLocale(argv[i], NULL);
        if (wargv[i] == NULL) {
            fprintf(stderr, "PymolX: 无法解析命令行参数\n");
            return 1;
        }
    }

    Py_SetProgramName(wargv[0]);
    Py_Initialize();
    /* 第三个参数为 0：不要把 .app 目录塞进 sys.path */
    PySys_SetArgvEx(argc, wargv, 0);
    /* 解释器本体在系统 Python 里，pymol 和 Qt 装在虚拟环境，启动前把该目录加进来 */
    snprintf(
        snippet,
        sizeof(snippet),
        "import site, sys\n"
        "site.addsitedir(r'%s')\n"
        "import pymol\n"
        "sys.exit(pymol.launch())\n",
        PYMOLX_SITE);
    rc = PyRun_SimpleString(snippet);
    if (rc != 0) {
        if (PyErr_ExceptionMatches(PyExc_SystemExit)) {
            code = pymolx_exit_code();
        } else {
            PyErr_Print();
            code = 1;
        }
    }
    Py_Finalize();
    return code;
}
