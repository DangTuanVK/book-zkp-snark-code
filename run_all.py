#!/usr/bin/env python3
"""Chạy toàn bộ test của snarklib, sinh lại ket_qua/hang_so.txt, ghi ket_qua/tests.txt.

    python run_all.py            # chạy test + sinh hằng số
    python run_all.py --no-const # chỉ chạy test

Thoát với mã 1 nếu có test lỗi hoặc hằng số lệch kỳ vọng.
"""

import importlib
import io
import json
import os
import platform
import sys
import time
import traceback

ROOT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(ROOT, "tests")
sys.path.insert(0, ROOT)
sys.path.insert(0, TESTS)


def test_modules():
    """Danh sách module test: lấy từ manifest.json, bổ sung file test_*.py còn thiếu."""
    names = []
    try:
        with open(os.path.join(ROOT, "manifest.json"), encoding="utf-8") as f:
            names = [t["module"] for t in json.load(f).get("tests", [])]
    except (OSError, ValueError, KeyError):
        pass
    for fn in sorted(os.listdir(TESTS)):
        if fn.startswith("test_") and fn.endswith(".py") and fn[:-3] not in names:
            names.append(fn[:-3])
    return names


def main():
    out = io.StringIO()

    def log(s=""):
        print(s)
        out.write(s + "\n")

    log("snarklib — kết quả chạy test")
    log("Python %s" % platform.python_version())
    log("")
    n_pass, n_fail, t0 = 0, 0, time.time()
    for mod_name in test_modules():
        mod = importlib.import_module(mod_name)
        funcs = sorted(k for k in dir(mod) if k.startswith("test_") and callable(getattr(mod, k)))
        for fn in funcs:
            t = time.time()
            try:
                getattr(mod, fn)()
                n_pass += 1
                log("ĐẠT  %-18s %-24s %6.2f s" % (mod_name, fn, time.time() - t))
            except Exception:  # noqa: BLE001
                n_fail += 1
                log("LỖI  %-18s %-24s" % (mod_name, fn))
                log(traceback.format_exc())
    log("")
    log("Tổng: %d đạt, %d lỗi, %.1f s" % (n_pass, n_fail, time.time() - t0))

    bad_const = []
    if "--no-const" not in sys.argv:
        from snarklib import hang_so
        path, _, bad_const = hang_so.write()
        log("Hằng số: đã ghi %s" % os.path.relpath(path, ROOT))
        for s in bad_const:
            log("  LỆCH KỲ VỌNG: " + s)

    os.makedirs(os.path.join(ROOT, "ket_qua"), exist_ok=True)
    with open(os.path.join(ROOT, "ket_qua", "tests.txt"), "w", encoding="utf-8") as f:
        f.write(out.getvalue())
    sys.exit(1 if (n_fail or bad_const) else 0)


if __name__ == "__main__":
    main()
