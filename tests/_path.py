"""Thêm thư mục gốc kho mã vào sys.path để import snarklib khi chạy trực tiếp."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def run_module(globs):
    """Chạy mọi hàm test_* trong module (khi gọi file trực tiếp)."""
    names = sorted(k for k in globs if k.startswith("test_") and callable(globs[k]))
    for k in names:
        globs[k]()
        print("ĐẠT", k)
