"""Bộ đếm phép toán nhóm và pairing dùng chung cho toàn sách.

Mọi lược đồ cam kết trong snarklib gọi các hàm bọc dưới đây thay vì dùng
trực tiếp toán tử của curve.Point, nhờ đó ta đếm được chi phí THẬT của từng
thuật toán (không phải công thức ước lượng):

    from snarklib import opcount as oc
    with oc.scope() as c:            # c là dict đếm riêng cho khối lệnh này
        C = oc.msm(srs, coeffs)      # một MSM cỡ len(coeffs)
        ok = oc.pairing_check([(A, B), (C, D)])
    c["msm"], c["msm_size"], c["pairing"]

Khóa đếm (mỗi khóa có hậu tố nhóm, ví dụ "smul_G1", "msm_G2"; khóa không hậu
tố là tổng mọi nhóm):
    add        số phép cộng/nhân đôi điểm cấp thấp (kể cả bên trong smul, msm)
    smul       số phép nhân vô hướng đứng riêng (gọi smul)
    msm        số lần gọi MSM;  msm_size: tổng cỡ các MSM
    msm_nontrivial: số vô hướng khác 0 và 1 trong các MSM (phép nhân thật)
    msm_sizes  danh sách cỡ từng MSM (để lập bảng)
    pairing    số pairing (= số vòng Miller);  final_exp: số lũy thừa cuối
Đường cong toy là pairing đối xứng (Type 1, G1 = G2); tham số group chỉ để
ghi sổ, giúp văn bản sách viết dạng Type 3.
"""

from contextlib import contextmanager

_STACK = [{}]          # ngăn xếp bộ đếm: mọi scope đang mở đều được cộng dồn


def _bump(key, n=1, group=None):
    for c in _STACK:
        c[key] = c.get(key, 0) + n
        if group is not None:
            k = "%s_%s" % (key, group)
            c[k] = c.get(k, 0) + n


def reset():
    """Xóa bộ đếm toàn cục (không ảnh hưởng các scope đang mở)."""
    _STACK[0].clear()


def get():
    """Bản sao bộ đếm toàn cục."""
    return dict(_STACK[0])


@contextmanager
def scope():
    """Mở một bộ đếm riêng; trả về dict được cập nhật trong khối with."""
    c = {"msm_sizes": []}
    _STACK.append(c)
    try:
        yield c
    finally:
        for i in range(len(_STACK) - 1, 0, -1):     # gỡ đúng đối tượng (so định danh)
            if _STACK[i] is c:
                del _STACK[i]
                break


# --- phép toán nhóm ---------------------------------------------------------

def add(P, Q, group="G1"):
    """P + Q (đếm 1 phép cộng)."""
    _bump("add", 1, group)
    return P + Q


def sub(P, Q, group="G1"):
    """P - Q (đếm 1 phép cộng; phủ định điểm coi như miễn phí)."""
    _bump("add", 1, group)
    return P - Q


def _double_and_add(P, k, group):
    """k·P bằng double-and-add, đếm từng phép cộng/nhân đôi."""
    R = P.E.infinity()
    if k == 0 or P.is_infinity():
        return R
    for bit in bin(k)[2:]:
        if not R.is_infinity():
            R = R + R
            _bump("add", 1, group)
        if bit == "1":
            if R.is_infinity():
                R = P
            else:
                R = R + P
                _bump("add", 1, group)
    return R


def smul(P, k, group="G1"):
    """Nhân vô hướng k·P (đếm 1 smul và các phép cộng bên trong)."""
    _bump("smul", 1, group)
    k = int(k)
    if k < 0:
        return _double_and_add(-P, -k, group)
    return _double_and_add(P, k, group)


def msm(points, scalars, group="G1"):
    """MSM Σ k_i·P_i, cài ngây thơ (từng nhân vô hướng rồi cộng dồn).

    Đếm 1 MSM cỡ len(points); vô hướng 0 bị bỏ qua, vô hướng 1 chỉ tốn phép cộng.
    Chương 10 thay phần thân bằng Pippenger nhưng giữ nguyên cách ghi sổ.
    """
    if len(points) != len(scalars):
        raise ValueError("số điểm và số vô hướng khác nhau")
    _bump("msm", 1, group)
    _bump("msm_size", len(points), group)
    for c in _STACK:
        c.setdefault("msm_sizes", []).append(len(points))
    acc = None
    for P, k in zip(points, scalars):
        k = int(k)
        if k == 0 or P.is_infinity():
            continue
        if k != 1:
            _bump("msm_nontrivial", 1, group)
        T = _double_and_add(P, k, group) if k != 1 else P
        if acc is None:
            acc = T
        else:
            acc = acc + T
            _bump("add", 1, group)
    if acc is None:
        if not points:
            raise ValueError("MSM rỗng")
        return points[0].E.infinity()
    return acc


# --- pairing ----------------------------------------------------------------

def pairing(P, Q):
    """e(P, Q) đầy đủ: 1 vòng Miller + 1 lũy thừa cuối."""
    from .pairing import tate
    _bump("pairing", 1)
    _bump("final_exp", 1)
    return tate(P, Q)


def multi_pairing(pairs):
    """Π e(P_i, Q_i) với MỘT lũy thừa cuối chung (cách verifier thực tế làm)."""
    from .pairing import miller
    from .field import GF2
    E = pairs[0][0].E
    K = GF2(E.q)
    f = K(1)
    _bump("pairing", len(pairs))
    _bump("final_exp", 1)
    for P, Q in pairs:
        if P.is_infinity() or Q.is_infinity():
            continue
        f = f * miller(P, Q)
    from .curve import R_TOY
    return f ** ((E.q * E.q - 1) // R_TOY)


def pairing_check(pairs):
    """Kiểm Π e(P_i, Q_i) = 1 (dạng chuẩn của mọi phương trình kiểm KZG)."""
    return multi_pairing(pairs) == 1


def summary(c, keys=("pairing", "msm", "msm_size", "smul", "add")):
    """Chuỗi tóm tắt các khóa chính của một bộ đếm."""
    return ", ".join("%s=%d" % (k, c.get(k, 0)) for k in keys)
