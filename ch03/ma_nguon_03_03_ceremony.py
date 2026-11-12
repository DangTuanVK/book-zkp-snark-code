"""Powers of Tau 3 người và các phép kiểm pairing (Mục 3.3, lst:ch3-ceremony-check).

SRS hiện tại: g1 = [[τ^i]_1, i = 0..D], g2 = [[1]_2, [τ]_2]. Người thứ j chọn
s_j bí mật, thay [τ^i]_1 bởi s_j^i·[τ^i]_1 và [τ]_2 bởi s_j·[τ]_2, công bố bằng
chứng ([s_j]_1, [s_j]_2). Verifier công khai kiểm từng đóng góp bằng pairing.
Toy: Type 1 (G2 = G1); số s_j cố định để output tái lập được.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from snarklib import opcount as oc
from snarklib.curve import toy_generator

D = 8                                          # bậc SRS
G1 = G2 = toy_generator()[2]


def bang(A, B, C, Dp):
    """e(A, B) = e(C, D)  <=>  e(A, B)·e(-C, D) = 1."""
    return oc.pairing_check([(A, B), (-C, Dp)])


def kiem_dong_gop(cu, moi, bc):
    """Trả về danh sách lỗi (rỗng = hợp lệ) của đóng góp cu -> moi, bằng chứng bc."""
    loi = []
    g1, (one2, tau2) = moi["g1"], moi["g2"]
    if g1[0] != G1 or one2 != G2 or bc["s1"].is_infinity():
        loi.append("K1: điểm sinh bị đổi hoặc s = 0")
    # K2: τ_mới = τ_cũ · s  (nối vào đóng góp trước, bằng [s]_2 đã công bố)
    if not bang(g1[1], one2, cu["g1"][1], bc["s2"]):
        loi.append("K2: không nối tiếp đóng góp trước")
    # K3: [s]_1 và [s]_2 cùng một s
    if not bang(bc["s1"], one2, G1, bc["s2"]):
        loi.append("K3: [s]_1, [s]_2 không khớp")
    # K4: [τ]_1 và [τ]_2 cùng một τ
    if not bang(g1[1], one2, G1, tau2):
        loi.append("K4: [tau]_1, [tau]_2 không khớp")
    # K5: e([τ^(i+1)]_1, [1]_2) = e([τ^i]_1, [τ]_2) với mọi i
    for i in range(1, len(g1) - 1):
        if not bang(g1[i + 1], one2, g1[i], tau2):
            loi.append("K5: sai lũy thừa tại i = %d" % i)
    return loi


def dong_gop(cu, s):
    """Đóng góp trung thực với bí mật s."""
    g1 = [oc.smul(P, pow(s, i, 449)) for i, P in enumerate(cu["g1"])]
    moi = {"g1": g1, "g2": [cu["g2"][0], oc.smul(cu["g2"][1], s, "G2")]}
    return moi, {"s1": oc.smul(G1, s), "s2": oc.smul(G2, s, "G2")}


def chay(ten, s_list, gian=None):
    srs = {"g1": [G1] * (D + 1), "g2": [G2, G2]}        # τ = 1 ban đầu
    print("== Kịch bản:", ten)
    for j, s in enumerate(s_list, 1):
        moi, bc = dong_gop(srs, s)
        if gian and j == 2:
            moi, bc = gian(srs, moi, bc)
        with oc.scope() as c:
            loi = kiem_dong_gop(srs, moi, bc)
        print("  người %d: %s  (%d pairing)" % (j, "HỢP LỆ" if not loi else "; ".join(loi), c["pairing"]))
        if not loi:
            srs = moi                    # đóng góp lỗi bị loại, người sau nối tiếp bản trước
    return srs


def thay_mot_phan_tu(cu, moi, bc):
    g1 = list(moi["g1"]); g1[5] = oc.smul(g1[5], 2)    # nhân đôi [τ^5]_1
    return {"g1": g1, "g2": moi["g2"]}, bc


def bo_lich_su(cu, moi, bc):
    t = 99                                              # dựng lại SRS từ τ' = 99 tự chọn
    g1 = [oc.smul(G1, pow(t, i, 449)) for i in range(D + 1)]
    return {"g1": g1, "g2": [G2, oc.smul(G2, t, "G2")]}, {"s1": oc.smul(G1, t), "s2": oc.smul(G2, t, "G2")}


if __name__ == "__main__":
    S = [17, 23, 31]
    tau = 17 * 23 * 31 % 449
    srs = chay("ba người trung thực", S)
    print("  tau cuối = 17·23·31 mod 449 = %d; [tau]_1 khớp: %s" % (tau, srs["g1"][1] == tau * G1))
    srs2 = chay("người 2 thay [tau^5]_1", S, thay_mot_phan_tu)
    srs3 = chay("người 2 bỏ lịch sử, dựng SRS từ tau' = 99", S, bo_lich_su)
    assert srs["g1"][1] == tau * G1
    assert srs2["g1"][1] == (17 * 31 % 449) * G1 and srs3["g1"][1] == (17 * 31 % 449) * G1
