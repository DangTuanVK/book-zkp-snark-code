"""Ch.11: sinh vector kiểm thử và wrapper cấu hình cho RTL (montmul, ntt_bfly).

Không phải mã in sách; rtl/run_all.sh gọi script này trước khi mô phỏng.
Sinh ra:
  rtl/montmul/montmul_cfg.sv     wrapper montmul_w16/w32/w64/w256 (hằng p, p' cố định)
  rtl/ntt_bfly/ntt_bfly_cfg.sv   wrapper ntt_bfly_f449 (W = 16) và ntt_bfly_gl64 (W = 64)
  rtl/tb/vec/<ten>.hex           mỗi dòng một vector, ghép {a, b, y} hoặc {u, v, w~, x, y}
  rtl/tb/vec/cfg_<ten>.svh       macro cấu hình cho testbench
Giá trị kỳ vọng tính bằng snarklib.montgomery (REDC) và snarklib.fft (NTT trên F_449).
"""
import os, sys, random
GOC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GOC)
from snarklib.montgomery import Montgomery
from snarklib.field import GF
from snarklib import fft

RTL = os.path.join(GOC, "rtl")
VEC = os.path.join(RTL, "tb", "vec")

# Bốn cấu hình montmul: (tên, W, p, ghi chú). R = 2^W, một từ (s = 1).
R_BLS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
CAU_HINH = [
    ("montmul_w16", 16, 449, "F_449 của sách, R = 2^16"),
    ("montmul_w32", 32, 2**31 - 2**27 + 1, "BabyBear 15·2^27 + 1 (31 bit)"),
    ("montmul_w64", 64, 2**64 - 2**32 + 1, "Goldilocks 2^64 - 2^32 + 1"),
    ("montmul_w256", 256, R_BLS, "r của BLS12-381 (255 bit)"),
]
SO_NGAU_NHIEN = {16: 0, 32: 20000, 64: 20000, 256: 5000}


def hx(v, w):
    return format(v, "0%dx" % (w // 4))


def ghi_cfg(ten, dut, w, nvec, kieu):
    with open(os.path.join(VEC, "cfg_%s.svh" % ten), "w") as f:
        f.write("// sinh bởi ch11/sinh_vector_rtl.py -- không sửa tay\n")
        f.write("`define TB_W %d\n`define TB_NVEC %d\n" % (w, nvec))
        f.write('`define TB_VEC "%s.hex"\n`define TB_DUT %s\n`define TB_NAME "%s"\n'
                % (ten, dut, ten))
        f.write("`define TB_%s\n" % kieu)


def montmul_vectors(w, p, M, rng):
    if w == 16:                     # vét cạn mọi cặp (ã, b̃) trong [0, 448]^2
        cap = [(a, b) for a in range(p) for b in range(p)]
    else:                           # giá trị biên + ngẫu nhiên
        bien = [0, 1, 2, p - 1, p - 2, M.R_mod_p, M.R2, (p - 1) // 2]
        cap = [(a, b) for a in bien for b in bien]
        cap += [(rng.randrange(p), rng.randrange(p)) for _ in range(SO_NGAU_NHIEN[w])]
    return [(a, b, M.mul(a, b)) for a, b in cap]


def bfly_ref(p, M, u, v, wm):
    """Butterfly: t = REDC(v·w~) = v·w mod p (w~ = w·R mod p); (u + t, u - t) mod p."""
    t = M.mul(v, wm)
    return (u + t) % p, (u - t) % p


def ntt_trace(p, n, omega, a):
    """Chạy NTT Cooley-Tukey lặp, ghi lại mọi butterfly (u, v, w) theo đúng thứ tự."""
    a = fft._bit_reverse([x % p for x in a])
    vet, length = [], 2
    while length <= n:
        wl = pow(omega, n // length, p)
        for start in range(0, n, length):
            t = 1
            for k in range(length // 2):
                i, j = start + k, start + k + length // 2
                u, v = a[i], a[j]
                vet.append((u, v, t))
                a[i], a[j] = (u + v * t) % p, (u - v * t) % p
                t = t * wl % p
        length *= 2
    return a, vet


def bfly_vectors(p, w, rng, n_ntt, omega, so_ngau_nhien, F=None):
    M = Montgomery(p, word=w)
    out = []
    for n in n_ntt:
        om = pow(omega, (p - 1) // n, p) if F is None else int(F.root_of_unity(n))
        assert pow(om, n, p) == 1 and pow(om, n // 2, p) != 1
        a = [rng.randrange(p) for _ in range(n)]
        kq, vet = ntt_trace(p, n, om, a)
        if F is not None:            # đối chiếu với snarklib.fft.ntt
            assert kq == [int(x) for x in fft.ntt(F, [F(x) for x in a])]
        out += vet
    out += [(rng.randrange(p), rng.randrange(p), rng.randrange(p)) for _ in range(so_ngau_nhien)]
    out += [(p - 1, p - 1, p - 1), (0, 0, 0), (0, p - 1, 1), (p - 1, 0, 5)]
    res = []
    for u, v, tw in out:
        wm = tw * M.R % p             # twiddle lưu sẵn ở dạng Montgomery
        x, y = bfly_ref(p, M, u, v, wm)
        assert x == (u + v * tw) % p and y == (u - v * tw) % p
        res.append((u, v, wm, x, y))
    return M, res


def main():
    os.makedirs(VEC, exist_ok=True)
    rng = random.Random(1102)
    tom_tat = []
    # --- montmul ---
    dong = ["// sinh bởi ch11/sinh_vector_rtl.py -- không sửa tay",
            "// Wrapper cố định p và p' = -p^{-1} mod 2^W cho từng cấu hình (dùng cho mô phỏng và DC).",
            "`timescale 1ns/1ps"]
    for ten, w, p, ghi_chu in CAU_HINH:
        M = Montgomery(p, word=w)
        assert M.k == w and M.s == 1
        vs = montmul_vectors(w, p, M, rng)
        with open(os.path.join(VEC, ten + ".hex"), "w") as f:
            for a, b, y in vs:
                f.write(hx(a, w) + hx(b, w) + hx(y, w) + "\n")
        ghi_cfg(ten, ten, w, len(vs), "MONTMUL")
        dong += ["// %s: p = %s" % (ten, ghi_chu),
                 "module %s (input logic clk, rst_n, in_valid," % ten,
                 "  input logic [%d:0] a, b, output logic out_valid, output logic [%d:0] y);" % (w - 1, w - 1),
                 "  montmul #(.W(%d), .P(%d'h%x), .PINV(%d'h%x)) u_mm (.*);" % (w, w, p, w, M.p_prime),
                 "endmodule"]
        tom_tat.append("%s: W=%d p=0x%x p'=0x%x R mod p=0x%x, %d vector"
                       % (ten, w, p, M.p_prime, M.R_mod_p, len(vs)))
    open(os.path.join(RTL, "montmul", "montmul_cfg.sv"), "w").write("\n".join(dong) + "\n")
    # --- butterfly ---
    F449 = GF(449)
    GL = 2**64 - 2**32 + 1
    # 7 sinh nhóm nhân của Goldilocks: kiểm 7^((p-1)/q) != 1 với mọi ước nguyên tố q của p-1
    assert all(pow(7, (GL - 1) // q, GL) != 1 for q in (2, 3, 5, 17, 257, 65537))
    assert (GL - 1) == 2**32 * 3 * 5 * 17 * 257 * 65537
    dong = ["// sinh bởi ch11/sinh_vector_rtl.py -- không sửa tay", "`timescale 1ns/1ps"]
    for ten, p, w, n_ntt, omega, sn, F in [
            ("ntt_bfly_f449", 449, 16, (8, 16, 32, 64), None, 50000, F449),
            ("ntt_bfly_gl64", GL, 64, (8, 16, 32, 64), 7, 20000, None)]:
        M, vs = bfly_vectors(p, w, rng, n_ntt, omega, sn, F)
        with open(os.path.join(VEC, ten + ".hex"), "w") as f:
            for row in vs:
                f.write("".join(hx(x, w) for x in row) + "\n")
        ghi_cfg(ten, ten, w, len(vs), "BFLY")
        dong += ["module %s (input logic clk, rst_n, in_valid," % ten,
                 "  input logic [%d:0] u, v, w_m, output logic out_valid, output logic [%d:0] x, y);" % (w - 1, w - 1),
                 "  ntt_bfly #(.W(%d), .P(%d'h%x), .PINV(%d'h%x)) u_bf (.*);" % (w, w, p, w, M.p_prime),
                 "endmodule"]
        tom_tat.append("%s: W=%d p=0x%x, %d vector (gồm vết NTT n=%s, đối chiếu snarklib.fft: %s)"
                       % (ten, w, p, len(vs), ",".join(map(str, n_ntt)), "có" if F else "không"))
    open(os.path.join(RTL, "ntt_bfly", "ntt_bfly_cfg.sv"), "w").write("\n".join(dong) + "\n")
    print("\n".join(tom_tat))


if __name__ == "__main__":
    main()
