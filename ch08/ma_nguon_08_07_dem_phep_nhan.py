"""Chương 8 — Đếm phép nhân trường của Prover: sum-check (Spartan, HyperPlonk) vs FFT của PLONK.

Mạch: c bản sao độc lập của mạch ví dụ (c = 1..64) -> R1CS 13c ràng buộc (đệm 2^4..2^10),
Plonkish 20c hàng (đệm n = 32c). Đếm THẬT bằng trường đếm (sumcheck.counting_field):
  - Spartan: phần IOP (Az, Bz, Cz; eq; hai sum-check; W~(r_y)), không gồm PCS;
  - HyperPlonk: zerocheck cổng + productcheck hoán vị (mức PIOP, không gồm PCS);
  - PLONK: 8 FFT cỡ 4n + 5 FFT cỡ 2n + 12 FFT cỡ n (gabizon2019plonk.txt:247-248),
    mỗi FFT đếm trên snarklib.fft.ntt; so với công thức ≈ 54·n·log2(n) (bảng 1 của bài).
Đây là ƯỚC TÍNH phần FFT (không gồm MSM, nhân điểm). F_449 chỉ có NTT tới 64 điểm nên
phép đếm chạy trên p = 998244353 = 119·2^23 + 1; số phép nhân không phụ thuộc p.
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.field import GF
from snarklib.transcript import Transcript
from snarklib.dsl import Circuit
from snarklib.fft import ntt
from snarklib import sumcheck as sc, spartan

P = 998244353
CF = sc.counting_field(GF(P))


def mach_nhan_ban(c):
    """c bản sao mạch phân tích (N_i công khai, a_i, b_i 5 bit, a_i·b_i = N_i)."""
    C = Circuit(GF(P))
    vals = {}
    for i in range(c):
        N, a, b = C.public("N%d" % i), C.private("a%d" % i), C.private("b%d" % i)
        C.to_bits(a, 5, prefix="a%d_" % i); C.to_bits(b, 5, prefix="b%d_" % i)
        C.assert_mul(a, b, N)
        vals.update({"N%d" % i: 391, "a%d" % i: 17, "b%d" % i: 23})
    return C, vals


def dem_fft(n):
    CF.counter["mul"] = 0
    ntt(CF, [CF(i) for i in range(n)])
    return CF.counter["mul"]


rows = []
for c in (1, 2, 4, 8, 16, 32, 64):
    C, vals = mach_nhan_ban(c)
    r1 = C.compile_r1cs(); x, w = C.r1cs_assignment(vals)
    inst = spartan.pad_r1cs(r1)
    tot, ph = spartan.prover_mul_count(CF, inst, [CF(v) for v in r1.z(x, w)], Transcript(b"dem", r=P))
    pk = C.compile_plonkish(); a, b, cc, xp = C.plonkish_assignment(vals)
    T = sc.plonkish_tables(CF, pk, a, b, cc, xp)
    n = len(T["a"]); lg = n.bit_length() - 1
    CF.counter["mul"] = CF.counter["inv"] = 0
    sc.hyperplonk_prove(CF, T, Transcript(b"dem-hp", r=P))
    hp_mul, hp_inv = CF.counter["mul"], CF.counter["inv"]
    fft = 8 * dem_fft(4 * n) + 5 * dem_fft(2 * n) + 12 * dem_fft(n)
    rows.append({"ban_sao": c, "r1cs_m": r1.m, "r1cs_dem": inst["size"], "spartan_iop": tot,
                 "spartan_sumcheck": ph["sc1"] + ph["sc2"], "plonkish_hang": pk.n, "n_plonk": n,
                 "plonk_fft": fft, "plonk_54nlogn": 54 * n * lg, "hyperplonk_piop": hp_mul,
                 "hyperplonk_inv": hp_inv})
    print(rows[-1])

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ket_qua", "ch08", "so_phep_nhan.csv")
with open(out, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
print("Tỉ số FFT PLONK / IOP Spartan:", ["%.1f" % (r["plonk_fft"] / r["spartan_iop"]) for r in rows])
print("Tỉ số FFT PLONK / PIOP HyperPlonk:", ["%.1f" % (r["plonk_fft"] / r["hyperplonk_piop"]) for r in rows])
s = rows[0]["r1cs_dem"].bit_length() - 1
assert rows[0]["spartan_sumcheck"] == sc.prover_cost(s, 4, 3, 2) + sc.prover_cost(s, 2, 2, 1)
assert all(r["r1cs_dem"] == 16 * r["ban_sao"] and r["n_plonk"] == 32 * r["ban_sao"] for r in rows)
assert all(r["plonk_fft"] > r["spartan_iop"] for r in rows)
