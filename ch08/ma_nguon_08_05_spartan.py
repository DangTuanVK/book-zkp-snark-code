"""Chương 8 — Mini-Spartan chứng minh mạch ví dụ N = 391 (lst:ch8-spartan).

R1CS 13 ràng buộc × 14 biến (snarklib.dsl) đệm thành 16 × 16 (s = 4 biến).
Cam kết: Ligero (Chương 4) trên bảng W = z với vị trí công khai (1, N) đặt 0.
MLE của ma trận do Verifier tự tính trong O(nnz) — KHÔNG dùng SPARK.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.transcript import Transcript
from snarklib.dsl import mach_phan_tich
from snarklib import spartan, ligero

circ, _ = mach_phan_tich()
r1cs = circ.compile_r1cs()
inst = spartan.pad_r1cs(r1cs)
print("R1CS: m = %d, n = %d -> đệm %d x %d (s = %d); nnz(A, B, C) = %s"
      % (r1cs.m, r1cs.n, inst["size"], inst["size"], inst["s"],
         [len(inst[k]) for k in ("A", "B", "C")]))
x, w = circ.r1cs_assignment({"N": 391, "a": 17, "b": 23})
z = r1cs.z(x, w)
print("z =", [int(v) for v in spartan.pad_z(inst, z)])
print("W (bảng được cam kết) =", [int(v) for v in spartan.witness_table(inst, spartan.pad_z(inst, z))])

T = 8                                         # số cột Ligero được mở
proof = spartan.prove(inst, z, Transcript(b"ch8-spartan"), pcs=ligero, t=T)
info = {}
ok = spartan.verify(inst, [391], proof, Transcript(b"ch8-spartan"), ligero, t=T, info=info)
print("τ =", info["tau"], "; r_x =", info["rx"], "; r_y =", info["ry"])
print("Sum-check #1: s_j(0..3) =", proof["sc1"]["polys"])
print("v_A, v_B, v_C =", info["vABC"], "; e_x =", info["ex"])
print("r_A, r_B, r_C =", info["rABC"], "; Sum-check #2: s_j(0..2) =", proof["sc2"]["polys"])
print("W~(r_y) =", proof["vW"], "; z~(r_y) = v_Z =", info["vZ"],
      "; A~,B~,C~(r_x,r_y) =", info["v123"], "; e_y =", info["ey"])
print("Verifier chấp nhận:", ok)
sz = spartan.proof_size(inst, ligero, t=T)
print("Số vòng sum-check: %d; số cột Ligero mở: %d (ma trận Ligero %d x %d)"
      % (sz["rounds"], sz["queries"], sz["m_ligero"], sz["k_ligero"]))
print("Proof: sum-check %d + %d, (v_A,v_B,v_C) 3, W~(r_y) 1, Ligero %d phần tử trường + %d băm; gốc 1 băm"
      % (sz["sumcheck1"], sz["sumcheck2"], sz["pcs_field"], sz["pcs_hash"]))
print("Tổng: %d phần tử trường + %d giá trị băm 32 byte (%d byte nếu phần tử 32 byte)"
      % (sz["field_total"], sz["hash_total"], 32 * (sz["field_total"] + sz["hash_total"])))

# Witness giả (28, 30) thỏa R1CS trên F_449 nên Spartan CŨNG chấp nhận (lỗi của mạch, Ch.12)
x2, w2 = circ.r1cs_assignment({"N": 391, "a": 28, "b": 30})
p2 = spartan.prove(inst, r1cs.z(x2, w2), Transcript(b"ch8-gia"), pcs=ligero, t=T)
ok2 = spartan.verify(inst, [391], p2, Transcript(b"ch8-gia"), ligero, t=T)
# Witness sai (17, 24): thử 4 nhãn, cho phép lọt tối đa 1 lần (~1/449 mỗi lần)
x3, w3 = circ.r1cs_assignment({"N": 391, "a": 17, "b": 24})
bad = sum(spartan.verify(inst, [391], spartan.prove(inst, r1cs.z(x3, w3), Transcript(b"s%d" % k),
          pcs=ligero, t=T), Transcript(b"s%d" % k), ligero, t=T) for k in range(4))
print("Nghiệm giả (28, 30) được chấp nhận:", ok2, "; witness sai (17, 24) lọt:", bad, "/ 4")
assert ok and ok2 and bad <= 1 and (r1cs.m, r1cs.n, inst["size"]) == (13, 14, 16)
