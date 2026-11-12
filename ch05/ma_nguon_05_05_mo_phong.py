"""Chương 5 -- Simulator của Groth16 (zero-knowledge): biết tau, không biết witness.

Sim chọn A, B ngẫu nhiên và giải C từ phương trình verify (Groth16, Mục 3.2).
Proof mô phỏng verify được cho cả N không phân tích được trong mạch -- đó là
lý do tau phải bị hủy sau setup.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, setup, verify, simulate

E, P, G = toy_generator()
c = mach_phan_tich()[0]
qap = r1cs_to_qap(c.compile_r1cs())
rng = random.Random(2016)
pk, vk, tau = setup(qap, G, G, rng)

for N in (391, 389):                          # 389 nguyên tố: không có witness
    pi = simulate(qap, tau, [N], G, G, rng)
    print("N = %d: proof mô phỏng A %s, B %s, C %s -> verify: %s"
          % (N, pi.A, pi.B, pi.C, verify(vk, [N], pi)))
    assert verify(vk, [N], pi)
