"""Chương 5 -- Một phần tử CRS thừa đủ để đổi đầu vào công khai (bản chất lỗi BCTV14).

Groth16 chỉ công bố [L_j(x)/γ]_1 cho cột công khai j <= ℓ. Nếu setup lỡ công bố
thêm K = [L_1(x)/δ]_1 cho cột N, từ một proof hợp lệ cho N ai cũng sửa được
    C' = C - (N' - N)·K
để có proof hợp lệ cho N' bất kỳ, không cần witness.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from snarklib.curve import toy_generator
from snarklib.dsl import mach_phan_tich
from snarklib.groth16 import r1cs_to_qap, setup, prove, verify

E, P, G = toy_generator()
c = mach_phan_tich()[0]
qap = r1cs_to_qap(c.compile_r1cs())
x, w = c.r1cs_assignment({"N": 391, "a": 17, "b": 23})
rng = random.Random(2016)
pk, vk, tau = setup(qap, G, G, rng)
proof = prove(pk, qap, [1] + x + w, rng)          # cùng proof với 05_03

# Setup "lỗi": công bố thêm phần tử của cột N (j = 1) chia cho δ
al, be, ga, de, xx = tau
L1 = be * qap.u[1](xx) + al * qap.v[1](xx) + qap.w[1](xx)
K = int(L1 / de) * G
print("phần tử thừa K = [L_1(x)/δ]_1 =", K)

# Không có K: giữ nguyên proof, chỉ đổi N -> bị bác
print("proof cũ với N' = 389:", verify(vk, [389], proof))

def doi_dau_vao(proof, N_cu, N_moi):
    return proof._replace(C=proof.C + int(qap.F(N_cu - N_moi)) * K)

gia = doi_dau_vao(proof, 391, 389)                # 389 không có witness
print("proof giả cho N' = 389: C' =", gia.C, "-> verify:", verify(vk, [389], gia))
dem = sum(verify(vk, [Nm], doi_dau_vao(proof, 391, Nm)) for Nm in range(449))
print("số N' trong [0, 448] có proof giả hợp lệ:", dem)
assert verify(vk, [389], gia) and not verify(vk, [389], proof) and dem == 449
