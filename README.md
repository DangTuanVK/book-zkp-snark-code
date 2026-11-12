# Mã nguồn sách "zk-SNARK chuyên sâu"

Kho mã đi kèm sách **zk-SNARK chuyên sâu -- IOP đa thức, cam kết đa thức và hệ
chứng minh thế hệ mới** (tập 2 của "Chứng minh mà không tiết lộ", TS. Đặng Minh
Tuấn). Mọi con số, bảng và output in trong sách đều được sinh từ mã trong kho
này; kết quả chạy lưu ở `ket_qua/`.

## Cấu trúc thư mục

```
snarklib/        thư viện Python thuần dùng chung cho toàn sách
  field.py       F_p tổng quát, Tonelli–Shanks, phần tử sinh, căn đơn vị; F_{q^2} = F_q[i]/(i^2+1)
  poly.py        đa thức một biến, chia có dư, nội suy Lagrange, Z_H, L_i
  fft.py         NTT/INTT radix-2, coset NTT
  curve.py       đường cong toy y^2 = x^3 + x trên F_5387, sinh G1 tất định
  pairing.py     pairing Tate rút gọn bậc nhúng 2 + ánh xạ méo (Type 1)
  transcript.py  Fiat–Shamir dựa trên SHA-256 có tách miền
  merkle.py      cây Merkle SHA-256
  mle.py         mở rộng đa tuyến tính, eq, fold
  r1cs.py, plonkish.py, ccs.py   cấu trúc ràng buộc + kiểm thỏa
  hang_so.py     tính và ghi ket_qua/hang_so.txt
tests/           test_*.py viết bằng assert thuần (không cần pytest)
chNN/            ma_nguon_NN_KK_<ten>.py -- mã in trong Chương NN (thêm dần)
rust/            lab_r1_groth16_marlin/, lab_r2_nova/, lab_r3_msm_ntt/ (arkworks)
rtl/             montmul/, ntt_bfly/, tb/, dc/ (SystemVerilog, VCS, Design Compiler)
ket_qua/         hang_so.txt, tests.txt, chNN/*.txt
manifest.json    danh mục module, test và (sau này) mã nguồn từng chương
run_all.py       chạy toàn bộ test, sinh lại hằng số
```

## Hằng số xuyên suốt

| Ký hiệu | Giá trị |
|---|---|
| r (trường vô hướng) | 449, r - 1 = 2^6 · 7, phần tử sinh g = 3 |
| ω_8, ω_16, ω_32, ω_64 | 122, 349, 221, 391 (ω_n = 3^((r-1)/n)) |
| q (trường cơ sở) | 5387, q ≡ 3 (mod 4) |
| E | y^2 = x^3 + x, #E = 5388 = 12 · 449, bậc nhúng 2 |
| G1 | (2658, 1349) = 12 · (2, 1603) |
| e(G1, G1) | 1780 + 1636·i ∈ F_{q^2} |
| k1, k2 (coset PLONK, n = 8..64) | 2, 3 |
| Ví dụ | N = 391 = 17 · 23; nghiệm giả (28, 30): 840 ≡ 391 (mod 449) |

Bảng đầy đủ kèm cách tính: `ket_qua/hang_so.txt` (sinh tự động, không sửa tay).

## Trình thông dịch

Python 3.11, **chỉ thư viện chuẩn** (không numpy) để bạn đọc theo được từng
phép tính. Trên máy tác giả dùng môi trường mamba `pqc`:

```bash
PY=python3   # Python ≥ 3.8; tác giả dùng 3.11
```

Python hệ thống 3.6 không được hỗ trợ.

## Cách chạy

```bash
$PY run_all.py              # chạy mọi test, ghi ket_qua/tests.txt và ket_qua/hang_so.txt
$PY run_all.py --no-const   # chỉ chạy test
$PY -m snarklib.hang_so     # chỉ in và ghi bảng hằng số
$PY tests/test_pairing.py   # chạy riêng một file test
```

`run_all.py` thoát với mã khác 0 nếu có test lỗi hoặc hằng số lệch giá trị kỳ vọng.

## Lưu ý an toàn

Các tham số ở đây (r = 449, q = 5387, bậc nhúng 2) chỉ để tính tay và minh họa.
Pairing Type 1 với bậc nhúng 2 hoàn toàn không an toàn ở kích thước thật (tấn
công MOV); hệ thực tế dùng các đường cong Type 3 như BN254, BLS12-381.
