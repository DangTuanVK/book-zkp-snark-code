# Bài thực hành Rust (R1–R3)

Ba crate binary đi kèm sách *zk-SNARK chuyên sâu*. Ở trạng thái hiện tại mỗi
`main.rs` chỉ là **bài chạy thử** (smoke test) xác nhận toolchain và thư viện
build, chạy được; bài lab đầy đủ sẽ viết ở đợt sau.

| Crate | Nội dung chạy thử | Dùng ở |
|---|---|---|
| `lab_r1_groth16_marlin` | Groth16 (ark 0.6) và Marlin (ark-marlin 0.3) cho mạch `a·b = c`, `17·23 = 391` | Ch.5 |
| `lab_r2_nova` | IVC Nova 5 bước của `z ← z² + 3` trên Pallas/Vesta, nén bằng Spartan + IPA | Ch.9 |
| `lab_r3_msm_ntt` | MSM trên G1 và NTT (Radix-2) trên Fr của BLS12-381, cỡ 2^10 | Ch.10 |

## Môi trường

Toolchain cài bằng `rustup` dưới `<THU_MUC_LAM_VIEC>/.rust/`
(không ghi vào `~/.rustup`, `~/.cargo`). Mọi lệnh cần nạp `env.sh` trước:

```bash
source rust/env.sh        # đặt RUSTUP_HOME, CARGO_HOME, PATH, CARGO_TARGET_DIR=rust/target
rustc --version           # rustc 1.99.1 (58cf8f9c3 2026-10-13)
```

`CARGO_TARGET_DIR` trỏ về `rust/target/` (dùng chung cho 3 crate, đã gitignore).

## Build và chạy

```bash
source rust/env.sh
cd rust/lab_r1_groth16_marlin && cargo run --release
cd ../lab_r2_nova             && cargo run --release
cd ../lab_r3_msm_ntt          && cargo run --release
```

Output thật của lần chạy gần nhất: `ket_qua/rust_smoke.txt`; cấu hình máy đo:
`ket_qua/may_do.txt`. Mọi số đo trong sách trích từ hai file này (hoặc bản cập
nhật của chúng).

## Phiên bản crate (theo `Cargo.lock`, được commit)

| Crate | Phụ thuộc chính | Phiên bản |
|---|---|---|
| R1 | `ark-groth16`, `ark-snark`, `ark-relations`, `ark-bls12-381`, `ark-std`, `ark-serialize` | 0.6.0 |
| R1 | `ark-marlin` | 0.3.0 |
| R1 | `ark-poly-commit`, `ark-relations`, `ark-bls12-381`, `ark-poly`, `ark-std` (đổi tên `*-v03`) | 0.3.0 |
| R1 | `blake2` (đổi tên `blake2-v09`) | 0.9.2 |
| R2 | `nova-snark` | 0.76.0 |
| R2 | `ff` · `bincode` · `halo2curves` (gián tiếp) | 0.13.1 · 2.0.1 · 0.9.0 |
| R3 | `ark-ec`, `ark-ff`, `ark-poly`, `ark-bls12-381`, `ark-std` (bật `parallel`) | 0.6.0 |

## Ghi chú kỹ thuật

- **Marlin**: `ark-marlin` mới nhất trên crates.io là 0.3.0, chỉ tương thích
  hệ arkworks 0.3. Crate R1 khai báo song song hai thế hệ arkworks (0.6 cho
  Groth16, 0.3 cho Marlin, đổi tên gói bằng `package = ...`); kiểu dữ liệu hai
  bên **không trộn được**, nên mạch được viết lại bằng API `ark-relations 0.3`
  trong module `marlin`. Cargo cảnh báo `ark-poly-commit 0.3.0` chứa mã sẽ bị
  từ chối ở bản Rust tương lai (future-incompat) — hiện vẫn build được.
- `ark-marlin 0.3` panic (`index out of bounds` trong `ark-poly-commit`) với
  mạch chỉ có **1 ràng buộc**; lặp ràng buộc `a·b = c` 2 lần là đủ
  (`SO_LAN_LAP = 2`). SRS phổ quát dựng với tham số (64, 64, 64).
- Dòng `PC::Check failed` trong output R1 do `ark-marlin` tự in khi kiểm chứng
  thất bại (trường hợp cố ý dùng `c = 392`), không phải lỗi.
- **Nova**: dùng chu trình Pallas/Vesta + IPA (minh bạch). Engine
  `Bn256EngineKZG` (HyperKZG) ở nova-snark 0.76 chặn `setup` ngẫu nhiên trong
  bản build thường; muốn dùng phải nạp file ptau
  (`PublicParams::setup_with_ptau_dir`) hoặc bật feature `test-utils` (không an toàn).
  `prove_step 0` gần như tức thời vì bước đầu chỉ khởi tạo trạng thái.
- Máy đo dùng chung, tải thay đổi; R3 lấy trung vị 20 lần đo. Ghi tải máy
  (`uptime`) kèm mỗi lần đo.
