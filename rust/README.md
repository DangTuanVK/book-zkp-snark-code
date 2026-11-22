# Bài thực hành Rust (R1–R3)

Ba crate binary đi kèm sách *zk-SNARK chuyên sâu*. Mỗi `main.rs` là bài lab đầy
đủ được dẫn trong chương tương ứng.

| Crate | Nội dung | Dùng ở | Output đã lưu |
|---|---|---|---|
| `lab_r1_groth16_marlin` | Groth16 (ark 0.6) và Marlin (ark-marlin 0.3) trên BLS12-381 cho mạch ví dụ, nhân rộng tới miền 2^10 .. 2^16 | Ch.5 | `ket_qua/ch05/lab_r1.txt`, `lab_r1.csv` |
| `lab_r2_nova` | IVC Nova trên Pallas/Vesta với hàm bước kiểm phân tích `N = a·b`, nén bằng Spartan + IPA | Ch.9 | `ket_qua/ch09/lab_r2.txt`, `lab_r2.csv` |
| `lab_r3_msm_ntt` | Đo MSM trên G1, NTT trên Fr của BLS12-381 và pairing (BLS12-381, BN254) | Ch.10 | `ket_qua/ch10/r3.txt`, `r3.csv` |

## Môi trường

Toolchain cài bằng `rustup` dưới thư mục `$RUST_ROOT` (mặc định `~/.rust`, đặt riêng trong `rust/env.local.sh`)
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

Số đo trong sách trích từ các file output ở bảng trên; cấu hình máy đo:
`ket_qua/may_do.txt`. `ket_qua/rust_smoke.txt` là output của lần chạy thử
toolchain ban đầu, giữ lại để tham khảo.

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
