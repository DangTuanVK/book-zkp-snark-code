# Nạp môi trường Rust của sách: source rust/env.sh
# RUST_ROOT là nơi chứa rustup + cargo (mặc định ~/.rust). Đặt giá trị riêng của
# máy trong rust/env.local.sh (không đưa vào git), ví dụ: export RUST_ROOT=/data/.rust
_RUST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$_RUST_DIR/env.local.sh" ] && . "$_RUST_DIR/env.local.sh"
RUST_ROOT=${RUST_ROOT:-$HOME/.rust}
export RUSTUP_HOME="$RUST_ROOT/rustup"
export CARGO_HOME="$RUST_ROOT/cargo"
export PATH="$CARGO_HOME/bin:$PATH"
# Thư mục build chung cho 3 crate (đã gitignore)
export CARGO_TARGET_DIR="$_RUST_DIR/target"
