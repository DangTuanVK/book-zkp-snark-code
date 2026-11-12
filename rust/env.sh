# Nạp môi trường Rust của sách: source rust/env.sh
export RUSTUP_HOME=<THU_MUC_LAM_VIEC>/.rust/rustup
export CARGO_HOME=<THU_MUC_LAM_VIEC>/.rust/cargo
export PATH="$CARGO_HOME/bin:$PATH"
# Thư mục build chung cho 3 crate (đã gitignore)
export CARGO_TARGET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/target"
