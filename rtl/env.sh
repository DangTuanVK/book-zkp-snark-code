# Môi trường công cụ Synopsys cho rtl/: được các script trong rtl/ tự nạp.
# Cấu hình riêng của từng máy (đường dẫn công cụ, máy chủ license) đặt trong
# rtl/env.local.sh (không đưa vào git), ví dụ:
#   export SNPSLMD_LICENSE_FILE=27000@may-chu-license
#   export VCS=/duong/dan/vcs/bin/vcs
#   export DC=/duong/dan/syn/bin/dc_shell
#   export PY=/duong/dan/python3.11
_RTL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$_RTL_DIR/env.local.sh" ] && . "$_RTL_DIR/env.local.sh"
: "${SNPSLMD_LICENSE_FILE:?Cần đặt SNPSLMD_LICENSE_FILE (máy chủ license Synopsys) trong môi trường hoặc rtl/env.local.sh}"
export SNPSLMD_LICENSE_FILE
VCS=${VCS:-vcs}
DC=${DC:-dc_shell}
PY=${PY:-python3}
# Thư mục thư viện lsi_10k: mặc định suy từ vị trí dc_shell (<syn>/bin/dc_shell -> <syn>/libraries/syn)
if [ -z "$SYN_LIB_DIR" ] && command -v "$DC" >/dev/null 2>&1; then
  SYN_LIB_DIR="$(dirname "$(dirname "$(command -v "$DC")")")/libraries/syn"
fi
export SYN_LIB_DIR
if [ -z "$VCS_HOME" ] && command -v "$VCS" >/dev/null 2>&1; then
  VCS_HOME="$(dirname "$(dirname "$(command -v "$VCS")")")"
fi
export VCS_HOME
