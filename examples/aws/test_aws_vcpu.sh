#!/usr/bin/env bash
# スポットの vCPU の残りの読み取り spot_vcpu_remaining と、1 本の上限 vcpu_share の検査（偽の aws を PATH に置く）。flake の check（aws-vcpu）が実行する。
#
# 使い方: examples/aws/test_aws_vcpu.sh [<aws_common.sh のパス>]
set -euo pipefail
common=${1:-$(dirname "$0")/aws_common.sh}
# shellcheck source=aws_common.sh
source "$common"

failures=0
fail() {
  echo "FAIL: $*" >&2
  failures=$((failures + 1))
}
expect() { # <名前> <期待> <実際>
  [[ $2 == "$3" ]] || fail "$1: 期待 '$2'、実際 '$3'"
}

# 偽の aws。クォータは 128、実行中のスポットは「コア数 スレッド数」の行で返す
fake=$(mktemp -d)
trap 'rm -rf "$fake"' EXIT
cat > "$fake/aws" <<'AWS'
#!/usr/bin/env bash
case "$*" in
  *get-service-quota*) echo "${FAKE_QUOTA:-128.0}" ;;
  *describe-instances*) printf '%s' "${FAKE_RUNNING:-}" ;;
  *) echo "偽の aws が知らない呼び出し: $*" >&2; exit 1 ;;
esac
AWS
chmod +x "$fake/aws"
PATH=$fake:$PATH

expect "実行中なし" "128 128 0" "$(spot_vcpu_remaining)"
expect "実行中 16 コア x 2 と 8 コア x 2" "80 128 48" "$(FAKE_RUNNING=$'16\t2\n8\t2\n' spot_vcpu_remaining)"

# 1 本の上限: 残りを起動するジョブの数で割る。新しい実行（入力 8 本）も起動し直しも同じ
expect "8 本は 128 を 16 ずつ" 16 "$(vcpu_share 128 8 16)"
expect "実行中ありの 3 本は 96 を 32 ずつ" 32 "$(vcpu_share 96 3 16)"
expect "1 本は残りの全部" 128 "$(vcpu_share 128 1 16)"
expect "割り切れない 7 本は切り捨て" 18 "$(vcpu_share 128 7 16)"
# 割った値が最小の機種に満たなくても、残りが収まるなら最小の機種までは許す
expect "割った値が最小に満たない" 16 "$(vcpu_share 64 8 16)"
# 残りが最小の機種に満たなければ、割った値のまま（全部の候補が外れて止まる）
expect "残りが最小に満たない" 1 "$(vcpu_share 8 8 16)"

((failures == 0)) || { echo "$failures 件失敗" >&2; exit 1; }
echo "aws_vcpu: ok"
