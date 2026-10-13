#!/usr/bin/env bash
# aws_user_data.sh（ADR-0155）が作る user-data の検査。flake の check（aws-user-data）が走らせる。
#
# 使い方: examples/vera/test_aws_user_data.sh [<aws_user_data.sh のパス>]
set -euo pipefail
# shellcheck source=aws_user_data.sh
source "${1:-$(dirname "$0")/aws_user_data.sh}"

failures=0
fail() { echo "FAIL: $*" >&2; failures=$((failures + 1)); }
contains() { grep -qxF -- "$2" <<< "$1" || fail "$3: 行がない: $2"; }
lacks() { ! grep -qF -- "$2" <<< "$1" || fail "$3: あってはならない: $2"; }
line_of() { grep -nxF -- "$2" <<< "$1" | head -1 | cut -d: -f1 || true; }
before() { # <text> <先の行> <後の行> <名前>
  local a b
  a=$(line_of "$1" "$2")
  b=$(line_of "$1" "$3")
  [[ -n $a && -n $b ]] && ((a < b)) || fail "$4: 「$2」が「$3」より前にない"
}

openmc=/nix/store/00000000000000000000000000000000-openmc-0.15.2
python=/nix/store/11111111111111111111111111111111-python3-env
awscli=/nix/store/22222222222222222222222222222222-awscli2-2.35.11
args=(us-east-2 xlc-openmc-test trial openmc-3A trials/aws/openmc-3A.toml)
tail_args=(840 300 lib0123456789abcd nix0123456789abcd "$openmc" "$python" "$awscli")

fresh=$(aws_user_data "${args[@]}" 0 "${tail_args[@]}")
restart=$(aws_user_data "${args[@]}" 1 "${tail_args[@]}")

for name in fresh restart; do
  text=${!name}
  bash -n <<< "$text" || fail "$name: bash の構文"
  [[ $(head -1 <<< "$text") == '#!/bin/bash' ]] || fail "$name: 1 行目が #!/bin/bash でない"
  contains "$text" 'set -euo pipefail' "$name"
  # 最長の時間の shutdown は、ほかの何よりも先に（失敗しうる命令の前に）予約する
  before "$text" 'set -euo pipefail' 'shutdown -h +840' "$name"
  before "$text" 'shutdown -h +840' 'trap finish EXIT' "$name"
  before "$text" 'trap finish EXIT' 'curl -fsSL https://install.determinate.systems/nix | sh -s -- install --no-confirm' "$name"
  contains "$text" 'export HOME=/root AWS_DEFAULT_REGION=us-east-2' "$name"
  contains "$text" 'BUCKET=xlc-openmc-test' "$name"
  contains "$text" 'PREFIX=s3://$BUCKET/runs/trial/openmc-3A' "$name"
  contains "$text" "AWS=$awscli/bin/aws" "$name"
  contains "$text" "\"\$NIX/nix-store\" --realise $awscli > /dev/null" "$name"
  contains "$text" '"$AWS" s3 sync --only-show-errors "s3://$BUCKET/store/nix/nix0123456789abcd" cache' "$name"
  contains "$text" "\"\$NIX/nix\" copy --no-check-sigs --from \"file://\$WORK/cache\" $openmc $python" "$name"
  contains "$text" '"$AWS" s3 sync --only-show-errors "s3://$BUCKET/store/openmc-library/lib0123456789abcd" examples/vera/openmc-library' "$name"
  contains "$text" "printf 'OPENMC_BIN=%s\\nPYTHON=%s\\n' $openmc/bin $python/bin/python3 > tools.env" "$name"
  contains "$text" '  if ((waited >= 300)); then' "$name"
  # done を書かずにジョブが消えたら、最長の時間まで待たずに失敗で終わる
  contains "$text" '  if ! pgrep -f -- " $OUT(/run)?( |$)" > /dev/null; then' "$name"
  # 同期は done と boot.log を除く。finish が同期、boot.log、done の順に置く（done があれば、ほかはそろっている）
  contains "$text" '  "$AWS" s3 sync --only-show-errors --delete --exclude done --exclude boot.log "${filters[@]}" "$WORK/$OUT" "$PREFIX"' "$name"
  before "$text" '      sync_out && { synced=1; break; }' '    "$AWS" s3 cp --only-show-errors /var/log/xlc-job.log "$PREFIX/boot.log" || true' "$name"
  before "$text" '    "$AWS" s3 cp --only-show-errors /var/log/xlc-job.log "$PREFIX/boot.log" || true' \
    '      "$AWS" s3 cp --only-show-errors "$WORK/$OUT/done" "$PREFIX/done" || true' "$name"
  # 同期に失敗したときと SIGTERM（中断）では done を置かない。失敗の終了では終了コードを done に書く
  contains "$text" '    if ((synced && code != 143)); then' "$name"
  contains "$text" '  if ((code != 0 && code != 143)); then' "$name"
  contains "$text" 'until [[ -f $OUT/done && $(cat "$OUT/done") =~ ^[0-9]+$ ]]; do' "$name"
  contains "$text" '    echo "$code" > "$WORK/$OUT/done"' "$name"
  contains "$text" "trap 'exit 143' TERM" "$name"
  # finish の中の失敗で止まると、電源が切れない
  before "$text" 'finish() {' '  set +e' "$name"
  lacks "$text" '"$AWS" s3 cp --only-show-errors "$OUT/done" "$PREFIX/done"' "$name"
  contains "$text" 'READY=0' "$name"
  contains "$text" '  ((READY)) || return 0' "$name"
  contains "$text" "  trap '' PIPE TERM" "$name"
  # 終わりは必ず電源を切る（finish が最後に shutdown -h now）
  contains "$text" '  shutdown -h now' "$name"
done

contains "$fresh" 'bash examples/vera/run_remote.sh trials/aws/openmc-3A.toml' fresh
# 作業ディレクトリが S3 の内容を含むまで同期しない（再開では S3 から戻したあと）
before "$fresh" 'READY=1' 'bash examples/vera/run_remote.sh trials/aws/openmc-3A.toml' fresh
before "$restart" '"$AWS" s3 sync --only-show-errors --exclude boot.log --exclude done "$PREFIX" "$OUT"' 'READY=1' restart
lacks "$fresh" '--restart' fresh
lacks "$fresh" '"$AWS" s3 sync --only-show-errors --exclude boot.log --exclude done "$PREFIX" "$OUT"' fresh
contains "$restart" 'bash examples/vera/run_remote.sh --restart trials/aws/openmc-3A.toml' restart
before "$restart" '"$AWS" s3 sync --only-show-errors --exclude boot.log --exclude done "$PREFIX" "$OUT"' \
  'bash examples/vera/run_remote.sh --restart trials/aws/openmc-3A.toml' restart

# sync_out の振る舞い: 生成した statepoints と sync_out を、aws の代わりに引数を出す関数で走らせる
functions=$(sed -n '/^statepoints() {$/,/^}$/p; /^sync_out() {$/,/^}$/p' <<< "$fresh")
sync_args() { # <done を置くか 0|1> <statepoint の世代>... : sync_out が aws に渡す引数を 1 行に 1 つ
  local put_done=$1 dir g
  shift
  dir=$(mktemp -d)
  mkdir -p "$dir/results/openmc-3A/run"
  for g in "$@"; do touch "$dir/results/openmc-3A/run/statepoint.$g.h5"; done
  ((put_done)) && echo 0 > "$dir/results/openmc-3A/done"
  (
    # shellcheck disable=SC2034 # eval した sync_out が使う
    WORK=$dir OUT=results/openmc-3A PREFIX=s3://b/runs/trial/openmc-3A READY=1
    # shellcheck disable=SC2034 # eval した sync_out が使う
    AWS=aws_stub
    aws_stub() { printf '%s\n' "$@"; }
    eval "$functions"
    sync_out
  )
  rm -rf "$dir"
}
# 走っている間は、最新（書いている途中でありうる）を除き、書き終えたものだけを含める
running=$(sync_args 0 020 040 060)
grep -qxF 'run/statepoint.*' <<< "$running" || fail "sync_out: 走っている間に statepoint を除いていない"
grep -qxF 'run/statepoint.020.h5' <<< "$running" && grep -qxF 'run/statepoint.040.h5' <<< "$running" ||
  fail "sync_out: 書き終えた statepoint を含めていない"
! grep -qxF 'run/statepoint.060.h5' <<< "$running" || fail "sync_out: 最新の statepoint を含めている"
# done のあとは除かない（最後の世代の statepoint を写す）
finished=$(sync_args 1 200)
! grep -qF 'run/statepoint' <<< "$finished" || fail "sync_out: done のあとに statepoint を除いている: $finished"
# statepoint がないときも走る
none=$(sync_args 0)
grep -qxF -- '--delete' <<< "$none" || fail "sync_out: statepoint がないときに aws s3 sync を呼ばない"

# 怪しい引数は user-data を作らずに止まる
rejects() { # <名前> <引数>...
  local name=$1
  shift
  if aws_user_data "$@" > /dev/null 2>&1; then fail "$name: 受け付けた"; fi
}
rejects "再開が 2" "${args[@]}" 2 "${tail_args[@]}"
rejects "最長の分が 0" "${args[@]}" 0 0 300 lib nix "$openmc" "$python" "$awscli"
rejects "store でないパス" "${args[@]}" 0 840 300 lib nix /usr/bin "$python" "$awscli"
rejects "空白を含む入力" us-east-2 xlc-openmc-test trial openmc-3A 'a b.toml' 0 "${tail_args[@]}"
rejects "引数が足りない" "${args[@]}" 0 840

if ((failures)); then
  echo "$failures 件の失敗" >&2
  exit 1
fi
echo "aws_user_data: ok"
