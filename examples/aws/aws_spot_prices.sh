#!/usr/bin/env bash
# リージョンと AZ ごとのスポット価格（直近の履歴の最新）を価格の昇順に出す。起動の前にどこが安いかを見る参考の情報
# （aws_run.sh は AZ を選べない。docs/environment.md の「外部の計算機（AWS）」）。
#
# 使い方: examples/aws/aws_spot_prices.sh [--type T] [--hours H] [--regions R1,R2] [--running]
#   --type     インスタンスタイプ（既定 c7a.4xlarge）
#   --hours    履歴を引く範囲（直近 H 時間。既定 1。AZ ごとに最新の価格を使う）
#   --regions  引くリージョン（カンマ区切り。既定は describe-regions で有効なものすべて）
#   --running  タグ XlcRun の走っている（pending、running）インスタンスの AZ を表に印で示す
# 組織の SCP が DescribeSpotPriceHistory を UnauthorizedOperation で拒むリージョンは、表の後に 1 行にまとめる。ほかのエラーは止まる。
# オンデマンドの価格は pricing API から取る（取れなければ比の列を省く）。
set -euo pipefail
here=$(dirname "$0")

usage() { sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; }
type=c7a.4xlarge
hours=1
regions=
running=0
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --type) type=${2:?--type に値が要ります}; shift ;;
    --hours) hours=${2:?--hours に値が要ります}; shift ;;
    --regions) regions=${2:?--regions に値が要ります}; shift ;;
    --running) running=1 ;;
    *) echo "知らない引数 $1" >&2; exit 2 ;;
  esac
  shift
done
[[ $type =~ ^[a-z0-9-]+\.[a-z0-9-]+$ ]] || { echo "--type はインスタンスタイプ（c7a.4xlarge など）: $type" >&2; exit 2; }
[[ $hours =~ ^[1-9][0-9]*$ ]] || { echo "--hours は正の整数: $hours" >&2; exit 2; }
[[ $regions =~ ^([a-z0-9-]+(,[a-z0-9-]+)*)?$ ]] || { echo "--regions はカンマ区切りのリージョン: $regions" >&2; exit 2; }
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"

if [[ -n $regions ]]; then
  IFS=, read -r -a region_list <<< "$regions"
else
  # 失敗で止まるよう、先に取る（< <(...) の失敗は set -e で止まらない）
  region_text=$(aws ec2 describe-regions --query 'Regions[].RegionName' --output text)
  mapfile -t region_list < <(tr '\t' '\n' <<< "$region_text")
fi
((${#region_list[@]})) && [[ -n ${region_list[0]} ]] || { echo "リージョンがありません" >&2; exit 1; }
mapfile -t region_list < <(printf '%s\n' "${region_list[@]}" | LC_ALL=C sort -u)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
start=$(date -u -d "$hours hours ago" +%FT%TZ)

# リージョン $1 を引き、$work/$1.{prices,denied,od,running} に書く。SCP の拒否以外の失敗は $work/$1.fail にエラーを残す
fetch_region() {
  local region=$1 err=$work/$1.err
  # start-time の時点で有効だった価格も返るので、価格が変わらなかった AZ も出る
  if ! aws ec2 describe-spot-price-history --region "$region" --instance-types "$type" --product-descriptions Linux/UNIX \
    --start-time "$start" --query 'SpotPriceHistory[].[AvailabilityZone,Timestamp,SpotPrice]' --output text \
    2> "$err" | sort -k1,1 -k2,2r | awk '!seen[$1]++ {print $1, $3}' > "$work/$region.prices"; then
    if grep -q 'service control policy' "$err"; then : > "$work/$region.denied"; else cp "$err" "$work/$region.fail"; fi
    return 0
  fi
  # 価格のないリージョン（タイプの提供がない）は、オンデマンドも実行中も引かない
  [[ -s $work/$region.prices ]] || return 0
  # オンデマンドの価格。pricing API は us-east-1 だけにある。取れなければ空（比の列を省く）
  aws pricing get-products --region us-east-1 --service-code AmazonEC2 --query 'PriceList[]' --output text \
    --filters "Type=TERM_MATCH,Field=instanceType,Value=$type" "Type=TERM_MATCH,Field=regionCode,Value=$region" \
    Type=TERM_MATCH,Field=operatingSystem,Value=Linux Type=TERM_MATCH,Field=tenancy,Value=Shared \
    Type=TERM_MATCH,Field=preInstalledSw,Value=NA Type=TERM_MATCH,Field=capacitystatus,Value=Used 2> "$work/$region.pricing-err" |
    sed 's/.*"OnDemand"//' | grep -o '"USD":"[0-9.]*"' | awk '!n++ {gsub(/[^0-9.]/, ""); print}' > "$work/$region.od" || : > "$work/$region.od"
  [[ -s $work/$region.od ]] || echo "$region: オンデマンド価格を取得できません（比の列を省く）: $(tr '\n' ' ' < "$work/$region.pricing-err")" >&2
  if ((running)); then
    # AZ とタグ XlcRun の値
    # shellcheck disable=SC2016 # JMESPath のバッククォート
    if ! aws ec2 describe-instances --region "$region" --filters Name=tag-key,Values=XlcRun Name=instance-state-name,Values=pending,running \
      --query 'Reservations[].Instances[].[Placement.AvailabilityZone,Tags[?Key==`XlcRun`]|[0].Value]' --output text \
      2> "$err" > "$work/$region.running"; then
      cp "$err" "$work/$region.fail"
    fi
  fi
}

pids=()
for region in "${region_list[@]}"; do
  fetch_region "$region" &
  pids+=($!)
done
failed=0
for i in "${!pids[@]}"; do
  wait "${pids[i]}" || { echo "${region_list[i]}: 問い合わせの処理が失敗しました" >&2; failed=1; }
done
for region in "${region_list[@]}"; do
  if [[ -f $work/$region.fail ]]; then
    echo "$region: $(cat "$work/$region.fail")" >&2
    failed=1
  fi
done
((!failed)) || exit 1

denied=()
empty=()
no_rows=0
echo "インスタンスタイプ $type、Linux/UNIX、直近 $hours 時間（$start 以降）の AZ ごとの最新のスポット価格"
printf '%-8s %-16s %10s %9s %12s  %s\n' 順位 AZ '価格$/時間' 最安比 オンデマンド比 印
rows=
for region in "${region_list[@]}"; do
  if [[ -f $work/$region.denied ]]; then denied+=("$region"); continue; fi
  if [[ ! -s $work/$region.prices ]]; then empty+=("$region"); continue; fi
  od=$(cat "$work/$region.od")
  while read -r az price; do
    [[ -n $az ]] || continue
    mark=
    if ((running)) && grep -q "^${az}[[:space:]]" "$work/$region.running"; then
      mark="実行中 XlcRun=$(awk -v a="$az" '$1 == a {print $2}' FS='\t' "$work/$region.running" | sort -u | paste -sd,)"
    fi
    rows+="$az $price ${od:--} $mark"$'\n'
  done < "$work/$region.prices"
done
if [[ -n $rows ]]; then
  LC_ALL=C sort -k2,2g <<< "${rows%$'\n'}" | awk '
    { az[NR] = $1; price[NR] = $2; od[NR] = $3; mark[NR] = ""; for (i = 4; i <= NF; i++) mark[NR] = mark[NR] (i > 4 ? " " : "") $i }
    END {
      for (i = 1; i <= NR; i++) {
        rel_od = (od[i] == "-" || od[i] + 0 <= 0) ? "-" : sprintf("%.0f%%", 100 * price[i] / od[i])
        m = mark[i]; if (i == 1) m = "最安" (m == "" ? "" : " ") m
        printf "%-8d %-16s %10.4f %8.0f%% %12s  %s\n", i, az[i], price[i], (price[1] > 0) ? 100 * price[i] / price[1] : 0, rel_od, m
      }
    }'
else
  echo "（表にできる価格がありません。--type の綴りと提供されるリージョンを確かめる）" >&2
  no_rows=1
fi
if ((${#denied[@]})); then
  echo "SCP で拒否: $(IFS=,; echo "${denied[*]}")"
fi
if ((${#empty[@]})); then
  echo "価格の履歴なし: $(IFS=,; echo "${empty[*]}")"
fi
exit "$no_rows"
