---
status: proposed
date: 2026-10-12
requires: []
supersedes: []
superseded-by:
issues: []
---

# ADR-0151: Run OpenMC reference solutions on AWS EC2 c7a spot instances

## 背景

OpenMC の参照解は 1 本が 100 万粒子 x 2300〜2500 世代で、r995（Ryzen 9 9950X、16 物理コア）の 16 スレッドで 1 世代 9〜10 秒、1 本 5〜6 時間かかる（examples/vera/verification.md の「参照解」）。いま要る 4 本（3A、3B、5B-2D、5C-2D）は約 350 物理コア時間（4 本 x 16 コア x 5.5 時間）。そのあと Problem 4（3x3 の 3D）と Problem 5（全炉心の 3D）で 1 本あたり 5〜20 倍の計算が数本〜十数本あり、M3 全体で 1000〜3000 物理コア時間の見込み。r995 だけでは 4 本が逐次で約 22 時間、M3 全体で 4〜12 日かかる。

持ち込むものは、OpenMC と Python API（nix の `.#openmc`、`.#openmc-python`）、核データ 591 MB、入力とスクリプト。結果は数十 MB。

価格はすべて 2026-10-12 に集計サイトで見たスナップショットで、各社の公式の価格表では確認していない。スポットは変動する。

## 決定の要因

- 費用は物理コア時間あたりで比べる（vCPU が SMT のスレッドか物理コアかで 2 倍違う）
- x86_64 の nix のクロージャをそのまま送れる（再ビルドしない）
- 中断のやり直しが小さい
- 準備と権限の手間が小さい

## 検討した案

| 案                          | 物理コア              | OD $/h | スポット $/h | 物理コア時間あたり（OD / スポット） | M3 2000 物理コア時間（OD / スポット） |
| --------------------------- | --------------------- | ------ | ------------ | ----------------------------------- | ------------------------------------- |
| AWS c7a.16xlarge, us-east-1 | 64（vCPU = 物理コア） | 3.284  | 0.989〜1.447 | $0.051 / $0.015〜0.023              | 約 $103 / 約 $31〜45                  |
| AWS c7a.16xlarge, 東京      | 64                    | 4.134  | 1.629        | $0.065 / $0.025                     | 約 $129 / 約 $51                      |
| GCP c4d-highcpu-32, 東京    | 16（32 vCPU は SMT）  | 1.631  | 0.470〜0.558 | $0.102 / $0.029〜0.035              | 約 $204 / 約 $59〜70                  |
| GCP c3d-standard-60, ソウル | 30                    | 3.498  | 1.011        | $0.117 / $0.034                     | 約 $233 / 約 $67                      |
| Azure D64as_v6, East US     | 32 か 64（未確認）    | 2.906  | 0.609        | 32 コアなら $0.091 / $0.019         | 32 コアなら 約 $182 / 約 $38          |
| Hetzner Cloud CCX43         | 8（16 vCPU）          | 0.4956 | なし         | $0.062                              | 約 $124                               |
| GitHub Actions 64 コア      | 64 vCPU               | 9.72   | なし         | 約 $0.15〜0.30                      | 約 $300〜600                          |

出典（2026-10-12 に取得）:

- AWS: [C7a 公式](https://aws.amazon.com/ec2/instance-types/c7a/)（「Every vCPU on a C7a instance is a physical CPU core.」）、価格は [devzero](https://www.devzero.io/instances/aws/c7a.16xlarge)、[Holori](https://calculator.holori.com/aws/ec2/c7a.16xlarge)、[Vantage](https://instances.vantage.sh/aws/ec2/c7a.16xlarge)
- GCP: [Holori c4d 東京](https://calculator.holori.com/gcp/vm/c4d-highcpu-32/asia-northeast1)、[devzero c4d-highcpu-32](https://www.devzero.io/instances/gcp/c4d-highcpu-32)、[Holori c3d ソウル](https://calculator.holori.com/gcp/vm/c3d-standard-60/asia-northeast3)。vCPU が SMT のスレッドであることは Google の文書で確認していない（一般知識）
- Azure: [Holori](https://calculator.holori.com/azure/vm/standard-d64as-v6)、[仕様](https://learn.microsoft.com/en-us/azure/virtual-machines/sizes/general-purpose/dasv6-series)（vCPU が物理コアかは書かれていない）
- Hetzner: [Spare Cores CCX43](https://sparecores.com/server/hcloud/ccx43)、値上げは [wz-it](https://wz-it.com/blog/hetzner-preiserhoehung-juni-2026-cpx-ccx-alternativen/)、[privatedevops](https://privatedevops.com/news/hetzner-june-2026-cloud-price-increase-what-to-do)、専用機は [prismix](https://prismix.dev/news/e099eaeb0874)
- GitHub: [Actions runner pricing](https://docs.github.com/en/enterprise-cloud@latest/billing/reference/actions-runner-pricing)
- Oracle: [Always Free Resources](https://docs.oracle.com/en-us/iaas/Content/FreeTier/resourceref.htm)
- 中古 EPYC: [eBay の 7702 の出品](https://www.ebay.com/p/4040735398) ほか
- 電気料金: [Selectra](https://selectra.jp/energy/electricity/price-per-kwh)
- 取り出しの無料枠: AWS re:Post の回答が引く料金表の記述（[例](https://www.repost.aws/ja/questions/QUzBhYwSqEQly68_-NdYvDVQ/aws-data-transfer-free-tier)）。公式の料金表は確認していない

### 借りない（r995 だけ）

2000 物理コア時間は r995 で約 125 時間（約 5 日）。電気代は 31 円/kWh（26〜41 円）、300 W の仮定で約 37.5 kWh、約 1,200 円。費用は最小。待ち時間が長いのと、ほかの計算（Slurm の待ち行列）と r995 を取り合う点が欠点。

### GCP c4d / c3d の Spot

単価は AWS と同程度だが、vCPU が SMT のスレッドなので物理コアあたりは約 2 倍。16 物理コアに 32 vCPU が要る。us-central1 の c3d-standard-60 のスポットは $0.4338/h（30 コア、約 $0.0145/物理コア時間。[Holori](https://calculator.holori.com/gcp/vm/c3d-standard-60)）で AWS と並ぶが、その地域のオンデマンドは取得していない。

### Azure

スポットは安いが、Dasv6 の vCPU が物理コアか確認できなかった。確認するまで選ばない。

### Hetzner

Cloud の CCX は 2026-06-15 の値上げで +113〜169%（CCX63 は €374.49 から €853.49）。CCX43 は $0.062/物理コア時間でスポットがなく、AWS のスポットの約 3〜4 倍。専用機（AX162-1 は €612.30/月とされるが新旧は未確認）は月額で、42 時間分の計算でも 1 か月分を払う。サーバーオークションは現在価格を取得できなかった。

### Oracle の Always Free

ARM 4 OCPU / 24 GB で無料。ただし 4 コアなので 1 本が 4 倍以上かかり、aarch64 の OpenMC を再ビルドする。2026 年に枠を半分にしたという記事があるが、公式で確認できなかった。

### GitHub Actions

64 コア Linux が $0.162/分（約 $9.72/h）で、物理コアあたり最も高い。

### 中古の EPYC 7702（64 コア Zen2）

CPU 単体の asking price が $465〜740、Supermicro AS-1014S-WTRT + 7702P の出品が $2,300。クラウドのスポット換算（約 $0.0155/物理コア時間）で $2,300 は約 148,000 物理コア時間分で、M3 の約 70 倍。電気代（上記、約 1,200 円）は購入費に比べて無視できる。

### OVHcloud、Scaleway、DigitalOcean、日本の大学・HPC の共用

32 vCPU 級の現在価格を取得できなかった、または調べていない。

## 決定と理由

OpenMC の参照解のうち、r995 で流すと待ち時間が長くなるものは、AWS EC2 の c7a のスポットで流す。理由は次のとおり。

- c7a は vCPU が物理コアで SMT がなく、物理コアあたりのスポット価格が検討した案の中で最も低い（us-east で約 $0.015〜0.023）
- x86_64 なので `nix copy --to ssh-ng://` でクロージャを再ビルドなしに送れる
- 中断は 1 本のやり直しで済む（約 6 時間、$2 以下）。c7a.4xlarge（16 vCPU）を 1 本に 1 台ずつ割り当てると、中断の影響が 1 本に閉じる。statepoint からの再開は、今のスクリプトが対応していない。途中の statepoint を残す（`statepoint.batches`）ことと、VM の終了で消えない保存先は、ツールの側で決める（未実装。feat/0073-openmc-library のブランチで作業中）。それまでは中断で 1 本をやり直す
- リージョンは単価の低い us-east-1 / us-east-2 を既定にする。バッチ計算なので日本からの遅延は影響しない。核データ 591 MB の送信は無料（受信側）

借りる前提は、r995 の待ち時間を縮める必要があるときだけ。借りなくてよい計算（短いもの、数分のビルド）は手元で流す。

### 帰結

- 外部の計算機で流した結果は、k と統計誤差の参照解として使う。時間は記録しない（1 世代あたりの時間は r995 の測定だけ）
- examples/*/verification.md の「計算」の節に、計算機の種類（c7a.4xlarge など）とコア数を書く。Zen4（c7a）と Zen5（r995）で k が統計誤差の中で一致することは、初回に r995 と同じ入力の 1 本で確かめてから使う（浮動小数点の演算順序で乱数の列は同一にならないので、差が k の標準誤差の 2 倍以内かで見る。基準はその時に決める）
- 費用の見積もり（us-east のスポット）: いまの 4 本が約 $5〜8（352 物理コア時間 x $0.015〜0.023）、M3 の残りが約 $25〜38（1648 物理コア時間）。単価は c7a.16xlarge の値で、c7a.4xlarge のスポット価格と空きは別に動く。4 本並列には 64 vCPU のスポットクォータが要る。取り出しは AWS の無料枠（月 100 GB まで）に収まる（結果は数十 MB）。EBS とスナップショットの少額費用がかかる
- 準備と手順は docs/environment.md の「外部の計算機（AWS）」。起動、回収のツールは未実装（feat/0073-openmc-library のブランチで作業中）

### 確認

確かめていないことがある。結論に影響するものを先に並べる。

- 新規アカウントのスポット vCPU クォータの既定（引き上げ申請が要るかもしれない）
- Zen4（c7a）と Zen5（r995）の OpenMC での 1 コアあたりの速さ。上の時間の見積もり（約 5.5 時間）は r995 の値で置いている
- Azure Dasv6 の vCPU が物理コアか
- GCP の vCPU が SMT のスレッドであること（Google の文書）
- 価格はすべて集計サイトのスナップショット（2026-01〜09 の混在）で、公式の価格表ではない。スポットは変動する
- Hetzner サーバーオークションの現在価格、OVHcloud、Scaleway、DigitalOcean、日本の大学・HPC
- 実際に借りて 1 本流したときの費用と所要時間（流したら、この ADR の帰結ではなく verification.md に書く）
