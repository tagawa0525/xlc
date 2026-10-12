---
status: proposed
date: 2026-10-12
requires: [ADR-0122, ADR-0123, ADR-0124, ADR-0130]
supersedes: []
superseded-by:
issues: ["#0095"]
---

# ADR-0148: 参照解と結果を Zenodo に DOI 付きで公開し、DVC の共有用リモートは別に決める

## 背景

xlc の OpenMC の参照解（VERA 1B〜1D、2A〜2Q、3A/3B、5-2D など）、連続エネルギーの核データのライブラリ、xlc の多群ライブラリと
結果は、手元の `references/`（DVC のキャッシュを含む。ADR-0123、ADR-0130）にだけある。参照解は 1 問題 16〜37 分の OpenMC
計算（ADR-0130）で、他の人が xlc や別のコードの検証に使うには、同じ計算をやり直すしかない。公開すれば役に立つ、という提案
（2026-10-12）を受けて、置き場を決める。

判断の軸:

- 長期に残り、引用できること（DOI、版）。論文や verification.md から版を指せること
- 無料で、規約が科学データの公開に向くこと
- 参照解の 1 問題は約 100 MB（statepoint）、ライブラリは 591 MB（`references/openmc/endfb80-565`、2026-10-12 の `du`）で、
  1 ファイルの上限と総量に収まること
- 既存の DVC の運用（ADR-0123）を壊さず、保守を増やさないこと

置き場の役割は 2 つある。(1) 公開と引用（読むだけの人向け）、(2) 作業用の共有（外部の計算機で作った結果の取り込みなど、
DVC の push と pull）。この ADR は (1) を決め、(2) は必要になるまで作らない方針を決める。

## 検討した案

各サービスの条件（2026-10-12 に公式文書を取得して確認。確認できなかった項目は「未確認」）:

| 項目             | Zenodo                                                                                            | Hugging Face Datasets                                                              | Cloudflare R2                                                          | GitHub Releases                            | Software Heritage                                          |
| ---------------- | ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------- | ------------------------------------------ | ---------------------------------------------------------- |
| 容量             | 1 レコード 50 GB、100 ファイル。一度限りで 200 GB まで増やせる（100 ファイル未満）                | 公開は「best-effort」。最初の数 GB を超えたら有用性が要る。PRO で 10 TB まで       | 無料枠は標準の保存 10 GB-month。バケットの容量は無制限（有料）         | リリースの総量に上限なし                   | 未確認（解説の要約では 300 MB 超のオブジェクトは入れない） |
| 1 ファイル       | 50 GB（レコード全体）                                                                             | 推奨 200 GB 未満、上限 500 GB                                                      | 5 TiB（マルチパート）                                                  | 2 GiB 未満、1 リリースに 1000 ファイルまで | 同上                                                       |
| 保存の期間       | 「リポジトリの寿命」。ホストの CERN の実験計画が少なくとも 20 年先まであり、保証は bit 単位の保存 | 期間の保証は文書にない。DOI をつけると削除、改名、公開範囲の変更に問い合わせが要る | 期間の保証は文書にない（課金が続く間）                                 | 期間の保証は文書にない                     | 保存が目的                                                 |
| DOI              | あり（DataCite）。版ごとに別のレコードと識別子                                                    | 「Generate DOI」で DataCite の DOI（版ごとに再発行）                               | なし                                                                   | なし（Zenodo の連携で付ける）              | SWHID（DOI ではない）                                      |
| 版の管理         | 版は別のレコード。新版は前の版のファイルを取り込める（重複を持たない）                            | git（LFS）。履歴が数千コミットを超えると体験が劣化する                             | なし（オブジェクトの上書き。S3 の版管理 API は未実装）                 | リリースのタグ                             | git の履歴                                                 |
| API              | REST。ファイルの追加は 1 レコード 100 ファイル以下                                                | `huggingface_hub`、git                                                             | S3 互換（一部の API は未実装）                                         | `gh`                                       | 「save code now」                                          |
| S3 互換          | なし                                                                                              | なし                                                                               | あり。出口（egress）の転送料は無料                                     | なし                                       | なし                                                       |
| 料金             | 無料                                                                                              | 公開は無料（上記の条件）                                                           | 無料枠 10 GB、書き込み 100 万回、読み出し 1000 万回 / 月を超えると課金 | 無料                                       | 無料                                                       |
| 規約と科学データ | 全分野の研究成果物が対象。公開するファイルには利用許諾（license）の指定が必須                     | 再利用される有用なデータが前提。利用が見込めないなら別の置き場を勧める             | `r2.dev` は開発用で流量が制限される。本番は独自ドメインが要る          | 制限なし                                   | 主にソースコード                                           |

出典:

- Zenodo の方針: <https://about.zenodo.org/policies/>（保存の期間、50 GB、利用許諾の必須、版は別レコード）
- Zenodo の上限と版: <https://help.zenodo.org/docs/deposit/manage-files/>（50 GB、100 ファイル）、
  <https://help.zenodo.org/docs/deposit/manage-versions/>（新版、ファイルの取り込み）、
  <https://help.zenodo.org/docs/deposit/about-records/>（DataCite の DOI、bit 単位の保存）。増量の
  <https://help.zenodo.org/docs/deposit/manage-files/quota-increase> は、検索結果の要約で 200 GB まで、一度限りと確認した（ページの
  直接取得は 404 で、本文は未確認）
- Hugging Face の容量: <https://huggingface.co/docs/hub/storage-limits>、DOI: <https://huggingface.co/docs/hub/doi>
- Cloudflare R2 の料金: <https://developers.cloudflare.com/r2/pricing/>、上限: <https://developers.cloudflare.com/r2/platform/limits/>、
  S3 互換: <https://developers.cloudflare.com/r2/api/s3/api/>、公開バケット: <https://developers.cloudflare.com/r2/buckets/public-buckets/>
- GitHub Releases: <https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases>
- DVC のリモート: <https://doc.dvc.org/command-reference/remote/add>（S3 互換、HTTP など。Hugging Face の型はない）、
  <https://doc.dvc.org/user-guide/integrations/huggingface>（`dvc get` と `dvc import` による取得だけで、push と pull のリモートではない）
- Software Heritage: 検索結果の要約（<https://arxiv.org/pdf/2310.10295> など）だけで、公式の仕様書は取得していない。
  300 MB 超を入れない、という記述は二次情報で、確認していない

公開と引用の置き場:

- 案 A（採用）: Zenodo。DOI と版があり、保存の期間の根拠（CERN）が明文化され、無料で、研究成果物を対象とする。制約は 1 レコード
  100 ファイルで、参照解は問題ごとのアーカイブにまとめる（ADR-0149）。版は別のレコードで、ファイルは入れた後に変えない前提
- 案 B: Hugging Face Datasets。DOI も git による版もあるが、期間の保証が文書になく、無料の公開は「有用性」の審査つきの best-effort。
  DVC の push と pull の型がない。機械学習の利用者向けの機能（Parquet のビューア）は、MOC の出力には使わない
- 案 C: Cloudflare R2 だけ。出口が無料で S3 互換だが、DOI がなく、期間の保証が文書になく、公開は独自ドメインが要る。引用の置き場に
  ならない
- 案 D: GitHub Releases。このリポジトリには GitHub のリモートがない（CLAUDE.md）。DOI もなく、2 GiB 未満の制約がある
- 案 E: Software Heritage。ソースコードのためのアーカイブで、データの置き場として公式に確認できていない（上記）

DVC のリモート（作業用の共有）:

- 案 F（採用）: 今は作らない。外部の計算機で結果を作る必要が生じたとき、S3 互換のリモートを足す。第一候補は R2（無料枠の
  10 GB、出口が無料）。DVC の S3 の型で R2 の `endpointurl` が動くことは未確認で、足すときに試す。足さなければ、結果は
  `rsync` などでキャッシュに入れて `dvc commit` できる
- 案 G: 今すぐ R2 を DVC のリモートにする。`references/openmc/balance-1B` だけで 20 GB あり（2026-10-12 の
  `du`）、無料枠の 10 GB を超える。公開する分だけを入れても、まだ読む人がいない段階で、アカウントとトークンの保守が先に増える
- 案 H: Hugging Face を DVC のリモートにする。DVC に型がない（上記）。git-lfs を介する方式は DVC の `dvc.lock` の md5 とは別の
  アドレスを持つ二重の管理になる

Zenodo と DVC の関係:

- 案 I（採用）: Zenodo は DVC のリモートにしない（Zenodo に DVC の型がない）。公開は、DVC の出力から公開用のアーカイブを作る
  一方向の書き出しで、データカード（ADR-0149）に `dvc.lock` の md5 を添える。読む人は Zenodo から取ってきて、md5 を照らせる
- 案 J: 公開物を DVC の outs にして、Zenodo への送信も DVC のステージにする。Zenodo の公開は取り消せず（レコードは消されない）、
  `dvc repro` の副作用として呼ぶのは不適切

## 決定と理由

案 A、F、I にする。

- 公開と引用は Zenodo に置く。節目（M3 の完了など、参照解か xlc のライブラリが変わったとき）に新しい版を作り、その版の DOI を
  examples/*/verification.md の「計算」の節に、数値を出したコミットと並べて書く（ADR-0122）
- 公開用の下書き（アーカイブ、メタデータ、ファイルの追加）は手元のツールの 1 コマンド（Zenodo の API を呼ぶ）で作り、DVC のステージにはしない。取り消せない公開（publish）の操作は人が行う。トークンは git の外に置く（#0095）
- DVC の共有用のリモートは、外部の計算機を使い始めるとき（ADR-0149 の出所の記録が計算機の種類を求める）に S3 互換で足す。
  足すときは新しい ADR で決め、R2 を第一候補として動作を確かめる

Hugging Face と R2 は、利用者が増えて Zenodo の取得（50 GB、100 ファイル、増量は一度限り）が足りなくなったときの、副の置き場として
見直す。

## 確認

置き場の条件は上の表の出典を 2026-10-12 に取得して確認した。確認できなかったこと: Zenodo の増量ページの本文、Zenodo の API の
流量制限、R2 の DVC での動作、新しい版に取り込んだ前の版のファイルが 100 ファイルと 50 GB の上限に数えられるか、R2 の無料枠にクレジットカードの登録が要るか、Software Heritage の公式の対象範囲、Cloudflare の
サービス別規約の本文（CDN の条項は、R2 に置いたファイルの配信を妨げないと同社のブログが述べているが、規約の本文は未取得）。
実装の確認は #0095 の完了の条件にある。

## 帰結

- 公開は #0095 で実装する。実装をマージするブランチで、この ADR の確認と帰結を直して accepted にする
- Zenodo のレコードは消せないので、最初の公開の前に、含めるファイルと利用許諾を ADR-0149 のとおりに確認する
- 参照解の版を変えるたびに新しい版のレコードを作る。ファイルは前の版から取り込めるので、変わらないアーカイブは重複しない
- DVC のリモートが無い間、外部の計算機の結果は取り込みの手順が人手になる（#0095 の完了後に必要になれば別の ADR）
