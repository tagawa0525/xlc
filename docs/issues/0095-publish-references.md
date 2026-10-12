---
status: open
labels: [infra, docs]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 参照解と結果の公開を実装する（Zenodo、DVC のリモート、データカード）

## 問題

OpenMC の参照解、連続エネルギーのライブラリ、xlc のライブラリと結果は手元の `references/` にだけある。公開すれば他の人が
検証に使える、というユーザーの提案（2026-10-12）を、ADR-0148（置き場）にproposed で書いた。公開する記録の単位、内容、出所、ライセンスは、別の ADR として次のブランチで書く。実装は別の作業である。

関係する ADR: ADR-0148（前付けの `issues` にこの issue がある）

## 完了の条件

実装の前に確認すること:

1. ENDF/B-VIII.0 から NJOY で処理した HDF5 の再配布の可否を、一次の出典（NNDC、LANL の配布の方針）で確認する。NNDC の
   ENDF/B-VIII.0 のページには条件の記述がなかった（2026-10-12）。確認できなければ、ライブラリは作り方だけを公開する
   （公開する記録の ADR の案）
2. リポジトリに LICENSE がなく Cargo.toml に license もないので、xlc のコードのライセンスを ADR にして決める
3. 公開するアーカイブを作って総量を測り、50 GB を超えるならレコードを分けるか、Zenodo の増量（200 GB まで、一度限り）を申請する。
   Zenodo の API の流量制限と、増量ページの本文（直接取得が 404 だった）を確認する

実装の手順:

1. アカウント: Zenodo のアカウントと ORCID。まず <https://sandbox.zenodo.org> で試し、公開の版（取り消せない）は確認の後に作る
2. トークン: Zenodo のトークンは git に入れない。置き場は `~/.config/xlc/zenodo-token`（権限 600）か環境変数とし、スクリプトは
   どちらも無ければエラーで止まる（既定値を持たない）
3. 道具: DVC の出力から公開用のアーカイブ（問題ごとの `.tar.zst`）と `README.md`、`provenance.json` を作り、Zenodo の API で
   下書きの作成、ファイルの追加、メタデータの設定までを行う。公開（publish）は人が押す。出所の項目が欠けたら止まる
4. 最初に公開する対象: 軽い順に (3) xlc のライブラリと結果、(2) 参照解の VERA 1B と 2A（1 問題の流れを試す）、そのあと残りの問題
5. 記録の更新: examples/*/verification.md の「計算」の節に、数値を出したコミットと並べて版の DOI を書く（ADR-0122、
   スキル verification-record）
6. DVC のリモート（外部の計算機を使い始めるとき）: S3 互換のリモートを足す。R2 を第一候補に、`endpointurl` つきの DVC の S3 の型が
   動くか、10 GB の無料枠に収まるかを確かめてから、新しい ADR で決める。外部の計算機の結果にも同じデータカードを付け、同じ種と
   入力で機種が違うときの k の差を測る
7. 実装をマージするブランチで、ADR-0148 と、公開する記録の ADR の確認と帰結を実装に合わせ、accepted にする

## 経緯

- 2026-10-12: ユーザーの提案で起こした。Zenodo、Hugging Face、Cloudflare R2、GitHub Releases、Software Heritage の条件を
  公式文書で比べた（ADR-0148）
