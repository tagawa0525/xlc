---
status: open
labels: [docs, dvc]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# topic-start の後片付けの節と説明の、dvc-relock の後半との食い違い

## 問題

ブランチ `docs/topic-start-drop-library-copy` で、スキル topic-start の「後片付け」の節から、作り直した
`examples/vera/vera.msgpack` を main へコピーする項を消した。そのブランチの `/code-review` の 1 周目が、残した項と
スキルのほかの箇所について次を挙げた。どれも手順のとおりに実行して失敗する再現を書けず、直さずにここに残す。

- 「後片付け」の残した項は、main で `scripts/dvc-checkout.py dvc.yaml` と `scripts/dvc-checkout.py examples/vera/dvc.yaml` を
  実行すると書き、スキル dvc-relock の「後半: マージのあと」を言い換えている。後半はブランチが足したか変えた例ごとの
  `dvc.yaml`（`examples/vera/p6/dvc.yaml` のような入れ子のものを含む）と `.dvc` ファイル（ADR-0222）も戻し、
  `nix develop -c` を前に付ける。残した項はそれより狭く、グローバルの CLAUDE.md の「ほかのスキルの規則は言い換えずに
  スキルの名前と節か項で指す」にも合わない
- 同じスキルの「main の作業ツリー（~/github/xlc）では dvc も計算も実行しない」（計算の前の rebase の節）が、マージの後に
  main で dvc-checkout.py を実行する「後片付け」の項と食い違って読める
- スキルの `description` の「生成物のコピー」「マージ前後の片付け」は、コピーの項を消した後の本文に合わない（本文は
  「生成物を戻す」で、キャッシュから戻す）
- コピーの項を消した理由は、コミットの本文にだけある（ユーザーの 2026-11-06 の判断）。ADR か issue に記録するかは
  決めていない

関係する ADR: なし

## 完了の条件

「後片付け」の項を dvc-relock の後半の節を指す形にし、上の食い違いと description を直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-06: `docs/topic-start-drop-library-copy` の `/code-review`（high）の 1 周目の指摘から起こした
