---
status: open
labels: [dvc, tooling]
opened: 2026-10-28
closed:
milestone:
blocked-by: []
resolution:
---

# 未実行のステージの出力で dvc checkout が終了コード 255 を返す

## 問題

引数なしの `dvc checkout` は、main の上でも終了コード 255 で終わる。`dvc.lock` に hash がまだ無いステージの出力（examples/depletion/p2 の
本番、examples/vera/endfb70/openmc の 1B と 2B の OpenMC、examples/depletion の openmc-mox1 など）が戻らないためで、/dvc-relock の
前半の手順 1 は、`Checkout failed for following targets:` に並んだ出力ごとに、その `dvc.lock` を `git diff main` で見て、ブランチが
変えた出力でないかを人が確かめるとしている。relock のたびに同じ目視が要り、2026-10-28 の ADR-0211 のマージでも親が 2 回の
コマンドで確かめた。

`dvc checkout --allow-missing` は、hash の無い出力を `WARNING: No file hash info found for '<出力>'. It won't be created.` にして
終了コード 0 で終わる（main 1ed9b10a8 の examples/vera/endfb70/openmc/dvc.yaml で確かめた。--allow-missing なしは 255 で同じ 5 件を並べた）。
ただし、hash はあるがキャッシュに無い出力（/dvc-relock が止まるべき場合）も同じく飛ばすかは確かめていない。

関係する ADR: ADR-0130

## 完了の条件

hash の無い出力（未実行）と、それ以外の失敗（キャッシュに無い、ブランチが変えた lock の出力）を機械で分けて、後者だけで止める
手段（`--allow-missing` の警告の読み分け、または lock に hash の無い出力を除いて checkout するスクリプト）を /dvc-relock の手順 1 に入れ、
目視の段を消す。

## 経緯

- 2026-10-28: /retro で、relock のたびに同じ目視が要ることから書いた
