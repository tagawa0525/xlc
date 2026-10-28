---
status: open
labels: [dvc, openmc, tooling]
opened: 2026-10-28
closed:
milestone:
blocked-by: []
resolution:
---

# references の devShell の store のパスが examples の references-shell.txt と食い違う

## 問題

main の 6634685a6（pyyaml を検証用の Python に足した。#0182）で、references の devShell の store のパスが
`/nix/store/qdaa3j64…-nix-shell` から `/nix/store/jrv2kvr1…-nix-shell` に変わった（`nix eval --raw
'git+file://…?rev=<コミット>#devShells.x86_64-linux.references.outPath'` を 6634685a6 の親と 6634685a6 で比べた、2026-10-28）。
追跡している 4 つの `references-shell.txt`（examples/depletion、examples/depletion/p2、examples/vera、examples/vera/endfb70/openmc）は
前のパスのまま。

`references-shell` のステージは `always_changed` で、出力の `references-shell.txt` は OpenMC と NJOY を使うステージの deps なので、
目標の上流にこのステージを含む `dvc repro`（例: `examples/depletion/dvc.yaml:compare-1C`）は、OpenMC の参照解（1C の燃焼は数時間）まで
再実行する。ADR-0211 のブランチ（perf/0181-free-gas-erfcx）は、このため 1C の比較のステージを流していない。

pyyaml を足しても OpenMC と NJOY の結果は変わらないはずだが、確かめていない。

関係する ADR: なし

## 完了の条件

ツールの結果を変えない devShell の変更で参照解を再実行しない扱いを決め（例: 参照解の deps を OpenMC と NJOY の store のパスだけにする、
または変更の後に `references-shell.txt` だけを `dvc commit` で付け直す）、4 つの例の `references-shell.txt` と dvc.lock を今の devShell に
合わせる。そのあとに、ADR-0211 の 1C の比較（`compare-1C`）を流して記録の 1C の表を更新する。

## 経緯

- 2026-10-28: ADR-0211 の 1C の燃焼の再測定の前の `dvc status` で見つけた
