---
status: closed
labels: [dvc, openmc, tooling]
opened: 2026-10-28
closed: 2026-10-28
milestone:
blocked-by: []
resolution: ADR-0212
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

`pythonPackages`（flake.nix の検証用の Python のパッケージ）は default の devShell の Python と `referencesPython`（references の
devShell）の両方が使うので、検証用のスクリプトのために足したパッケージも references の devShell の store のパスを変える。

マージの前に見つからなかった。#0182 のブランチは /git-merge と /dvc-relock を通ったが、/dvc-relock の手順 3 は `references-shell` を
「`always_changed` なので `dvc status` に常に出る。対象から外す」とだけ書き、出力の `references-shell.txt` の中身が変わったかを
確かめる段が無い。#0182 の依頼文と CLAUDE.md の「計算」にも、flake.nix の変更が references の devShell の store のパスを変えうる
ことは書いていない。

同じ種類の問題（結果を変えない deps の変更で、重いステージが changed になる）: #0145（deps のコメントだけの変更。`dvc commit` で
付け直したのはその場の判断で、再発を防ぐ仕組みは .py の docstring とコメントを除く ADR-0210 のほかは決めていない）。ADR-0211 のブランチは、
`bin/xs` のバイトだけが変わる refactor（eaecf6144）を、ライブラリの作り直しの 2 本を避けるために戻した（8926b94e4）。

関係する ADR: ADR-0212

## 完了の条件

ツールの結果を変えない devShell の変更で参照解を再実行しない扱いを決め（例: 参照解の deps を OpenMC と NJOY の store のパスだけにする、
または変更の後に `references-shell.txt` だけを `dvc commit` で付け直す）、4 つの例の `references-shell.txt` と dvc.lock を今の devShell に
合わせる。そのあとに、ADR-0211 の 1C の比較（`compare-1C`）を流して記録の 1C の表を更新する。

flake.nix か flake.lock を変えたブランチで、マージの前に references の devShell の store のパスと 4 つの `references-shell.txt` の
食い違いを見つける仕組み（/dvc-relock の段か、flake.nix と flake.lock を対象にした pre-commit の hook）を足す。

## 経緯

- 2026-10-28: ADR-0211 の 1C の燃焼の再測定の前の `dvc status` で見つけた
- 2026-10-28: /retro で、#0182 のマージで見つからなかった理由（/dvc-relock が `references-shell` を対象から外すだけで中身を比べない）を足した
- 2026-10-28: ADR-0212 で、参照のステージの Python のパッケージの並びを検証用から分けた（pyyaml は default と python-tests の devShell
  だけ）。references の devShell は pyyaml の前と同じ `qdaa3j64…` に戻り、4 つの `references-shell.txt` と lock は付け直さずに合った。
  pre-commit の hook `references-shell` が、flake.nix、flake.lock、docs/patches/ の変更で食い違うコミットを止める。`compare-1C` を流し、記録の 1C の表を更新した（9382e2a29 で実行、lock は 3a273a1bc）。
  「pyyaml を足しても OpenMC と NJOY の結果は変わらないはず」は確かめていないまま、確かめずに済む形にした
