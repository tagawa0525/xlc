---
status: open
labels: [openmc, reference]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# The statepoint pruner does not run during the OpenMC run

## 問題

examples/openmc_resume.py の `run_model` は、走行中に 30 秒ごとに `prune(work)` で最新の `KEEP`（2）個のほかの statepoint を消す
つもりだが、走行中は 1 つも消えていない。`openmc.Model.run(cwd=work)` は `change_directory(cwd)` でプロセスの作業ディレクトリを
`work` に移す（OpenMC 0.15.2 の openmc/model/model.py の 705 行）ので、別のスレッドの `prune(work)` が相対パスの `work`
（`run_remote.sh` は `results/<名前>/run` を渡す）を `work/work` として探し、何も見つけない（候補。コードから読んだ）。走行の後の
`prune(work)` は作業ディレクトリが戻ってから走るので効き、最後に 2 個が残る（試走の「最新の 2 個が残る」はこれ）。

観測（2026-10-13、xlc のコミット 6360a566、ADR-0155 の確認、c7a.large、examples/vera/aws-check/openmc-3A.toml の間隔 20）: S3 に
1 分ごとに `--delete` で写した `run/` に、走行中の statepoint が 020 から 180 まですべて残り、完了の後の整理で 200 だけになった。

影響: 本番（間隔 100、2500 世代）では、走行中に約 25 個の statepoint が VM の EBS（40 GB）と S3 の `runs/` にたまる（1 個の大きさは
本番の粒子数で確かめていない）。完了すれば消えるので結果には影響しない。

## 完了の条件

走行中に statepoint が最新の `KEEP` 個に保たれる（例: `run_model` で `work` を絶対パスにする）ことを、テストで確かめる。

## 経緯

- 2026-10-13: ADR-0155 の再開の確認で、S3 の `run/` に古い statepoint が残ることから見つけた
