---
status: open
labels: [aws, depletion]
opened: 2026-10-15
closed:
milestone:
blocked-by: []
resolution:
---

# Share the operator setup between openmc_depletion.py and the resumable run

## 問題

`examples/depletion/openmc_depletion_resume.py`（ADR-0174）の `solve_step` は、`TimedOperator` と `CECMIntegrator` の引数
（`normalization_mode`、`fission_yield_mode`、`fission_yield_opts`、`reaction_rate_mode`、`solver`、`timestep_units`）を
`openmc_depletion.py` の `solve` からコピーしている。`solve` だけを変えると、区間に分けた AWS の実行が参照解の `openmc-1C` と違う設定で
流れ、比較が黙って崩れる。両方を検査するテストはない。

`openmc_depletion.py` の `solve` を変えると `openmc-1C` の重い段階（1.5 時間以上）が古くなるので、ADR-0174 ではコピーした。
共通の関数（例: オペレーターと積分器を作る `make_operator`）を `openmc_depletion.py` に置き、両方が呼ぶ形が構造的な解決で、
段階を 1 回再実行する費用が要る。

関係する ADR: ADR-0174

## 完了の条件

共通の関数にして `openmc-1C` を再実行する（dvc.lock を更新する）か、コピーが一致することを検査するテストを足して、閉じる。

## 経緯

- 2026-10-15: ADR-0174 の実装のレビューから
