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

`openmc_depletion.py` の `solve` を変えると `openmc-1C` の重いステージ（1.5 時間以上）が古くなるので、ADR-0174 ではコピーした。
共通の関数（例: オペレーターと積分器を作る `make_operator`）を `openmc_depletion.py` に置き、両方が呼ぶ形が構造的な解決で、
ステージを 1 回再実行する費用が要る。

関係する ADR: ADR-0174、ADR-0178

## 完了の条件

共通の関数にして `openmc-1C` を再実行する（dvc.lock を更新する）か、コピーが一致することを検査するテストを足して、閉じる。

## 経緯

- 2026-10-15: ADR-0174 の実装のレビューから
- 2026-10-15: 集合体の燃焼（ADR-0178）の `openmc_depletion_assembly.py` の `operator` が同じ引数の 3 つ目のコピーになった（`openmc_depletion.py` の `TimedOperator` と `YIELD_CUTOFF` は import している）。共通の関数にするときは 3 つとも置き換える
- 2026-10-15: ADR-0174 の確認 (c)（区間に分けた AWS の実機の実行が手元の分けない実行と統計の範囲で一致するか）は済み。結果は examples/depletion/verification.md の「AWS での燃焼の参照解の確認」。中断からの再開の実機の確認は残る
- 2026-10-27: MOX のピンセルの燃焼（ADR-0206）の `openmc_depletion_mox.py` の `solve` が、同じ引数の 4 つ目のコピーになった
  （`TimedOperator` と `YIELD_CUTOFF` は import している）。共通の関数にするときは 4 つとも置き換える
