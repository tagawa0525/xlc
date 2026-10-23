---
status: open
labels: [lattice, depletion]
opened: 2026-10-15
closed:
milestone: M4
blocked-by: []
resolution:
---

# 格子計算の燃焼は燃料のリングの物質しか燃やせず、Pyrex、IFBA、WABA の吸収体を燃焼させられない

## 問題

集合体の燃焼の検証（ADR-0178）の xlc の入力のうち、`[depletion] materials` に Pyrex（2F）、IFBA の被覆（2M、2N）、WABA の吸収体（2N）を
挙げたものは、lattice が燃料のリングの物質だけを燃焼の領域にできないので止まる（1 ステップの試行、2026-10-15、examples/depletion/p2/lattice-<ジョブ>.toml
から燃焼度の点を 0.1 だけにしたもの）:

- 2F: `cell 63: the depleted material is not shielded fuel`（crates/lattice/src/depletion.rs の `depletion_regions` の後。Pyrex は自己遮蔽する燃料でない）
- 2M、2N: `cell 20: region 60 is outside the fuel rings`（IFBA の被覆は燃料のリングの外の領域）

OpenMC の参照解はこれらの吸収体も燃焼させる（openmc_depletion_assembly.py）ので、xlc の領域に吸収体を足すか、決めて比べる範囲を絞るまで、2F、2M、2N の
xlc と OpenMC の比較ができない。2A、2A-density、2A-rods（燃料だけ）と、2O、2P（Gd 燃料のリング。ADR-0182）は走る。

関係する ADR: ADR-0178、ADR-0179

## 完了の条件

燃料でない物質（B-10 などの吸収体、Pyrex の Si と O）を燃焼の領域にでき（反応率は自己遮蔽しない微視的断面積、または吸収体の自己遮蔽を決める）、
2F、2M、2N が走る。IFBA の 10 μm の被覆の中の flux の取り方（薄い領域の MOC の平均）を、OpenMC の燃焼の B-10 の数密度で確かめる。

## 経緯

- 2026-10-15: 集合体の燃焼の入力（ADR-0178）の xlc の試行で分かった
- 2026-10-15: 吸収体を燃焼の領域にできるようにした（ブランチ fix/0139-deplete-absorbers。決め方は ADR-0179）。
  古い p2 のライブラリ（ADR-0177 の前）で 2F と 2M を 1 ステップ（0.1 MWd/kg）流し、走ることと、Pyrex（5 領域）の B-10 が −1.3%、
  IFBA の被覆（19 領域）の B-10 が −2.0〜−2.2% 減ることを確かめた（記録に書かない動作確認。Slurm のジョブ 405）
- 2026-10-21: main（f8f1bd45。ADR-0177 の後の p2 のライブラリ、90 核種）へ rebase し、2F と 2M を同じく 1 ステップ（0.1 MWd/kg）流した。
  走り、Pyrex（5 領域）の B-10 が −1.33〜−1.36%、IFBA の被覆（19 領域）の B-10 が −2.01〜−2.19% 減る（古いライブラリとほぼ同じ。
  記録に書かない動作確認。Slurm のジョブ 412、2F と 2M で合わせて 12 分）。2N も同じく流し、IFBA の被覆（16 領域）の B-10 が
  −1.78〜−2.28%、WABA（4 領域）の B-10 が −1.39〜−1.40% 減る（ジョブ 413）
- 残り（この issue を閉じる条件）: IFBA の 10 μm の被覆の中の中性子束の取り方を、OpenMC の燃焼の B-10 の数密度で確かめる。比べるのは
  p2 の比較のステージ（examples/depletion/p2/dvc.yaml の `compare@2M`、`compare@2N`。compare_depletion_assembly.py）で、AWS の
  m4-assembly の OpenMC の参照解（openmc-2M、openmc-2N）がそろってから流す。2026-10-21 の時点で、8 本のうち完了は 1 本
- 残り（決めること）: 燃焼させる物質の温度が複数のときの輸送する核種の集合（ADR-0179 では決めず、温度が複数なら止める）。決めるときに
  調べること: p2 のライブラリがどの核種をどの温度で持つか。OpenMC の参照解が物質ごとの温度で何を輸送しているか（`nuclides_with_data`
  と温度の扱い）。燃料と吸収体で別の集合（領域の温度ごとに、その温度のデータがある核種）にする案と、すべての温度の共通部分にする案
  （吸収体の温度のデータがない燃料の核種が燃料でも輸送されなくなり得る）の比較。今の入力はどれも 1 つの温度なので、今は止まらない
- 関係して確かめていないこと: 吸収体を燃焼させる重さ（#0137）と、中心のリングを埋める吸収体の棒の UFG の範囲（ADR-0179 の「帰結」）
