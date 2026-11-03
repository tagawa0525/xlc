---
status: open
labels: [lattice, shielding, performance]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# 間引かない格子の 2A の燃焼で自己遮蔽が機械のメモリを使い切る

## 問題

examples/depletion/dvc.yaml の `lattice-2A`（`lattice-2A.toml`。1 ステップ、`ufg-tolerance` なし、`few-group` あり）を、main f86d7fda8 の
dvc.lock の `bin/lattice`（md5 36667823）で流すと、最初の輸送の自己遮蔽（`shielding:dancoff` の close の後、約 100 秒）で、カーネルの
OOM killer がプロセスを止める（Slurm 644、2026-10-31、exit 137）。

- `journalctl -k`（2026-10-31 20:15:04 JST）: `global_oom`、`Killed process 711567 (lattice) total-vm:94595084kB, anon-rss:88861200kB`。
  Slurm の制限ではない（r995 の `TaskPlugin = (null)`、`JobAcctGatherParams = (null)`。ReqMem 86,000M は要求の値で、超えても止めない。
  2026-11-03 に `scontrol show config` で確かめた）。r995 のメモリは 91 GB
- Slurm 644 の MaxRSS 88,276,164K
- 同じ入力の前の実行（Slurm 457、2026-10-23、ADR-0213 と ADR-0214 の前の bin、16 スレッド）の MaxRSS は 12,959,824K（ジョブの全体。
  `lattice-2A` を含む）

2A の一生の燃焼（examples/depletion/p2/lattice-2A.toml）と history/ の入力（Slurm 582、617、621、629、630、638 など。MaxRSS 20〜26 GB）は
`ufg-tolerance = 1e-4` の間引いた格子で、この入力とは違う。

88 GB に増えたのは bin の版の違いではない（候補が 1 つに絞れるので、bin の二分は流していない）:

- 間の bin（621c10a69 の 22cf3181 と c67424147 の 36667823）の lattice の差は history.rs、residual.rs、main.rs の摂動の解の `h.micro`
  （分岐と履歴の経路）と format のテーブルで、自己遮蔽（shielding.rs、ufg/）と燃焼の経路に変更はない（`git diff 621c10a69 c67424147`）
- この入力は 2026-10-23 の後に流していない（dvc.lock の出力のログの時刻は 2026-10-22T17:46Z）。その間に ADR-0213（自由ガスの
  散乱源の区間を、使用可能メモリ − 2 GiB の 8 割まで持つ）と ADR-0214（スレッドの数の既定を論理 CPU の 32 に）が入った

原因の候補（確かめていない）: ADR-0213 の持つ量の全体は、自己遮蔽の始めの使用可能メモリ（r995 の排他のジョブで 80 GB ほど）から
（80 − 2）× 0.8 ≈ 62 GB で、間引かない格子（ADR-0213 の「検討した案」で 1 回の減速計算 4.8 GB）のピンを 32 スレッドで同時に解くと、
持つ区間がそこまで増える。ADR-0213 の「帰結」は、ほかの作業の配列（格子点、飛程の表、衝突確率）が残りの 2 割と 2 GiB に収まる
必要があると書く。間引かない格子の 32 本の減速計算の作業の配列がそれを超えたことが候補で、持った量と作業の配列の大きさは測っていない
（`UFG held memory` は debug のログ）。

関係する ADR: なし（本文の ADR は accepted で、この issue をフロントマターに挙げていない）

## 完了の条件

- 間引かない格子の 2A（この入力）の 1 ステップの、持った量（`UFG held memory`）と、持たない（`ufg-hold-memory = 0`）ときの MaxRSS を
  測り、上の候補を確かめる
- 既定の持つ量が、ほかの作業の配列と足して機械のメモリを超えない形にする（ADR-0213 を置き換える ADR）か、この入力の扱いを決めて、
  examples/depletion の `lattice-2A` が流れる

## 経緯

- 2026-10-31: bin の変更で changed のステージの流し直し（Slurm 644）で見つけた。流し直しはこの issue を閉じるまで止めている
- 2026-11-03: 前の 20〜25 GB のジョブが別の入力（間引いた格子）で、間の bin の差が自己遮蔽に触れないことを確かめ、bin の二分を流さずに
  ADR-0213 の持つ量を候補にした
