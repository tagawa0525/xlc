---
status: accepted
date: 2026-10-22
requires: [ADR-0159, ADR-0130]
supersedes: []
superseded-by:
amends: [ADR-0159]
issues: []
---

# ADR-0184: 重いステージは、目標のステージを指定して 1 本ずつ Slurm に流す

## 背景

ADR-0159 は ADR-0123 から、DVC のステージを実行する形を引き継いだ。「実行するのは `dvc repro examples/<問題>/dvc.yaml`
（すべては `dvc repro -P`）で、Slurm で流すときは `dvc repro` 全体を流す」。CLAUDE.md の「計算」の「手でステージを選ばない」も
同じ形である。

一方、docs/plans/m4-depletion.md の「流すだけの手順」は、8 本の lattice を各 2〜25 時間かけて、
`dvc repro examples/depletion/p2/dvc.yaml:lattice@2A` のように 1 本ずつ Slurm に流す。全体を 1 つのジョブにすると、
排他の時間が 1 本のジョブに連なって待ち行列を塞ぎ、途中で止まると見積もり（`-t`）の外れも大きく、ジョブ単位の
時間も記録できない。また全体は、比較のステージが参照解（OpenMC の結果）を要るので、取り出しの前は失敗する。

2026-10-22 に、ブランチ chore/0145-relock-after-history-rewrite の /code-review が、「lattice@2A をステージ名で流した」ことを
規則違反として指摘した。規則の文言と計画の手順が食い違っていた。

規則が守りたいのは、`-s`（単独のステージ）で上流を飛ばし、古い上流の出力から作った出力を `dvc.lock` に記録しないことである
（ADR-0130、/dvc-relock）。

判断の軸:

- 古い上流の出力から下流を作らないこと
- 1 本ずつ見積もりと時間を記録できること
- 参照解がまだない例でも、重いステージを流せること

## 検討した案

- 案 A: `dvc repro examples/<問題>/dvc.yaml` 全体を 1 つの Slurm ジョブで流す（ADR-0159 の文面のまま）。上流の漏れはないが、
  8 本が 1 つのジョブに連なり、本数ごとの見積もりと時間が取れない。参照解を待つ比較のステージで失敗する
- 案 B（採用）: 重いステージは目標のステージを指定して（`dvc repro <yaml>:<ステージ>`）1 本ずつ流す。DVC は目標の
  上流のうち changed のものも実行する（`-s` を付けない限り）ので、古い上流の出力から作る誤りは起きない。ジョブごとに
  見積もりと時間を記録できる
- 案 C: `-s` で 1 ステージだけを流す。上流を実行しないので、上流が changed のときに古い出力から作った md5 を
  `dvc.lock` に書く。数秒で終わり、上流が changed でないと分かっているステージ（/dvc-relock）に限れば足りる

## 決定と理由

案 B にする。ADR-0159 の実行の形を次のとおり置き換える（ADR-0159 のほかの決定は変わらない）。

- 例の全部は `dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）。1 分を超えるときの Slurm も、全体を 1 つのジョブにする
- 重いステージを 1 本ずつ Slurm に流すときは `dvc repro <yaml>:<ステージ>`（`-s` を付けない）。目標の上流のうち
  changed のステージは DVC が実行するので、ジョブの時間にそれが入る
- 比較のステージなどが参照解の待ちで失敗する例では、全体を流せないので、目標のステージを指定する
- `-s` は次の 2 つの例外だけにする
  - /dvc-relock の rebase のあとに数秒で終わり、上流が changed でないと分かっているステージ
  - 数秒で終わる下流のステージを、上流の出力を `dvc.lock` の記録のまま使って再実行するとき。上流の changed が deps
    （`bin/` の実行ファイルなど）だけで、作業ツリーの上流の outs が `dvc.lock` の記録と同じときに限る。確かめ方は、
    `dvc status --with-deps --json <yaml>:<ステージ>` の出力で、目標の上流のステージがどれも `changed deps` だけを持ち、ほかの種類
    （`changed outs`、`changed command`、`changed checksum`）を持たないこと。目標のステージ自身の `changed deps` は再実行の理由
    なので数えず、`always_changed` のステージ（`references-shell`。いつも `always changed` を出す）も、出力を
    `dvc.lock` から変えないので数えない。上流は changed のまま `dvc.lock` に残る

規則が守りたいのは、上流の出力のうち `dvc.lock` の記録と違うもの（`-s` で上流を飛ばしたときに作業ツリーに残った古い出力）から
作った md5 を `dvc.lock` に書かないことである。2 つ目の例外の上流は、今の deps（`bin/`）で流せば変わりうる出力ではあるが、
outs は `dvc.lock` の記録と同じなので、下流が `dvc.lock` に書く deps の md5 は上流の `dvc.lock` の outs の md5 と一致し、
`dvc.lock` は上流の記録とつながったまま残る。記録の出どころ（実行したコミットと分岐点。ADR-0130）は、上流と下流で別の
コミットになるので、記録には両方を書く。後に `-s` を付けずに `dvc repro` すると、上流が今の deps で流れ、下流も作り直される。

## 確認

確かめるものがない約束事である。DVC が目標の上流の changed を実行することは、ADR-0159 の確認 (d) と、2026-10-22 の
lattice@2A のジョブ 421（docs/plans/m4-depletion.md）の実行で見ている。

2 つ目の `-s` の例外の前提は、2026-11-03 に次で見ている:

- `dvc status --with-deps --json examples/depletion/dvc.yaml:compare-mox1`（main a0fd48c0e の作業ツリー）は、上流の `lattice-mox1` に
  `{"changed deps": {"bin/lattice": "modified"}}` だけを、`references-shell` に `always changed` を出し、ほかの種類を出さない
- `-s` で再実行したコミット a0cb82d30 の examples/depletion/dvc.lock で、`compare-1C@rings1` の deps の
  `lattice-1C-rings1.msgpack` の md5（8e76694f）は再実行の前後で変わらず、`compare-1C@rings4` の `lattice-1C-rings4.msgpack`（3b565010）も同じ

## 帰結

- この ADR は ADR-0159 の実行の形の一文だけを置き換える。ADR 全体の置き換えではないので、ADR-0159 のフロントマターは
  変えない（ADR-0134 と同じ扱い）。ADR-0123 は ADR-0159 が置き換え済みである
- CLAUDE.md の「計算」と docs/plans/m4-depletion.md の「流すだけの手順」は、この ADR に合わせた
- 2 つ目の `-s` の例外（「## 変更」の 2026-11-03）は、CLAUDE.md の「計算」と .claude/skills/dvc-relock/SKILL.md の 3 に合わせた

## 変更

- 2026-11-03: 決定と理由の `-s` の項に、上流の changed が deps だけで outs が `dvc.lock` と同じとき、数秒の下流のステージに
  `-s` を許す例外を足した。前は「`-s` は、/dvc-relock の rebase のあとに数秒で終わり、上流が changed でないと分かっているステージの
  例外だけにする」。きっかけは 2026-11-03 の MOX-1 の記録（main のマージ 54ab434ff、コミット a0cb82d30）で、compare_depletion.py に
  核種を足して `compare-1C@rings1`、`compare-1C@rings4`、`compare-mox1`（数秒）を再実行するとき、上流の `lattice-1C@rings1`、
  `lattice-1C@rings4`、`lattice-mox1` が main の `bin/` と食い違って changed だったが、作業ツリーの lattice の出力は
  `dvc.lock` の記録（上流を実行したブランチの分岐点は main の 3dda20789。examples/depletion/verification.md）どおりで、`-s` で
  再実行した（比較の分岐点は main の e04aa5ecc）。例外を許すことはユーザーが決めた（ユーザー 2026-11-03）。
  比べて採らなかった案:
  - 規則どおり `-s` を付けずに上流も流す: 新しい `bin/` は ADR-0226 で燃焼の履歴の形式（`cells`）が変わった後の版で、
    上流の履歴が作り直され、記録のほかの数値（ADR-0226 の前の履歴から出したもの）と出どころが混ざる
  - 上流を `dvc commit` で付け直す: `bin/` のバイトの差は ADR-0212 の (b) の付け直してよい条件に当たらない（同じ計算かを
    機械で示せない）
