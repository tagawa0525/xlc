---
status: accepted
date: 2026-10-22
requires: [ADR-0159, ADR-0130]
supersedes: []
superseded-by:
amends: [ADR-0159]
issues: ["#0220"]
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
    なので数えない。`always_changed` のステージ（`references-shell`）がいつも出す `always changed` の項目は、deps と outs の
    差を表さないので数えないが、同じステージの `changed outs` は数える。`references-shell.txt` は devShell が変わっても
    `dvc status` に出ないので、`scripts/check-references-shell.sh` で確かめる（/dvc-relock の 3 と同じ）。`python-ast` は、
    下流の cmd が出力（`.ast`）でなく元の .py を実行するので、changed なら先に実行して changed でなくする（ADR-0210）。上流は
    changed のまま `dvc.lock` に残り、記録には上流の出どころ（上流を実行したコミットと分岐点）を下流と別に書く

2 つ目の例外は、背景と案 C が防ぐ危険の一部を、条件を付けて受け入れる。上流の deps が changed なので、上流の outs は
今の deps（`bin/` など）から見れば古く、今の deps で流せば変わりうる。一方、outs は `dvc.lock` の記録と同じなので、下流が
`dvc.lock` に書く deps の md5 は上流の `dvc.lock` の outs の md5 と一致し、`dvc.lock` の上流と下流の記録はつながったまま残る
（どの出力から下流を作ったかを `dvc.lock` が正しく示す）。案 C が `dvc.lock` に書いてしまうと書いた「古い出力から作った md5」の
うち、`dvc.lock` と食い違う md5 は書かれず、書かれるのは上流の deps から見た古さだけである。その古さは、記録の出どころ
（実行したコミットと分岐点。ADR-0130）を上流と下流で別に書くことで読めるようにする。後に `-s` を付けずに `dvc repro` すると、
上流が今の deps で流れ、下流も作り直される。数秒の下流に限るのは、上流を今の deps で再実行したときに下流も数秒で作り直せる
ためである。

## 確認

1 つ目から 3 つ目の項と `-s` の 1 つ目の例外は、確かめるものがない約束事である。DVC が目標の上流の changed を実行することは、ADR-0159 の確認 (d) と、2026-10-22 の
lattice@2A のジョブ 421（docs/plans/m4-depletion.md）の実行で見ている。

2 つ目の `-s` の例外の前提は、2026-11-03 に次で見ている（`-s` で実行した時点の作業ツリーの `dvc status` は残っていない）:

- 再実行の後の main 16167a058 の作業ツリーで、`dvc status --with-deps --json examples/depletion/dvc.yaml:compare-mox1` は、上流の
  `lattice-mox1` に `{"changed deps": {"bin/lattice": "modified"}}` だけを、`references-shell` に `always changed` を出し、ほかの種類を
  出さない。上流の lattice の outs は、再実行の前と同じ `dvc.lock` の記録のままである（次の項）
- 再実行の結果を書いた examples/depletion/dvc.lock のコミット e242fb780（比較を実行したコミットは 5381ca6d。examples/depletion/verification.md の
  「計算」）で、3 本の比較のステージの deps のうち変わったのは compare_depletion.py の正規形だけで、上流の出力の md5
  （`lattice-1C-rings1.msgpack` 8e76694f、`lattice-1C-rings4.msgpack` 3b565010、`lattice-mox1-rings1.msgpack` 6fe00e52 と OpenMC の履歴）は
  前後で変わらない

## 帰結

- この ADR は ADR-0159 の実行の形の一文だけを置き換える。ADR 全体の置き換えではないので、ADR-0159 のフロントマターは
  変えない（ADR-0134 と同じ扱い）。ADR-0123 は ADR-0159 が置き換え済みである
- CLAUDE.md の「計算」と docs/plans/m4-depletion.md の「流すだけの手順」は、この ADR に合わせた
- 2 つ目の `-s` の例外（「## 変更」の 2026-11-03）は、CLAUDE.md の「計算」と .claude/skills/dvc-relock/SKILL.md の 3 に合わせた

## 変更

- 2026-11-03: 決定と理由の `-s` の項に、上流の changed が deps だけで outs が `dvc.lock` と同じとき、数秒の下流のステージに
  `-s` を許す例外を足した。前は「`-s` は、/dvc-relock の rebase のあとに数秒で終わり、上流が changed でないと分かっているステージの
  例外だけにする」。きっかけは 2026-11-03 の MOX-1 の記録（main のマージ d2cf5859d、コミット e242fb780）で、compare_depletion.py に
  核種を足して `compare-1C@rings1`、`compare-1C@rings4`、`compare-mox1`（数秒）を再実行するとき、上流の `lattice-1C@rings1`、
  `lattice-1C@rings4`、`lattice-mox1` が main の `bin/` と食い違って changed だったが、作業ツリーの lattice の出力は
  `dvc.lock` の記録（上流を実行したブランチの分岐点は main の dcc9036a2。examples/depletion/verification.md）どおりで、`-s` で
  再実行した（比較を実行したコミットは 5381ca6d、記録の分岐点は main の c0009a046）。例外を許すことはユーザーが決めた（ユーザー 2026-11-03）。
  比べて採らなかった案:
  - 規則どおり `-s` を付けずに上流も流す: 新しい `bin/` は ADR-0226 で燃焼の履歴の形式（`cells`）が変わった後の版で、
    上流の履歴が作り直され、記録のほかの数値（ADR-0226 の前の履歴から出したもの）と出どころが混ざる
  - 上流を `dvc commit` で付け直す: `bin/` のバイトの差は ADR-0212 の (b) の付け直してよい条件に当たらない（同じ計算かを
    機械で示せない）
  残り（openmc-reference の `-s`、確かめ方の機械の判定）は #0220
