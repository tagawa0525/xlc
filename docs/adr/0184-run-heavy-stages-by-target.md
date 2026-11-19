---
type: ADR
status: accepted
date: 2026-10-22
requires: [ADR-0159, ADR-0130]
supersedes: []
superseded-by:
amends: [ADR-0159]
issues: ["#0220", "#0288", "#0268"]
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
  - /dvc-relock の rebase のあとに再実行するステージ。再実行するかは、その結果が要るかと時間の見積もりで、流す側が決める。
    `scripts/dvc-relock.py --run <yaml>:<ステージ>` に挙げたものだけを、2 つ目の例外と同じ確かめ方（`scripts/dvc-repro.py -s`）で実行する
  - 下流のステージを、上流の出力を `dvc.lock` の記録のまま使って実行するとき（ステージの重さによらない）。上流の changed が deps
    （`bin/` の実行ファイルなど）だけで、作業ツリーの上流の outs が `dvc.lock` の記録と同じときに限る。確かめ方は、
    `dvc status --with-deps --json <yaml>:<ステージ>` の出力で、目標の上流のステージがどれも `changed deps` だけを持ち、ほかの種類
    （`changed outs`、`changed command`、`changed checksum`）を持たないこと。目標のステージ自身の `changed deps` は再実行の理由
    なので数えない。`always_changed` のステージ（`references-shell`）がいつも出す `always changed` の項目は、deps と outs の
    差を表さないので数えないが、同じステージの `changed outs` は数える。`references-shell.txt` は devShell が変わっても
    `dvc status` に出ないので、`scripts/check-references-shell.sh` で確かめる（`scripts/dvc-relock.py` も同じものを呼ぶ）。正規形のステージ
    （`python-ast`、`toml-values`、`binary-forms`）は、下流の cmd が出力（正規形）でなく元の .py、TOML、`bin/` を読むので、changed なら
    先に実行して changed でなくする（ADR-0210、ADR-0254、ADR-0255）。例のステージを目標にした `-s` は
    `scripts/dvc-repro.py -s <yaml>:<ステージ>` で実行し、このラッパーがこれらを確かめてから目標だけを実行する（`dvc status` で
    分かる部分は scripts/check-skip-upstream.py の判定、changed の正規形のステージは先に実行、`references-shell.txt` は
    `scripts/check-references-shell.sh`、参照解の目標は `--allow-reference`。ADR-0258）。重いステージは、このラッパーを Slurm に
    投げる。上流は changed のまま `dvc.lock` に残り、記録には上流の出どころ（上流を実行したコミットと分岐点）を下流と別に書く

2 つ目の例外は、背景と案 C が防ぐ危険の一部を、条件を付けて受け入れる。上流の deps が changed なので、上流の outs は
今の deps（`bin/` など）から見れば古く、今の deps で流せば変わりうる。一方、outs は `dvc.lock` の記録と同じなので、下流が
`dvc.lock` に書く deps の md5 は上流の `dvc.lock` の outs の md5 と一致し、`dvc.lock` の上流と下流の記録はつながったまま残る
（どの出力から下流を作ったかを `dvc.lock` が正しく示す）。案 C が `dvc.lock` に書いてしまうと書いた「古い出力から作った md5」の
うち、`dvc.lock` と食い違う md5 は書かれず、書かれるのは上流の deps から見た古さだけである。その古さは、記録の出どころ
（実行したコミットと分岐点。ADR-0130）を上流と下流で別に書くことで読めるようにする。後に `-s` を付けずに `dvc repro` すると、
上流が今の deps で流れ、上流の出力が変われば下流も作り直される（同じバイトなら下流は飛び、記録は `-s` で作ったもののまま残る）。重い下流を
`-s` で作るときは、この作り直しの費用を見積もりに入れる。

## 確認

1 つ目から 3 つ目の項は、確かめるものがない約束事である。`-s` の 1 つ目の例外は、`scripts/test_dvc_relock.py`（挙げたステージだけを
`scripts/dvc-repro.py -s` で実行し、挙げなければ計算のステージを実行しないこと）で確かめる。DVC が目標の上流の changed を実行することは、ADR-0159 の確認 (d) と、2026-10-22 の
lattice@2A のジョブ 421（docs/plans/m4-depletion.md）の実行で見ている。

2 つ目の `-s` の例外の前提は、2026-11-03 に次で見ている（`-s` で実行した時点の作業ツリーの `dvc status` は残っていない）:

- 再実行の後の main a0fd48c0e の作業ツリーで、`dvc status --with-deps --json examples/depletion/dvc.yaml:compare-mox1` は、上流の
  `lattice-mox1` に `{"changed deps": {"bin/lattice": "modified"}}` だけを、`references-shell` に `always changed` を出し、ほかの種類を
  出さない。上流の lattice の outs は、再実行の前と同じ `dvc.lock` の記録のままである（次の項）
- 再実行の結果を書いた examples/depletion/dvc.lock のコミット a0cb82d30（比較を実行したコミットは 374fa19e。examples/depletion/verification.md の
  「計算」）で、3 本の比較のステージの deps のうち変わったのは compare_depletion.py の正規形だけで、上流の出力の md5
  （`lattice-1C-rings1.msgpack` 8e76694f、`lattice-1C-rings4.msgpack` 3b565010、`lattice-mox1-rings1.msgpack` 6fe00e52 と OpenMC の履歴）は
  前後で変わらない

## 帰結

- この ADR は ADR-0159 の実行の形の一文だけを置き換える。ADR 全体の置き換えではないので、ADR-0159 のフロントマターは
  変えない（ADR-0134 と同じ扱い）。ADR-0123 は ADR-0159 が置き換え済みである
- CLAUDE.md の「計算」と docs/plans/m4-depletion.md の「流すだけの手順」は、この ADR に合わせた
- 2 つ目の `-s` の例外（「## 変更」の 2026-11-03）は、CLAUDE.md の「計算」と .claude/skills/dvc-relock/SKILL.md の 3 に合わせた。
  .claude/skills/openmc-reference/SKILL.md の `-s` は外し、確かめ方のうち `dvc status` で分かる部分は scripts/check-skip-upstream.py が
  判定する（#0220。changed の `python-ast`、`toml-values`、`binary-forms` も止める。`references-shell.txt` は判定の外）

## 変更

- 2026-11-20: 決定と理由の `-s` の 1 つ目の例外を、時間の線で決めない形にし、確かめ方を 2 つ目の例外と同じにした。
  `references-shell.txt` の確かめ方の注の指す先を、スキルの手順から `scripts/dvc-relock.py` に替えた。前は「/dvc-relock の
  rebase のあとに数秒で終わり、上流が changed でないと分かっているステージ」「（/dvc-relock の 3 と同じ）」。「数秒」の線は
  過去のセッションが置いたもので、ユーザーの発言に出どころが無い（#0300 の 1）。rebase の後の付け直しは `scripts/dvc-relock.py`
  になり、呼ぶ側が名前で挙げたステージだけを `scripts/dvc-repro.py -s` の経路で実行するので、上流の確かめ方は 2 つ目の例外と
  同じになる（ADR-0258 の決定と理由の「rebase の後の付け直し」と、案 T4）
- 2026-11-20: 決定と理由の `-s` の 2 つ目の例外から、数秒で終わる下流に限る条件を外し、確かめ方を `scripts/dvc-repro.py -s` に
  移した。最後の段落の、数秒に限る理由を、重い下流の作り直しの費用の注意に替えた。前は「数秒で終わる下流のステージを、上流の
  出力を `dvc.lock` の記録のまま使って再実行するとき」「`dvc status` で分かる部分は scripts/check-skip-upstream.py が判定する」
  「数秒の下流に限るのは、上流を今の deps で再実行したときに下流も数秒で作り直せるためである」。lock の上流と下流のつながりは
  上流の outs が lock のとおりかだけで決まり、下流の時間に依らない。数秒の限定はユーザーが求めたものでなく、時間のかかる計算
  ほど上流を作り直さずに流せる経路が要る（ADR-0258 の背景の「上流を飛ばす実行を数秒のステージに限る規則の経緯」と、決定と理由の
  「実行」）。1 つ目の例外の「数秒で終わり」は、ADR-0258 の移行の段 3 で直す
- 2026-11-19: 決定と理由の `-s` の 2 つ目の例外の、先に実行する正規形のステージに `toml-values` と `binary-forms` を足し、判定のスクリプトを
  指した。前は「`python-ast` は、下流の cmd が出力（`.ast`）でなく元の .py を実行するので、changed なら先に実行して changed でなくする
  （ADR-0210）」。この ADR の後に ADR-0254（TOML の値の正規形）と ADR-0255（実行ファイルの正規形）が同じ形の deps を足し、下流が元の
  ファイルを読むので、changed のまま `-s` で下流を流すと `dvc.lock` に古い正規形の md5 が入る（#0220 のレビューの指摘。python-ast は
  examples/vera/openmc_model.py に注釈を足して、判定が止めることを確かめた）
- 2026-11-03: 決定と理由の `-s` の項に、上流の changed が deps だけで outs が `dvc.lock` と同じとき、数秒の下流のステージに
  `-s` を許す例外を足した。前は「`-s` は、/dvc-relock の rebase のあとに数秒で終わり、上流が changed でないと分かっているステージの
  例外だけにする」。きっかけは 2026-11-03 の MOX-1 の記録（main のマージ 54ab434ff、コミット a0cb82d30）で、compare_depletion.py に
  核種を足して `compare-1C@rings1`、`compare-1C@rings4`、`compare-mox1`（数秒）を再実行するとき、上流の `lattice-1C@rings1`、
  `lattice-1C@rings4`、`lattice-mox1` が main の `bin/` と食い違って changed だったが、作業ツリーの lattice の出力は
  `dvc.lock` の記録（上流を実行したブランチの分岐点は main の 3dda20789。examples/depletion/verification.md）どおりで、`-s` で
  再実行した（比較を実行したコミットは 374fa19e、記録の分岐点は main の e04aa5ecc）。例外を許すことはユーザーが決めた（ユーザー 2026-11-03）。残り（openmc-reference の `-s`、確かめ方の機械の判定）は #0220 に書いた。
  比べて採らなかった案:
  - 規則どおり `-s` を付けずに上流も流す: 新しい `bin/` は ADR-0226 で燃焼の履歴の形式（`cells`）が変わった後の版で、
    上流の履歴が作り直され、記録のほかの数値（ADR-0226 の前の履歴から出したもの）と出どころが混ざる
  - 上流を `dvc commit` で付け直す: `bin/` のバイトの差は ADR-0212 の (b) の付け直してよい条件に当たらない（同じ計算かを
    機械で示せない）
