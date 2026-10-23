---
status: accepted
date: 2026-10-23
requires: [ADR-0097, ADR-0122, ADR-0130, ADR-0159, ADR-0173, ADR-0181, ADR-0184, ADR-0187, ADR-0193, ADR-0194, ADR-0195]
supersedes: []
superseded-by:
issues: ["#0165"]
---

# ADR-0196: 格子計算の燃焼と分岐を、点と状態ごとのチェックポイントから再開する

## 背景

格子計算の燃焼（[ADR-0173](0173-deplete-lattice-pins-by-ce-cm.md)）は 40〜70 点で 35 分〜3 時間（examples/depletion/verification.md の
「時間」。2A 34 分 33 秒、2M 2 時間 58 分 12 秒）、燃焼度の軸を持つ分岐のテーブル（[ADR-0187](0187-tabulate-burnup-branches-from-history.md)）は数十〜数百の
状態で 1 状態 25〜78 秒かかる。履歴（[ADR-0181](0181-write-depletion-pin-powers-into-history.md)）もテーブル（[ADR-0097](0097-tabulate-branches-on-a-full-grid.md)）も計算の終わりに 1 度だけ
書くので、途中で止めると全部を失う。r995 の Slurm は排他のジョブを 1 本ずつ流し、この決定を始めた時点では割り込みがなく
（PreemptMode は OFF）、急ぐ短い計算を前に出すには、走っているジョブを `scancel` で止めて投げ直すしかなかった。ユーザーの要望
（2026-10-23）: 「リスタート計算機能がほしい。長い計算を途中でぶった切って、急ぐ計算の優先度を上げられるため」。同じ日のうちに、
`urgent` のパーティションが走っているジョブを止めて（SUSPEND。`scontrol show config` の PreemptMode は GANG,SUSPEND）先に走る
ようになり（nixfiles の slurm-run の `--urgent`）、急ぐ計算を前に出すのに `scancel` は要らなくなった。チェックポイントは、ジョブを
取り消したとき、ノードやプロセスが止まったとき（停電、OOM killer、時間切れ）に、書き終えた点から続けるために残る。

`dvc repro` はステージを実行する前に outs を消す。燃焼と分岐のステージ（examples/depletion/p2/dvc.yaml の `lattice` など）は
`bin/lattice <入力>` を実行し、出力と標準出力、ログを outs にする（[ADR-0130](0130-drive-reference-results-with-dvc.md)、[ADR-0159](0159-depend-stages-on-built-binaries.md)）。

判断の軸:

- 再開して完成した結果が、一気に流した結果とビットで一致すること（ADR-0122 の「コミットから出し直せる」を崩さない）
- 入力、読むファイル、コードを変えたあとに、古い途中の結果から続けないこと
- 止めたジョブを同じコマンド（`dvc repro`、`slurm-run.sh` の同じ引数）で投げ直せば再開し、dvc.yaml、`dvc.lock`、入力を変えないこと
- 書く費用が、点（輸送 2 回、2A で約 50 秒）や状態（25〜78 秒）に比べて小さいこと

## 検討した案

書く単位と書き方:

- 案 W1（採用）: 点（状態）を解き終えるたびに、それまでの全体を一時ファイルに書き、`fsync` してから rename で置き換える。殺されても
  ノードが止まっても、前の版か新しい版のどちらかが残る
- 案 W2: 点ごとに追記する（点ごとのファイル、または追記できる形式）。書く量は 1 点分で済むが、読み手と形式が増え、途中で切れた
  追記を見分ける必要がある。W1 の書く費用は「確認」の測定で足りる大きさなので採らない
- 案 W3: 数点おきに書く。失う量が増え、間隔に根拠のない既定値が要る

燃焼の再開点（CE/CM の 1 ステップは、始まりと中点の輸送 2 回。点はステップの始まり）:

- 案 P1（採用）: ステップを解き終えた時点で、点 0..j と、次の点の組成 n_{j+1}（そのステップの修正子で進めたもの）を書く。再開は
  n_{j+1} の輸送から始め、輸送の回数（`solves`）も一気に流したものと同じになる。失うのは解いている途中の 1 ステップまで
- 案 P2: 点 j の組成から再開し、点 j の輸送を解き直す。履歴の点だけで足りるが、再開のたびに輸送が 1 回余計になる

置き場（DVC との整合）:

- 案 L1（採用）: 出力と別のファイル（入力の拡張子を `.checkpoint` にしたパス。`lattice-2A.toml` なら `lattice-2A.checkpoint`）。
  outs でないので `dvc repro` は消さず、git は `*.checkpoint` を無視する。計算を終えて出力を書いたら消す。出力のパスには完成した
  ものだけが書かれ、下流（`[history]` の分岐）が未完の履歴を読むことがない
- 案 L2: 出力そのものを途中で書き、outs を `persist: true` にする（DVC 3.67.1）。dvc.yaml の変更で stage の定義が変わり、途中の履歴が
  出力のパスに残るので、未完と完成を見分ける印と、下流での確認が要る

再開の指定:

- 案 R1: 入力の項目（`[depletion]` の `resume = true` など）。入力は DVC の deps なので、項目を変えるとステージが changed になり、
  投げ直しのコマンドで再開できない。入力を指紋に入れると、項目を変えただけで指紋が変わる
- 案 R2: コマンドの引数（`lattice --resume <入力>`）。dvc.yaml の cmd を変えなければ `dvc repro` で使えず、cmd にいつも付けるなら R3 と同じ
- 案 R3（採用）: チェックポイントがあれば、指紋が一致するときだけ続け、再開した点（状態）を info のログに出す。一致しなければ、
  チェックポイントのパスと、消せば始めから解くことを示してエラーで止める。「入力やコードを変えたのに古い途中の結果から続ける」事故は、
  指紋が入力、読むファイル、実行ファイル、スレッドの数を含むので起きない。再開しても結果は変わらない（判断の軸の 1 つ目。「確認」）ので、
  再開を入力で選ぶ理由がない

指紋:

- 案 F1（採用）: SHA-256（RustCrypto の `sha2`）で、実行ファイル（`/proc/self/exe`）、入力の本文、読むファイル（ライブラリか巨視的
  断面積、チェーン、燃焼の履歴、媒質）、rayon のスレッドの数（出力がスレッドの数によらないかは確かめていない。#0060）を、名前と長さを
  前置して続けて hash する
- 案 F2: ファイルの大きさと更新時刻。中身が変わっても同じになりうる。`dvc checkout` や worktree の作り直しで時刻が変わる
- 案 F3: 標準の `DefaultHasher`。算法が版で決まっていない（std の文書）ので、内容の同一性の指紋に向かない

## 決定と理由

W1、P1、L1、R3、F1 にする。

- 燃焼（`depletion::run`）は、CE/CM のステップを解き終えるたびに、それまでの履歴（点、時間の内訳、輸送の回数）と次の点の組成を
  チェックポイントに書く。再開は書いた履歴に続けて、次の点の組成からステップを進める
- 分岐（`lattice` の状態のループ。`checkpoint::solve_states`）は、状態を解き終えるたびに、それまでの状態の k、出力ごとの少数群定数、
  媒質をチェックポイントに書く（燃焼の最後の点と同じく、最後の状態の後は書かない）。再開は書いた状態を標準出力に出し直し、続きの状態から解く。状態が 1 つの問題（分岐と `[history]` の
  ない入力）は書かない
- MOC は前の輸送の解から反復を始める（ADR-0195）ので、チェックポイントはその解（`moc::Start`）も持つ。燃焼は、書く時点の直前の
  輸送（ステップの中点の輸送）の収束した解を、分岐は軸ごとの解（`branches::Starts`）と次に解く状態を書き、再開はそこから始める。
  持たないと、再開した最初の輸送が平らな初期値から始まり、結果が収束の判定の程度で一気に流したものと違う（#0165）
- チェックポイントの形式は `xlc-lattice-checkpoint`（フィールド名付きの map の MessagePack）で、`fingerprint` と、燃焼なら
  `history`（`xlc-depletion-history` のバイト列）と `next` と `start`、分岐なら `states` と `starts` を持つ
- 履歴の `seconds` と `solves` は、再開の前後の実行の和にする。ログの `total_seconds` は最後の実行だけの値

- 燃焼の履歴の各点は、`[depletion]` の `few-group`（ADR-0173）を書けば、その点の輸送で均質化した少数群定数（`pins = true` なら
  形状関数と角の不連続因子も）を、分岐のテーブルの状態と同じ型（`xlc-table` の状態 1 つ）で持つ。チェックポイントの履歴もそれを持つ。
  書く費用はこれを含めて測った（「確認」）

## 確認

- `depletion::tests::each_finished_step_saves_the_points_so_far_and_the_next_composition`: ステップごとに、それまでの点、輸送の回数、
  次の点の組成を書く
- `depletion::tests::a_run_resumed_from_a_checkpoint_matches_the_straight_run_bit_for_bit`: 1C のピンセルの 3 ステップを、1〜3 回目の
  書き込みの後で止め、チェックポイントのファイルから再開した履歴が、時間の内訳を除いて一気に流したものとバイトで一致する
  （少数群定数を含む）
- `depletion::tests::a_finished_history_cannot_be_resumed`、`a_resumed_history_needs_a_point`、`a_resumed_history_of_another_shape_is_an_error`
- `checkpoint::tests::states_resume_after_the_saved_ones_and_match_the_straight_run`（状態の結果が始める解によるふりの解き方で、
  再開した状態が一気に流したものと同じ解から始まる）、`saved_states_and_starts_read_back_bit_for_bit`、
  `saved_starts_must_follow_the_saved_states`、`the_last_state_is_not_saved`、
  `more_saved_states_than_the_run_has_are_an_error`、
  `a_checkpoint_of_another_fingerprint_is_an_error`、`a_checkpoint_of_the_other_kind_is_an_error`、
  `without_a_checkpoint_the_run_starts_from_the_beginning`、`fingerprints_tell_apart_where_the_bytes_belong`、
  `restored_states_that_do_not_fit_the_run_stop_before_solving`、`a_checkpoint_needs_exactly_one_kind`、
  `an_unreadable_checkpoint_names_its_file`
- 2A の燃焼度の点を最初の 10 点に切った入力（examples/depletion/p2/lattice-2A-10.toml。各点の少数群定数と形状関数も書く）を、5 点目のステップの後で `scancel` して投げ直した
  履歴が、一気に流した履歴と `seconds` を除いてバイトで一致し、11 点の k∞ と数密度、ピン出力が本番の 2A の履歴の最初の 11 点と
  ビットで一致した。書く費用はステップごとに 5.7〜30 ミリ秒（2.75〜14.9 MB）で、1 ステップは約 56 秒。履歴の 1 点は約 1.35 MB で、
  少数群定数と形状関数は 7.6 KB（examples/depletion/verification.md の「時間」、xlc 6ba9ca28、2026-10-23、Slurm 495、496、497）。
  この試行は ADR-0195 の前のコード（MOC を平らな初期値から始め、チェックポイントは `start` を持たない）で流した。`start` を持たせた後の
  再開の一致は、上の単体テストで確かめ、2A では再実行していない
- 2A の `start` は 36.2 MB（MessagePack。代表の領域 3,617 × 69 群のスカラー中性子束 2.3 MB、高次のモーメント 20.2 MB、代表のトラックの端の
  角度中性子束 13.7 MB）で、1 つ目のステップの後のチェックポイントは 39.0 MB だった（examples/depletion/verification.md の「時間」、
  xlc 8c6e4632、2026-10-23、排他でない Slurm 508）
- 分岐の状態のループは、本物の格子計算で止めて再開する試行をしていない（状態の結果の読み書きと、続きから解く順は上のテスト）

## 帰結

- W1 は書くたびにそれまでの全体を書くので、書く量の和は点の数の 2 乗で増える。70 点の 2M の最後の 1 回は約 95 MB（1 点 1.35 MB）で、
  10 点の測定（14.9 MB で 30 ミリ秒）から 1 秒以下と見込む（測っていない）。2M の 1 ステップ（自己遮蔽 67 秒の輸送 2 回）に比べて
  小さいので、追記（W2）にしない。見込みが外れたら W2 を比べ直す
- 細かい群の媒質（ADR-0128）を出力する分岐は、1 状態の媒質が群の数の 2 乗の散乱行列を持ち、書く量がさらに大きい。測っていない。
  そうした分岐を何百状態も流すときは、書く時間（ログの `checkpoint` の行）を見て W2 と比べ直す
- 履歴のファイルの md5 は、一気に流しても実行ごとに違う（`seconds` が時間の内訳）。再開した出力と一気に流した出力の一致は、
  `seconds` を除いて比べる。dvc.lock の出力の md5 も、再開の有無によらず実行ごとに変わる（今までと同じ）
- 再開したジョブの経過時間とログの `total_seconds` は最後の実行の分だけで、時間の記録に使わない。始めから測るときはチェック
  ポイントを消す（.claude/skills/verification-record）
- 燃焼のチェックポイントは、履歴に加えて MOC の解 1 つ（2A で 36.2 MB）を毎回書く。分岐は軸ごとの解のうち異なるものだけ（多くて軸の数。
  同じ解は 1 度だけ書き、軸からは番号で指す。2A の格子で 1 つ 36.2 MB）を毎回書く。1 状態 25〜78 秒に比べて書く時間は小さいと見込む（1 つ目のステップの 39.0 MB は 40 ミリ秒。排他でない
  ジョブの値で、記録の時間には使わない）
- 反復を始める解や、前の輸送から持ち越す量を足すとき（ADR-0194 の Dancoff 係数のように、入力が同じなら解き直しても同じになるもの
  を除く）は、チェックポイントにも足す。再開のテストがビットの不一致として検出する
- 燃焼の履歴の `few-group` を本番の入力（examples/depletion/p2）に書くかと、公称の状態のテーブルを履歴から輸送なしで作るかは、
  この ADR で決めない（入力を変えると DVC のステージが changed になる）
- 急ぐ計算を前に出すのは `--urgent`（SUSPEND）で、`scancel` しない（nixfiles の slurm-run のスキル）。チェックポイントは取り消しと停止のため
