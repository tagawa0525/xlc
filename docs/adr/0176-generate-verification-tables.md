---
status: proposed
date: 2026-10-15
requires: [ADR-0117, ADR-0122, ADR-0130, ADR-0143, ADR-0159]
supersedes: []
superseded-by:
issues: ["#0130"]
---

# ADR-0176: 検証の記録の数値の表を、比較の段階が書く JSON から描き、記録と一致することを検査する

## 背景

検証の記録（examples/*/verification.md）の数値の表は、比較のツールの出力（`compare-*.out`、文章の形式）を読んで手で打っている。
2026-10-14 の 3A の記録では、手計算の転記の誤りが 2 つ見つかった（OpenMC と 2σ 以内の層の数、ピンの最大の差の倍率）。
表は 13 の記録に 91 個（`grep -c '^| -'` の合計。examples/vera/verification.md に 42）あり、どれも人間（サブエージェントを含む）が
読んで打つ。ユーザーの指摘（2026-10-15）は、結果の表などを自動で処理する仕組みを整えたほうが楽、というもの。

決めること:

- 比較のツールが、表に出す量を機械可読の形でも書くか
- 表を描くツールと、記録の中の置き場の目印
- 表の「計算: <コミット>、<日付>」（CLAUDE.md。ADR-0117、ADR-0122）を自動で出すか
- 記録と描いた表が一致することの検査をどこに組み込むか
- ツールの置き場

判断の軸:

- 正確さ: 転記の誤りが構造的に起きないこと。手で直した表がそのまま通らないこと
- 再現性: 記録の数値が、checkout すれば出る結果（ADR-0122）と結びついたままであること。rebase で測り直さないこと（commit-refs remap）
- 保守: ツールと目印が増える分の読む量。標準ライブラリで済むこと
- 検査が実際に走ること: 記録が手で崩れたときに、コミットとマージの前で気づくこと

## 検討した案

### 表の作り方

- 案 A: 比較のツールが JSON を段階の outs に書き、別の 1 本のツールが記録の目印の間をそこから描いた表で置き換える（採る）
- 案 B: 比較のツールが Markdown の断片（表）を直接書き、記録がそれを取り込む。ツールが表の形（列、桁、符号、全角の幅）を
  知ることになり、同じ量を別の表（問題の組み合わせや列の並びが違うもの）に使い回せない。複数の問題の行を 1 つの表にまとめる
  のは、問題ごとに走る段階の外の仕事なので、結局まとめるツールが要る。数値と文章の分離ができず、`--check` で比べるのも文字列になる
- 案 C: 記録を Markdown でなく生成物にする（テンプレートと JSON から全体を作る）。記録は文章（原因の候補、確かめていないこと、
  再現のコマンド）が大半で、表は一部である。文章をテンプレートに移すと、書き直し（ADR-0117 の「今の姿」）のたびにテンプレートと
  記録の 2 つを読むことになり、読みやすさと保守の両方で悪い。git の差分でも記録が読めなくなる
- 案 D: 手で打つまま、転記の確認を `/code-review` に任せる。誤りの検出が人間の注意に依存し、昨日の 2 件を見つけたのも偶然に近い。採らない

### JSON の形

- 案 A: 表に出す量だけを、割合は小数、k の差は pcm で持つ（採る）。ツールは書式（桁数、符号、単位）だけを持つ
- 案 B: 書式化した文字列を JSON に持つ。表の形が比較のツールに戻り、案 B（上）と同じ問題になる
- 案 C: 層ごとの値など、すべての中間値を持つ。`.out` が既に持つ。表に出さない量は、必要になるまで足さない

### 計算の欄（コミットと日付）

表の下の「計算: <コミット>、<日付>」（実際の記録では、本文に「数値は 60cf97d0（分岐点 4b644566）、2026-10-14」の形）。

- 案 A: 段階が `git rev-parse HEAD` を呼んで JSON に書く。JSON の md5 が、段階を実行したコミットごとに変わる。rebase で hash が変わると、
  JSON の中の hash を `commit-refs remap` で書き換える必要があり（対象が JSON にも広がる）、書き換えると `dvc.lock` の md5 と合わなくなる。
  ADR-0122 は rebase のあとに測り直さず付け替えると決めており、段階を再実行しない限り JSON は変えられない。また、1 つの表の数値は
  複数の段階が出す（xlc の炉心は core の段階、OpenMC は実行 m3、比較は compare の段階）ので、比較の段階の HEAD だけでは足りない
- 案 B: ツールが `git log -1 --format=%h -- <dvc.lock>` で、その JSON の出力を最後に変えた `dvc.lock` のコミットを取る。取れるのは
  `dvc.lock` のコミット（実行したコミットの次）であって、実行したコミットではない（topic-start の「ジョブが終わったら」の通り、
  記録に書くのは実行したときの HEAD で、`dvc.lock` はその次のコミット）。relock（rebase のあとに数秒の段階を再実行する）でも変わり、
  分岐点の復元の記録（「分岐点 4b644566」）も出せない
- 案 C: 今までどおり手で書く。コミットは `commit-refs remap` が付け替える（ADR-0122）。表の数値は検査され、欄は書く人が書く
  （採る）。表の数値と欄の対応が崩れる余地は残る

### 検査の組み込み

- 案 A: `.pre-commit-config.yaml` の hook と、crates/docs/tests/consistency.rs の `cargo test`（`run-checks.sh` が実行する）の両方（採る）。
  ADR-0143 の `no-inline-tests` と同じ形。`run-checks.sh` は pre-commit の設定のうち design-records の hook しか実行せず、`cargo test` は
  実行する。hook だけだと /git-merge の検査（`run-checks.sh --base`）を通る道がない
- 案 B: DVC のステージにする（JSON と記録を deps、検査結果を outs）。記録は文章なので `dvc.lock` が毎回変わり、`dvc repro` を実行しないと
  検査が走らない。記録が手で崩れても、`dvc repro` を実行するまで誰も気づかない。採らない
- 案 C: hook だけ。`cargo test` と二重にならないが、上のとおり /git-merge の道に乗らない

### JSON の置き場

- 案 A: DVC のキャッシュに置く（最初の試作。ほかの出力と同じく git で追跡しない）。チェック（`--check`、`cargo test`）が作業ツリーの
  DVC の出力に依存し、`dvc checkout` していない worktree と、nix の sandbox（CI）で落ちる。落ちる前提を守る手順が要り、
  チェックが目的（記録が崩れたときの検出）より環境の準備を求める。捨てた
- 案 B: `cache: false` で git で追跡する（採る）。`dvc.lock` は md5 を持ったままなので、手で直した JSON は照合で落ちる

### 置き場

- 案 A: `examples/verification_tables.py`（採る）。記録は examples/<問題>/ にあり、`examples/openmc_*.py` のように、複数の例が使う
  Python のツールが既にある。テストは同じ階層の `test_verification_tables.py`（`test_compare_p3.py` と同じ流儀）
- 案 B: `scripts/`。現在の中身は hook から呼ぶシェルと AWS のツールで、記録の数値に結びついた Python は無い。hook が呼ぶツールは
  `scripts/check-*.sh` の形だが、このツールは `--check` のほかに置き換えも行う（検査だけのツールではない）

## 決定と理由

案 A（JSON + 1 本のツール + 手で書く計算の欄 + hook と `cargo test`）。

1. **比較のツール**: `compare_p3.py --json <出力先>` が、標準出力に出す数と同じ値を JSON にも書く。表に出す量だけを持つ（割合は小数、
   k の差は pcm、位置は [行, 列]、層は 1 始まり）。形は `{"format": "xlc-compare-p3", "problem", "k": {...}, "axial_vs_keno": {...},
   "axial_vs_openmc": {...}, "pins_vs_openmc": {...}, "layer_pins": {"levels": [...]}, ...}`（examples/vera/compare_p3.py の docstring）。
   既存の `compare-3A.out` と `compare-3B.out` は変えない（md5 が同じ）。OpenMC の参照解は DVC のステージの出力でない
   （AWS の VM の実行 m3 から取り出した `references/` のファイル）ので、それを deps にする段階 `compare-<問題>-openmc` を足し、
   `.out` と JSON を outs にする（JSON は `cache: false` で git で追跡する。「JSON の置き場」）。今後の比較のツール（compare_p5.py、compare_keno.py など）も、表の元にする量を `--json` で書く
2. **表を描くツール**: `examples/verification_tables.py <記録>`。記録の `<!-- table: <名前> -->` から `<!-- /table -->` までを、
   その名前の表（空行 + 表 + 空行）で置き換える。名前と JSON と描く関数の対応はツールの `TABLES`。目印の中は丸ごと置き換わるので、
   脚注などの文章は目印の外に書く。未知の名前、閉じていない目印、入れ子、同じ名前の 2 回の使用は
   エラー。列は prettier と同じ幅（全角を 2）でそろえ、負の符号は記録と同じ U+2212。標準ライブラリだけ
3. **`--check`**: 置き換えず、描いた表と記録が違えば差分を出して 1 で終わる。JSON が無い（段階を実行していない。`dvc repro` の指示を出す）、JSON の md5 が
   同じ場所の `dvc.lock` の値と違う（再実行して relock していない、手で直した）ときも 1 で終わる。黙って読み飛ばさない
4. **計算の欄は手で書く**。案 A は rebase と相性が悪く、案 B は取れる値が違う（上）。欄は `commit-refs remap` が付け替え、
   表の数値は `--check` が検査する。残る穴（表の数値と欄の hash の対応）は、JSON の md5 が `dvc.lock` と一致するかを検査することと、
   記録に書く手順（verification-record）で抑える
5. **検査**: `.pre-commit-config.yaml` の hook `verification-tables`（`examples/verification_tables.py`、`*/verification.md`、`*/dvc.lock` の変更で走る）と、
   `crates/docs/tests/consistency.rs` の `record_tables_match_the_compare_outputs`（`cargo test` に含まれるので `run-checks.sh` でも走る）
6. **試作の範囲**: Problem 3 の 2 つの表（`p3-reference-results`、`p3-openmc-comparison`）。残りは #0130

## 確認

試作（ブランチ feat/verification-tables、分岐点は main の 1e0c5941、2026-10-15）で確認した。

- 描いた 2 つの表は、既存の手で打った表と文字まで同じ（誤りの修正は 0 件）。段階は 4 つ（compare-3A、compare-3B、compare-3A-openmc、
  compare-3B-openmc）を `dvc repro -s` で再実行した（秒。コード、炉心の出力、OpenMC の参照解は変えていない）。`compare-3A.out` と
  `compare-3B.out` の md5 は変わらない
- テスト: `examples/vera/test_compare_p3.py`（`axial_report` が返す値、`--json` の出力が標準出力と同じ k を持つこと）、
  `examples/test_verification_tables.py`（表の行、列の幅、置き換えの冪等、目印の誤り、`dvc.lock` との不一致、JSON が無いこと、
  手で直した表が `--check` で落ちること）
- 本文の数値を JSON と突き合わせた（2σ 以内の層の数 3A 1 と 3B 4、ピンの最大の差の倍率 3A 5.0 と 3B 6.3、軸方向のオフセットの差の倍率 22.9 と 35.9）。
  表の外の文章の数値は検査しない（#0130）

## 帰結

- 表の数値の転記の誤りは、描く側では起きない。手で直した表は `--check` で落ちる
- 表の外の文章（本文の中の数値）は、手で書くまま。同じ量を 2 か所に書くと食い違いうる。本文の数値を JSON から埋める方式は、必要に
  なるまで作らない（#0130）
- JSON は DVC の出力（`cache: false`）のまま git で追跡する。`--check` と `cargo test` は、worktree でも nix の sandbox でも
  `dvc checkout` なしに走る。JSON が手で直されると `dvc.lock` の md5 と合わず落ちる。JSON は数 KB の要約なので、リポジトリは
  ほとんど増えない（試作の 2 つで 14 KB ずつ）
- Python のテストは `run-checks.sh` の pytest が走らない配置（リポジトリの直下に `test_*.py` が無い）なので、手で
  `python3 -m unittest discover -s examples -p 'test_verification_tables.py'` を実行する（`test_compare_p3.py` と同じ）
- 比較のツールが `--json` を持つたびに `p3_inputs.py` のような生成元と TABLES を足す。表の形の知識は `verification_tables.py` に集まる
- 試作の 2 つの表の元の `compare-<問題>-openmc` は、OpenMC の参照解の 4 つのファイルを deps にする。参照解を作り直したら（実行 m3 でなくなったら）
  deps のパスを直し、段階を再実行する
