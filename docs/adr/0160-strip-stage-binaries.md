---
status: accepted
date: 2026-10-14
requires: [ADR-0159]
supersedes: []
superseded-by:
issues: []
---

# ADR-0160: DVC の build の段階は、bin/ に置く実行ファイルを strip する

## 背景

ADR-0159 は、計算の段階の deps をビルドした実行ファイル（`bin/xs`、`bin/lattice`、`bin/core`）にし、出力を変えない変更では
実行ファイルのバイトが変わらず段階が走らないことを狙った。worktree によって変わる RUNPATH は `NIX_NO_SELF_RPATH` で
取り除き、strip は RUNPATH を消さないので退けた（ADR-0159 の案 P2）。

テストを tests.rs に移すブランチ（#0089。refactor/0089-move-tests、main の 703e611c の上）で `build` を走らせると、3 つの
実行ファイルの md5 が main の `dvc.lock` と違った。計算に使う節は同じで、違うのは記号表だけだった（「確認」）。記号表は
実行のときに読まれず、計算の結果に関わらない。

判断の軸は ADR-0159 と同じで、段階が走るかどうかを、段階が実際に行う計算で決めることである。

## 検討した案

- 案 A（採用）: `scripts/build-binaries.sh` の `target/dvc` のビルドに `--config 'profile.release.strip="symbols"'` を付け、
  リンクのときに記号表（`.symtab`、`.strtab`）を落とす。記号の名前の違いで md5 が変わらない。設定はスクリプト（`build` の
  deps）にあり、`Cargo.toml` の release の設定と手元の `target/release` は変えない。`bin/` の実行ファイルでは perf などが
  関数の名前を出せなくなるので、そのときは `target/release` を使う。ADR-0159 が strip を退けたのは RUNPATH の解決として
  で、ここで扱う記号の名前は別の原因なので、その判断と食い違わない
- 案 A2: ビルドのあと、binutils の `strip -o` で写す。記号表を落とす効果は同じ（この 2 つは別のバイトになるが、どちらも
  #0089 の前後で一致した。「確認」）。ただし、PATH の先頭の `strip` の版でバイトが変わりうるので、devShell の binutils で
  あることを確かめる手間が増える。案 A は cargo とリンカ（`build` の deps と ADR-0159 の範囲）だけで済む
- 案 B: 実行ファイルはそのまま置き、md5 の代わりに、割り当てられる節（ロードされる節）だけの digest を書く段階を挟んで
  deps にする。記号表を残せるが、判定を自前の道具で持つことになり（ADR-0123 の「デファクトの道具を使う」に反する）、
  どの節を含めるかを自分で決めて保守することになる。ADR-0159 の案 D と同じく、間に 1 つファイルが増える
- 案 C: strip せず、テストのモジュールの有無で md5 が変わるのを受け入れる。#0089 のテストの移動で、最新だった xs の
  63 段階（ADR-0159 の移行で揃えたもの）が古くなり、ADR-0159 の目的（出力を変えない変更で流し直さない）が、その例で
  果たせない

## 決定と理由

案 A にする。ADR-0159 は置き換えず、`build` のビルドの設定だけを補う。`build` は、ADR-0159 の手順のビルドに
`profile.release.strip="symbols"` を足し、記号を落とした実行ファイルを `bin/` に写す。worktree のパスなどの確認は、写す
実行ファイルそのものに対して行う。

記号表は実行のときに読まれないので、strip しても計算は変わらない。#0089 の観察では、記号表を落とすとテストの移動の前後で
3 つともバイト一致するので、ADR-0159 の目的がその例で成り立つ。

この変更で `bin/` の md5 が変わるので、計算の段階の deps が変わる。ADR-0159 の移行と同じく、この変更の直前の main
（703e611c）で最新だった段階は `dvc commit` で取り込み、古かった段階は古いまま残す。最新だったかは、切り替えの前の
`dvc status` に出なかったことで決め、取り込んだ段階は、切り替えのあとの `dvc status` で deps の差が `bin/xs` だけだった
ものに限った。記号を落としても計算が変わらないことは、記号表が実行のときに読まれないことによる（案 A と案 A2 で
落とした版は、strip しない版と同じ命令を持つと見ているが、節ごとには比べていない）。

## 確認

main の 703e611c と refactor/0089-move-tests の dd0aef10（703e611c の上）を、それぞれ `scripts/build-binaries.sh`（ADR-0159
の版）でビルドして比べた（2026-10-14）。

- strip しない md5: main は xs `8376dbec…`、lattice `bfe45232…`、core `2f204012…`（main の `dvc.lock` と同じ）、#0089 は
  xs `b0cc2562…`、lattice `40f1cda8…`、core `e426a664…`
- 割り当てられる節（`readelf -S` の各節を `objcopy -O binary --only-section` で取り出したもの）は、3 つとも全節がバイト
  一致した。`.symtab` と `.strtab` は大きさが同じか近く（xs の `.strtab` は 0xb8b1e と 0xb8b24）、中身が違った
- `nm` の記号の名前の差は、xs 6、lattice 40、core 10 行で、どれも `anon.<hash>.<番号>.llvm.<数>` の `.llvm.<数>` の添字
  だけだった（例: xs の `anon.009487…29.llvm.9964443520753796892` と `….llvm.10936925308106898528`）。LLVM が
  コード生成の単位をまたぐ内部の記号に付ける添字で、テストのモジュールを消すと変わる、と見ているが、決まり方は確かめて
  いない
- binutils の `strip -o` で落とした md5（案 A2）は、main と #0089 で 3 つとも一致した: xs `731fd73e…`、lattice `3e17fcf2…`、
  core `51134de4…`
- `--config 'profile.release.strip="symbols"'` でビルドした md5（案 A）も、main と #0089 で 3 つとも一致した: xs `e77257b3…`、
  lattice `85a3c566…`、core `fe2d8761…`。`.symtab` はなかった。このブランチの `dvc repro -s build` が書いた `bin/` も同じ md5
- 移行: 703e611c で最新だった段階のうち、`bin/` の md5 が変わって古くなったのは 63（examples/vera の xs 系 62 と
  examples/reflector の xs。deps の差は `bin/xs` だけ）で、`dvc commit` で取り込んだ（案 A の `bin/xs`）。`dvc.lock` の差分は `bin/xs` の md5 と
  大きさだけだった。lattice と core を呼ぶ段階は前から古く、`dvc.lock` を変えていない。そのあとの `dvc status` は、
  変更の前と同じ 127 項目（`references-shell` と、前から古い 126 の段階）

## 帰結

- #0089 のテストの移動のように、記号の名前だけを変える変更では、`bin/` の md5 が変わらない（観察は #0089 の 1 例。添字の
  決まり方を確かめていないので、ほかの変更で割り当てられる節まで変わるかは分からない）。#0089 は、rebase の
  あとに `build` を走らせ直し、直下の `dvc.lock` を更新するだけで済む
- `bin/` の実行ファイルには記号がない。perf などで関数の名前が要るときは `target/release` を使う（docs/environment.md）。
  パニックのメッセージの位置（ファイルと行）は記号表でなく割り当てられる節にあるので残る。release はデバッグ情報を
  持たないので、`RUST_BACKTRACE` のバックトレースの関数の名前も `bin/` では出なくなる。段階の中のパニックを調べるときは、
  同じ入力を `target/release` の実行ファイルで走らせ直す
- 記号を落としても `.comment`（コンパイラとリンカの版）と `.shstrtab` は残る。リンカなどの版が変わると、計算が同じでも `bin/` の
  md5 は変わる（ADR-0159 の帰結のとおり）
- 行の位置を変える整理は、割り当てられる節を変えるので、strip しても段階を古くする（#0110）
