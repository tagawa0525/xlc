---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# xs が MF5 の LF=12（Madland-Nix）と MF1 の LNU=1（多項式の ν）を読めず、ENDF/B-VII.0 と VII.1 の燃焼のライブラリを作れない

## 問題

- 版ごとのピンセルの燃焼（#0303）の VII.0 と VII.1 のライブラリ（ブランチ feat/0303-pin-depletion-by-endf-version の
  examples/depletion/endfb70/xs.toml、endfb71/xs.toml）の 65 核種を、分離共鳴の SLBW を読めるようにした xs（#0305）の `xs nuclide`
  で 1 つずつ実行すると、次の 5 つが読み込みで止まる（2026-11-20、手元の試行。どの核種も 30 秒で打ち切り、それまでに終わった
  52、打ち切った 73 は読み込みでは止まらなかった。打ち切った核種が後の処理で止まるかは確かめていない）
  - VII.0 と VII.1 の Am-241、Am-243: 核分裂スペクトル（MF5 MT18）が LF=12（Madland-Nix）で、「LF=12 is not supported」
    （crates/xs/src/endf/mf5.rs は LF=1、5、7、9 だけを読む）
  - VII.0 の Np-239: ν（MF1）が LNU=1（多項式）で、「LNU=1 is not supported」
- VIII.0 と VIII.1 の同じ 65 核種がこれらの形式を持たないかは、`xs range` が MF2 の分離共鳴の上端までしか読まないので、この調べ方では
  確かめていない（VIII.0 は本番のライブラリ `depletion.msgpack` が作れている）

関係する ADR: なし

## 完了の条件

- [x] MF5 の LF=12 を ENDF-102 の式で読み、NJOY と比べる（式と群への振り分けは GROUPR と同じで、NJOY の誤りに当たらず、形の決定は
  ENDF-102 の定義のとおりなので ADR にしない。比較は examples/depletion/verification.md の「版ごとのライブラリの Madland-Nix の核分裂
  スペクトルと多項式の ν」）
- [x] MF1 の LNU=1 を読む（MT452 の定数と 1 次の多項式。4 つの版の全評価の LNU=1 はどれも MT452 の NC=2。同じ節）
- [ ] VII.0 と VII.1 の 65 核種の `xs nuclide` が、打ち切らずに最後まで終わる

## 経緯

- 2026-11-20: SLBW を読めるようにした後に VII.0 と VII.1 の 65 核種を `xs nuclide` に通して見つけ、起こした（#0305）
- 2026-11-20: LF=12 と LNU=1 を読めるようにし、VII.0 の Am-241、Am-243、Np-239 と VII.1 の Am-241、Am-243 の `xs nuclide` が終了コード 0 で
  終わることを確かめた（ブランチ fix/0307-read-madland-nix-and-polynomial-nu）。65 核種の全部は、マージの後に確かめる
- 2026-11-20: LF=12 の累積の閉じた式は、E'/T_M が小さいところで 1 程度の項の差になり、1e-16 ほどの絶対誤差が残る（1e-2〜1 eV の群の
  割合が 1.6e-7 ずれ、1e-4 eV より下の群では負になりうる。レビューの指摘）。小さいところを、g を √u E1(u) の積分に書いた 8 点の
  Gauss-Legendre の二重の積分で求める形を試したが、Am-241（VII.0）の `compare_njoy_gendf` が閉じた式だけの 4.1 秒から 36.7 秒になった
  ので捨て、同じ積分の級数の始めの 3 項にした（2.4 秒。閾値の前後の相対誤差は、閉じた式で 3e-9 以下、級数で 1.1e-11 以下）
