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

- [ ] MF5 の LF=12 を ENDF-102 の式で読み、NJOY と比べる（形を決める ADR）
- [ ] MF1 の LNU=1 を読む
- [ ] VII.0 と VII.1 の 65 核種の `xs nuclide` が、打ち切らずに最後まで終わる

## 経緯

- 2026-11-20: SLBW を読めるようにした後に VII.0 と VII.1 の 65 核種を `xs nuclide` に通して見つけ、起こした（#0305）
