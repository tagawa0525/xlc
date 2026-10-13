---
status: open
labels: [lattice, core]
opened: 2026-10-09
closed:
milestone: M3
blocked-by: ["#0093"]
resolution:
---

# 軸方向の反射体の定数

## 問題

ADR-0128 で、軸方向の反射体の定数を、単一集合体から書き出した均質な燃料の媒質（`xlc-medium`）と、層の断面の中で
体積で混ぜた物質（`mix`）を集合体の高さ全体に並べた平板の格子計算（`axial = true`）で作り、少数群定数の z の辺
（6 つの辺と `height`）を core の z の面で径方向と同じ規則で使うようにした。VERA Problem 3A の平板の例
（examples/reflector の `axial-3A-*.toml`）で、core が平板の MOC の k とノード平均の中性子束を再現する
（examples/reflector/verification.md）。残るのは次のとおり。

- 上の層（プレナム、上の端栓とギャップ）のテーブルの制御棒の離散の軸（ADR-0128 の「制御棒の軸」）。テーブルの形式の
  軸（`rod`、値の名前）と core の区間ごとの選択、棒を切る高さ（炉心の形の `[rods]` の `top`）は #0063 で入れた
  （ADR-0127）。残るのは、平板の格子計算で上の層の制御棒の値の定数を作ること（格子計算の制御棒の分岐 #0093 の書き方に
  そろえる）と、core が核分裂のない層の部分挿入のノードを扱うこと（z の1ノードの問題を ANM にする。いまはエラー）
- Problem 3A と 3B の core の解（examples/vera/p3。燃料は 2D の格子、支持格子は ADR-0140、軸方向は平板）を、同じ核データの
  OpenMC（3D）と比べる（ADR-0128 の確認）。KENO の値との比較は済み（examples/vera/verification.md の「Problem 3 の炉心の
  計算」）。OpenMC の本番は #0073 の経緯のとおり外部の計算機で流す。燃料の両端の層の軸方向の出力の差（3A で −6.1% と −3.2%）の原因は確かめていない
- 上の層の吸収体を体積で混ぜる近似（ADR-0128）の大きさが分かっていない。Problem 5 の全引き抜きで、上の層の
  制御棒の値を棒なしに替えた k の差を測り、KENO との差に比べて無視できなければ、棒のセルを中性子束で均質化した媒質を
  混ぜる形を ADR にする
- 平板の列の幅は、格子計算の CMFD の収束（#0069）で 1 cm 程度に制限している

関係する ADR: [ADR-0128](../adr/0128-make-axial-reflector-constants-in-a-full-height-slab.md)、[ADR-0134](../adr/0134-use-current-weighted-d-where-koebke-has-no-root.md)、[ADR-0140](../adr/0140-place-spacer-grids-as-axial-layers-of-framed-lattices.md)、[ADR-0141](../adr/0141-use-medium-b1-diffusion-in-slab-fuel-nodes.md)

## 完了の条件

上の層のテーブルに制御棒の離散の軸を入れ（#0093 の分岐の書き方で）、核分裂のない層の部分挿入のノードを core が扱い、Problem 3 を OpenMC（3D）と KENO の値と比べる
（ADR-0095 と ADR-0128 の確認）。

## 経緯

- 2026-10-09: ADR-0095 の I 型を実装したときに残した
- 2026-10-11: 入力と出力の形、core の z の面の扱い、層を ADR-0128（proposed）に書いた。単一集合体から書き出す均質な
  燃料の媒質（`xlc-medium`）、集合体の高さ全体の平板、少数群定数の z の辺（6つの辺と `height`）、z の面の f* と真空の比。
  上の層の制御棒の離散の軸（全引き抜きでも棒はプレナムの層にある。制御棒の部分挿入の ADR-0127 の
  「反射体のテーブルの制御棒の軸」と #0063）を含む
- 2026-10-11: ADR-0128 を実装した（媒質、`mix`、`axial = true`、少数群定数の z の辺、core の z の面）。VERA 3A の平板で、
  定数を作った体系を core が再現することを確かめ、ADR-0128 を accepted にした。上の層の制御棒の軸は #0063 に合わせて
  入れるので、#0063 を待つ（`blocked-by`）。平板を 0.25 cm の列で解くと格子計算の CMFD が収束しなかった（#0069）
- 2026-10-12: #0063 でテーブルの制御棒の軸と core の選択を入れたので、残りは上の層の値の定数（#0093 の格子計算の分岐を
  待つ）と核分裂のない層の部分挿入のノードになった。`blocked-by` を #0093 にした
- 2026-10-12: Problem 3A と 3B を core で解き、KENO の値と比べた（k の差は +84 pcm と +193 pcm）。平板の燃料のノードを
  炉心の 49 層にそろえると、上端の燃料のノードで中性子流の重みの D が負になったので、媒質だけのノードを媒質の B1 の D に
  した（ADR-0141）
- 2026-10-13: examples/reflector の平板を ADR-0141 で再実行した。f* の方法の燃料のノードの差の最大は燃料 16 ノードで
  8.5% から 4.2% に小さくなり、32 ノードで 3.0% から 5.5% に大きくなった（理由は確かめていない。
  examples/reflector/verification.md）
