---
status: open
labels: [verification]
opened: 2026-11-04
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0226 の後の再実行で、記録の表の外の本文が古い値のまま残る

## 問題

ブランチ chore/rerun-after-adr-0226（2026-11-03〜04）は、bin の変更で changed のステージを今の bin で再実行し、比較の JSON から描く
表と「計算」の欄を今の実行に合わせた。表の外の本文（手で書いた表、結論、節の冒頭の文）は、依頼の範囲の外として書き直していない。
次の箇所が、再実行の前の値のまま残り、描き直した表や dvc.lock の出力と食い違う。

- examples/vera/verification.md: 冒頭の「ADR-0209 の後に作り直していない」の文。「Problem 3 の炉心の計算」の出どころ（1f650c35、
  ジョブ 475、`vera.msgpack` 05132f1d）と手で書いた KENO の表と分析。P3 の「OpenMC との比較」の本文（xlc − OpenMC の +38、+41 pcm
  など。出どころは今の実行に替えたので、本文と食い違う）。格子の「## OpenMC との比較」の表と出どころ（ジョブ 465、470）
- examples/vera/verification-endfb70.md: 「結論」の数値（VII.0 の +80〜+245 pcm など）
- examples/depletion/verification.md: 冒頭の、1C と 2A だけが ADR-0209 の既定という文。MOX-1 の冒頭の −219 pcm（../mox/ の比較は今
  −212 pcm）。2F〜2A-rods の問題ごとの k の差の文。2A の履歴を前の md5 7efea664 で「今の履歴」と引く節。「core のミクロ燃焼」、
  「サイクルのステップ」、「制御棒の分岐と値ごとの格子計算の結合の一致」の節は #0213 の残りに書いた
- examples/mox、examples/m0、examples/m1、examples/core、examples/reflector、examples/colorset、examples/split の verification.md の
  k と表、冒頭の ADR-0209 の文

関係する ADR: なし

## 完了の条件

表の外の本文を、表を描き直したステージの JSON（JSON の無いステージは dvc.lock の出力）の値に合わせて書き直し、上の一覧の箇所が
今の実行と食い違わなくなる。

## 経緯

- 2026-11-04: ADR-0226 の後の再実行（ブランチ chore/rerun-after-adr-0226）の報告とレビューで、残りとして挙げた
