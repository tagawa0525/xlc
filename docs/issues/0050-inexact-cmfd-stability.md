---
status: open
labels: [lattice, moc, verification]
opened: 2026-10-10
closed:
milestone:
blocked-by: []
resolution:
---

# 緩めた CMFD の精度が、収束の遅い体系で MOC の反復を増やさないかの確認

## 問題

ADR-0110 は、CMFD の精度を MOC のスイープによる変化に比例して緩める（`moc::inexact_cmfd`）。確かめたのは VERA 2B の
集合体と C5G7 だけで、どちらも MOC の反復の回数は変わらなかった。ADR-0110 の「結果」にも、収束の遅い体系で回数が
増えたら緩め方を見直す、とある。

先行例を軽く調べた（2026-10-10、web 検索のみ。論文の本文は読んでいない）。

- CMFD の許容誤差を MOC の残差に連動させる方式は、見つからなかった。OpenMC の CMFD は固定の許容誤差（絶対 1e-10、
  相対 1e-5）で、STREAM は線形ソルバの呼び出しを経験的な固定回数（5 回）にしている。調べた範囲では、固定値が普通
- 内側の精度を外側の残差に比例させる考え方は、不正確な Newton 法の強制項（Eisenstat と Walker、1996）として一般的で、
  NDA を JFNK で解く研究（Knoll、Park、Smith、2011 など）の系統にある。ただし CMFD に当てた例は、見つからなかった
- 内側の収束と安定性について、報告が逆向きに分かれる。consistent CMFD は内側反復の収束が甘いと不安定になる
  （pCMFD の導入の動機）。一方、熱流動との連成では、低次元の拡散方程式を厳しく収束させると不安定になる、との報告もある。
  どちらが当たるかは体系によると思われる
- CMFD は、光学的に厚い格子や高散乱の問題で発散や劣化が起きる（lpCMFD、2D/1D の Fourier 解析）。緩めた CMFD が
  この領域で反復を増やすかは、文献からは分からない

関係する ADR: なし

出典: OpenMC cmfd.xml、pCMFD（ANS）、Low-Order Multiphysics Coupling Techniques（OSTI 1785622）、Practical Considerations
for the Adoption of Anderson Acceleration（OSTI 2580831）、Knoll・Park・Smith（OSTI 1016204）、Willert ら（OSTI 22382121）、
PHYSOR 2020 lpCMFD、2D/1D CMFD の Fourier 解析（Ann. Nucl. Energy）。

## 完了の条件

- 光学的に厚い、または高散乱の体系（メッシュを粗くした格子、反射体の厚い体系など）で、設定の精度で毎回解く場合と
  緩める場合の MOC の反復の回数と k を比べ、examples/*/verification.md に記録する
- 反復が増える、または収束しない体系があれば、緩め方を直す案を ADR にする（ADR-0110 は accepted なので、新しい ADR で
  置き換える）。なければ、このまま閉じる

## 経緯

- 2026-10-10: ユーザーが ADR-0110 の方式に前例があるかを尋ね、軽く調べた。同じ方式の先行例は見つからなかった。
  ADR-0110 は accepted で書き換えられないので、調査の結果と未確認の点をここに書いた
