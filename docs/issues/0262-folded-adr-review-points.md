---
status: open
labels: [docs]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# ブランチの上の ADR を畳む規則と issue の書き直しの手順のレビューで残した点

## 問題

docs/adr/README.md の「運用」の畳む項と rejected の項、docs/issues/README.md の「運用」の書き直しの前の項に残る点:

- /git-merge の「proposed の ADR」の節の「試して、採らないと決めた」は、同じブランチの別の ADR で置き換える場合を分けず、
  rejected にして残すと読める。畳む規則はリポジトリの README にだけあり、マージの手順の側からは読めない。README の rejected の項は
  その項を、置き換えないときに当たるものとして指しており、スキルの項より狭い条件で引いている
- 畳む項は、消した全文の在り処に「消す前のコミットとパス」を許す。ブランチのコミットは rebase で変わる（`remap.commitRefs` の
  付け替えに頼る）が、blob は変わらない
- 畳む規則そのものと、比べた案（rejected の ADR を別のブランチに分けてマージする）は ADR に書いていない
- 書き直しの前に探す文書を ADR、検証の記録、計画、ほかの issue に限っている。入力の TOML の注釈、コードのコメント、CLAUDE.md も
  issue を引きうる（リポジトリの全体を `git grep` で探す形もある）
- 同じ未マージのブランチで作って書き直した issue には、書き直す前の版のある main のコミットが無い
- docs/adr/README.md の「書き方」の rejected の項の括弧が、「運用」の畳む項を短く言い換えている
- /adr の status の `withdrawn`（決める前に取り下げた）と、畳んで消す形の関係が書かれていない

関係する ADR: なし

## 完了の条件

上の点ごとに、README か /git-merge（nixfiles）を直すか、直さない理由を書いてから閉じる。

## 経緯

- 2026-11-07: 畳む規則と書き直しの前の手順を足したブランチのローカルのレビューで、直さなかった点から起こした
