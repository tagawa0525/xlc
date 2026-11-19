---
status: open
labels: [docs, terms]
opened: 2026-11-19
closed:
milestone:
blocked-by: []
resolution:
---

# xlc の語の辞書を置き、nominal と branch の訳語を業界の語に直す

## 問題

燃焼のテーブルの nominal と branch を「公称」「分岐」と訳してきたが、「公称」は炉物理の業界の語でなく、Gd の集合体の
燃焼点の報告で何を指すか伝わらなかった。「分岐」も、git の分岐点やプログラムの条件分岐と同じ語で、専門用語と一般の語の
区別が付かない。語を業界で通じる語に揃え、崩れた語の辞書で止める。語は xlc でしか使わないので、辞書はリポジトリの
`.ja-terms/prh.yml`（claude-harness の ADR-0096）に置き、共通の辞書にある xlc でしか使わない語の規則（claude-harness の #316 の
23 規則）も移す（ユーザー 2026-11-19）。

直す語（ユーザー 2026-11-19。英語の元の語で示す。訳した旧い語は `.ja-terms/prh.yml` の規則にある）: nominal state →
基準条件、off-nominal state → ブランチ条件、nominal depletion → ベース計算、nominal table → ベース計算のテーブル、
branch calculation → ブランチ計算、branch table → ブランチのテーブル、history branch → 履歴ブランチ、instantaneous
branch → 瞬時ブランチ。同じ対応の表の案のうち、nominal boron → 基準ホウ素濃度、no branches → ブランチなし、
control-rod（温度、密度、ホウ素）branch → 制御棒（温度、密度、ホウ素）のブランチも、ユーザーに見せて異論が出なかった
（2026-11-19）。branch states（`[branches]` の軸の直積で、基準条件を含む）は、基準条件を除くブランチ条件と分けて
「ブランチ計算の状態」とした（レビューで、`burnup` と `off-nominal-burnup` の説明が同じ語になると分かったため）。

機械的な置き換えはできないので、辞書は前後の語で専門用語と決まる形だけを止める。残りは次のとおり:

- 辞書が止めない形（「公称の点」「公称の組成」「公称の <数> ppm」「履歴からの分岐」「1 状態の分岐」、単独の「分岐」など）は、
  文脈で直す語が変わる。行ごとに決めて直す
- 例の入力（TOML）は、ほかの入力のコピーやスクリプトが書いて手で直さないものが多く、元と一緒に直さないと食い違う。DVC の
  ステージの deps（入力の TOML と python-ast の対象の .py）は、直すとステージが changed になる。入力の TOML、それを書く
  スクリプト、deps の .py と、その出力の文を確かめるテストの 555 ファイルは `.textlintignore` で検査から外してあり、元ごとに
  直し、deps はステージごとに付け直しの条件（ADR-0212 の (b)）か再実行を決める
- 記録の見出しを新しい語に直したので、外した入力の注釈が引く見出し（例: examples/depletion/p9-gd/branch-2O-all.toml が
  ADR-0201 の自己遮蔽の節を引く旧い見出し、examples/depletion/branches/shielding-*.toml）は今の見出しと食い違う。入力を直すときに
  見出しの引用も直す
- claude-harness の共通の辞書から 23 規則を消す（claude-harness の #316。この辞書がマージされた後）

関係する ADR: なし

## 完了の条件

- [ ] 辞書が止めない形の行を、直すか一般の語と決める
- [ ] 例の入力とそれに結びつく 555 ファイルを直して `.textlintignore` から外す
- [ ] `ja-terms check` が、追跡しているファイル（docs/issues/0142 を除く）で 0 件になる

## 経緯

- 2026-11-19: リポジトリの辞書を置き、辞書が止める形を例の入力の外で直した。例の入力とそれに結びつく 555 ファイルは検査から外した
