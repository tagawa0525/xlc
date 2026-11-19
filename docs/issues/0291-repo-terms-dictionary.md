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

- [x] 辞書が止めない形の行を、直すか一般の語と決める
- [x] 例の入力とそれに結びつく 555 ファイルを直して `.textlintignore` から外す
- [ ] `ja-terms check` が、追跡しているファイル（docs/issues/0142 を除く）で 0 件になる

## 経緯

- 2026-11-19: リポジトリの辞書を置き、辞書が止める形を例の入力の外で直した。例の入力とそれに結びつく 555 ファイルは検査から外した
- 2026-11-19: 辞書が止めない「公称」「分岐」の残りの .md の 814 行を決めた（直したものと keep。ブランチ docs/0291-remaining-terms-md）。決められない 16 行は `ask` のまま残し、親が決める
- 2026-11-19: 例の入力とそれに結びつく 555 ファイル、crates の注釈、残りの .py を直し、`.textlintignore` から外した
  （ブランチ chore/0291-terms-in-example-inputs）。入力を書くスクリプトの文字列から直し、コピーは元と同じ行で直した。直した TOML は
  どれも直す前と同じ値に読める。.md の担当が決められなかった ADR-0215、ADR-0225、ADR-0250 と examples/depletion/verification.md
  の行も直した。docs/issues/0198 の「分岐が減る」はプログラムの条件分岐、ADR-0200 の「変更」の前の文の引用は辞書が止めない
  ので、どちらもそのまま。docs/issues/0279 と 0289 の 3 行は、測定の担当が直しているファイルなので残した
- 2026-11-19: 同じブランチの DVC のステージ。直したファイルを deps に持つステージは 192（dvc.lock の deps のパスと `git diff --name-only main`
  の照合）。OpenMC の参照解の openmc-1C@rings1、@rings4、openmc-mox1 は ADR-0212 の (b) で付け直し、数秒のステージ（比較 13、
  OpenMC の extract 3、python-ast 6）は再実行した（出力の md5 は extract、hzp-compare とも変わらず、比較の JSON は入力の md5 だけが
  変わった）。残りの 167 は、lock の md5 と今のファイルが違う deps に bin/lattice（110）、bin/core（53）、bin/xs（4）があり、main で
  既に changed なので、(b) の条件 2 に当たらず changed のまま残した。比較の JSON の入力の md5 の食い違いを止める検査は #0295
- 2026-11-19: 同じブランチで、生成器を再実行すると値が変わる生成物があった: issue-0285 の production-2A-60-*.toml の 4 本（元の
  p9-2A/branch-2A.toml の燃料の温度の点が後で 4 点になった）、issue-0180/reuse と issue-0279 の生成物（読む
  p9-product/lattice-2A.msgpack が後で変わった）、issue-0045 と issue-0129 の生成物（元の入力が後で変わった）。追跡している入力が
  記録の出どころなので、これらは再実行せず、注釈の行が生成器と一致することだけを確かめた
