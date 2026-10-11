---
status: accepted
date: 2026-10-11
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0121: OpenMC の参照解の出力は生成物として git で追跡しない

## 背景

OpenMC の参照解の出力（VERA 2B〜2I のピン出力 `examples/vera/p2/openmc-*.msgpack`、MOX-1 と MOX-2 の燃料の実効断面積
`examples/mox/openmc-*.msgpack`）を git で追跡していた（MOX は .gitignore の例外にして）。どれも `openmc_model.py` が OpenMC の
作業ディレクトリに書く `pins.msgpack`、`fuel-xs.msgpack` の写しで、作り直すと OpenMC の計算に 1 問題 16〜34 分かかる
（examples/vera/verification.md、examples/mox/verification.md）。同じ例の xlc のライブラリ（`vera.msgpack`）、xlc の
ピン出力、OpenMC の作業ディレクトリ（references/openmc/）は追跡していない。

- 生成物と入力を、追跡するかどうかで区別できること
- 中身を git の差分で読めないバイナリを、作り直すたびに履歴に足さないこと
- 比べる道具と手順を変えずに済むこと

## 検討した案

- 案 A: msgpack のまま追跡する（これまで）。作り直さずに比べられるが、生成物のうちこれだけが追跡され、差分を読めない
  バイナリが作り直すたびに履歴に残る
- 案 B: テキスト（TOML や CSV）にして追跡する。差分は読めるが、`openmc_model.py` と比べる道具に書き出しと読み込みを
  もう 1 つ足すことになり、生成物を追跡することは変わらない
- 案 C: 追跡をやめ、比べる道具を references/openmc/ の作業ディレクトリから直接読む形に直す。写しがなくなるが、
  道具の引数と検証の記録の手順を書き直すことになり、作業ディレクトリの名前（vera-2B と work-2F など）もそろっていない
- 案 D（採用）: 追跡をやめ、置き場所は変えない（.gitignore で無視する）

## 決定と理由

案 D にする。OpenMC の参照解の出力は、ほかの生成物と同じく git で追跡しない。置き場所は今のまま
（VERA は [ADR-0116](0116-place-vera-examples-by-problem.md) の `examples/vera/p<n>/openmc-<問題>.msgpack`、MOX は
`examples/mox/openmc-<問題>.msgpack`）とし、比べる道具と手順を変えない。

これまでに追跡したものは、履歴からも消す（ユーザーの判断）。生成物を履歴に残す理由はなく、書き換えの前に GitHub の
非公開リポジトリへ全ブランチを push して保険にする。書き換えで変わるコミットの番号は、記録の中で付け直す。

作り直しの時間は、追跡しない理由を覆すほどではない。参照解は OpenMC の作業ディレクトリ（references/openmc/、追跡しない）
にも残っており、手元で失っても手順（docs/environment.md、examples/mox/verification.md）で作り直せる。

## 確認

約束事なので、確かめるテストはない。`git ls-files '*.msgpack'` が空で、`git check-ignore` が参照解の出力を無視する
ことを確かめた。

## 帰結

- worktree では、xlc のライブラリと同じく、main の参照解の出力をコピーして使う
- 追跡をやめるマージで、main の作業ツリーからも参照解の出力が消える。マージの前に写しを取り、マージの後に戻す。
  VERA は references/openmc/ の作業ディレクトリの `pins.msgpack` にも残る。MOX は作業ディレクトリに残っていないので、
  失ったら examples/mox/verification.md の手順で作り直す
