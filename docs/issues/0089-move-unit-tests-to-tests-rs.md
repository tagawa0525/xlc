---
status: closed
labels: [refactor]
opened: 2026-10-12
closed: 2026-10-14
milestone:
blocked-by: []
resolution: コミット e0323140
---

# 単体テストを tests.rs に移し、本番のファイルのテストのコードを検査で止める

## 問題

main（1842033d、2026-10-12）で、`crates/*/src` の 64 個のファイルがテストをファイルの中に持つ。テストだけの項目が
本番のファイルに 4 個ある（crates/format/src/fewgroup.rs の `sample`、`sample_pins`、crates/numeric/src/pins/node.rs の
`modes_at`、`flux`）。置き方と、移したあとに成り立つべきことは ADR-0143 で決めた。確認の方法と検査の置き場所は、移す実装のブランチで決める。

関係する ADR: [ADR-0143](../adr/0143-put-unit-tests-in-tests-rs.md)

## 完了の条件

- `run-checks.sh` で走る検査が、本番のファイルの中のテストのコードを見つけて失敗し、移したあとの main で通る
- 64 個のファイルのテストと 4 個の項目が `tests.rs` にあり、移す前後で release の実行ファイル（位置の情報を除く）と
  実行されるテストの集まりが変わらないことを確かめ、ADR-0143 の「確認」に実際の確認の方法と結果を書いて、
  `status: accepted` と `date` をその日にしてある
  （docs/adr/README.md の「運用」）
- build の段階が作る strip した実行ファイルの md5 が、移す前の main の `dvc.lock` と一致する（ADR-0159、ADR-0160）。一致すれば計算の段階は古くならず、`dvc repro -P` は流さない
- proposed の ADR と開いた issue に書いたテストの場所のパスが、`tests.rs` に直してある

## 経緯

- 2026-10-12: コードのレビューで、大きなファイルの半分以上がテストと分かった。ADR-0143 で置き方を決めた
- 2026-10-13: refactor/0089-move-tests で移した（crate ごとの移動は xs 32、lattice 20、numeric 5、format 4、docs 2、core 1 のファイル）。
  確認の方法と結果は ADR-0143 の「確認」に書いた。検査の置き場所は、`.pre-commit-config.yaml` の hook と cargo test（consistency.rs）の
  両方にした。proposed の ADR と開いた issue には、テストの場所をパスで書いたものが crates/core/src/rods/tests.rs の 1 件
  （すでに `tests.rs`）だけで、直すものはなかった
