---
status: open
labels: [tooling]
opened: 2026-10-31
closed:
milestone:
blocked-by: []
resolution:
---

# dvc checkout の hook が dvc を実行しない形を止める

## 問題

`.claude/hooks/block-raw-dvc-checkout.py`（ADR-0223）は、引用の外の `$(…)`、バッククォート、`(` で単純なコマンドを区切るので、
置換の後の語が `dvc checkout` で始まると、dvc を実行しないコマンドも止める。feat/dvc-checkout-hook の 1 周目の /code-review
（2026-10-31）が、`echo $(date) dvc checkout`、``echo `date` dvc checkout``、`a=(dvc checkout)` で `runs_raw_dvc_checkout` が
True を返すことを示した。また、出力に触れない `dvc checkout --help` と `-h` も止め、`scripts/dvc-checkout.py` が dvc checkout に
渡す引数を読む手段がない（`nix develop -c dvc checkout --help` も止まる）。

2 周目（同日）は、bash の値を取るオプションのうち、テストの無いもの（`-O`、`--rcfile`、`--init-file`、`-c -- '…'`）を
挙げた。

どれも普段の書き方ではないので、レビューの基準（直すのは普段の書き方の誤り）で直さなかった。

関係する ADR: ADR-0223

## 完了の条件

置換を語の一部として読み（中のコマンドは並びに入れ、外のコマンドを区切らない）、`--help` と `-h` だけの `dvc checkout` を
通す。テストに上の例と、bash の値を取るオプションの例を足す。

## 経緯

- 2026-10-31: feat/dvc-checkout-hook の 1 周目と 2 周目のレビューの指摘から書いた
