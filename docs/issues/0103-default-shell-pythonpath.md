---
status: open
labels: [verification]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# default の devShell が外の shell の PYTHONPATH を外さない

## 問題

references の devShell は shellHook で `PYTHONPATH` を外す（ADR-0124）が、default の devShell は外さない。外の shell に
`PYTHONPATH` が残っていると、検証用の Python（`withPackages`）も DVC のステージの python3 も、そのパッケージを site-packages
より先に読む。ステージが検証用でない版のパッケージで黙って走り、数値が変わりうる。

調べた #0101 で見つけた。2026-10-12 のジョブ 247 で、dvc のラッパーの素の python3 が msgpack を読めた理由の候補だが、
確かめていない。

関係する ADR: なし

## 完了の条件

default の devShell で外の `PYTHONPATH` を外すか、外さない理由を決めて記録する。flake の check `dvc-stage-python`
（#0101）に、`PYTHONPATH` を置いた場合を足すかも決める。

## 経緯

- 2026-10-13: #0101 のレビューで、閉じる issue の中にだけ書いていたので分けた（fix/0101-dvc-stage-python）
