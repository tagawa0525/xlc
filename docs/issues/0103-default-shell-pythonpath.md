---
type: Issue
status: closed
labels: [verification]
opened: 2026-10-13
closed: 2026-11-19
milestone:
blocked-by: []
resolution: コミット 05e9a2460
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
- 2026-11-19: 閉じた。default の devShell の shellHook で `PYTHONPATH` を外した（`PYTHONPATH=/nonexistent` で入ると、前は残り、今は外れる。
  検証用の Python は msgpack と yaml を読める）。flake の check `dvc-stage-python` には足さない: check は nix のサンドボックスで外の環境を
  持たず、外の `PYTHONPATH` は shell に入るときの問題で、shellHook が扱う
