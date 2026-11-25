---
type: Issue
status: open
labels: [verification, vera]
opened: 2026-11-25
closed:
milestone:
blocked-by: []
resolution:
---

# VERA の問題 4、5 を ENDF/B-VII.0 の検証用ライブラリで解いて KENO-VI と比べるかを決める

## 問題

ADR-0264（proposed）は、本番の核データを ENDF/B-VII.1 にし、仕様書の KENO-VI（ENDF/B-VII.0）と比べる VII.0 の検証用ライブラリ
（examples/vera/endfb70）とその比べ方（ADR-0175）は変えないとした。その検証用のステージに、問題 4、5 が無い。

- VII.0 の検証用のステージは examples/vera/endfb70/dvc.yaml（ライブラリと 1B〜1D、2A〜2K、2Q）と examples/vera/endfb70/p3/dvc.yaml（3A、3B）
  だけで、ADR-0175 の決定と理由が KENO-VI と比べる問題に挙げるのもこの組
- 仕様書は、問題 4（4-2D を含む）と 5（5-2D を含む）と REF1-2D にも KENO-VI の参照解を持つ（ADR-0175 の背景）。REF1-2D は本番にも
  ステージが無い（examples/reflector/verification.md は REF1-2D の比較を #0043 に残す）
- 本番の 5-2D（examples/vera/p5/dvc.yaml）は、本番のライブラリ（今は VIII.0）の xlc の解を、KENO-VI（VII.0）と OpenMC（5A-2D。同じ VIII.0）と
  比べている（examples/vera/verification.md の「Problem 5-2D の炉心の計算」）。KENO-VI との差には核データの版の差が入り、本番を VII.1 に
  した後も入る
- 問題 4 は本番にもステージが無い。問題 4、5 の制御棒の価値を KENO-VI と比べる作業は #0098 にあり、どの版のライブラリで比べるかは
  書いていない

決めることは、問題 4、5（と REF1-2D）を VII.0 の検証用ライブラリでも解いて KENO-VI と比べるか（ステージを足すか）、本番の版の解を同じ版の
OpenMC と比べるだけにするか（5A-2D の今の OpenMC の参照解は VIII.0 で、本番を VII.1 にすると版が合わなくなる）。ADR-0175 の決定と理由は
3D の全炉心と問題 4 の OpenMC を作らないとしているので、それらを後者にするなら ADR-0175 の「## 変更」が要る。

関係する ADR: ADR-0175、ADR-0264

## 完了の条件

- [ ] 問題 4、5（と REF1-2D）の KENO-VI との比較に使うライブラリを決めて ADR に書く（ADR-0175 の「## 変更」か新しい ADR）。VII.0 で解くと決めたら、
  examples/vera/endfb70 にステージを足す

## 経緯

- 2026-11-25: ADR-0264 を案 C（本番は VII.1、KENO-VI との比較は VII.0 の検証用ライブラリ）に書き直したとき、検証用のステージに
  問題 4、5 が無いことを分けて起こした
