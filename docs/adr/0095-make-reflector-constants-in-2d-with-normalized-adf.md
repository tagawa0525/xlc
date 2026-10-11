---
status: accepted
date: 2026-10-09
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0095: 径方向の反射体の定数は、不均質な燃料の集合体とバッフル・水を並べた 2D の格子計算で、面の向きの種類ごとに作り、反射体の ADF は燃料の ADF で割って持つ

## 背景

反射体の定数は、燃料の集合体と反射体を並べた格子計算から作る（[ADR-0084](0084-make-reflector-constants-with-fuel.md)）。炉心の外の構造は、バッフルと周りの
水までを反射体の定数に入れる（[ADR-0078](0078-put-baffle-in-reflector.md)）。ADF は GET で、core と同じ NEM で各ノードを解いて作る（[ADR-0094](0094-define-adf-by-one-node-nem.md)）。
残っているのは、格子の形（1D か 2D か、燃料の集合体の数、炉心の外周の階段の角）、反射体のノードの範囲と外側の境界、
燃料側と反射体側の ADF の使い方、軸方向の反射体（ノズル、炉心板、プレナム、上下の水）の作り方である。

- 炉心の外周の集合体の出力（VERA REF1-2D と 5-2D）
- 定数が、隣の燃料の種類にどれだけ左右されないか
- 格子計算の回数と、形状の手間（いまの格子は長方形のセルを持つ。[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）

## 検討した案

格子の形:

- 1D の平板（均質化した燃料、ギャップ、バッフル、水）: Scarabee はこの形で、燃料の集合体1体分の均質な燃料と反射体を
  S_N で解く（https://scarabee.readthedocs.io/en/latest/_modules/scarabee/reseau/reflector.html ）。CASMO も 1D の反射体の計算を持つ
  （J. Rhodes, K. Smith, D. Lee, "CASMO-5 Development and Applications," PHYSOR 2006、黒体の境界を反射体の計算のために
  加えたとある）。軽いが、燃料の不均質と、外周のピンの近くの熱中性子束の形が入らない。BWR の下部反射体では、CASMO5 の
  1D の均質化が熱群の吸収断面積を 26% 大きく、高速群を 62% 小さく見積もった（J. Loberg et al., Nucl. Technol. 177(1)
  (2012)、https://www.ans.org/pubs/journals/nt/article-13323 。要旨だけ読んだ）
- 不均質な燃料の集合体1体と、集合体に接する平板の反射体の 2D（Polaris の REFL の形状。SCALE 6.3 の Polaris の
  マニュアル、3.2.4.2 節と slab。燃料の集合体との並べ方は手元の写しで確かめていない）: 燃料の不均質が入る。平らな面
  だけを表す
- 2D の小さな問題を面の向きの種類ごとに作る: 燃料と反射体が1面で接する I 型は集合体2つ分、2面で接する L 型と、角だけで
  接する角型は 2×2 の集合体の問題とする（nTRACER-RENUS。J. H. Jeong, Y. S. Ban, H. G. Joo, "Establishment and Preliminary
  Verification of the nTRACER-RENUS Core Analysis System," Trans. KNS Spring Meeting (2017)、2.3 節）
- 多くの集合体（REF1-2D のような 3×3、または全炉心の 2D）を解いて、そこから反射体の定数を切り出す（DRAGON の例:
  N. Bejaoui, G. Marleau, "DRAGON Simulation of Radial PWR Reflector," CNS (2012)、要旨だけ読んだ。CASMO-5M は全炉心を
  解ける。MPACT は全炉心の計算から節点の定数を出す。MPACT 4.4 Theory Manual、13 章）: 環境は実際に近いが、定数が
  その炉心の装荷に依存する

反射体の ADF の使い方:

- 反射体の問題の燃料側と反射体側の ADF を、そのまま炉心で使う: 炉心で隣にある燃料は、反射体の問題の燃料と種類も状態も
  違うので、燃料側の ADF が合わない
- 反射体側の ADF を、同じ問題の燃料側の ADF で割って持ち（f* = f_反射体 / f_燃料）、炉心では隣の燃料の単一集合体の ADF を
  掛けて使う（採用）。面の連続の条件は ADF の比だけを決める（面ごとに自由度が1つ）ので、比を反射体の問題から、燃料側の
  大きさを実際に隣にある燃料から取る。Smith の方法として Scarabee が実装し（K. S. Smith, "Nodal diffusion methods and
  lattice physics data in LWR analyses: Understanding numerous subtle details," Prog. Nucl. Energy 101, 360–369 (2017)
  を引く。この論文は読んでいない）、nTRACER-RENUS も同じ式を使う（上の論文の式 (4)、(5)）

軸方向の反射体:

- 3D の格子計算: xlc の格子計算は 2D である
- 均質化した燃料（単一集合体の計算から、ライブラリの群で中性子束の重みで均質化したもの）と、軸方向の層を並べた 1D の
  平板を、格子計算で解く（長方形のセルを並べ、横の辺を鏡面反射にする）（採用）。SIMULATE-3 も、径方向と軸方向の
  反射体に不連続の考え方を使い、境界のアルベドを入力しない（NRC の SCE-9001 の安全評価、
  https://www.nrc.gov/docs/ML1332/ML13326A883.pdf ）。作り方の詳細は確かめていない

## 決定と理由

径方向の反射体の定数は、不均質な燃料の集合体とバッフル・水を並べた 2D の格子計算で、I 型、L 型、角型の面の向きごとに
作る。径方向の反射体のノードが要るのは REF1-2D、5-2D、Problem 5 である（Problem 3 と 4 は径方向の外周が鏡面反射で、
径方向の反射体を持たない。仕様書の Problem 3、4 の節）。REF1-2D は東側の1列だけで I 型で足りるが、5-2D と Problem 5 の
炉心の外周は階段の形で、L 型と角型の反射体のノードがあり、バッフルも角で折れる。多くの集合体を解いて切り出す案より、
定数が装荷に依存しない。

- I 型は、燃料の集合体1体と反射体のノード1つ（集合体のピッチ 21.50 cm の幅）。L 型は燃料3体と反射体1つ、角型は燃料1体と
  反射体3つの 2×2 とし、角型の定数は燃料と角だけで接するノードから取る。バッフルは VERA の仕様書の 1.13 節（燃料との
  隙間 0.19 cm、厚さ 2.85 cm、SS304）のとおり階段に沿って置く
- 燃料側の外周の辺と、反射体の横の辺は鏡面反射、反射体の外側の辺は真空（入射 0）とする。core でも反射体のノードの
  外側を入射の部分中性子流が 0 の境界とし、格子計算と炉心計算の境界をそろえる。core の NEM にはこの境界がまだなく
  （いまは正味の流れ 0 と中性子束 0 だけ）、加える必要がある（[ADR-0094](0094-define-adf-by-one-node-nem.md)。[ADR-0099](0099-couple-core-nodes-by-get-with-vacuum-boundary.md) で加え、[ADR-0101](0101-impose-transport-current-ratio-at-vacuum-faces.md) で格子計算の J/φ を課す形にした）。反射体のノードより外の水とバレルは入れない（[ADR-0078](0078-put-baffle-in-reflector.md)。比べ方の偏りは #0016）
- 反射体のノードは、ギャップ、バッフル、水をまとめて中性子束と体積で重みづけして均質化する。拡散係数は
  [ADR-0104](0104-weight-transport-cross-sections-by-currents.md) のとおり、反射体のノード（核分裂がない）では Koebke の等価理論の D、燃料のノードでは中性子流で
  重みづけした outflow の D とする。反射体のノードは core と1ノードの問題の両方で ANM で解く（[ADR-0102](0102-solve-non-fissile-nodes-with-anm.md)）
- 反射体の ADF は、[ADR-0103](0103-define-adf-by-one-node-kernels.md) の方法で燃料のノードと反射体のノードの両方について求め、燃料と接する面ごとに f*= f_反射体 /
  f_燃料 を持つ。core では、その面の反射体側の ADF を f* × （隣の燃料の、単一集合体の格子の ADF）とし、燃料側は燃料自身の
  単一集合体の ADF を使う。反射体のノードどうしの面の ADF は 1 とする（同じ物質どうしの面で、比べる参照の差が小さいと
  見込む。5-2D で確かめる）。外側の面（入射 0 の境界）の ADF は、ADR-0103 の1ノードの解から作り、core の境界の条件で
  均質な面の中性子束に掛けて使う（入射 0 の条件は面の中性子束を含むので、等価性にはその面の ADF が要る）
- 格子計算は、入力でノード（格子のセルの長方形の範囲。格子を重ならずに覆う）を並べ、ノードごとに1つの少数群定数の
  ファイルを書く（`few-group` の `nodes`）。1つのファイルに複数のノードを入れる形にしないのは、core がいまの単一の
  少数群定数の形をそのまま読め、分岐のテーブル（[ADR-0097](0097-tabulate-branches-on-a-full-grid.md)）でも燃料と反射体は別のテーブルになるためである。
  ノードの辺の不均質な中性子束と正味の流れは、MOC の粗メッシュ（格子のセル）の面ごとの値を辺にわたって平均する。
  面の中性子束は、領域の体積と同じく、トラックから推定した面の長さで割る。1ノードの問題の内側の辺は、core の隣の
  ノードと同じく、隣のノードの幅と参照の横方向の漏れで2次式を作る。反射体のノードの印（`reflector`）を付けたノードは、
  反射体でないノードと接する辺の ADF を f*で持ち、少数群定数の `normalized_adf` にその辺を記す。core はその辺の ADF を
  f* × 隣のノードの同じ面の ADF とする
- 燃料の集合体は、その位置に実際にある種類（VERA 5 では外周の集合体の濃縮度と Pyrex の配置）の代表を1つ選ぶ。f* で
  割るので、燃料の種類への依存は小さいと見込み、5-2D で確かめる

軸方向の反射体（下: 炉心板、ボトムノズル、端栓とギャップ、上: プレナム、端栓、トップノズル、炉心板、水）は、均質化した
燃料と、層ごとの均質な物質を並べた 1D の平板で作る。層の物質は VERA の仕様書の混ぜた組成（ノズルと炉心板は仕様書の
表 5 と混合の表）とし、プレナムのようにピンのある層も層ごとに均質化する。仕様書の Problem 3 の模型に合わせて、炉心板
の外に 50 cm の水を置き、外側を真空（入射 0）とする。層の境の ADF は ADR-0103 の方法で軸方向に作り、燃料と接する面は径方向と同じく
燃料側の ADF で割って持つ。燃料の層どうしの軸方向の面の ADF は 1 とする（燃料の格子は軸方向に一様）。

## 確認

- 格子のノードごとの均質化（`crates/lattice/src/fewgroup.rs`）: 1つの物質の無限媒質を2つのノードに分けると物質のまま
  で ADF は 1（`nodes_of_an_infinite_medium_are_the_medium`）、燃料と水の帯で、ノードの領域だけの反応率の保存、燃料の
  ノードの中性子流の重みの D、隣のノードの辺を持つ1ノードの問題の ADF、真空の辺の J/φ、反射体のノードの f*と
  `normalized_adf`、反射体のノードの x の2つの面の不連続因子が等しいこと（`strip_nodes_make_get_constants_and_normalized_adfs`）、
  ノードの誤り（`invalid_nodes_are_errors`）。MOC の面の中性子束と正味の流れ（`faces_of_an_infinite_medium_see_the_mean_flux`、
  `faces_agree_with_sides_and_mirror_lattices`）、入力（`parses_few_group_nodes`、`invalid_few_group_nodes_are_errors`）、
  少数群定数の形式（`normalized_sides_are_checked`）、core の f* の面（`normalized_adfs_take_the_neighbor_adf`、
  `normalized_adfs_need_a_plain_neighbor`）、隣のノードの辺の1ノードの問題（`neighbor_edges_shape_the_transverse_leakage`、
  core の `get_adfs_reproduce_inner_faces_with_transverse_leakage`）
- I 型の格子（examples/reflector/verification.md の「I 型の反射体の GET の等価性と f*の方法の誤差」、8dbf48a、
  2026-10-09。自己遮蔽なし）: core の NEM（燃料）と ANM（反射体）は、格子計算の MOC の k を +1.5 pcm、ノード平均の中性子束を
  10⁻⁴ 以内で再現した。残る差は MOC の鏡面反射の上の辺のノードごとの正味の流れ（#0032）による。燃料のノードを単一集合体の
  定数にした f* の方法では k が −129 pcm、ノード平均の中性子束が最大 1.4% 違った
- REF1-2D を、core の NEM（I 型の反射体）と、同じ多群定数の xlc の MOC（REF1-2D の 3×3 の形状をそのまま解く）、
  OpenMC（2D）、KENO の値（仕様書の付録 I。k = 0.993677 ± 0.000021）で比べる。判定の閾値（案）: xlc の MOC に対して
  k ±50 pcm、集合体出力の差の最大 1%、外周のピン出力（ピン出力の再構成のあと）の相対の差の RMS 2%（ユーザーの決定が要る）。
  まだしていない（core の装荷図、#0031）
- 5-2D で、L 型と角型の反射体の効きと、反射体のノードどうしの ADF を 1 とする近似と、代表の燃料を選ぶ近似の効き（別の燃料で
  作った f* との差）を測る（#0029）
- 軸方向は Problem 3 を OpenMC（3D）と KENO の値と比べる（#0030）

## 結果

- I 型だけを実装した。L 型と角型の格子、core の反射体のノードどうしの面の ADF を 1 とする扱いと、反射体の定数を炉心の
  位置に合わせて回す扱いは #0029、軸方向の反射体は #0030、core に種類の違うノードを並べる入力（装荷図）は #0031
- 反射体のノードの均質な解は、NEM でも、ANM で中性子束の重みの D でも、真空の面で負になった。ANM と Koebke の等価理論の
  D で正になり、ADF を作れる（ADR-0102、ADR-0104。examples/reflector/verification.md）
- f* の方法の誤差（I 型で −129 pcm）は、燃料のノードの定数の差から来る。REF1-2D で、燃料が反射体の片側にしかない体系の
  誤差を測る
- 反射体の問題の格子は外周が真空になるので、Dancoff 係数の計算での真空の境界の扱いが先に要る（#0008）。I 型の例は
  自己遮蔽なしで解いた
- 反射体の定数はホウ素濃度と減速材の温度（水の密度）に強く依存するので、分岐のテーブルに入れる（[ADR-0097](0097-tabulate-branches-on-a-full-grid.md)）
- 5-2D の KENO の値はバレルなどを含むので、比べ方の偏りが入る（#0016）
- 均質化した燃料で作る軸方向の反射体は、燃料の上下の端のピンの不均質を表さない。Problem 3 の上下の端の出力の差が大きい
  ときに見直す
- 細部は、xlc の格子計算（MOC）の結果との比較を通じて段階的に詰め、新しい ADR か issue に記録する（2026-10-08 の
  ユーザーの方針）
