---
status: accepted
date: 2026-10-07
supersedes: [ADR-0023]
superseded-by:
issues: ["#0002", "#0003"]
---

# ADR-0068: 共鳴の自己遮蔽は、ピンごとの等価なセルの超微細群（UFG）の減速計算で直接求める

## 背景

自己遮蔽を等価原理だけで扱う方式（[ADR-0023](0023-shield-by-bondarenko-and-equivalence.md)）は、UFG の重みのライブラリ（[ADR-0072](0072-weight-shielded-data-by-flux-calculator.md)）と組み合わせると U-238 の共鳴吸収を
大きく見積もった。

- 精度（評価済みデータだけから決まり、経験的な補正を持たないこと）
- 計算時間が実用的であること

## 検討した案

試して採らなかった方式は、どちらも方法そのものの限界で、文献でも同じ向きと大きさの誤差が報告されている（Gibson の MIT
博士論文 2016 の Table 9.72 と 6.7〜6.8）。

- 等価原理だけで σ0 を求める方式（[ADR-0023](0023-shield-by-bondarenko-and-equivalence.md)。ピンごと・群ごとの Dancoff 係数と Carlvik の2項有理近似）: UFG の重みの
  ライブラリと組み合わせると、U-238 の共鳴吸収を +4〜6% 多く見積もった（VERA 1B で −1887 pcm。NR の重みのライブラリの
  ときは誤差が打ち消し合って +232 pcm だった。docs/plans/m2-assembly.md、32fca1e と 42f67de、2026-10-07）。燃料の中の
  散乱源を平坦とすることと、有理近似による誤差である
- ESSM（埋め込み自己遮蔽法。SCALE/Polaris と同じ。ブランチ `feat/essm` に比較用に残す。main にはマージしない）:
  NR 近似の固定源の MOC Ω·∇ψ + Σ_i (Σa,i + Σp,i) ψ = Σ_i Σp,i / 4π の、ピンの燃料の平均の中性子束 φ から
  Σb = Σa φ / (1 − φ) とし、σ0,R = (Σb − N_R σp,R + Σ_{j≠R} N_j σa,j) / N_R でライブラリを引く（均質な媒質では厳密）。
  VERA 1B では U-238 の σ0 が約 52 b になり（Carlvik は約 46 b）、k が −3080 pcm だった（docs/plans/m2-assembly.md、
  a764ab0、2026-10-07）。脱出を群で平均した吸収（燃料は光学的に薄い）で評価するが、実際の共鳴吸収の大半は燃料が
  ほぼ黒体になる共鳴のピークで起きるため、脱出を大きく見積もる
- サブグループ法: ESSM の後に一時は進む方向としたが、λ（IR）や確率テーブルのフィットのような経験的な補正を持つ
- ピンごとの UFG の減速計算（MPACT の WSA/CPM、nTRACER の PSM、SCALE の CENTRM と同じ系統。採用）

## 決定と理由

共鳴の自己遮蔽は、ピンごとの超微細群（UFG）の減速計算で直接求める。調べた方式のうち、評価済みデータだけから決まり、
λ（IR）や確率テーブルのフィットのような経験的な補正を持たず、計算時間も実用的なもの。

- 集合体の黒体の燃料の MOC でピンごとの Dancoff 係数を求め、Dancoff 係数が一致するよう外径を決めた円筒の等価なセルを
  作る（燃料、ギャップ、被覆管、減速材のリング）。外径はピンごとに1つで、UFG の範囲の群ごとにその群の全断面積で
  合わせた外径の平均とする。Dancoff 係数は群ごとに変わるが（VERA 1B で 0.363〜0.388）、その違いはほぼ減速材の全断面積
  によるもので、合わせた外径の違いは 0.1% だった（0.6936〜0.6942 cm。正方形のセルと面積の等しい円の 0.711 cm より
  小さい。docs/plans/m2-assembly.md、a011cb8、2026-10-07）
- そのセルを、多リングの衝突確率法と、全核種の弾性散乱による減速（重心系で等方、標的は静止）で解く。セルの中の全核種
  （燃料、被覆管の Zr など）の共鳴自身の散乱、共鳴の干渉、リングごとの分布と温度が入る
- リングごとの中性子束で群定数に縮約する

## 確認

`fuel_rings_are_shielded_by_the_cell`、`equivalent_radius_recovers_the_cell`、`resonance_is_shielded_more_inside_the_fuel`、
`pins_next_to_water_are_shielded_less`。VERA 1B〜1D を同じ核データの OpenMC と比べる（コード同士の検証。
examples/vera/README.md）。等価なセルの誤差は、集合体全体を UFG の MOC で直接解く参照のツールで、案内管や集合体の端の
そばのピンで測る。OpenMC とは、中性子束の重みの群定数と直接比べる（反応率から逆算した σ0 には、MOC の角度の扱いの
誤差が混ざる）。

## 結果

- 範囲は [ADR-0069](0069-cap-cell-slowing-down-at-20-kev.md)、散乱核の近似は [ADR-0070](0070-treat-targets-at-rest.md)、MOC に合わせる SPH 因子は [ADR-0071](0071-correct-with-sph-factors.md)、ライブラリに加えるポイントワイズの断面積は [ADR-0073](0073-store-pointwise-data-in-library.md)、
  範囲の外は [ADR-0074](0074-shield-outside-range-by-equivalence.md)
- 名前の「UFG」は中身（ポイントワイズ）と合わない（#0003）
- lattice の格子の細かさは [ADR-0004](0004-keep-solver-settings-per-stage.md)、#0002
