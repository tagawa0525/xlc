---
type: Issue
status: closed
labels: [xs, endf]
opened: 2026-11-20
closed: 2026-11-21
milestone:
blocked-by: []
resolution: コミット aa4a2ae19
---

# 対数の補間則が、ln をとる側の正でない値で NaN を黙って返す

## 問題

- ENDF の補間則の式（`Interpolation::interpolate`、crates/xs/src/endf/tab1.rs）は、ln をとる補間則（INT=3 lin-log、
  INT=4 log-lin、INT=5 log-log）で、ln をとる側の y が負か、符号の違う 2 点か、ln をとる側の x が正でないとき、NaN を
  返す。NaN はどの比較でも偽なので、後の検査（ADR-0267 の 0 K の負の値の検査など）もすり抜ける。ENDF/B-VII.0 の Cf-253 の
  群定数の NaN（#0320）は、非分離共鳴の遮蔽した捕獲が 8500 eV で正、10 keV で負の区間を log-log で補間した値が直接の
  原因（別のセッションの調べ。このブランチでは式を読んだだけで、Cf-253 を流して確かめていない）
- 読み込みの TAB1 は `Tab1::check` が負の y と正でない x を止めていたが、計算の途中で作る表（非分離共鳴の自己遮蔽の表、
  LSSF=0 の非分離共鳴の範囲の 0 K の和）と、TAB1 でない値の並び（MF4 と出射分布の入射エネルギーの間）は、この検査を
  通らずに式に渡っていた
- ENDF-102（references/papers/endf-102-bnl-224854-2023-inre.txt の 0.5 節の表 16）は INT=3〜5 を「ln(x) か ln(y) が線形」と
  定義するだけで、正でない値の扱いを書かない。0 の値は、INT=6 の説明が「σ1 か σ2 が 0 なら結果は明らかに 0」と書き、
  MF7 の説明が「S = 0.0 はとても小さい数を表す」と書く。NJOY の terp1（references/NJOY2016/src/endf.f90）は、log-log で
  y1 = 0 なら 0 を返し、y2 = 0 なら exp(−∞) で 0 になる。負の値は扱わない

関係する ADR: なし

## 対応

- 定義域は `Interpolation::check_x`（ln x をとる補間則で両端の x が正）と `check_interval`（それに加えて、ln y をとる補間則で
  両端の y が負でない）の 1 か所に置いた。NaN の端は定義域の外とする。y の端点が 0 の区間は、ln y → −∞ の極限として 0 のまま残す（上の ENDF-102 と
  NJOY と同じ。根拠は式の注釈）
- 式に値を渡す経路（`interpolate` を呼ぶ 11 か所）と、値を作る所の検査:
  - `Tab1` の補間（`eval`、`ClampedTab1`、MF4 の表の μ、MF6 の `limits`、`linearize`）: 補間則つきの `Tab1` を作る所は 4 つで、
    どれも `Tab1::check` を通す。読み込みの TAB1（`Reader::tab1`）、MF6 の LAW=1 の表、しきい値に合わせて作り直した表
    （`pointwise` の threshold shift）、LSSF=0 の非分離共鳴の範囲の平均とバックグラウンドの和（`replace_with_unresolved`。負の和を
    log-log で線形化すると NaN だった）。ほかの `Tab1` は `Tab1::linear`（lin-lin で、定義域の制限が無い）
  - MF4 の入射エネルギーの間（Legendre 係数と表の密度）: 入射エネルギーを `check_x` で確かめる（`check_energy_laws`）。係数の
    ln y の補間則と負の密度は前から読み込みで止めている
  - MF5 と MF6 の出射分布の入射エネルギーの間（crates/xs/src/multigroup/outgoing.rs の群ごとの確率、単位基底の範囲の両端と
    密度）: `check_incident_domain`（mf5.rs）が、入射エネルギー、両端の分布の最小の密度、単位基底では出射エネルギーの範囲の
    両端を確かめる。
    MF5 の `parse_outgoing` と、MF6 の LAW=1 と LAW=7 の `OutgoingDistributions` を作る所から呼ぶ
  - 非分離共鳴の自己遮蔽の補正（`UnresolvedShielding::adjust`）: `UnresolvedShielding::new` が、σ0 ごとと無限希釈の表の成分ごとに
    確かめ、σ0 の番号、成分、格子点、値を示す
- `interpolate` は f64 を返すまま。中でパニックにする形と `Result` にする形は採らなかった（下の「検討して捨てた形」）

## 走査の結果

`survey_log_interpolation`（crates/xs/examples/survey_log_interpolation.rs。xlc 02be59700、2026-11-20）で、4 つの版の中性子と
熱中性子散乱のファイルの全部を数えた。レコードは行の形で見つける（例の先頭の注釈）。入れ子の TAB2（MF6 の LAW=7）の部分
レコードの値と、単位基底の出射エネルギーの範囲の両端は数えていない。この 2 つは、読み込みの走査（下）が `check_incident_domain`
を通して確かめた。

| 版     | 中性子のファイル | ln をとる区間のある TAB1 / TAB2 | y < 0 | x ≤ 0 | y = 0                                                       |
| ------ | ---------------- | ------------------------------- | ----- | ----- | ----------------------------------------------------------- |
| VII.0  | 393              | 16446 / 10                      | 0     | 0     | 1（S-34 の MF3 MT1。4.8e5 eV の不連続点で、補間しない区間） |
| VII.1  | 423              | 16184 / 6                       | 0     | 0     | 0                                                           |
| VIII.0 | 557              | 9447 / 5                        | 0     | 0     | 0                                                           |
| VIII.1 | 558              | 7007 / 5                        | 0     | 0     | 0                                                           |

- 熱中性子散乱のファイル（20、21、34、114 個）は、MF7 MT4 の S(α, β) の log-lin の表だけに y = 0 があり（ln をとる区間のある
  TAB1 のうち、VII.0 で 1088、VII.1 で 1105、VIII.0 で 4114、VIII.1 で 13994）、負の y と正でない x は無い。S(α, β) は
  `interpolate` を通さず、ln S の表に直して扱う（crates/xs/src/multigroup/sab.rs）
- 読み込みで新しく止まる評価は無い。読み込みの走査（read-survey のステージ）をこのブランチで流し直し、4 つの版の報告は main と
  同じだった（xlc e7c82cb18、2026-11-21）
- 新しく止まりうるのは、非分離共鳴の自己遮蔽の表に負の値ができる評価で、読み込みの走査では数えられない（`xs nuclide` の走査
  （#0292 の process-survey）で数える。このブランチでは流していない）。前はその表の補間が NaN になり、群定数に黙って入っていた。
  負の値のでき方は #0319、Cf-253 は #0320
- LSSF=0 の 0 K の和が負の評価は、main でも格子点の負の値を ADR-0267 の検査が止めていた（このブランチの合成のテストを main の
  コードで流すと、ADR-0267 のエラーで止まった）ので、止まる評価は増えない。止まる理由の文が、補間則の定義域の文に変わる
- 評価の誤りとして読み替える仕組み（ADR-0272 の `corrections`、ADR-0267 の `negative-to-zero`）は、読み込みで止まる評価が
  無いので要らない。計算の途中の表の負の値は評価の値の誤りでないので、これらでは読み替えられない（#0319 が扱う）
- 既存の核種の出力: VIII.0 の U-235、U-238、Zr-90 の 600 K の 1 点の `xs nuclide` の出力は、main（3d17e4e6a）のビルドと
  バイトまで同じだった（xlc 937899d6d、2026-11-21）

## 完了の条件

- [x] 合成の入力のテストで、ln をとる補間則の定義域の外の区間を持つ表がエラーで止まり、`interpolate` が NaN を返さない
- [x] 4 つの版の全ファイルで、ln をとる補間則の区間に正でない値を持つ表を数え、既定をエラーにして止まる評価を挙げる
- [x] 既存の核種の出力が変わらない（`xs nuclide` の出力のバイトが main のビルドと同じ）
- [x] 読み込みの走査（read-survey）で新しく止まる評価が無いかを確かめる

## 検討して捨てた形

- `interpolate` の中で、定義域の外ならパニックで止める: 本番のコードの panic と assert はリポジトリの lint（clippy.toml の
  `disallowed-macros`、Cargo.toml の `[workspace.lints.clippy]` の `panic`）が禁じ、エラーを返す形を求める
- `interpolate` を `Result` にする: 呼ぶ所が 11 か所あり、f64 を返す関数（MF4 の密度、`ClampedTab1::eval`、出射分布の
  `accumulate`）の中から呼ぶので、全部の呼び出しの形が変わる。どの経路も値を作る所で検査でき、そこでエラーにすればどの表の
  どの区間かを示せるので、値を作る所の検査にした

## 経緯

- 2026-11-20: Cf-253 の NaN（#0320）の直接の原因として、別のセッションの担当が見つけた
- 2026-11-21: 4 つの版を走査し、値を作る所で定義域を確かめるようにして閉じた
