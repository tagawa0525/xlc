//! 円筒のセル（同心円のリング、外周は白色境界）の衝突確率法。
//!
//! 弦の高さ y で積分する。ある弦の上で、区間 a に一様な等方線源があるとき、その先の区間 b で初めて衝突する数は
//! Bickley 関数 Ki3 の2階差分 D = Ki3(τab) − Ki3(τab+τa) − Ki3(τab+τb) + Ki3(τab+τa+τb) で表せる（τab は a と b の
//! 間の光学的距離）。弦の向きを逆にした区間の対は、弦が中心について対称なので、鏡映した対と同じ値になる。
//! そこで向きを一つに固定し、高さを [0, R] に限って2倍する:
//!
//! - 区間の対: V_iΣ_iP_ij = 2∫dy Σ_{a∈i が b∈j より前} D
//! - 同じ区間: 2∫dy [τa − Ki3(0) + Ki3(τa)]
//! - 外への脱出: 2∫dy [Ki3(τ後) − Ki3(τ後+τa)]
//!
//! 一つの区間からの3つの和は τa に打ち消し合う。そこで脱出はこの和の残りとして求め（下の積分表示と表の値の
//! どちらを使っても厳密になる）、体積も同じ求積で 2∫dy Σ_a ℓa として求めておけば（MOC で追跡の体積を使うのと
//! 同じ）、中性子の保存と相反則は y の求積の精度によらず丸め誤差の範囲で成り立つ。
//!
//! 真空のリング（Σ = 0）を扱うため、衝突の代わりに飛程（単位の線源密度あたりの、リングの中を飛ぶ長さ）
//! A_ij = V_iP_ij/Σ_j を求める。それには上の差分を Σ で割る。τ が小さいと差分が打ち消し合い、Ki3 の表の誤差が
//! 1/τ² で増えるので、そのときは差分を τ で割った形の積分表示（[`trapezoid`]）で直接求める。Σ = 0 はその極限になる。

use std::f64::consts::PI;
use std::sync::OnceLock;
use xlc_numeric::quadrature::{NodesError, gauss_legendre};

/// y の求積で、小区間（リングの間の区間を、外の境界の特異点に合わせて分けたもの）ごとに使う Gauss–Legendre の点の数
const POINTS: usize = 16;

/// 円筒のセルの計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum CpError {
    #[error("ring radii must be positive and increasing: {0:?}")]
    Radii(Vec<f64>),
    #[error("the y quadrature needs at least one point per ring")]
    Points,
    #[error(transparent)]
    Nodes(#[from] NodesError),
    #[error("{given} values for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("ring {0}: total cross sections must be finite and not negative")]
    Total(usize),
    #[error("ring {0}: sources must be finite")]
    Source(usize),
    #[error("every ring is void")]
    Void,
}

/// 弦の上の区間（リングの番号と長さ）と、その弦の y の求積の重み
struct Chord {
    weight: f64,
    segments: Vec<(usize, f64)>,
}

/// 同心円のリングからなる円筒のセル。
pub struct Cylinder {
    radii: Vec<f64>,
    chords: Vec<Chord>,
    volumes: Vec<f64>,
}

impl Cylinder {
    /// `radii` はリングの外半径 (cm) を内側から並べたもの。
    pub fn new(radii: Vec<f64>) -> Result<Self, CpError> {
        Self::with_points(radii, POINTS)
    }

    /// 小区間ごとに `points` 点の求積で弦をとる。
    pub fn with_points(radii: Vec<f64>, points: usize) -> Result<Self, CpError> {
        let increasing = radii.windows(2).all(|w| w[0] < w[1]);
        if radii.is_empty() || !increasing || !radii.iter().all(|r| r.is_finite() && *r > 0.0) {
            return Err(CpError::Radii(radii));
        }
        if points == 0 {
            return Err(CpError::Points);
        }
        let nodes = gauss_legendre(points)?;
        let mut chords = Vec::new();
        let mut inner = 0.0;
        for (k, &outer) in radii.iter().enumerate() {
            let mut push = |y: f64, weight: f64| {
                chords.push(Chord {
                    weight,
                    segments: segments(&radii, k, y),
                })
            };
            // 外のリングの境界 r_{k+1} にも弦の長さの √(r_{k+1} − y) の特異点があり、それが区間の上端に近いと
            // 求積が収束しにくい。上端から特異点までの距離 g = r_{k+1} − r_k を幅とする小区間をとり、その下は幅を
            // 倍々にする。どの小区間でも、特異点までの距離が幅を下回らない
            let gap = radii.get(k + 1).map_or(f64::INFINITY, |r| r - outer);
            // 上端の小区間は y = r_k − d t² として、弦の長さの √(r_k − y) の振る舞いを t について滑らかにする
            let top = gap.min(outer - inner);
            for &(x, w) in &nodes {
                let t = 0.5 * (x + 1.0);
                push(outer - top * t * t, w * top * t);
            }
            let (mut upper, mut width) = (outer - top, top);
            while top < outer - inner && upper > inner {
                width *= 2.0;
                let lower = (upper - width).max(inner);
                let half = 0.5 * (upper - lower);
                for &(x, w) in &nodes {
                    push(lower + half * (x + 1.0), w * half);
                }
                upper = lower;
            }
            inner = outer;
        }
        let mut volumes = vec![0.0; radii.len()];
        for chord in &chords {
            for &(ring, length) in &chord.segments {
                volumes[ring] += 2.0 * chord.weight * length;
            }
        }
        Ok(Self {
            radii,
            chords,
            volumes,
        })
    }

    /// リングの体積（単位の高さあたりの面積, cm²）。衝突確率と同じ求積で求めたもので、π(r_k² − r_{k−1}²) とは
    /// 求積の誤差だけ違う。
    pub fn volumes(&self) -> Vec<f64> {
        self.volumes.clone()
    }

    /// 白色境界の衝突確率 P[i][j]（リング i で生まれた中性子が、リング j で初めて衝突する確率）。
    /// 真空のリングの列は 0 になる。
    pub fn probabilities(&self, totals: &[f64]) -> Result<Vec<Vec<f64>>, CpError> {
        let white = self.tracks(totals)?;
        let volumes = &self.volumes;
        Ok(white
            .iter()
            .zip(volumes)
            .map(|(row, v)| row.iter().zip(totals).map(|(w, t)| w * t / v).collect())
            .collect())
    }

    /// リングごとの等方線源 `sources` (1/cm³/s) に対するスカラー中性子束。真空のリングの中性子束も求める。
    pub fn fixed_source(&self, totals: &[f64], sources: &[f64]) -> Result<Vec<f64>, CpError> {
        self.check(sources.len())?;
        if let Some(i) = sources.iter().position(|q| !q.is_finite()) {
            return Err(CpError::Source(i));
        }
        let white = self.tracks(totals)?;
        let volumes = &self.volumes;
        Ok((0..self.radii.len())
            .map(|j| {
                let tracks: f64 = sources.iter().zip(&white).map(|(q, row)| q * row[j]).sum();
                tracks / volumes[j]
            })
            .collect())
    }

    fn check(&self, given: usize) -> Result<(), CpError> {
        let rings = self.radii.len();
        if given != rings {
            return Err(CpError::Rings { given, rings });
        }
        Ok(())
    }

    /// 白色境界の飛程 W[i][j]（リング i の単位の線源密度あたりの、リング j の中を飛ぶ長さの合計）。対称になる。
    pub fn tracks(&self, totals: &[f64]) -> Result<Vec<Vec<f64>>, CpError> {
        self.check(totals.len())?;
        if let Some(i) = totals.iter().position(|t| !t.is_finite() || *t < 0.0) {
            return Err(CpError::Total(i));
        }
        let (mut tracks, escapes) = self.vacuum(totals);
        // 外周から入る中性子（余弦分布）が初めて衝突する確率 P_sj = 4V_jΣ_jP_js/S、飛程は 4e_j/S
        let surface = 2.0 * PI * self.radii[self.radii.len() - 1];
        let entering: Vec<f64> = escapes.iter().map(|e| 4.0 * e / surface).collect();
        // 外周から入って、どこかで衝突する確率 1 − P_ss（P_ss から引くと、透明に近いセルで桁が落ちる）
        let collides: f64 = entering.iter().zip(totals).map(|(l, t)| l * t).sum();
        if collides == 0.0 {
            return Err(CpError::Void);
        }
        for (row, e) in tracks.iter_mut().zip(&escapes) {
            for (w, l) in row.iter_mut().zip(&entering) {
                *w += e * l / collides;
            }
        }
        Ok(tracks)
    }

    /// 外周が真空のときの飛程 A[i][j] と、単位の線源密度あたりの脱出の数 e[i] = V_iP_is。
    fn vacuum(&self, totals: &[f64]) -> (Vec<Vec<f64>>, Vec<f64>) {
        let n = self.radii.len();
        let mut tracks = vec![vec![0.0; n]; n];
        let mut escapes = vec![0.0; n];
        // 弦の上の区間の境界の対 (p, q) の間の光学的距離 τ_p + … + τ_{q−1} の Ki3 を [p * longest + q] に置く。
        // 光学的距離は区間 p から順に足す（弦の始めから足した和の差にすると、厚い区間の後で丸めが大きくなる）。
        // `between`（区間 a+1 から b−1 まで順に足したもの）は (a+1, b) の和と同じ順に足している
        let longest = self
            .chords
            .iter()
            .map(|c| c.segments.len())
            .max()
            .unwrap_or(0)
            + 1;
        let mut memo = vec![0.0; longest * longest];
        let mut memo2 = vec![f64::NAN; longest * longest];
        // 弦ごとに作り直さず再利用する
        let mut path: Vec<Path> = Vec::with_capacity(longest);
        let mut collided: Vec<f64> = Vec::with_capacity(longest);
        let mut taus: Vec<f64> = Vec::with_capacity(longest);
        for chord in &self.chords {
            let weight = 2.0 * chord.weight;
            path.clear();
            path.extend(chord.segments.iter().map(|&(ring, length)| Path {
                ring,
                length,
                total: totals[ring],
                tau: totals[ring] * length,
            }));
            let last = path.len() - 1;
            // Ki2 は真空の区間のある弦でだけ使う
            if path.iter().any(|p| p.tau == 0.0) {
                memo2[..(path.len() + 1) * longest].fill(f64::NAN);
            }
            taus.clear();
            taus.extend(path.iter().map(|p| p.tau));
            // 区間の対 (a, b) は a + b ≤ L − 1 なので、使う境界の対 (p, q) は q ≤ L + 1 − p（下のループの4項）。
            // 光学的距離は q について増えるだけなので、表の範囲を超えたら先はどれも 0
            for p in 0..=path.len() {
                let row = &mut memo[p * longest..(p + 1) * longest];
                let end = (path.len() + 1 - p).min(path.len());
                let mut optical = 0.0;
                for q in p..=end {
                    if q > p {
                        optical += taus[q - 1];
                    }
                    if optical >= NEGLIGIBLE {
                        row[q..=end].fill(0.0);
                        break;
                    }
                    row[q] = ki3(optical);
                }
            }
            let k3 = |p: usize, q: usize| -> f64 { memo[p * longest + q] };
            // 真空と組になる対の Ki2 も、境界の対ごとに1回だけ求める（区間 b が真空のときの引数 (a, b) は区間 a から
            // 足すので、between + τa と足す順が違い、最後の桁が違いうる）
            let mut k2 = |p: usize, q: usize| -> f64 {
                let slot = &mut memo2[p * longest + q];
                if slot.is_nan() {
                    let optical: f64 = taus[p..q].iter().sum();
                    *slot = small(ki2, optical);
                }
                *slot
            };
            // 区間ごとの、その区間から前に飛んで衝突する数 Σ_b Σb A_ab（自己の項を含む）
            collided.clear();
            collided.resize(path.len(), 0.0);
            for (a, from) in path.iter().enumerate() {
                let own = from.own();
                tracks[from.ring][from.ring] += weight * own;
                collided[a] += own * from.total;
                // 弦が中心について対称なので、区間の対 (a, b) と鏡映した対 (L−1−b, L−1−a) は同じ値になる。
                // a + b ≤ L − 1 の対だけを求め、a + b < L − 1 の対は鏡映の分も A[i][j] の側に足す
                let mut between = 0.0;
                for (b, to) in path.iter().enumerate().take(last + 1 - a).skip(a + 1) {
                    // 光学的距離は弦に沿って増えるだけなので、表の範囲を超えたら先の対はどれも無視できる（[`NEGLIGIBLE`]。
                    // Ki3 の差分の最初の項 (a+1, b) の引数も同じ値）
                    if between >= NEGLIGIBLE {
                        break;
                    }
                    // Ki3 の差分の4項は、境界の対 (a+1, b)、(a, b)、(a+1, b+1)、(a, b+1) の値
                    // 片方が真空なら、D/(τa τb) は真空の側の τ → 0 の極限で [Ki2(τab) − Ki2(τab + τ)]/τ（τ は他方）。
                    // τab + τ は、a が真空なら境界の対 (a+1, b+1)、b が真空なら (a, b) の光学的距離
                    let value = if from.direct(to) {
                        (k3(a + 1, b) - k3(a, b) - k3(a + 1, b + 1) + k3(a, b + 1))
                            / (from.total * to.total)
                    } else if from.tau == 0.0 && to.tau >= VOID_TABLE {
                        from.length * to.length * (k2(a + 1, b) - k2(a + 1, b + 1)) / to.tau
                    } else if to.tau == 0.0 && from.tau >= VOID_TABLE {
                        from.length * to.length * (k2(a + 1, b) - k2(a, b)) / from.tau
                    } else {
                        from.thin(to, between)
                    };
                    collided[a] += value * to.total;
                    if a + b < last {
                        tracks[from.ring][to.ring] += 2.0 * weight * value;
                        collided[last - b] += value * from.total;
                    } else {
                        tracks[from.ring][to.ring] += weight * value;
                    }
                    between += to.tau;
                }
            }
            // 衝突しなければ外へ出る（ℓa = Σ_b Σb A_ab + e_a。表と積分表示のどちらの値でも保存則が厳密になる）
            for (from, c) in path.iter().zip(&collided) {
                escapes[from.ring] += weight * (from.length - c);
            }
        }
        // 鏡映の分を片側に足したので、平均で対称にする
        for i in 1..n {
            let (above, below) = tracks.split_at_mut(i);
            for (j, row) in above.iter_mut().enumerate() {
                let mean = 0.5 * (below[0][j] + row[i]);
                below[0][j] = mean;
                row[i] = mean;
            }
        }
        (tracks, escapes)
    }
}

/// 弦の上の区間
struct Path {
    ring: usize,
    length: f64,
    total: f64,
    tau: f64,
}

/// Ki3 の差分の大きさの目安（τ を 1 で頭打ちにした積）がこれより小さいと、差分を積分表示で直接求める。
/// Ki3 の表の誤差（2e-15 程度）が、差分の値の 1e-8 程度に収まる境目
const CANCELLATION: f64 = 1e-6;

/// 真空の区間と、光学的厚さがこれ以上の区間の対は、Ki2 の差分を表で求める（`Cylinder::vacuum`）。Ki2 の表の誤差は
/// 1e-3 以下で 4e-10、それより上で 1e-12 程度で、真空に隣り合う対の Ki2(0) は表の点の値なので、差分の値の
/// 1e-8 程度に収まる
const VOID_TABLE: f64 = 1e-2;

/// 光学的距離がこれ以上なら Ki3、Ki2 を 0 とする（表の上端。Ki2(48) ≈ 2e-22 で、ほかの項（τ や Ki3(0) = π/4 の程度）
/// に足しても倍精度の丸めより小さい。表の外を台形則で積分し直す手間を省く）
const NEGLIGIBLE: f64 = MAX;

impl Path {
    /// 単位の線源密度あたりの、同じ区間の中を飛ぶ長さ: [τ − Ki3(0) + Ki3(τ)]/Σ²
    fn own(&self) -> f64 {
        if self.tau == 0.0 {
            // τ → 0 の極限: ℓ² Ki1(0)/2
            return self.length.powi(2) * PI / 4.0;
        }
        if self.tau >= NEGLIGIBLE {
            return (self.tau - ki3(0.0)) / (self.total * self.total);
        }
        if self.tau.min(1.0).powi(2) < CANCELLATION {
            // [τ − Ki3(0) + Ki3(τ)]/τ² = ∫ φ2(τ cosh t)/cosh t dt
            // 指数の因子がない（B = 0）
            self.length.powi(2) * trapezoid(coarse_stride(0.0), |c| phi2(self.tau * c) / c)
        } else {
            (self.tau - ki3(0.0) + ki3(self.tau)) / (self.total * self.total)
        }
    }

    /// 区間の対の D/(Σa Σb) を、Ki3 の差分そのままで求めてよいか（打ち消し合いが小さい。真空の区間を含めば積は 0 で
    /// 偽）。偽なら、真空と厚い区間の対は Ki2 の差分（`Cylinder::vacuum`）、ほかは積分表示（[`Path::thin`]）で求める
    fn direct(&self, to: &Path) -> bool {
        self.tau.min(1.0) * to.tau.min(1.0) >= CANCELLATION
    }

    /// Ki3 の差分では打ち消し合い（[`Path::direct`] が偽）、真空と光学的に厚い区間の対でもない、薄い区間の対の、
    /// 単位の線源密度あたりの、光学的距離 `between` だけ先の区間 `to` の中を飛ぶ長さ D/(Σa Σb)
    fn thin(&self, to: &Path, between: f64) -> f64 {
        // D/(τa τb) = ∫ exp(−τab cosh t) φ1(τa cosh t) φ1(τb cosh t)/cosh t dt
        let integral = trapezoid(coarse_stride(between), |c| {
            (-between * c).exp() * phi1(self.tau * c) * phi1(to.tau * c) / c
        });
        self.length * to.length * integral
    }
}

/// 引数が [`NEGLIGIBLE`] 以上なら 0 とした Ki_n
fn small(ki: fn(f64) -> f64, x: f64) -> f64 {
    if x >= NEGLIGIBLE { 0.0 } else { ki(x) }
}

/// φ1(z) = (1 − exp(−z))/z
fn phi1(z: f64) -> f64 {
    if z == 0.0 { 1.0 } else { -(-z).exp_m1() / z }
}

/// φ2(z) = (z − 1 + exp(−z))/z²。z が小さいと分子が打ち消し合うので、級数で求める（打ち切りの誤差は z⁵/5040）
fn phi2(z: f64) -> f64 {
    if z < 1e-2 {
        0.5 - z * (1.0 / 6.0 - z * (1.0 / 24.0 - z * (1.0 / 120.0 - z / 720.0)))
    } else {
        (z + (-z).exp_m1()) / (z * z)
    }
}

/// 高さ y の弦が横切る区間を、弦に沿って順に並べる。y はリング k の中（r_{k−1} ≤ y < r_k）にある。
fn segments(radii: &[f64], k: usize, y: f64) -> Vec<(usize, f64)> {
    // 弦の半分の長さ √(r² − y²) の差を、桁落ちしないよう (r_m² − r_{m−1}²)/(h_m + h_{m−1}) で求める
    let half = |r: f64| ((r - y) * (r + y)).max(0.0).sqrt();
    let mut outer: Vec<(usize, f64)> = Vec::with_capacity(radii.len() - k);
    for m in k + 1..radii.len() {
        let (r0, r1) = (radii[m - 1], radii[m]);
        outer.push((m, (r1 - r0) * (r1 + r0) / (half(r1) + half(r0))));
    }
    let mut chord: Vec<(usize, f64)> = outer.iter().rev().copied().collect();
    chord.push((k, 2.0 * half(radii[k])));
    chord.extend(outer);
    chord
}

/// 表の区間: [0, 1) は細かく（Ki の 3 階以上の導関数が 0 で発散するため）、[1, MAX) は粗く
const FINE: f64 = 8192.0;
const COARSE: f64 = 1024.0;
const MAX: f64 = 48.0;

/// Ki_n（n = 2, 3）の表。点ごとに、値と、刻みを掛けた導関数 −Ki_{n−1} × 刻み（補間の式の係数）を並べる
struct Bickley {
    fine: [Vec<[f64; 2]>; 2],
    coarse: [Vec<[f64; 2]>; 2],
}

fn table() -> &'static Bickley {
    static TABLE: OnceLock<Bickley> = OnceLock::new();
    TABLE.get_or_init(|| {
        let nodes = |n: i32, count: f64, start: f64, step: f64| -> Vec<[f64; 2]> {
            (0..=count as usize)
                .map(|i| {
                    let x = start + i as f64 * step;
                    [bickley(n, x), -bickley(n - 1, x) * step]
                })
                .collect()
        };
        let both = |count: f64, start: f64, step: f64| [2, 3].map(|n| nodes(n, count, start, step));
        Bickley {
            fine: both(FINE, 0.0, 1.0 / FINE),
            coarse: both((MAX - 1.0) * COARSE, 1.0, 1.0 / COARSE),
        }
    })
}

/// Bickley 関数 Ki3(x) = ∫_0^{π/2} cos²θ exp(−x/cosθ) dθ
#[inline]
pub fn ki3(x: f64) -> f64 {
    interpolate(x, 3)
}

/// Ki2(x)
fn ki2(x: f64) -> f64 {
    interpolate(x, 2)
}

/// Ki_n(x)（n = 2, 3）を、表の値と導関数 −Ki_{n−1} の3次 Hermite 補間で求める
#[inline]
fn interpolate(x: f64, n: usize) -> f64 {
    let table = table();
    // 刻みは2の冪なので、刻みで割る代わりに点の数を掛けても同じ値になる
    let (nodes, start, per) = if x < 1.0 {
        (&table.fine[n - 2], 0.0, FINE)
    } else if x < MAX {
        (&table.coarse[n - 2], 1.0, COARSE)
    } else {
        return bickley(n as i32, x);
    };
    let s = (x - start) * per;
    let i = (s as usize).min(nodes.len() - 2);
    let t = s - i as f64;
    let ([f0, d0], [f1, d1]) = (nodes[i], nodes[i + 1]);
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * f0
        + (t3 - 2.0 * t2 + t) * d0
        + (-2.0 * t3 + 3.0 * t2) * f1
        + (t3 - t2) * d1
}

/// Ki_n(x) = ∫_0^∞ exp(−x cosh t)/cosh^n t dt
fn bickley(n: i32, x: f64) -> f64 {
    trapezoid(1, |c| (-x * c).exp() / c.powi(n))
}

/// cosθ = 1/cosh t とした、t についての積分 ∫_0^∞ f(cosh t) dt の台形則。f は t について減る関数とする。
/// 被積分関数は t について偶で、|Im t| < π/2 で解析的なので、台形則の誤差は exp(−π²/h) で減る
/// （x が大きいと exp(−x cosh t) の幅が狭まるが、表の範囲 x < 48 では相対誤差 1e-13 程度に収まる）。
///
/// 刻みは 0.125 × `stride`（[`coarse_stride`]）。
fn trapezoid(stride: usize, f: impl Fn(f64) -> f64) -> f64 {
    const STEP: f64 = 0.125;
    // 1/cosh t は exp(−t) で減るので、t = 40 で 1e-17 を下回る
    const NODES: usize = 320;
    static COSH: OnceLock<Vec<f64>> = OnceLock::new();
    let cosh = COSH.get_or_init(|| (0..=NODES).map(|k| (k as f64 * STEP).cosh()).collect());
    let mut sum = 0.5 * f(1.0);
    for &c in cosh.iter().step_by(stride).skip(1) {
        let term = f(c);
        sum += term;
        // すべてアンダーフローしたときも止まるよう、等号を含める
        if term <= 1e-17 * sum {
            break;
        }
    }
    sum * STEP * stride as f64
}

/// 被積分関数 exp(−B cosh t) g(cosh t)/cosh t（g は帯の上で有界。φ1、φ2 など）の積分を台形則で求める刻みの倍率。
/// 帯 |Im t| < π/2 の端では exp(−B cosh t) が 1 の程度なのに、積分は e^−B の程度なので、刻み h の相対誤差は
/// およそ e^B exp(−π²/h)（g が積分を小さくする分だけ相対では大きくなる）になる。B < 4 なら刻み 0.25 で 3e-14 以下
/// （B = 3.99、τ = 3 の φ1 で 1.5e-14。打ち消し合いの目安 [`CANCELLATION`] で許す差分の誤差 1e-8 よりずっと小さい）。
/// それより大きいと刻み 0.125 のまま
fn coarse_stride(b: f64) -> usize {
    if b < 4.0 { 2 } else { 1 }
}

#[cfg(test)]
mod tests;
