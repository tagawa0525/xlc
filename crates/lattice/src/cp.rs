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
//! 一つの区間からの3つの和は、Ki3 の値によらず τa に打ち消し合うので、中性子の保存と相反則は y の求積の精度に
//! よらず丸め誤差の範囲で成り立つ。
//!
//! 真空のリング（Σ = 0）を扱うため、衝突の代わりに飛程（単位の線源密度あたりの、リングの中を飛ぶ長さ）
//! A_ij = V_iP_ij/Σ_j を求める。真空のリングに入る項は、Σ → 0 の極限として Ki2 で表す。

use std::f64::consts::PI;
use std::sync::OnceLock;

/// y の求積で、リングの間の区間ごとに使う Gauss–Legendre の点の数
const POINTS: usize = 16;

/// 円筒のセルの計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum CpError {
    #[error("ring radii must be positive and increasing: {0:?}")]
    Radii(Vec<f64>),
    #[error("the y quadrature needs at least one point per ring")]
    Points,
    #[error("Gauss-Legendre nodes for {0} points did not converge")]
    Nodes(usize),
    #[error("{given} values for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("ring {0}: total cross sections must be finite and not negative")]
    Total(usize),
    #[error("ring {0}: sources must be finite, and zero in void rings")]
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
}

impl Cylinder {
    /// `radii` はリングの外半径 (cm) を内側から並べたもの。
    pub fn new(radii: Vec<f64>) -> Result<Self, CpError> {
        Self::with_points(radii, POINTS)
    }

    /// リングの間の区間ごとに `points` 点の求積で弦をとる。
    pub fn with_points(radii: Vec<f64>, points: usize) -> Result<Self, CpError> {
        let increasing = radii.windows(2).all(|w| w[0] < w[1]);
        if radii.is_empty() || !increasing || !radii.iter().all(|r| r.is_finite() && *r > 0.0) {
            return Err(CpError::Radii(radii));
        }
        if points == 0 {
            return Err(CpError::Points);
        }
        let nodes = gauss_legendre(points)?;
        let mut chords = Vec::with_capacity(radii.len() * points);
        let mut inner = 0.0;
        for (k, &outer) in radii.iter().enumerate() {
            // y = r_k − (r_k − r_{k−1}) t² で、弦の長さの √(r_k − y) の振る舞いを t について滑らかにする
            let width = outer - inner;
            for &(x, w) in &nodes {
                let t = 0.5 * (x + 1.0);
                let y = outer - width * t * t;
                chords.push(Chord {
                    weight: w * width * t,
                    segments: segments(&radii, k, y),
                });
            }
            inner = outer;
        }
        Ok(Self { radii, chords })
    }

    /// リングの体積（単位の高さあたりの面積, cm²）。
    pub fn volumes(&self) -> Vec<f64> {
        let mut inner = 0.0;
        self.radii
            .iter()
            .map(|&r| {
                let v = PI * (r * r - inner * inner);
                inner = r;
                v
            })
            .collect()
    }

    /// 白色境界の衝突確率 P[i][j]（リング i で生まれた中性子が、リング j で初めて衝突する確率）。
    /// 真空のリングの列は 0 になる。
    pub fn probabilities(&self, totals: &[f64]) -> Result<Vec<Vec<f64>>, CpError> {
        let white = self.white(totals)?;
        let volumes = self.volumes();
        Ok(white
            .iter()
            .zip(&volumes)
            .map(|(row, v)| row.iter().zip(totals).map(|(w, t)| w * t / v).collect())
            .collect())
    }

    /// リングごとの等方線源 `sources` (1/cm³/s) に対するスカラー中性子束。真空のリングの中性子束も求める。
    pub fn fixed_source(&self, totals: &[f64], sources: &[f64]) -> Result<Vec<f64>, CpError> {
        self.check(sources.len())?;
        for (i, (q, t)) in sources.iter().zip(totals).enumerate() {
            // 真空のリングどうしの飛程は求めていないので、真空に線源を置けない
            if !q.is_finite() || (*t == 0.0 && *q != 0.0) {
                return Err(CpError::Source(i));
            }
        }
        let white = self.white(totals)?;
        let volumes = self.volumes();
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
    /// 真空のリングどうしの値は求めない（0 のまま）。
    fn white(&self, totals: &[f64]) -> Result<Vec<Vec<f64>>, CpError> {
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
        for chord in &self.chords {
            let weight = 2.0 * chord.weight;
            let path: Vec<Path> = chord
                .segments
                .iter()
                .map(|&(ring, length)| Path {
                    ring,
                    length,
                    total: totals[ring],
                    tau: totals[ring] * length,
                })
                .collect();
            let mut rest: f64 = path.iter().map(|p| p.tau).sum();
            for (a, from) in path.iter().enumerate() {
                rest -= from.tau;
                escapes[from.ring] += weight * from.escape(rest.max(0.0));
                if !from.is_void() {
                    tracks[from.ring][from.ring] += weight * from.own();
                }
                let mut between = 0.0;
                for to in &path[a + 1..] {
                    tracks[from.ring][to.ring] += weight * from.to(to, between);
                    between += to.tau;
                }
            }
        }
        // 弦が中心について対称なので、A[i][j] と A[j][i] は丸め誤差を除いて等しい
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

impl Path {
    fn is_void(&self) -> bool {
        self.total == 0.0
    }

    /// 単位の線源密度あたりの、同じ区間の中を飛ぶ長さ: [τ − Ki3(0) + Ki3(τ)]/Σ²
    fn own(&self) -> f64 {
        (self.tau - ki3(0.0) + ki3(self.tau)) / (self.total * self.total)
    }

    /// 単位の線源密度あたりの、光学的距離 `between` だけ先の区間 `to` の中を飛ぶ長さ
    fn to(&self, to: &Path, between: f64) -> f64 {
        match (self.is_void(), to.is_void()) {
            (false, false) => {
                let d = ki3(between) - ki3(between + self.tau) - ki3(between + to.tau)
                    + ki3(between + self.tau + to.tau);
                d / (self.total * to.total)
            }
            // Σ → 0 の極限
            (false, true) => to.length * (ki2(between) - ki2(between + self.tau)) / self.total,
            (true, false) => self.length * (ki2(between) - ki2(between + to.tau)) / to.total,
            // 真空どうしは使わない（真空に線源を置かない）
            (true, true) => 0.0,
        }
    }

    /// 単位の線源密度あたりの、光学的距離 `rest` の先で外へ出る数
    fn escape(&self, rest: f64) -> f64 {
        if self.is_void() {
            self.length * ki2(rest)
        } else {
            (ki3(rest) - ki3(rest + self.tau)) / self.total
        }
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

/// [−1, 1] の n 点の Gauss–Legendre の点と重み（Newton 法で Legendre 多項式の零点を求める）。
fn gauss_legendre(n: usize) -> Result<Vec<(f64, f64)>, CpError> {
    let mut nodes = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut converged = false;
        for _ in 0..100 {
            let (p, dp) = legendre(n, x);
            let step = p / dp;
            x -= step;
            if step.abs() <= 1e-15 {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(CpError::Nodes(n));
        }
        let (_, dp) = legendre(n, x);
        nodes.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    Ok(nodes)
}

/// Legendre 多項式 P_n(x) とその導関数
fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, x);
    for k in 2..=n {
        let k = k as f64;
        (p0, p1) = (p1, ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k);
    }
    if n == 0 {
        return (1.0, 0.0);
    }
    (p1, n as f64 * (x * p1 - p0) / (x * x - 1.0))
}

/// 表の区間: [0, 1) は細かく（Ki の 3 階以上の導関数が 0 で発散するため）、[1, MAX) は粗く
const FINE: f64 = 8192.0;
const COARSE: f64 = 256.0;
const MAX: f64 = 48.0;

/// 表の点での (Ki1, Ki2, Ki3)
struct Bickley {
    fine: Vec<[f64; 3]>,
    coarse: Vec<[f64; 3]>,
}

fn table() -> &'static Bickley {
    static TABLE: OnceLock<Bickley> = OnceLock::new();
    TABLE.get_or_init(|| {
        let nodes = |count: f64, start: f64, step: f64| -> Vec<[f64; 3]> {
            (0..=count as usize)
                .map(|i| bickley(start + i as f64 * step))
                .collect()
        };
        Bickley {
            fine: nodes(FINE, 0.0, 1.0 / FINE),
            coarse: nodes((MAX - 1.0) * COARSE, 1.0, 1.0 / COARSE),
        }
    })
}

/// Bickley 関数 Ki_n(x) = ∫_0^{π/2} cos^{n−1}θ exp(−x/cosθ) dθ を、表の値と導関数 −Ki_{n−1} の3次 Hermite 補間で求める
fn interpolate(x: f64, n: usize) -> f64 {
    let table = table();
    let (nodes, start, step) = if x < 1.0 {
        (&table.fine, 0.0, 1.0 / FINE)
    } else if x < MAX {
        (&table.coarse, 1.0, 1.0 / COARSE)
    } else {
        return bickley(x)[n - 1];
    };
    let s = (x - start) / step;
    let i = (s as usize).min(nodes.len() - 2);
    let t = s - i as f64;
    let (a, b) = (&nodes[i], &nodes[i + 1]);
    let (f0, f1) = (a[n - 1], b[n - 1]);
    let (d0, d1) = (-a[n - 2] * step, -b[n - 2] * step);
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * f0
        + (t3 - 2.0 * t2 + t) * d0
        + (-2.0 * t3 + 3.0 * t2) * f1
        + (t3 - t2) * d1
}

/// Ki3(x)
pub fn ki3(x: f64) -> f64 {
    interpolate(x, 3)
}

/// Ki2(x)
fn ki2(x: f64) -> f64 {
    interpolate(x, 2)
}

/// (Ki1, Ki2, Ki3)(x) を、cosθ = 1/cosh t として Ki_n(x) = ∫_0^∞ exp(−x cosh t)/cosh^n t dt の台形則で求める。
/// 被積分関数は t について偶で、|Im t| < π/2 で解析的なので、台形則の誤差は exp(−π²/h) で減る。
fn bickley(x: f64) -> [f64; 3] {
    const STEP: f64 = 0.125;
    let mut sum = [0.0; 3];
    let mut k = 0;
    loop {
        let t = k as f64 * STEP;
        let c = t.cosh();
        let f1 = (-x * c).exp() / c;
        let factor = if k == 0 { 0.5 } else { 1.0 };
        let terms = [f1, f1 / c, f1 / (c * c)];
        for (s, f) in sum.iter_mut().zip(terms) {
            *s += factor * f;
        }
        if f1 < 1e-18 * sum[0] {
            break;
        }
        k += 1;
    }
    sum.map(|s| s * STEP)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn bickley_function_has_known_values() {
        // Ki3(0) = π/4、Ki3 は単調に減る
        assert!(close(ki3(0.0), PI / 4.0, 1e-12));
        // Ki3(x) = ∫_0^{π/2} cos²θ exp(−x/cosθ) dθ（SciPy の適応求積で求めた値）
        assert!(close(ki3(1.0), 0.237_845_082_192_855, 1e-8), "{}", ki3(1.0));
        assert!(close(ki3(2.0), 0.076_963_590_311_658, 1e-8), "{}", ki3(2.0));
        assert!(ki3(2.0) < ki3(1.0) && ki3(20.0) > 0.0);
    }

    /// 燃料・ギャップ・被覆管・水のセル
    const RADII: [f64; 4] = [0.4096, 0.418, 0.475, 0.711];

    #[test]
    fn probabilities_are_conserved_and_reciprocal() {
        let totals = [0.6, 0.0, 0.3, 1.2];
        let p = Cylinder::new(RADII.to_vec())
            .unwrap()
            .probabilities(&totals)
            .unwrap();
        let volumes = Cylinder::new(RADII.to_vec()).unwrap().volumes();
        for i in [0, 2, 3] {
            // 白色境界では、どこかで必ず衝突する
            let sum: f64 = (0..4).map(|j| p[i][j]).sum();
            assert!(close(sum, 1.0, 1e-8), "row {i}: {sum}");
            for j in [0, 2, 3] {
                let (a, b) = (
                    volumes[i] * totals[i] * p[i][j],
                    volumes[j] * totals[j] * p[j][i],
                );
                assert!(close(a, b, 1e-8), "{i},{j}: {a} vs {b}");
            }
            // 真空では衝突しない
            assert_eq!(p[i][1], 0.0);
        }
    }

    #[test]
    fn homogeneous_cell_is_an_infinite_medium() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let flux = cell.fixed_source(&[0.8; 4], &[2.4; 4]).unwrap();
        for f in &flux {
            assert!(close(*f, 3.0, 1e-8), "{f}");
        }
    }

    #[test]
    fn void_ring_carries_the_surrounding_flux() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let flux = cell
            .fixed_source(&[0.8, 0.0, 0.8, 0.8], &[2.4, 0.0, 2.4, 2.4])
            .unwrap();
        for i in [0, 2, 3] {
            assert!(close(flux[i], 3.0, 1e-8), "{i}: {}", flux[i]);
        }
    }

    #[test]
    fn rejects_bad_radii() {
        assert!(Cylinder::new(vec![]).is_err());
        assert!(Cylinder::new(vec![0.5, 0.4]).is_err());
        assert!(Cylinder::new(vec![0.0, 0.4]).is_err());
    }
}
