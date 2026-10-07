//! R-matrix limited 形式（LRF=7）の Reich-Moore の断面積（0 K）。
//!
//! 捕獲は消去されたチャネルとして共鳴エネルギーの虚部に入れ（E_λ − iΓγ/2）、開いた中性子のチャネルの
//! R 行列から X = (I − iK)⁻¹K（K = P^½ R P^½）を求める。衝突行列は U = e^{−i(φ+φ')} (I + 2iX)。
//! 透過因子は剛体球のもので、シフト因子は使わない。閾値より下の非弾性散乱のチャネルは閉じているので除く
//! （NJOY2016 samm と同じ）。
//!
//! 捕獲の幅がない共鳴のエネルギーちょうどでは R 行列が発散する（極）が、X は有限の極限を持つ。
//! M = I − iK とすると X = i(I − M⁻¹) で、極の共鳴の項を K = K' + w wᵀ/d と分けると、d → 0 の極限は
//! M⁻¹ = M'⁻¹ − (M'⁻¹w)(wᵀM'⁻¹) / (wᵀM'⁻¹w)（Sherman-Morrison の公式）。NJOY はこの点で発散する。

use std::f64::consts::PI;

use num_complex::Complex64;

use super::factors::{penetrability, phase_shift};
use crate::endf::{OrbitalL, ParticlePair, RMatrixLimited};
use crate::physics::wave_number_factor;

/// R-matrix limited の断面積の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum RmlError {
    #[error("particle pair {pair}: {what} is not supported")]
    Pair { pair: usize, what: String },
    #[error("{0} elastic particle pairs (MT=2); exactly one is needed")]
    ElasticPairs(usize),
    #[error("more than one particle pair for MT{0}")]
    DuplicateMt(u32),
    #[error("spin group {group}: more than one capture channel")]
    CaptureChannels { group: usize },
    #[error(
        "spin group {group}: the resonance at {er} eV has a width in channel {channel}, \
         which is closed at that energy"
    )]
    ClosedWidth {
        group: usize,
        er: f64,
        channel: usize,
    },
    #[error("spin group {group}: several resonances without capture width at {energy} eV")]
    Poles { group: usize, energy: f64 },
    #[error("spin group {group}: singular level matrix at {energy} eV")]
    Singular { group: usize, energy: f64 },
}

/// 断面積の計算に使う形に前処理した R-matrix limited のパラメータ。
#[derive(Debug, Clone)]
pub struct RMatrixLimitedXs {
    /// 非弾性散乱の準位の MT（粒子対の順）
    levels: Vec<u32>,
    groups: Vec<Group>,
}

/// 粒子対の役割。
#[derive(Debug, Clone, Copy, PartialEq)]
enum Role {
    Capture,
    Elastic,
    /// `levels` の添字
    Level(usize),
}

#[derive(Debug, Clone)]
struct Group {
    /// スピンの統計因子 g_J
    g: f64,
    /// 捕獲の幅を持つ共鳴があるか（なければ捕獲は厳密に 0）
    captures: bool,
    /// 中性子のチャネル（ファイルの順。捕獲のチャネルは除く）
    channels: Vec<Channel>,
    resonances: Vec<Level>,
}

#[derive(Debug, Clone)]
struct Channel {
    role: Role,
    l: OrbitalL,
    /// 閾値（実験室系, eV）
    threshold: f64,
    /// 波数 k = factor × sqrt(E − threshold)（10^-12 cm 単位）
    factor: f64,
    /// 透過因子に使う半径 (10^-12 cm)
    apt: f64,
    /// 位相シフトに使う半径 (10^-12 cm)
    ape: f64,
}

impl Channel {
    fn wave_number(&self, e: f64) -> f64 {
        self.factor * (e - self.threshold).sqrt()
    }
}

#[derive(Debug, Clone)]
struct Level {
    er: f64,
    /// Γγ / 2
    half_capture: f64,
    /// チャネルごとの換算幅の振幅 γ（`Group::channels` の順）
    amplitudes: Vec<f64>,
}

impl RMatrixLimitedXs {
    pub fn new(data: &RMatrixLimited) -> Result<Self, RmlError> {
        let mut levels = Vec::new();
        let roles = data
            .pairs
            .iter()
            .enumerate()
            .map(|(i, pair)| role(i + 1, pair, &mut levels))
            .collect::<Result<Vec<_>, _>>()?;
        let elastic: Vec<&ParticlePair> = (data.pairs.iter())
            .zip(&roles)
            .filter(|(_, r)| **r == Role::Elastic)
            .map(|(p, _)| p)
            .collect();
        let [elastic] = elastic[..] else {
            return Err(RmlError::ElasticPairs(elastic.len()));
        };
        // 実験室系と重心系のエネルギーの比 A/(A+1)
        let arat = elastic.mb / (elastic.mb + elastic.ma);
        let spins = (2.0 * elastic.ia.abs() + 1.0) * (2.0 * elastic.ib.abs() + 1.0);

        let groups = data
            .spin_groups
            .iter()
            .enumerate()
            .map(|(group, sg)| {
                let mut capture = None;
                let mut channels = Vec::new();
                for (k, c) in sg.channels.iter().enumerate() {
                    let pair = &data.pairs[c.pair];
                    match roles[c.pair] {
                        Role::Capture if capture.is_some() => {
                            return Err(RmlError::CaptureChannels { group: group + 1 });
                        }
                        Role::Capture => capture = Some(k),
                        role => {
                            let reduced = pair.ma * pair.mb / (pair.ma + pair.mb);
                            channels.push((
                                k,
                                Channel {
                                    role,
                                    l: c.l,
                                    threshold: -pair.q / arat,
                                    factor: wave_number_factor() * (reduced * arat).sqrt(),
                                    apt: c.apt,
                                    ape: c.ape,
                                },
                            ));
                        }
                    }
                }
                let resonances: Vec<Level> = sg
                    .resonances
                    .iter()
                    .map(|r| {
                        let amplitudes = channels
                            .iter()
                            .map(|(k, c)| amplitude(r.er, r.widths[*k], c, group + 1, *k + 1))
                            .collect::<Result<_, _>>()?;
                        Ok(Level {
                            er: r.er,
                            half_capture: capture.map_or(0.0, |k| 0.5 * r.widths[k].abs()),
                            amplitudes,
                        })
                    })
                    .collect::<Result<_, _>>()?;
                Ok(Group {
                    g: (2.0 * sg.aj.abs() + 1.0) / spins,
                    captures: resonances.iter().any(|r| r.half_capture > 0.0),
                    channels: channels.into_iter().map(|(_, c)| c).collect(),
                    resonances,
                })
            })
            .collect::<Result<_, _>>()?;
        Ok(Self { levels, groups })
    }

    /// `cross_sections` の値の反応: 全断面積（MT1）、弾性散乱（MT2）、捕獲（MT102）、準位ごとの非弾性散乱。
    pub fn reactions(&self) -> Vec<u32> {
        [1, 2, 102]
            .into_iter()
            .chain(self.levels.iter().copied())
            .collect()
    }

    /// 共鳴の半値幅（ピークのまわりの格子点に使う）。
    pub fn peaks(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.groups.iter().flat_map(|g| {
            g.resonances.iter().map(|r| {
                let neutron: f64 = (r.amplitudes.iter().zip(&g.channels))
                    .map(|(&a, c)| {
                        let e = (r.er - c.threshold).abs();
                        2.0 * a * a * penetrability(c.l, c.wave_number(c.threshold + e) * c.apt)
                    })
                    .sum();
                (r.er, r.half_capture + 0.5 * neutron)
            })
        })
    }

    /// エネルギー e (eV) での断面積 (barn)。`reactions` の順に並べる。
    pub fn cross_sections(&self, e: f64) -> Result<Vec<f64>, RmlError> {
        let mut sigma = vec![0.0; 3 + self.levels.len()];
        for (group, g) in self.groups.iter().enumerate() {
            let open: Vec<&Channel> = g.channels.iter().filter(|c| e > c.threshold).collect();
            let Some(entrance) = open.iter().find(|c| c.role == Role::Elastic) else {
                continue;
            };
            let n = open.len();
            let open_index: Vec<usize> = (g.channels.iter().enumerate())
                .filter(|(_, c)| e > c.threshold)
                .map(|(i, _)| i)
                .collect();
            let root_p: Vec<f64> = open
                .iter()
                .map(|c| penetrability(c.l, c.wave_number(e) * c.apt).sqrt())
                .collect();
            let phi: Vec<f64> = open
                .iter()
                .map(|c| phase_shift(c.l, c.wave_number(e) * c.ape))
                .collect();

            // K = P^½ R P^½、R = Σ γγ / (E_λ − E − iΓγ/2)。極の共鳴は w = P^½ γ として分ける
            let mut k = vec![Complex64::ZERO; n * n];
            let mut pole: Option<Vec<f64>> = None;
            for r in &g.resonances {
                let d = r.er - e;
                let denominator = d * d + r.half_capture * r.half_capture;
                if denominator == 0.0 {
                    let w: Vec<f64> = (0..n)
                        .map(|a| r.amplitudes[open_index[a]] * root_p[a])
                        .collect();
                    // 開いたチャネルに幅がない共鳴は寄与しない
                    if w.iter().all(|&v| v == 0.0) {
                        continue;
                    }
                    if pole.is_some() {
                        return Err(RmlError::Poles {
                            group: group + 1,
                            energy: e,
                        });
                    }
                    pole = Some(w);
                    continue;
                }
                let shape = Complex64::new(d / denominator, r.half_capture / denominator);
                for a in 0..n {
                    for b in 0..n {
                        let gamma = r.amplitudes[open_index[a]] * r.amplitudes[open_index[b]];
                        k[a * n + b] += shape * (gamma * root_p[a] * root_p[b]);
                    }
                }
            }
            // (I − iK) X = K
            let mut m: Vec<Complex64> = k.iter().map(|&v| -Complex64::I * v).collect();
            for a in 0..n {
                m[a * n + a] += 1.0;
            }
            let singular = || RmlError::Singular {
                group: group + 1,
                energy: e,
            };
            let x = match pole {
                None => solve(m, k, n).ok_or_else(singular)?,
                Some(w) => pole_limit(m, &w, n).ok_or_else(singular)?,
            };

            let kn = entrance.wave_number(e);
            let pifac = g.g * PI / (kn * kn);
            // この spin group の全断面積・弾性散乱・非弾性散乱
            let mut part = vec![0.0; sigma.len()];
            for (a, ca) in open.iter().enumerate() {
                if ca.role != Role::Elastic {
                    continue;
                }
                let (s2, c2) = (2.0 * phi[a]).sin_cos();
                let xaa = x[a * n + a];
                // 1 − Re U_aa と Im U_aa（1 − cos 2φ を 2 sin²φ にして桁落ちを避ける）
                let one_minus_re =
                    2.0 * phi[a].sin().powi(2) + 2.0 * xaa.im * c2 - 2.0 * xaa.re * s2;
                let im = -s2 * (1.0 - 2.0 * xaa.im) + 2.0 * xaa.re * c2;
                part[0] += 2.0 * pifac * one_minus_re;
                part[1] += pifac * (one_minus_re * one_minus_re + im * im);
                for (b, cb) in open.iter().enumerate() {
                    // 位相の因子は絶対値 1 なので、a ≠ b では |U_ab|² = 4|X_ab|²
                    let transfer = 4.0 * pifac * x[a * n + b].norm_sqr();
                    match cb.role {
                        Role::Elastic if b != a => part[1] += transfer,
                        Role::Level(i) => part[3 + i] += transfer,
                        _ => {}
                    }
                }
            }
            // 捕獲は残りの吸収。捕獲の幅がなければ、差を取らずに厳密に 0 とする（丸めの誤差を残さない）
            if g.captures {
                part[2] = part[0] - part[1] - part[3..].iter().sum::<f64>();
            }
            for (s, p) in sigma.iter_mut().zip(part) {
                *s += p;
            }
        }
        Ok(sigma)
    }
}

/// 粒子対の役割を決める（扱えない粒子対はエラー）。非弾性散乱の準位は `levels` に MT を加える。
fn role(index: usize, pair: &ParticlePair, levels: &mut Vec<u32>) -> Result<Role, RmlError> {
    let unsupported = |what: String| Err(RmlError::Pair { pair: index, what });
    if pair.za != 0.0 {
        return unsupported(format!("a charged particle (ZA={})", pair.za));
    }
    if pair.shf != 0 {
        return unsupported(format!("the shift factor (SHF={})", pair.shf));
    }
    match (pair.ma, pair.mt) {
        // 光子と捕獲。透過因子は計算しない（PNT=0 は MT=102 なら計算しないこと）
        (0.0, 102) if pair.pnt <= 0 => Ok(Role::Capture),
        (0.0, 102) => unsupported(format!("penetrability for capture (PNT={})", pair.pnt)),
        (1.0, mt @ (2 | 51..=90)) if pair.pnt >= 0 => {
            if mt == 2 {
                if pair.q != 0.0 {
                    return unsupported(format!("elastic scattering with Q={}", pair.q));
                }
                Ok(Role::Elastic)
            } else if levels.contains(&mt) {
                Err(RmlError::DuplicateMt(mt))
            } else {
                levels.push(mt);
                Ok(Role::Level(levels.len() - 1))
            }
        }
        (1.0, 2 | 51..=90) => {
            unsupported(format!("neutrons without penetrability (PNT={})", pair.pnt))
        }
        (ma, mt) => unsupported(format!("MA={ma} MT={mt}")),
    }
}

/// 共鳴エネルギーでの幅 Γ から換算幅の振幅 γ = sign(Γ) √(|Γ| / 2P) を求める（ENDF の約束）。
/// 透過因子はチャネルのエネルギー |E_λ − 閾値| で求める（負の共鳴は |E_λ|）。
fn amplitude(
    er: f64,
    width: f64,
    c: &Channel,
    group: usize,
    channel: usize,
) -> Result<f64, RmlError> {
    if width == 0.0 {
        return Ok(0.0);
    }
    let closed = || RmlError::ClosedWidth { group, er, channel };
    // 非弾性散乱のチャネルは、閾値より下では閉じていて幅を定義できない
    if c.role != Role::Elastic && er <= c.threshold {
        return Err(closed());
    }
    let e = (er - c.threshold).abs();
    let p = penetrability(c.l, c.wave_number(c.threshold + e) * c.apt);
    if p <= 0.0 {
        return Err(closed());
    }
    Ok(width.signum() * (width.abs() / (2.0 * p)).sqrt())
}

/// 極での X = i(I − M⁻¹) の極限。`m` は極の共鳴を除いた M' = I − iK'、`w` は極の共鳴の P^½ γ。
fn pole_limit(m: Vec<Complex64>, w: &[f64], n: usize) -> Option<Vec<Complex64>> {
    let identity: Vec<Complex64> = (0..n * n)
        .map(|i| {
            if i % (n + 1) == 0 {
                Complex64::ONE
            } else {
                Complex64::ZERO
            }
        })
        .collect();
    let mut inverse = solve(m, identity, n)?;
    // u = M'⁻¹ w（M' は複素対称なので wᵀ M'⁻¹ = uᵀ）
    let u: Vec<Complex64> = (0..n)
        .map(|a| (0..n).fold(Complex64::ZERO, |acc, b| acc + inverse[a * n + b] * w[b]))
        .collect();
    let s = (0..n).fold(Complex64::ZERO, |acc, a| acc + u[a] * w[a]);
    // 開いたチャネルに幅がなければ、極の共鳴は寄与しない
    if s.norm_sqr() > 0.0 {
        let scale = s.inv();
        for a in 0..n {
            for b in 0..n {
                inverse[a * n + b] -= u[a] * u[b] * scale;
            }
        }
    }
    Some(
        (0..n * n)
            .map(|i| {
                // i(δ − M⁻¹)
                let identity = if i % (n + 1) == 0 {
                    Complex64::ONE
                } else {
                    Complex64::ZERO
                };
                Complex64::I * (identity - inverse[i])
            })
            .collect(),
    )
}

/// 複素 n×n 行列の方程式 M X = B を部分ピボットつき Gauss の消去法で解く（行列は行優先）。
/// 特異なら `None`。
fn solve(mut m: Vec<Complex64>, mut b: Vec<Complex64>, n: usize) -> Option<Vec<Complex64>> {
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            m[i * n + col]
                .norm_sqr()
                .total_cmp(&m[j * n + col].norm_sqr())
        })?;
        if m[pivot * n + col].norm_sqr() == 0.0 {
            return None;
        }
        for j in 0..n {
            m.swap(col * n + j, pivot * n + j);
            b.swap(col * n + j, pivot * n + j);
        }
        let p = m[col * n + col].inv();
        for j in 0..n {
            m[col * n + j] *= p;
            b[col * n + j] *= p;
        }
        for row in 0..n {
            if row != col {
                let f = m[row * n + col];
                for j in 0..n {
                    m[row * n + j] = m[row * n + j] - f * m[col * n + j];
                    b[row * n + j] = b[row * n + j] - f * b[col * n + j];
                }
            }
        }
    }
    Some(b)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::endf::{
        OrbitalL, ParticlePair, RMatrixLimited, RmlChannel, RmlResonance, SpinGroup,
    };
    use crate::resonance::reich_moore_tests::{AP, AWRI, single_level_exact};

    fn pair(ma: f64, mb: f64, q: f64, mt: u32) -> ParticlePair {
        ParticlePair {
            ma,
            mb,
            za: 0.0,
            zb: 92.0,
            ia: if ma == 0.0 { 1.0 } else { 0.5 },
            ib: 0.0,
            q,
            pnt: 0,
            shf: 0,
            mt,
        }
    }

    fn channel(pair: usize, l: OrbitalL, sch: f64) -> RmlChannel {
        RmlChannel {
            pair,
            l,
            sch,
            bnd: 0.0,
            ape: AP,
            apt: AP,
        }
    }

    /// 捕獲と弾性散乱の粒子対に、`levels` の Q 値の非弾性散乱の粒子対（MT51 から）を加える。
    fn pairs(levels: &[f64]) -> Vec<ParticlePair> {
        let mut pairs = vec![pair(0.0, AWRI + 1.0, 0.0, 102), pair(1.0, AWRI, 0.0, 2)];
        pairs.extend(
            levels
                .iter()
                .zip(51..)
                .map(|(&q, mt)| pair(1.0, AWRI, q, mt)),
        );
        pairs
    }

    /// Reich-Moore の `single_level` と同じ、l=0 の単一準位。
    fn single_level(gn: f64, gg: f64) -> RMatrixLimited {
        RMatrixLimited {
            pairs: pairs(&[]),
            spin_groups: vec![SpinGroup {
                aj: 0.5,
                pj: 0.0,
                channels: vec![channel(0, OrbitalL::S, 0.0), channel(1, OrbitalL::S, 0.5)],
                resonances: vec![RmlResonance {
                    er: 6.674,
                    widths: vec![gg, gn],
                }],
            }],
        }
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
            "{what}: actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn single_level_matches_breit_wigner() {
        let xs = RMatrixLimitedXs::new(&single_level(1.5e-3, 2.3e-2)).unwrap();
        assert_eq!(xs.reactions(), vec![1, 2, 102]);
        for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
            let expected = single_level_exact(e, 1.5e-3, 2.3e-2, 0.0);
            let [total, elastic, capture] = xs.cross_sections(e).unwrap()[..] else {
                panic!("expected total, elastic and capture");
            };
            assert_close(total, expected.total, &format!("total at {e}"));
            assert_close(elastic, expected.elastic, &format!("elastic at {e}"));
            assert_close(capture, expected.capture, &format!("capture at {e}"));
        }
    }

    /// 弾性散乱と、Q = -1 keV の準位への非弾性散乱のチャネルを持つ共鳴（捕獲なし）。
    pub(crate) fn inelastic_level(gg: f64) -> RMatrixLimited {
        RMatrixLimited {
            pairs: pairs(&[-1.0e3]),
            spin_groups: vec![SpinGroup {
                aj: 0.5,
                pj: 0.0,
                channels: vec![
                    channel(0, OrbitalL::S, 0.0),
                    channel(1, OrbitalL::S, 0.5),
                    channel(2, OrbitalL::S, 0.5),
                ],
                resonances: vec![RmlResonance {
                    er: 2.0e3,
                    widths: vec![gg, 3.0, -2.0],
                }],
            }],
        }
    }

    #[test]
    fn inelastic_channel_opens_at_threshold() {
        let xs = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
        assert_eq!(xs.reactions(), vec![1, 2, 102, 51]);
        // 閾値（実験室系）は |Q| (A+1)/A
        let threshold = 1.0e3 * (AWRI + 1.0) / AWRI;
        let below = xs.cross_sections(threshold * 0.999).unwrap();
        assert_eq!(below[3], 0.0);
        for e in [threshold * 1.001, 1.5e3, 1.999e3, 2.0e3, 2.001e3, 3.0e3] {
            let [total, elastic, capture, level] = xs.cross_sections(e).unwrap()[..] else {
                panic!("expected four reactions");
            };
            assert!(level > 0.0, "level at {e}");
            // 捕獲がなければ、全断面積は弾性散乱と非弾性散乱の和（ユニタリ性）
            assert!(capture.abs() <= 1e-9 * total, "capture {capture} at {e}");
            assert_close(elastic + level, total, &format!("unitarity at {e}"));
        }
    }

    #[test]
    fn capture_takes_the_rest_of_the_absorption() {
        let xs = RMatrixLimitedXs::new(&inelastic_level(0.5)).unwrap();
        let [total, elastic, capture, level] = xs.cross_sections(2.0e3).unwrap()[..] else {
            panic!("expected four reactions");
        };
        assert!(capture > 0.0);
        assert_close(elastic + capture + level, total, "sum of parts");
    }

    #[test]
    fn pole_of_a_resonance_without_capture_is_the_limit() {
        // 捕獲がないと共鳴エネルギーは R 行列の極になるが、断面積は連続。両側の値の平均（傾きの分が
        // 打ち消し合う）と比べる
        let xs = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
        let at = xs.cross_sections(2.0e3).unwrap();
        let below = xs.cross_sections(2.0e3 * (1.0 - 1.0e-9)).unwrap();
        let above = xs.cross_sections(2.0e3 * (1.0 + 1.0e-9)).unwrap();
        for ((a, b), c) in at.iter().zip(&below).zip(&above) {
            assert_close(*a, 0.5 * (b + c), "pole");
        }
    }

    #[test]
    fn resonance_without_widths_does_not_count_as_a_pole() {
        // 幅がすべて 0 の共鳴は何も寄与しないので、同じエネルギーに極があっても問題ない
        let mut data = inelastic_level(0.0);
        data.spin_groups[0].resonances.push(RmlResonance {
            er: 2.0e3,
            widths: vec![0.0, 0.0, 0.0],
        });
        let xs = RMatrixLimitedXs::new(&data).unwrap();
        let single = RMatrixLimitedXs::new(&inelastic_level(0.0)).unwrap();
        assert_eq!(
            xs.cross_sections(2.0e3).unwrap(),
            single.cross_sections(2.0e3).unwrap()
        );
    }

    #[test]
    fn rejects_width_of_a_channel_closed_at_the_resonance() {
        let mut data = inelastic_level(0.5);
        // 共鳴が閾値より下にあると、非弾性散乱の幅は定義できない
        data.spin_groups[0].resonances[0].er = 500.0;
        let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
        assert!(error.contains("closed"), "{error}");
    }

    #[test]
    fn rejects_shift_factor() {
        let mut data = single_level(1.5e-3, 2.3e-2);
        data.pairs[1].shf = 1;
        let error = RMatrixLimitedXs::new(&data).unwrap_err().to_string();
        assert!(error.contains("SHF=1"), "{error}");
    }

    #[test]
    fn rejects_charged_pairs() {
        let mut data = single_level(1.5e-3, 2.3e-2);
        data.pairs[1].za = 1.0;
        assert!(RMatrixLimitedXs::new(&data).is_err());
    }
}
