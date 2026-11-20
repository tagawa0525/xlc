//! R-matrix limited 形式（LRF=7）の Reich-Moore の断面積（0 K）。
//!
//! 捕獲は消去されたチャネルとして共鳴エネルギーの虚部に入れ（E_λ − iΓγ/2）、開いた中性子と核分裂のチャネルの
//! R 行列から X = P^½ (I − R L)⁻¹ R P^½（ENDF-102 の式 D.42、D.47。L = (S − B) + iP）を求める。K = P^½ R P^½ と
//! Λ = L/P（対角）で書くと X = (I − KΛ)⁻¹K で、シフト因子がなければ Λ = i なので X = (I − iK)⁻¹K。
//! 衝突行列は U = e^{−i(φ+φ')} (I + 2iX)。R には、チャネルのバックグラウンドの R 行列（SAMMY の対数の形、
//! ENDF-102 の式 D.77）を対角に足す。
//! 中性子の透過因子とシフト因子は剛体球のもの（シフト因子は SHF=1 の粒子対だけ）。荷電粒子の放出のチャネル（MT600〜849）の透過因子は
//! Coulomb の波動関数のもの（ENDF-102 の式 D.84。`coulomb`）で、出口のチャネルにだけなるので位相シフトを使わない。
//! 核分裂のチャネルは透過因子を 1 とし（ENDF-102 の D.1.7）、
//! 幅は Γ = 2γ² になる（LRF=3 の Reich-Moore の Γfa、Γfb と同じチャネル）。位相シフトは 0 とする（核分裂は入射の
//! チャネルにならず、断面積は |U_ab|² = 4|X_ab|² だけで位相によらない）。
//! 閾値より下の非弾性散乱のチャネルは閉じているので除く（NJOY2016 samm と同じ）。
//!
//! 捕獲の幅がない共鳴のエネルギーちょうどでは R 行列が発散する（極）が、X は有限の極限を持つ。
//! M = I − KΛ とすると X = (M⁻¹ − I)Λ⁻¹ で、極の共鳴の項を K = K' + w wᵀ/d と分けると M = M' − w (Λw)ᵀ/d で、
//! d → 0 の極限は M⁻¹ = M'⁻¹ − (M'⁻¹w)((Λw)ᵀM'⁻¹) / ((Λw)ᵀM'⁻¹w)（Sherman-Morrison の公式）。Λ = i では
//! X = i(I − M⁻¹) になる。NJOY はこの点で発散する。

use std::f64::consts::PI;

use num_complex::Complex64;

use super::coulomb;
use super::factors::{penetrability, phase_shift, shift};
use super::linear::{Singular, pole_limit, solve};
use crate::endf::{LogBackground, OrbitalL, ParticlePair, RMatrixLimited};
use crate::physics::{AMASSN, AMU, CLIGHT, EV, HBAR, wave_number_factor};

/// NJOY2016 samm の coulx が η の大きいときの近似の式に切り替える η（`eta.ge.five`）。
const NJOY_LARGE_ETA: f64 = 5.0;

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
    #[error(
        "spin group {group}: the resonance at {er} eV has a width in channel {channel}, \
         whose penetrability underflows to 0 at that energy"
    )]
    Underflow {
        group: usize,
        er: f64,
        channel: usize,
    },
    #[error(transparent)]
    Coulomb(#[from] coulomb::NotConverged),
    #[error("spin group {group}: several resonances without capture width at {energy} eV")]
    Poles { group: usize, energy: f64 },
    #[error("spin group {group}: singular level matrix at {energy} eV")]
    Singular { group: usize, energy: f64 },
    #[error("spin group {group}: a background R-matrix on the capture channel")]
    CaptureBackground { group: usize },
    #[error(
        "spin group {group}: {energy} eV is outside the background R-matrix range \
         ({ed}, {eu}) eV of channel {channel}"
    )]
    BackgroundRange {
        group: usize,
        channel: usize,
        energy: f64,
        ed: f64,
        eu: f64,
    },
}

/// 断面積の計算に使う形に前処理した R-matrix limited のパラメータ。
#[derive(Debug, Clone)]
pub struct RMatrixLimitedXs {
    /// 出口のチャネルの反応の MT（非弾性散乱の準位と荷電粒子の放出。粒子対の順）
    exits: Vec<u32>,
    /// 核分裂の粒子対があるか（あれば `cross_sections` の最後の値が核分裂）
    fission: bool,
    groups: Vec<Group>,
}

/// 粒子対の役割。
#[derive(Debug, Clone, Copy, PartialEq)]
enum Role {
    Capture,
    Elastic,
    /// `exits` の添字
    Exit(usize),
    Fission,
}

#[derive(Debug, Clone)]
struct Group {
    /// スピンの統計因子 g_J
    g: f64,
    /// 捕獲の幅を持つ共鳴があるか（なければ捕獲は厳密に 0）
    captures: bool,
    /// 中性子と核分裂のチャネル（ファイルの順。捕獲のチャネルは除く）
    channels: Vec<Channel>,
    resonances: Vec<Level>,
}

#[derive(Debug, Clone)]
struct Channel {
    role: Role,
    l: OrbitalL,
    /// 透過因子と位相シフトを計算するか（しなければ P = 1、φ = 0）
    penetrable: bool,
    /// Coulomb のパラメータと波数の積 ηk（10^12 cm⁻¹。荷電粒子の対でなければ 0 で、透過因子は剛体球のもの）
    eta_k: f64,
    /// 閾値（実験室系, eV）
    threshold: f64,
    /// 波数 k = factor × sqrt(E − threshold)（10^-12 cm 単位）
    factor: f64,
    /// 透過因子に使う半径 (10^-12 cm)
    apt: f64,
    /// 位相シフトに使う半径 (10^-12 cm)
    ape: f64,
    /// シフト因子を計算するときの境界条件 B（SHF=1 の粒子対だけ。なければ S − B = 0）
    boundary: Option<f64>,
    background: Option<LogBackground>,
    /// spin group の中のチャネルの番号（ファイルの順、1 から。エラーの文に使う）
    index: usize,
}

impl Channel {
    fn wave_number(&self, e: f64) -> f64 {
        self.factor * (e - self.threshold).sqrt()
    }

    fn penetrability(&self, e: f64) -> Result<f64, RmlError> {
        if !self.penetrable {
            return Ok(1.0);
        }
        let k = self.wave_number(e);
        if self.eta_k > 0.0 {
            Ok(coulomb::penetrability(
                self.l.number(),
                self.eta_k / k,
                k * self.apt,
            )?)
        } else {
            Ok(penetrability(self.l, k * self.apt))
        }
    }

    /// 剛体球の位相シフト。入射のチャネル（中性子の弾性散乱）にだけ使う。
    fn phase_shift(&self, e: f64) -> f64 {
        if self.penetrable {
            phase_shift(self.l, self.wave_number(e) * self.ape)
        } else {
            0.0
        }
    }

    /// S − B（シフト因子を計算しないチャネルは `None`）。
    fn shift_minus_boundary(&self, e: f64) -> Option<f64> {
        self.boundary
            .map(|b| shift(self.l, self.wave_number(e) * self.apt) - b)
    }

    /// バックグラウンドの R 行列（ENDF-102 の式 D.77）。範囲 (ED, EU) の外では対数が定義できないのでエラー。
    fn background(&self, e: f64, group: usize) -> Result<f64, RmlError> {
        let Some(b) = &self.background else {
            return Ok(0.0);
        };
        if !(b.ed < e && e < b.eu) {
            return Err(RmlError::BackgroundRange {
                group,
                channel: self.index,
                energy: e,
                ed: b.ed,
                eu: b.eu,
            });
        }
        Ok(b.r0 + (b.r1 + b.r2 * e) * e
            - b.s1 * (b.eu - b.ed)
            - (b.s0 + b.s1 * e) * ((b.eu - e) / (e - b.ed)).ln())
    }
}

#[derive(Debug, Clone)]
struct Level {
    er: f64,
    /// Γγ / 2
    half_capture: f64,
    /// ピークの半値幅 Γγ/2 + Σ γ²P（P はチャネルのエネルギー |E_λ − 閾値| で求める。`peaks` が使う）
    half_width: f64,
    /// チャネルごとの換算幅の振幅 γ（`Group::channels` の順）
    amplitudes: Vec<f64>,
}

impl RMatrixLimitedXs {
    pub fn new(data: &RMatrixLimited) -> Result<Self, RmlError> {
        let mut exits = Vec::new();
        let roles = data
            .pairs
            .iter()
            .enumerate()
            .map(|(i, pair)| role(i + 1, pair, &mut exits))
            .collect::<Result<Vec<_>, _>>()?;
        let fission = match roles.iter().filter(|r| **r == Role::Fission).count() {
            0 => false,
            1 => true,
            _ => return Err(RmlError::DuplicateMt(18)),
        };
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
                        Role::Capture if c.background.is_some() => {
                            return Err(RmlError::CaptureBackground { group: group + 1 });
                        }
                        Role::Capture => capture = Some(k),
                        role => {
                            let reduced = pair.ma * pair.mb / (pair.ma + pair.mb);
                            channels.push((
                                k,
                                Channel {
                                    role,
                                    l: c.l,
                                    penetrable: role != Role::Fission,
                                    eta_k: coulomb_eta_k(pair, reduced),
                                    threshold: -pair.q / arat,
                                    factor: wave_number_factor() * (reduced * arat).sqrt(),
                                    apt: c.apt,
                                    ape: c.ape,
                                    boundary: (pair.shf == 1).then_some(c.bnd),
                                    background: c.background,
                                    index: k + 1,
                                },
                            ));
                        }
                    }
                }
                let resonances: Vec<Level> = sg
                    .resonances
                    .iter()
                    .map(|r| {
                        let amplitudes: Vec<f64> = channels
                            .iter()
                            .map(|(k, c)| amplitude(r.er, r.widths[*k], c, group + 1, *k + 1))
                            .collect::<Result<_, _>>()?;
                        let half_capture = capture.map_or(0.0, |k| 0.5 * r.widths[k].abs());
                        let mut neutron = 0.0;
                        for (&a, (_, c)) in amplitudes.iter().zip(&channels) {
                            let e = c.threshold + (r.er - c.threshold).abs();
                            neutron += 2.0 * a * a * c.penetrability(e)?;
                        }
                        Ok::<_, RmlError>(Level {
                            er: r.er,
                            half_capture,
                            half_width: half_capture + 0.5 * neutron,
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
        Ok(Self {
            exits,
            fission,
            groups,
        })
    }

    /// `cross_sections` の値の反応: 全断面積（MT1）、弾性散乱（MT2）、捕獲（MT102）、出口のチャネル（準位ごとの非弾性散乱と
    /// 荷電粒子の放出）、核分裂（MT18。核分裂の粒子対があるときだけ）。
    pub fn reactions(&self) -> Vec<u32> {
        [1, 2, 102]
            .into_iter()
            .chain(self.exits.iter().copied())
            .chain(self.fission.then_some(18))
            .collect()
    }

    /// 荷電粒子のチャネルの幅のうち、共鳴エネルギーで Coulomb の障壁の下（ρ が転回点 η + √(η² + l(l+1)) より小さい）にあり、
    /// η が `NJOY_LARGE_ETA` 以上のものの反応の MT（幅ごとに 1 つ）。NJOY2016 RECONR（samm の coulx と coulfg）は、そこで透過因子を
    /// 近似の式（JWKB、η の大きいときの漸近級数）で求め、1e-35 以下の透過因子を 0 にする（共鳴の幅の換算では 1 にする）ので、
    /// その反応は NJOY と違う（ADR-0273）。
    pub fn njoy_approximated_widths(&self) -> Vec<u32> {
        let mut mts = Vec::new();
        for g in &self.groups {
            for r in &g.resonances {
                for (&a, c) in r.amplitudes.iter().zip(&g.channels) {
                    let Role::Exit(i) = c.role else { continue };
                    if a == 0.0 || c.eta_k == 0.0 {
                        continue;
                    }
                    let k = c.wave_number(r.er);
                    let (eta, rho) = (c.eta_k / k, k * c.apt);
                    let l = f64::from(c.l.number());
                    if eta >= NJOY_LARGE_ETA && rho < eta + (eta * eta + l * (l + 1.0)).sqrt() {
                        mts.push(self.exits[i]);
                    }
                }
            }
        }
        mts
    }

    /// 共鳴のピークのエネルギーと半値幅（ピークのまわりの格子点に使う）。シフト因子を計算するチャネルがあれば、
    /// ピークは準位のずれ Δ = −Σ γ²(S − B)（S は共鳴エネルギーで求める）だけ E_λ から動く。
    pub fn peaks(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.groups.iter().flat_map(|g| {
            g.resonances.iter().map(|r| {
                let mut level_shift = 0.0;
                for (&a, c) in r.amplitudes.iter().zip(&g.channels) {
                    let e = c.threshold + (r.er - c.threshold).abs();
                    if let Some(d) = c.shift_minus_boundary(e) {
                        level_shift -= a * a * d;
                    }
                }
                (r.er + level_shift, r.half_width)
            })
        })
    }

    /// エネルギー e (eV) での断面積 (barn)。`reactions` の順に並べる。
    pub fn cross_sections(&self, e: f64) -> Result<Vec<f64>, RmlError> {
        let fission_index = 3 + self.exits.len();
        let mut sigma = vec![0.0; fission_index + usize::from(self.fission)];
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
                .map(|c| c.penetrability(e).map(f64::sqrt))
                .collect::<Result<_, _>>()?;

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
            // バックグラウンドの R 行列は対角だけ: K_aa += P_a R_bkg
            for (a, c) in open.iter().enumerate() {
                if c.background.is_some() {
                    k[a * n + a] += root_p[a] * root_p[a] * c.background(e, group + 1)?;
                }
            }
            // Λ = i + (S − B)/P。シフト因子を計算するチャネルがなければ Λ = i のまま (I − iK) X = K を解く
            let lambda: Option<Vec<Complex64>> = if open.iter().any(|c| c.boundary.is_some()) {
                Some(
                    open.iter()
                        .zip(&root_p)
                        .map(|(c, rp)| {
                            let ratio = c.shift_minus_boundary(e).map_or(0.0, |d| d / (rp * rp));
                            Complex64::new(ratio, 1.0)
                        })
                        .collect(),
                )
            } else {
                None
            };
            // (I − KΛ) X = K
            let mut m: Vec<Complex64> = match &lambda {
                None => k.iter().map(|&v| -Complex64::I * v).collect(),
                Some(lambda) => (0..n * n).map(|i| -k[i] * lambda[i % n]).collect(),
            };
            for a in 0..n {
                m[a * n + a] += 1.0;
            }
            let singular = || RmlError::Singular {
                group: group + 1,
                energy: e,
            };
            let x = match pole {
                None => {
                    solve(&mut m, &mut k, n).map_err(|Singular| singular())?;
                    k
                }
                Some(w) => match &lambda {
                    None => x_at_pole(m, &w, n),
                    Some(lambda) => x_at_shifted_pole(m, &w, lambda, n),
                }
                .map_err(|Singular| singular())?,
            };

            let kn = entrance.wave_number(e);
            let pifac = g.g * PI / (kn * kn);
            // この spin group の全断面積・弾性散乱・非弾性散乱・核分裂
            let mut part = vec![0.0; sigma.len()];
            for (a, ca) in open.iter().enumerate() {
                if ca.role != Role::Elastic {
                    continue;
                }
                let phi = ca.phase_shift(e);
                let (s2, c2) = (2.0 * phi).sin_cos();
                let xaa = x[a * n + a];
                // 1 − Re U_aa と Im U_aa（1 − cos 2φ を 2 sin²φ にして桁落ちを避ける）
                let one_minus_re = 2.0 * phi.sin().powi(2) + 2.0 * xaa.im * c2 - 2.0 * xaa.re * s2;
                let im = -s2 * (1.0 - 2.0 * xaa.im) + 2.0 * xaa.re * c2;
                part[0] += 2.0 * pifac * one_minus_re;
                part[1] += pifac * (one_minus_re * one_minus_re + im * im);
                for (b, cb) in open.iter().enumerate() {
                    // 位相の因子は絶対値 1 なので、a ≠ b では |U_ab|² = 4|X_ab|²
                    let transfer = 4.0 * pifac * x[a * n + b].norm_sqr();
                    match cb.role {
                        Role::Elastic if b != a => part[1] += transfer,
                        Role::Exit(i) => part[3 + i] += transfer,
                        Role::Fission => part[fission_index] += transfer,
                        Role::Elastic | Role::Capture => {}
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

/// 粒子対の役割を決める（扱えない粒子対はエラー）。出口のチャネル（非弾性散乱の準位と荷電粒子の放出）は `exits` に MT を加える。
fn role(index: usize, pair: &ParticlePair, exits: &mut Vec<u32>) -> Result<Role, RmlError> {
    let unsupported = |what: String| Err(RmlError::Pair { pair: index, what });
    // SHF=1 は計算し、0 と -1 は計算しない（-1 は ENDF-102 の 2005 年版の「計算しない」。2023 年版は 0 と書く）
    if !matches!(pair.shf, -1..=1) {
        return unsupported(format!("the shift factor flag SHF={}", pair.shf));
    }
    let role = role_of(index, pair, exits)?;
    // シフト因子は透過因子と同じ剛体球の量なので、透過因子を計算しない捕獲と核分裂には定義がない
    if pair.shf == 1 && matches!(role, Role::Capture | Role::Fission) {
        return unsupported(format!("the shift factor (SHF=1) for MT={}", pair.mt));
    }
    // 荷電粒子のシフト因子は Coulomb の波動関数のもので、剛体球の式は使えない（4 つの版の評価に SHF=1 の荷電粒子の対は無い。ADR-0273）
    if pair.shf == 1 && pair.za != 0.0 {
        return unsupported(format!(
            "the shift factor (SHF=1) for a charged particle (ZA={})",
            pair.za
        ));
    }
    Ok(role)
}

fn role_of(index: usize, pair: &ParticlePair, exits: &mut Vec<u32>) -> Result<Role, RmlError> {
    let unsupported = |what: String| Err(RmlError::Pair { pair: index, what });
    let mut exit = |mt: u32| {
        if exits.contains(&mt) {
            Err(RmlError::DuplicateMt(mt))
        } else {
            exits.push(mt);
            Ok(Role::Exit(exits.len() - 1))
        }
    };
    let neutron = pair.za == 0.0;
    match (pair.ma, pair.mt) {
        // 光子と捕獲。透過因子は計算しない（PNT=0 は MT=102 なら計算しないこと）
        (0.0, 102) if pair.pnt <= 0 => Ok(Role::Capture),
        (0.0, 102) => unsupported(format!("penetrability for capture (PNT={})", pair.pnt)),
        (1.0, 2) if neutron && pair.pnt >= 0 => {
            if pair.q != 0.0 {
                return unsupported(format!("elastic scattering with Q={}", pair.q));
            }
            Ok(Role::Elastic)
        }
        (1.0, mt @ 51..=90) if neutron && pair.pnt >= 0 => exit(mt),
        (1.0, 2 | 51..=90) if neutron => {
            unsupported(format!("neutrons without penetrability (PNT={})", pair.pnt))
        }
        // 荷電粒子の放出（p、d、t、³He、α の準位の MT600〜849。ENDF-102 の 3.4.6）。質量と電荷から Coulomb の透過因子を求める
        (ma, mt @ 600..=849) if pair.za > 0.0 && ma > 0.0 && pair.pnt >= 0 => exit(mt),
        (_, 600..=849) if pair.za > 0.0 && pair.pnt < 0 => unsupported(format!(
            "a charged particle without penetrability (PNT={})",
            pair.pnt
        )),
        // 核分裂の片割れの粒子は決まらないので MA を見ない。核分裂の透過因子は 1（ENDF-102 の D.1.7）なので、PNT=-1 と、
        // 既定の PNT=0 を受ける。PNT=0 の既定を 2.2.1.6 は「fission (MT=19)」を -1、ほかを +1 と書くが、D.1.7 と NJOY2016 samm の
        // ppdefs（MT=18 の lpent=0）に合わせる。PNT=1（剛体球の透過因子）は核分裂に定義がないのでエラーにする
        (_, 18) if pair.pnt <= 0 => Ok(Role::Fission),
        (_, 18) => unsupported(format!("penetrability for fission (PNT={})", pair.pnt)),
        (ma, mt) => unsupported(format!("MA={ma} ZA={} MT={mt}", pair.za)),
    }
}

/// Coulomb のパラメータ η = z Z e² μ / (ħ² k)（ENDF-102 の式 D.80）に k を掛けた値 z Z α μ c / ħ（10^12 cm⁻¹）。`reduced` は
/// 換算質量（中性子の質量の単位）。片方の電荷が 0 なら 0。
fn coulomb_eta_k(pair: &ParticlePair, reduced: f64) -> f64 {
    // 微細構造定数 α = e²/(ħc)。e² を EV（erg/eV）で書くと EV² c / (10^16 ħ)（NJOY2016 phys の finstri の逆数）
    let alpha = EV * EV * CLIGHT / (1.0e16 * HBAR);
    pair.za * pair.zb * alpha * reduced * AMASSN * AMU * CLIGHT / HBAR * 1.0e-12
}

/// 共鳴エネルギーでの幅 Γ から換算幅の振幅 γ = sign(Γ) √(|Γ| / 2P) を求める（ENDF の約束）。
/// 透過因子はチャネルのエネルギー |E_λ − 閾値| で求める（負の共鳴は |E_λ|）。透過因子を計算しないチャネルは
/// P = 1 で、閉じることがない。
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
    if !c.penetrable {
        return Ok(width.signum() * (width.abs() / 2.0).sqrt());
    }
    let closed = || RmlError::ClosedWidth { group, er, channel };
    // 非弾性散乱のチャネルは、閾値より下では閉じていて幅を定義できない
    if c.role != Role::Elastic && er <= c.threshold {
        return Err(closed());
    }
    let e = (er - c.threshold).abs();
    let p = c.penetrability(c.threshold + e)?;
    // 荷電粒子のチャネルは、障壁の深い下で透過因子が倍精度の 0 になる
    if p <= 0.0 && e > 0.0 && c.eta_k > 0.0 {
        return Err(RmlError::Underflow { group, er, channel });
    }
    if p <= 0.0 {
        return Err(closed());
    }
    Ok(width.signum() * (width.abs() / (2.0 * p)).sqrt())
}

/// 極での X = i(I − M⁻¹) の極限。`m` は極の共鳴を除いた M' = I − iK'、`w` は極の共鳴の P^½ γ。
fn x_at_pole(mut m: Vec<Complex64>, w: &[f64], n: usize) -> Result<Vec<Complex64>, Singular> {
    let mut inverse: Vec<Complex64> = (0..n * n)
        .map(|i| {
            if i % (n + 1) == 0 {
                Complex64::ONE
            } else {
                Complex64::ZERO
            }
        })
        .collect();
    solve(&mut m, &mut inverse, n)?;
    pole_limit(&mut inverse, w, n)?;
    Ok((0..n * n)
        .map(|i| {
            // i(δ − M⁻¹)
            let identity = if i % (n + 1) == 0 {
                Complex64::ONE
            } else {
                Complex64::ZERO
            };
            Complex64::I * (identity - inverse[i])
        })
        .collect())
}

/// シフト因子があるときの極での X = (M⁻¹ − I)Λ⁻¹ の極限。`m` は極の共鳴を除いた M' = I − K'Λ、`w` は極の共鳴の
/// P^½ γ。M' は対称でないので、`pole_limit`（M' が複素対称で左右のベクトルが同じ形）は使えず、
/// M⁻¹ → M'⁻¹ − (M'⁻¹w)(vᵀM'⁻¹) / (vᵀM'⁻¹w)、v = Λw をここで求める。
fn x_at_shifted_pole(
    mut m: Vec<Complex64>,
    w: &[f64],
    lambda: &[Complex64],
    n: usize,
) -> Result<Vec<Complex64>, Singular> {
    let identity = |i: usize| {
        if i.is_multiple_of(n + 1) {
            Complex64::ONE
        } else {
            Complex64::ZERO
        }
    };
    let mut inverse: Vec<Complex64> = (0..n * n).map(identity).collect();
    solve(&mut m, &mut inverse, n)?;
    let u: Vec<Complex64> = (0..n)
        .map(|a| (0..n).map(|b| inverse[a * n + b] * w[b]).sum())
        .collect();
    let t: Vec<Complex64> = (0..n)
        .map(|b| (0..n).map(|a| lambda[a] * w[a] * inverse[a * n + b]).sum())
        .collect();
    let s: Complex64 = (0..n).map(|a| t[a] * w[a]).sum();
    if s == Complex64::ZERO || s.is_nan() {
        return Err(Singular);
    }
    let scale = s.inv();
    Ok((0..n * n)
        .map(|i| {
            let (a, b) = (i / n, i % n);
            (inverse[i] - u[a] * t[b] * scale - identity(i)) / lambda[b]
        })
        .collect())
}

#[cfg(test)]
pub(crate) mod tests;
