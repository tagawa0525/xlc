//! MF4（出射中性子の角度分布）。等方（LTT=0）、Legendre 係数（LTT=1）、表（LTT=2）と、低いエネルギーで
//! Legendre 係数・高いエネルギーで表を使う混合（LTT=3）を扱う。

use super::Error;
use super::mf6::Frame;
use super::tab1::{Interpolation, Tab1};
use super::tape::{Reader, Section};
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_8, legendre};

/// 入射エネルギーごとの Legendre 係数 a_l（l = 1..NL）。f(μ) = Σ (2l+1)/2 a_l P_l(μ)、a_0 = 1。
#[derive(Debug, Clone, PartialEq)]
pub struct LegendreTable {
    /// 入射エネルギー間の補間則の範囲（NBT）と補間則
    pub breaks: Vec<usize>,
    pub laws: Vec<Interpolation>,
    /// 入射エネルギー (eV)
    pub energies: Vec<f64>,
    /// [入射エネルギー][l − 1]
    pub coefficients: Vec<Vec<f64>>,
}

/// 入射エネルギーごとの散乱角の余弦 μ の確率密度 f(μ) の表。
#[derive(Debug, Clone, PartialEq)]
pub struct AngleTables {
    pub breaks: Vec<usize>,
    pub laws: Vec<Interpolation>,
    pub energies: Vec<f64>,
    pub tables: Vec<Tab1>,
}

/// 角度分布の形。
#[derive(Debug, Clone, PartialEq)]
pub enum Angular {
    /// LTT=0（LI=1）: 等方
    Isotropic,
    /// LTT=1
    Legendre(LegendreTable),
    /// LTT=2
    Tabulated(AngleTables),
    /// LTT=3: `legendre` の最後の入射エネルギーまでは Legendre 係数、その上は表
    Mixed {
        legendre: LegendreTable,
        tabulated: AngleTables,
    },
}

/// MF4 の1反応。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf4 {
    pub frame: Frame,
    pub angular: Angular,
}

impl Mf4 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let mut r = section.reader();
        let head = r.cont()?;
        let ltt = head.l2;
        // LVT = 1 なら、Legendre 係数の変換行列の LIST が続く（ENDF/B-VIII.0 では使われていない）
        if head.l1 != 0 {
            return Err(Error::new(
                r.line(),
                "MF4 transformation matrices (LVT = 1) are not supported",
            ));
        }
        let cont = r.cont()?;
        let frame = match cont.l2 {
            1 => Frame::Laboratory,
            2 => Frame::CenterOfMass,
            lct => return Err(Error::new(r.line(), format!("LCT={lct} is not supported"))),
        };
        let angular = match (ltt, cont.l1) {
            (0, 1) => Angular::Isotropic,
            (1, 0) => Angular::Legendre(legendre_table(&mut r)?),
            (2, 0) => Angular::Tabulated(angle_tables(&mut r)?),
            (3, 0) => {
                let legendre = legendre_table(&mut r)?;
                let tabulated = angle_tables(&mut r)?;
                // 表の部分は Legendre の部分の最後のエネルギーから始まる
                if legendre.energies.last() != tabulated.energies.first() {
                    return Err(Error::new(
                        r.line(),
                        "the tabulated part must start where the Legendre part ends",
                    ));
                }
                Angular::Mixed {
                    legendre,
                    tabulated,
                }
            }
            (ltt, li) => {
                return Err(Error::new(
                    r.line(),
                    format!("LTT={ltt} with LI={li} is not supported"),
                ));
            }
        };
        Ok(Self { frame, angular })
    }

    /// 入射エネルギー `e` での角度分布の Legendre 係数 a_0..a_order（a_0 = 1。`frame` の系の μ について）。
    pub fn moments(&self, e: f64, order: usize) -> Vec<f64> {
        self.at(e).moments(order)
    }

    /// 入射エネルギー `e` での確率密度 f(μ)（`frame` の系の μ について。∫ f dμ = 1）。
    pub fn density(&self, e: f64, mu: f64) -> f64 {
        self.at(e).density(mu)
    }

    /// 入射エネルギー `e` での角度分布。入射エネルギーの間は、Legendre 係数、または同じ μ での表の値を、
    /// ENDF の補間則で補間する。範囲の外は端の分布とし、同じエネルギーが続くところでは後の分布を使う。
    pub fn at(&self, e: f64) -> Distribution<'_> {
        match &self.angular {
            Angular::Isotropic => Distribution::Isotropic,
            Angular::Legendre(table) => table.at(e),
            Angular::Tabulated(tables) => tables.at(e),
            Angular::Mixed {
                legendre,
                tabulated,
            } => {
                if legendre.energies.last().is_some_and(|&top| e <= top) {
                    legendre.at(e)
                } else {
                    tabulated.at(e)
                }
            }
        }
    }

    /// 分布を与えている入射エネルギー (eV)。分布はここで折れ曲がるので、入射エネルギーの積分の区切りにする
    pub fn energies(&self) -> Vec<f64> {
        match &self.angular {
            Angular::Isotropic => Vec::new(),
            Angular::Legendre(table) => table.energies.clone(),
            Angular::Tabulated(tables) => tables.energies.clone(),
            Angular::Mixed {
                legendre,
                tabulated,
            } => [&legendre.energies[..], &tabulated.energies[..]].concat(),
        }
    }
}

/// 1つの入射エネルギーでの角度分布。
#[derive(Debug, Clone, PartialEq)]
pub enum Distribution<'a> {
    Isotropic,
    /// Legendre 係数 a_1..a_NL
    Legendre(Vec<f64>),
    /// 2つの入射エネルギーの表の間（同じ表なら `fraction` は使わない）。f は同じ μ の値を `law` で補間する
    Tabulated {
        lower: (f64, &'a Tab1),
        upper: (f64, &'a Tab1),
        law: Interpolation,
        energy: f64,
    },
}

impl Distribution<'_> {
    /// 確率密度 f(μ)
    pub fn density(&self, mu: f64) -> f64 {
        match self {
            Self::Isotropic => 0.5,
            Self::Legendre(a) => {
                let p = legendre(mu, a.len());
                0.5 + a
                    .iter()
                    .zip(&p[1..])
                    .enumerate()
                    .map(|(k, (a, p))| (2 * k + 3) as f64 / 2.0 * a * p)
                    .sum::<f64>()
            }
            Self::Tabulated {
                lower,
                upper,
                law,
                energy,
            } => {
                let f = |t: &Tab1| table_value(t, mu);
                if std::ptr::eq(lower.1, upper.1) {
                    return f(lower.1);
                }
                law.interpolate(*energy, (lower.0, f(lower.1)), (upper.0, f(upper.1)))
            }
        }
    }

    /// Legendre 係数 a_0..a_order（a_0 = 1）。表の分布は、補間した密度を両方の表の μ の点で区切って積分し、
    /// ∫ f dμ で規格化する（評価済みデータの表は丸めで規格化がわずかにずれる）
    pub fn moments(&self, order: usize) -> Vec<f64> {
        let mut a = vec![0.0; order + 1];
        a[0] = 1.0;
        match self {
            Self::Isotropic => {}
            Self::Legendre(coefficients) => {
                for (a, c) in a[1..].iter_mut().zip(coefficients) {
                    *a = *c;
                }
            }
            Self::Tabulated { .. } => {
                let cuts = self.kinks();
                let mut norm = 0.0;
                // 対数の補間則では密度が多項式でないので、区間をさらに分ける
                for w in cuts.windows(2) {
                    let width = (w[1] - w[0]) / TABLE_PIECES as f64;
                    let half = 0.5 * width;
                    for piece in 0..TABLE_PIECES {
                        let mid = w[0] + width * (piece as f64 + 0.5);
                        for &(node, weight) in &GAUSS_LEGENDRE_8 {
                            let mu = mid + half * node;
                            let f = self.density(mu);
                            norm += weight * half * f;
                            for (a, p) in a[1..].iter_mut().zip(&legendre(mu, order)[1..]) {
                                *a += weight * half * f * p;
                            }
                        }
                    }
                }
                if norm > 0.0 {
                    for v in &mut a[1..] {
                        *v /= norm;
                    }
                }
            }
        }
        a
    }

    /// 密度が折れ曲がる μ（表の分布の両方の表の点。昇順）。それ以外の分布では μ = −1, 1 だけ
    pub fn kinks(&self) -> Vec<f64> {
        let Self::Tabulated { lower, upper, .. } = self else {
            return vec![-1.0, 1.0];
        };
        let mut cuts: Vec<f64> = lower.1.x.iter().chain(&upper.1.x).copied().collect();
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        cuts
    }

    /// Legendre 係数の分布の次数（それ以外は 0）
    pub fn order(&self) -> usize {
        match self {
            Self::Legendre(a) => a.len(),
            _ => 0,
        }
    }
}

/// 表の分布のモーメントで、表の点の間をさらに等分する数
const TABLE_PIECES: usize = 4;

/// 表の μ での値（表は [−1, 1] を覆うことを読むときに確かめてある。点が2つ未満の表は 0）
fn table_value(t: &Tab1, mu: f64) -> f64 {
    let n = t.x.len().min(t.y.len());
    if n < 2 {
        return 0.0;
    }
    let k = t.x.partition_point(|&x| x <= mu).max(1).min(n - 1);
    t.law_of_interval(k - 1)
        .interpolate(mu, (t.x[k - 1], t.y[k - 1]), (t.x[k], t.y[k]))
}

/// 入射エネルギー `e` を挟む点の番号（範囲の外と、ちょうど点の上では同じ番号を2つ）と、その間の補間則。
/// 同じエネルギーが続くところでは後の点を使う。
fn bracket(
    breaks: &[usize],
    laws: &[Interpolation],
    energies: &[f64],
    e: f64,
) -> (usize, usize, Interpolation) {
    let n = energies.len();
    let k = energies.partition_point(|&x| x <= e);
    if k == 0 {
        return (0, 0, Interpolation::LinLin);
    }
    if k == n || energies[k - 1] == e {
        return (k - 1, k - 1, Interpolation::LinLin);
    }
    // 区間 (k−1, k) の補間則（TAB2 の読み込みで最後の NBT が点の数に等しいことを確かめてある）
    let law = laws[breaks.partition_point(|&b| b <= k)];
    (k - 1, k, law)
}

impl LegendreTable {
    fn at(&self, e: f64) -> Distribution<'_> {
        let (i, j, law) = bracket(&self.breaks, &self.laws, &self.energies, e);
        let order = self.coefficients[i].len().max(self.coefficients[j].len());
        let get = |k: usize, l: usize| self.coefficients[k].get(l).copied().unwrap_or(0.0);
        let a = (0..order)
            .map(|l| {
                if i == j {
                    get(i, l)
                } else {
                    law.interpolate(
                        e,
                        (self.energies[i], get(i, l)),
                        (self.energies[j], get(j, l)),
                    )
                }
            })
            .collect();
        Distribution::Legendre(a)
    }
}

impl AngleTables {
    fn at(&self, e: f64) -> Distribution<'_> {
        let (i, j, law) = bracket(&self.breaks, &self.laws, &self.energies, e);
        Distribution::Tabulated {
            lower: (self.energies[i], &self.tables[i]),
            upper: (self.energies[j], &self.tables[j]),
            law,
            energy: e,
        }
    }
}

/// TAB2 と、それに続く入射エネルギーごとの LIST（Legendre 係数）を読む。
fn legendre_table(r: &mut Reader<'_>) -> Result<LegendreTable, Error> {
    let tab2 = r.tab2()?;
    let laws = laws(&tab2.codes, r.line())?;
    // 係数は負にもなるので、ln を取る補間則は使えない
    if let Some(law) = laws
        .iter()
        .find(|law| matches!(law, Interpolation::LogLin | Interpolation::LogLog))
    {
        return Err(Error::new(
            r.line(),
            format!("{law:?} interpolation of Legendre coefficients is not defined"),
        ));
    }
    let mut energies = Vec::with_capacity(tab2.nz);
    let mut coefficients = Vec::with_capacity(tab2.nz);
    for _ in 0..tab2.nz {
        let list = r.list()?;
        energies.push(list.c2);
        coefficients.push(list.values);
    }
    increasing(&energies, r.line())?;
    Ok(LegendreTable {
        breaks: tab2.breaks,
        laws,
        energies,
        coefficients,
    })
}

/// TAB2 と、それに続く入射エネルギーごとの TAB1（f(μ)）を読む。μ は負にもなるので、ln μ を使う補間則
/// （INT=3, 5）は定義できず、エラーにする。表は μ = −1〜1 を覆い、f ≥ 0 でなければならない。
fn angle_tables(r: &mut Reader<'_>) -> Result<AngleTables, Error> {
    let tab2 = r.tab2()?;
    let laws = laws(&tab2.codes, r.line())?;
    let mut energies = Vec::with_capacity(tab2.nz);
    let mut tables = Vec::with_capacity(tab2.nz);
    for _ in 0..tab2.nz {
        let record = r.tab1()?;
        if let Some(law) = record
            .table
            .laws
            .iter()
            .find(|law| matches!(law, Interpolation::LinLog | Interpolation::LogLog))
        {
            return Err(Error::new(
                r.line(),
                format!("{law:?} interpolation in mu is not defined for negative mu"),
            ));
        }
        let t = &record.table;
        if t.x.first() != Some(&-1.0) || t.x.last() != Some(&1.0) {
            return Err(Error::new(
                r.line(),
                "angular tables must cover mu from -1 to 1",
            ));
        }
        if let Some(bad) = t.y.iter().find(|&&f| !(f >= 0.0 && f.is_finite())) {
            return Err(Error::new(
                r.line(),
                format!("negative or invalid density {bad}"),
            ));
        }
        energies.push(record.c2);
        tables.push(record.table);
    }
    increasing(&energies, r.line())?;
    Ok(AngleTables {
        breaks: tab2.breaks,
        laws,
        energies,
        tables,
    })
}

fn laws(codes: &[i64], line: usize) -> Result<Vec<Interpolation>, Error> {
    codes
        .iter()
        .map(|&c| {
            Interpolation::from_code(c)
                .ok_or_else(|| Error::new(line, format!("invalid interpolation code {c}")))
        })
        .collect()
}

fn increasing(energies: &[f64], line: usize) -> Result<(), Error> {
    if energies.is_empty() {
        return Err(Error::new(line, "no incident energies"));
    }
    // 同じエネルギーが続くのは不連続（後の分布をそのエネルギーから使う）
    if let Some(w) = energies.windows(2).find(|w| w[0] > w[1]) {
        return Err(Error::new(
            line,
            format!("incident energies must not decrease: {w:?}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
