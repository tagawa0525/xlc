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
        if head.n2 != 0 {
            return Err(Error::new(
                r.line(),
                "MF4 transformation matrices (NM > 0) are not supported",
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
            (3, 0) => Angular::Mixed {
                legendre: legendre_table(&mut r)?,
                tabulated: angle_tables(&mut r)?,
            },
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

    /// Legendre 係数 a_0..a_order（a_0 = 1）。表の分布は、補間した密度を両方の表の μ の点で区切って積分する
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
            Self::Tabulated { lower, upper, .. } => {
                let mut cuts: Vec<f64> = lower.1.x.iter().chain(&upper.1.x).copied().collect();
                cuts.sort_by(f64::total_cmp);
                cuts.dedup();
                // 対数の補間則では密度が多項式でないので、区間をさらに分ける
                for w in cuts.windows(2) {
                    let width = (w[1] - w[0]) / TABLE_PIECES as f64;
                    let half = 0.5 * width;
                    for piece in 0..TABLE_PIECES {
                        let mid = w[0] + width * (piece as f64 + 0.5);
                        for &(node, weight) in &GAUSS_LEGENDRE_8 {
                            let mu = mid + half * node;
                            let f = self.density(mu);
                            for (a, p) in a[1..].iter_mut().zip(&legendre(mu, order)[1..]) {
                                *a += weight * half * f * p;
                            }
                        }
                    }
                }
            }
        }
        a
    }
}

/// 表の分布のモーメントで、表の点の間をさらに等分する数
const TABLE_PIECES: usize = 4;

/// 表の μ での値（表は [−1, 1] を覆うことを読むときに確かめてある）
fn table_value(t: &Tab1, mu: f64) -> f64 {
    let k = t.x.partition_point(|&x| x <= mu).clamp(1, t.x.len() - 1);
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
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

    fn r(f: [&str; 6]) -> String {
        rec(f, 825, 4, 2)
    }

    fn parse(text: &str) -> Result<Mf4, Error> {
        let tape = Tape::parse(text).unwrap();
        Mf4::parse(&tape.section(825, 4, 2).unwrap())
    }

    fn close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len(), "{actual:?}");
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() < 1e-12, "{actual:?} vs {expected:?}");
        }
    }

    /// Legendre 係数（LTT=1）の部分: 入射エネルギー 1e5 と 2e5 eV で a_1 = 0.1 と 0.3、a_2 = 0.02 と 0.04
    fn legendre_part() -> String {
        [
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "1.0e5", "0", "0", "2", "0"]),
            r(["0.1", "0.02", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "2", "0"]),
            r(["0.3", "0.04", "", "", "", ""]),
        ]
        .concat()
    }

    /// 表（LTT=2）の部分: 2e5 eV で f(μ) = 0.5 + 0.3 μ、4e5 eV で等方
    fn tabulated_part() -> String {
        [
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.2", "1.0", "0.8", "", ""]),
            r(["0.0", "4.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
        ]
        .concat()
    }

    fn head(ltt: &str, li: &str, lct: &str) -> String {
        [
            r(["8016.0", "15.85751", "0", ltt, "0", "0"]),
            r(["0.0", "15.85751", li, lct, "0", "0"]),
        ]
        .concat()
    }

    #[test]
    fn isotropic_distribution() {
        let mf4 = parse(&head("0", "1", "1")).unwrap();
        assert_eq!(mf4.frame, Frame::Laboratory);
        assert_eq!(mf4.angular, Angular::Isotropic);
        close(&mf4.moments(1.0e6, 2), &[1.0, 0.0, 0.0]);
    }

    #[test]
    fn legendre_coefficients_are_interpolated_between_energies() {
        let mf4 = parse(&(head("1", "0", "2") + &legendre_part())).unwrap();
        assert_eq!(mf4.frame, Frame::CenterOfMass);
        close(&mf4.moments(1.5e5, 3), &[1.0, 0.2, 0.03, 0.0]);
        // 次数を打ち切る
        close(&mf4.moments(1.0e5, 1), &[1.0, 0.1]);
    }

    #[test]
    fn tabulated_distribution_gives_its_moments() {
        // f(μ) = 0.5 + 0.3 μ の a_1 = ∫ f μ dμ = 0.2、a_2 = 0
        let mf4 = parse(&(head("2", "0", "2") + &tabulated_part())).unwrap();
        close(&mf4.moments(2.0e5, 2), &[1.0, 0.2, 0.0]);
        close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
    }

    #[test]
    fn mixed_distribution_switches_at_the_last_legendre_energy() {
        let text = head("3", "0", "2") + &legendre_part() + &tabulated_part();
        let mf4 = parse(&text).unwrap();
        close(&mf4.moments(1.5e5, 2), &[1.0, 0.2, 0.03]);
        close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
    }

    #[test]
    fn log_lin_tables_are_integrated_exactly() {
        // ln f が μ に線形（INT=4）: f(μ) = 0.2 e^{b(μ+1)}、b = ln 3 / 2（f(1) = 0.6）の
        // a_1 = ∫ μ f dμ = 0.2 e^b [e^{bμ}(μ/b − 1/b²)]_{−1}^{1}
        let text = head("2", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "1"]),
                r(["1", "2", "", "", "", ""]),
                r(["0.0", "1.0e5", "0", "0", "1", "2"]),
                r(["2", "4", "", "", "", ""]),
                r(["-1.0", "0.2", "1.0", "0.6", "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        let b = 3.0_f64.ln() / 2.0;
        let primitive = |mu: f64| (b * mu).exp() * (mu / b - 1.0 / (b * b));
        // ∫ f dμ = 0.2 (e^{2b} − 1)/b で割って規格化する
        let norm = 0.2 * ((2.0 * b).exp() - 1.0) / b;
        let expected = 0.2 * b.exp() * (primitive(1.0) - primitive(-1.0)) / norm;
        let a = mf4.moments(1.0e5, 1);
        assert!((a[1] - expected).abs() < 1e-12, "{a:?} vs {expected}");
    }

    #[test]
    fn repeated_incident_energy_is_a_discontinuity() {
        // 先頭の入射エネルギーが2回ある（Pb-204 など）: 2つ目の分布から使う
        let text = head("1", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "3"]),
                r(["3", "2", "", "", "", ""]),
                r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
                r(["0.5", "", "", "", "", ""]),
                r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
                r(["0.0", "", "", "", "", ""]),
                r(["0.0", "2.0e5", "0", "0", "1", "0"]),
                r(["0.2", "", "", "", "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        // 1e-5 eV の a_1 = 0 と 2e5 eV の 0.2 の間の線形補間
        let expected = 0.2 * (1.0e5 - 1.0e-5) / (2.0e5 - 1.0e-5);
        close(&mf4.moments(1.0e5, 1), &[1.0, expected]);
    }

    #[test]
    fn density_of_each_form() {
        let iso = parse(&head("0", "1", "1")).unwrap();
        assert_eq!(iso.density(1.0e6, 0.3), 0.5);
        // Legendre: f(μ) = 1/2 + 3/2 a_1 μ + 5/2 a_2 P_2(μ)、1.5e5 eV で a_1 = 0.2、a_2 = 0.03
        let legendre = parse(&(head("1", "0", "2") + &legendre_part())).unwrap();
        let mu = 0.4_f64;
        let expected = 0.5 + 1.5 * 0.2 * mu + 2.5 * 0.03 * 0.5 * (3.0 * mu * mu - 1.0);
        assert!((legendre.density(1.5e5, mu) - expected).abs() < 1e-12);
        // 表: 2e5 eV の 0.5 + 0.3 μ と 4e5 eV の 0.5 を 3e5 eV で線形に補間
        let tabulated = parse(&(head("2", "0", "2") + &tabulated_part())).unwrap();
        assert!((tabulated.density(3.0e5, mu) - (0.5 + 0.15 * mu)).abs() < 1e-12);
    }

    #[test]
    fn repeated_first_energy_uses_the_later_distribution_there() {
        let text = head("1", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "3"]),
                r(["3", "2", "", "", "", ""]),
                r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
                r(["0.5", "", "", "", "", ""]),
                r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
                r(["0.0", "", "", "", "", ""]),
                r(["0.0", "2.0e5", "0", "0", "1", "0"]),
                r(["0.2", "", "", "", "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        close(&mf4.moments(1.0e-5, 1), &[1.0, 0.0]);
    }

    #[test]
    fn tabulated_moments_follow_the_interpolated_density() {
        // 入射エネルギーの間を INT=4（ln f が E に線形）で補間する: モーメントは補間した密度のモーメント
        let text = head("2", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "2"]),
                r(["2", "4", "", "", "", ""]),
                r(["0.0", "2.0e5", "0", "0", "1", "2"]),
                r(["2", "2", "", "", "", ""]),
                r(["-1.0", "0.2", "1.0", "0.8", "", ""]),
                r(["0.0", "4.0e5", "0", "0", "1", "2"]),
                r(["2", "2", "", "", "", ""]),
                r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        let e = 3.0e5;
        let nodes = xlc_numeric::quadrature::gauss_legendre(64).unwrap();
        let integral = |g: &dyn Fn(f64) -> f64| -> f64 {
            nodes
                .iter()
                .map(|&(x, w)| w * g(x) * mf4.density(e, x))
                .sum()
        };
        // 補間した密度は規格化されていないので、∫ f dμ で割る
        let expected = integral(&|x| x) / integral(&|_| 1.0);
        let a = mf4.moments(e, 1);
        assert!((a[1] - expected).abs() < 1e-12, "{a:?} vs {expected}");
    }

    #[test]
    fn tables_must_cover_all_angles() {
        let text = head("2", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "1"]),
                r(["1", "2", "", "", "", ""]),
                r(["0.0", "2.0e5", "0", "0", "1", "2"]),
                r(["2", "2", "", "", "", ""]),
                r(["-0.5", "0.5", "1.0", "0.5", "", ""]),
            ]
            .concat();
        assert!(parse(&text).is_err());
    }

    #[test]
    fn transformation_matrix_is_not_supported() {
        // LVT = 1（HEAD の L1。Legendre 係数の変換行列が続く）は使われていないので扱わない
        let text = [
            r(["8016.0", "15.85751", "1", "0", "0", "0"]),
            r(["0.0", "15.85751", "1", "1", "0", "0"]),
        ]
        .concat();
        assert!(parse(&text).is_err());
    }

    #[test]
    fn tabulated_moments_are_normalized() {
        // f(μ) = 0.98 (0.5 + 0.3 μ)（∫ f = 0.98）でも a_1 = 0.2
        let text = head("2", "0", "1")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "1"]),
                r(["1", "2", "", "", "", ""]),
                r(["0.0", "2.0e5", "0", "0", "1", "2"]),
                r(["2", "2", "", "", "", ""]),
                r(["-1.0", "0.196", "1.0", "0.784", "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        close(&mf4.moments(2.0e5, 1), &[1.0, 0.2]);
    }

    #[test]
    fn mixed_parts_must_meet() {
        // Legendre の部分が 2e5 eV で終わり、表の部分が 4e5 eV から始まる
        let tabulated = [
            r(["0.0", "0.0", "0", "0", "1", "1"]),
            r(["1", "2", "", "", "", ""]),
            r(["0.0", "4.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
        ]
        .concat();
        let text = head("3", "0", "2") + &legendre_part() + &tabulated;
        assert!(parse(&text).is_err());
    }
}
