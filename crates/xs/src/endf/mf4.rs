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
    /// 入射エネルギーの間は、係数（表の分布ではそのモーメント）を ENDF の補間則で補間する。範囲の外は端の値とする。
    pub fn moments(&self, e: f64, order: usize) -> Vec<f64> {
        let mut a = vec![0.0; order + 1];
        a[0] = 1.0;
        let higher = match &self.angular {
            Angular::Isotropic => return a,
            Angular::Legendre(table) => table.moments(e, order),
            Angular::Tabulated(tables) => tables.moments(e, order),
            Angular::Mixed {
                legendre,
                tabulated,
            } => {
                if legendre.energies.last().is_some_and(|&top| e <= top) {
                    legendre.moments(e, order)
                } else {
                    tabulated.moments(e, order)
                }
            }
        };
        a[1..].copy_from_slice(&higher);
        a
    }
}

impl LegendreTable {
    /// a_1..a_order（与えられた次数より上は 0）
    fn moments(&self, e: f64, order: usize) -> Vec<f64> {
        let at = |k: usize| -> Vec<f64> {
            (1..=order)
                .map(|l| self.coefficients[k].get(l - 1).copied().unwrap_or(0.0))
                .collect()
        };
        interpolate(&self.breaks, &self.laws, &self.energies, e, at)
    }
}

impl AngleTables {
    /// a_l = ∫ f(μ) P_l(μ) dμ（l = 1..order）。f は区間ごとに線形（INT=2）として積分する
    fn moments(&self, e: f64, order: usize) -> Vec<f64> {
        let at = |k: usize| -> Vec<f64> {
            let t = &self.tables[k];
            let mut a = vec![0.0; order];
            for i in 1..t.x.len() {
                let (x0, x1) = (t.x[i - 1], t.x[i]);
                let (mid, half) = (0.5 * (x0 + x1), 0.5 * (x1 - x0));
                for &(node, weight) in &GAUSS_LEGENDRE_8 {
                    let mu = mid + half * node;
                    let f = t.y[i - 1] + (t.y[i] - t.y[i - 1]) * (mu - x0) / (x1 - x0);
                    for (a, p) in a.iter_mut().zip(legendre(mu, order).into_iter().skip(1)) {
                        *a += weight * half * f * p;
                    }
                }
            }
            a
        };
        interpolate(&self.breaks, &self.laws, &self.energies, e, at)
    }
}

/// 入射エネルギー `e` での値を、入射エネルギーごとの値 `at(k)` から補間則で求める（範囲の外は端の値）。
fn interpolate(
    breaks: &[usize],
    laws: &[Interpolation],
    energies: &[f64],
    e: f64,
    at: impl Fn(usize) -> Vec<f64>,
) -> Vec<f64> {
    let n = energies.len();
    if e <= energies[0] {
        return at(0);
    }
    if e >= energies[n - 1] {
        return at(n - 1);
    }
    let k = energies.partition_point(|&x| x <= e);
    let law = breaks
        .iter()
        .position(|&b| k < b)
        .map_or(Interpolation::LinLin, |i| laws[i]);
    let (lo, hi) = (at(k - 1), at(k));
    lo.iter()
        .zip(&hi)
        .map(|(&y0, &y1)| law.interpolate(e, (energies[k - 1], y0), (energies[k], y1)))
        .collect()
}

/// TAB2 と、それに続く入射エネルギーごとの LIST（Legendre 係数）を読む。
fn legendre_table(r: &mut Reader<'_>) -> Result<LegendreTable, Error> {
    let tab2 = r.tab2()?;
    let laws = laws(&tab2.codes, r.line())?;
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

/// TAB2 と、それに続く入射エネルギーごとの TAB1（f(μ)）を読む。f は線形補間（INT=2）だけを扱う。
fn angle_tables(r: &mut Reader<'_>) -> Result<AngleTables, Error> {
    let tab2 = r.tab2()?;
    let laws = laws(&tab2.codes, r.line())?;
    let mut energies = Vec::with_capacity(tab2.nz);
    let mut tables = Vec::with_capacity(tab2.nz);
    for _ in 0..tab2.nz {
        let record = r.tab1()?;
        if record
            .table
            .laws
            .iter()
            .any(|&law| law != Interpolation::LinLin)
        {
            return Err(Error::new(
                r.line(),
                format!(
                    "angular tables must be lin-lin in mu, got {:?}",
                    record.table.laws
                ),
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
    if let Some(w) = energies.windows(2).find(|w| w[0] >= w[1]) {
        return Err(Error::new(
            line,
            format!("incident energies must be increasing: {w:?}"),
        ));
    }
    if energies.is_empty() {
        return Err(Error::new(line, "no incident energies"));
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
        // ln f が μ に線形（INT=4）: f(μ) = c e^μ（c = 1/(e − 1/e)）の a_1 = c ∫ μ e^μ dμ = 2c/e
        let c = 1.0 / (1.0_f64.exp() - (-1.0_f64).exp());
        let (lo, hi) = (
            format!("{}", c / 1.0_f64.exp()),
            format!("{}", c * 1.0_f64.exp()),
        );
        let text = head("2", "0", "2")
            + &[
                r(["0.0", "0.0", "0", "0", "1", "1"]),
                r(["1", "2", "", "", "", ""]),
                r(["0.0", "1.0e5", "0", "0", "1", "2"]),
                r(["2", "4", "", "", "", ""]),
                r(["-1.0", &lo, "1.0", &hi, "", ""]),
            ]
            .concat();
        let mf4 = parse(&text).unwrap();
        let a = mf4.moments(1.0e5, 1);
        assert!((a[1] - 2.0 * c / 1.0_f64.exp()).abs() < 1e-12, "{a:?}");
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
        close(&mf4.moments(1.0e5, 1), &[1.0, 0.1]);
    }

    #[test]
    fn transformation_matrix_is_not_supported() {
        // NM > 0（Legendre 係数の変換行列）は使われていないので扱わない
        let text = [
            r(["8016.0", "15.85751", "0", "1", "0", "1"]),
            r(["0.0", "15.85751", "0", "2", "0", "1"]),
        ]
        .concat();
        assert!(parse(&text).is_err());
    }
}
