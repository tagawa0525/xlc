//! MF5（出射中性子のエネルギー分布）。表形式（LF=1）と、θ(E) で決まる形のうち
//! 一般の蒸発スペクトル（LF=5）、単純な核分裂スペクトル（LF=7）、蒸発スペクトル（LF=9）、
//! エネルギーによる Watt の分布（LF=11）、Madland-Nix の核分裂スペクトル（LF=12）を扱う。
//! それ以外の形式は、対象核種で必要になるまで追加しない。

use super::Error;
use super::tab1::{ClampedTab1, Interpolation, Tab1};
use super::tape::{Reader, Section};

/// 入射エネルギーごとに与えられた、出射エネルギーの確率密度の表。
#[derive(Debug, Clone, PartialEq)]
pub struct OutgoingDistributions {
    /// 入射エネルギー間の補間則の範囲（NBT）
    pub breaks: Vec<usize>,
    /// 入射エネルギー間の補間則
    pub laws: Vec<IncidentLaw>,
    /// 入射エネルギー (eV)
    pub incident: Vec<f64>,
    /// 各入射エネルギーでの出射エネルギーの確率密度 (1/eV)
    pub spectra: Vec<Tab1>,
}

/// 部分分布の形。LF=5、7、9、11 の解析式の分布は、出射エネルギー E' を 0 ≤ E' ≤ E - U に制限する。
#[derive(Debug, Clone, PartialEq)]
pub enum EnergyLaw {
    /// LF=1: 入射エネルギーごとの表
    Tabulated(OutgoingDistributions),
    /// LF=5: 一般の蒸発スペクトル f(E') ∝ g(E'/θ)。g(x) は表（`shape`）。θ はエネルギーによらない値だけを扱う
    /// （ENDF/B-VIII.0 の遅発中性子は θ = 1）
    GeneralEvaporation {
        restriction: f64,
        theta: f64,
        shape: Tab1,
    },
    /// LF=7: 単純な核分裂スペクトル f(E') ∝ √E' exp(-E'/θ(E))
    Maxwellian {
        restriction: f64,
        theta: ClampedTab1,
    },
    /// LF=9: 蒸発スペクトル f(E') ∝ E' exp(-E'/θ(E))
    Evaporation {
        restriction: f64,
        theta: ClampedTab1,
    },
    /// LF=11: エネルギーによる Watt の分布 f(E') ∝ exp(-E'/a(E)) sinh(√(b(E) E'))（ENDF-102 の 5.1.1.5）
    Watt {
        restriction: f64,
        /// a(E) (eV)
        a: ClampedTab1,
        /// b(E) (1/eV)
        b: ClampedTab1,
    },
    /// LF=12: Madland-Nix の核分裂スペクトル f(E') = (g(E', E_F(L)) + g(E', E_F(H)))/2（ENDF-102 の 5.1.1.6）。
    /// E' は 0〜∞ で、U で切らない
    MadlandNix {
        /// 軽い核分裂片の核子あたりの平均運動エネルギー E_F(L) (eV)
        light: f64,
        /// 重い核分裂片の E_F(H) (eV)
        heavy: f64,
        /// 核分裂片の最大の温度 T_M(E) (eV)
        temperature: ClampedTab1,
    },
}

/// MF5 の1反応。部分分布の確率 p_k(E) と、その分布の組。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf5 {
    pub subsections: Vec<(Tab1, EnergyLaw)>,
}

impl Mf5 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let mut r = section.reader();
        let head = r.cont()?;
        let nk = usize::try_from(head.n1)
            .map_err(|_| Error::new(r.line(), format!("invalid NK {}", head.n1)))?;
        let mut subsections = Vec::new();
        for _ in 0..nk {
            let start = r.line() + 1;
            let probability = r.tab1()?;
            // θ(E) で決まる分布では、p(E) の TAB1 の C1 が U で、θ(E) の TAB1 が続く（LF=5 はさらに g(x) の TAB1）。
            // LF=11 は a(E) と b(E) の TAB1 が続く（ENDF-102 の 5.2.5）
            let restriction = probability.c1;
            let law = match probability.l2 {
                1 => EnergyLaw::Tabulated(parse_outgoing(&mut r)?),
                5 => EnergyLaw::GeneralEvaporation {
                    restriction,
                    theta: constant_theta(&mut r)?,
                    shape: general_evaporation_shape(&mut r)?,
                },
                7 => EnergyLaw::Maxwellian {
                    restriction,
                    theta: positive_theta(&mut r)?,
                },
                9 => EnergyLaw::Evaporation {
                    restriction,
                    theta: positive_theta(&mut r)?,
                },
                11 => watt(restriction, &mut r)?,
                // ENDF-102 の 5.2.6 で p(E) の C1 は 0.0。E' は 0〜∞（5.1.1.6）なので、U を書いた評価は意味が決まらない
                12 if restriction != 0.0 => {
                    return Err(Error::new(
                        start,
                        format!("LF=12 has no restriction U, got {restriction}"),
                    ));
                }
                12 => madland_nix(&mut r)?,
                lf => return Err(Error::new(start, format!("LF={lf} is not supported"))),
            };
            subsections.push((probability.table, law));
        }
        Ok(Self { subsections })
    }
}

/// 解析式の分布の θ(E) を読む。θ は分布の温度なので正でなければならない。
fn positive_theta(r: &mut Reader<'_>) -> Result<ClampedTab1, Error> {
    let theta = r.tab1()?.table;
    clamped_positive(theta, r.line(), "theta")
}

/// LF=11 の c = a b/4 の上限。群への振り分け（multigroup::outgoing の `watt_cumulative`）は、x = E'/a = 1 で閉じた式と級数を
/// 切り替え、その前後の精度をこの範囲で確かめてある。閉じた式は a b/4 が大きいと桁落ちし（a b/4 = 25 の x = 1 で相対誤差
/// 1.3e-8、100 で桁が合わない）、級数は 700 ほどで溢れる。4 つの版の LF=11（ENDF/B-VII.0、VII.1 の U-233）は 0.71 以下
const WATT_MAX_C: f64 = 5.0;

/// Watt の分布の a(E)、b(E) を読む。a は分布の温度なので正。b = 0 では sinh が 0 になり分布が無く、負の b は sinh の引数が虚数になる。
fn watt(restriction: f64, r: &mut Reader<'_>) -> Result<EnergyLaw, Error> {
    let a = clamped_positive(r.tab1()?.table, r.line(), "a")?;
    let b = clamped_positive(r.tab1()?.table, r.line(), "b")?;
    // 表の間の補間の値は端の点の間にあり、表の外は端の値なので、a b/4 は両方の表の最大の積で抑えられる
    let largest = |t: &ClampedTab1| t.table().y.iter().copied().fold(0.0, f64::max);
    let c = 0.25 * largest(&a) * largest(&b);
    if c > WATT_MAX_C {
        return Err(Error::new(
            r.line(),
            format!("LF=11 with a b/4 up to {c} is not supported (accurate up to {WATT_MAX_C})"),
        ));
    }
    Ok(EnergyLaw::Watt { restriction, a, b })
}

/// Madland-Nix のスペクトルの E_F(L)、E_F(H)（TAB1 の C1、C2）と T_M(E) を読む。どれも正のエネルギーでなければならない。
fn madland_nix(r: &mut Reader<'_>) -> Result<EnergyLaw, Error> {
    let record = r.tab1()?;
    for fragment in [record.c1, record.c2] {
        if !(fragment > 0.0 && fragment.is_finite()) {
            return Err(Error::new(
                r.line(),
                format!("E_F must be positive, got {fragment}"),
            ));
        }
    }
    Ok(EnergyLaw::MadlandNix {
        light: record.c1,
        heavy: record.c2,
        temperature: clamped_positive(record.table, r.line(), "T_M")?,
    })
}

/// 分布の正の量（θ、T_M、Watt の a と b）の表を、範囲の外で端の値を使う表にする。`name` はエラーの文に出す量の名前。
fn clamped_positive(table: Tab1, line: usize, name: &str) -> Result<ClampedTab1, Error> {
    if let Some(bad) = table.y.iter().find(|&&t| t <= 0.0 || t.is_nan()) {
        return Err(Error::new(
            line,
            format!("{name} must be positive, got {bad}"),
        ));
    }
    ClampedTab1::new(table).ok_or_else(|| Error::new(line, format!("{name} table is empty")))
}

/// 一般の蒸発スペクトルの θ(E) を読む。エネルギーによる θ は、群の割合の折れ目の入射エネルギーが θ の逆関数で決まる
/// （積分の区間の切れ目に要る）ので、必要になるまで扱わない。
fn constant_theta(r: &mut Reader<'_>) -> Result<f64, Error> {
    let theta = positive_theta(r)?;
    let values = &theta.table().y;
    if values.iter().any(|&t| t != values[0]) {
        return Err(Error::new(
            r.line(),
            "LF=5 with an energy-dependent theta is not supported",
        ));
    }
    Ok(values[0])
}

/// 一般の蒸発スペクトルの g(x) を読む。g は確率密度の形なので負であってはならない。x = 0 で正の値から始めるのは、
/// E - U > 0 のどの入射エネルギーでも、群の割合の分母（E'/θ までの g の積分）が正であるため。
fn general_evaporation_shape(r: &mut Reader<'_>) -> Result<Tab1, Error> {
    let shape = r.tab1()?.table;
    if let Some(bad) = shape.y.iter().find(|&&g| g < 0.0 || g.is_nan()) {
        return Err(Error::new(
            r.line(),
            format!("g(x) must not be negative, got {bad}"),
        ));
    }
    if shape.x.first() != Some(&0.0) || shape.y[0] <= 0.0 {
        return Err(Error::new(
            r.line(),
            "g(x) must start at x = 0 with a positive value",
        ));
    }
    Ok(shape)
}

/// TAB2 と、それに続く NE 個の TAB1（C2 が入射エネルギー）を読む。
/// MF5 の LF=1 と MF6 の出射分布で共通の形。
pub fn parse_outgoing(r: &mut Reader<'_>) -> Result<OutgoingDistributions, Error> {
    let tab2 = r.tab2()?;
    let line = r.line();
    let mut incident = Vec::new();
    let mut spectra = Vec::new();
    for _ in 0..tab2.nz {
        let record = r.tab1()?;
        incident.push(record.c2);
        spectra.push(record.table);
    }
    let counts: Vec<usize> = spectra.iter().map(|t| t.x.len()).collect();
    let laws = incident_laws(&tab2.codes, &tab2.breaks, &counts, line)?;
    check_unit_base(&laws, &tab2.breaks, &spectra, line)?;
    check_incident_domain(&laws, &tab2.breaks, &incident, &spectra, line)?;
    Ok(OutgoingDistributions {
        breaks: tab2.breaks,
        laws,
        incident,
        spectra,
    })
}

/// 入射エネルギー間の出射分布の補間則（ENDF の INT）。
///
/// INT=1〜5 はふつうの補間、21〜25 は単位基底（出射エネルギーの範囲を 0〜1 に正規化して補間する）。
/// 11〜15 は対応点の補間で、隣り合う分布の点の数が違うと定義できないので、NJOY と同じく単位基底に
/// 切り替える。点の数が同じ対応点の補間は扱わない。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IncidentLaw {
    /// 補間の式（単位基底では正規化した座標で使う）
    pub law: Interpolation,
    pub unit_base: bool,
}

/// INT の値の方式（十の位）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Scheme {
    Cartesian,
    CorrespondingPoint,
    UnitBase,
}

/// 入射エネルギー間の補間則を、出射分布の点の数を見て決める。`counts` は各入射エネルギーの分布の点の数。
pub(super) fn incident_laws(
    codes: &[i64],
    breaks: &[usize],
    counts: &[usize],
    line: usize,
) -> Result<Vec<IncidentLaw>, Error> {
    let unsupported = |code: i64| {
        Error::new(
            line,
            format!("incident-energy interpolation INT={code} is not supported"),
        )
    };
    codes
        .iter()
        .enumerate()
        .map(|(range, &code)| {
            let scheme = match code / 10 {
                0 => Scheme::Cartesian,
                1 => Scheme::CorrespondingPoint,
                2 => Scheme::UnitBase,
                _ => return Err(unsupported(code)),
            };
            let law = Interpolation::from_code(code % 10).ok_or_else(|| unsupported(code))?;
            // この補間則の範囲の区間（点 k と k+1 の間）
            let first = if range == 0 { 0 } else { breaks[range - 1] - 1 };
            let intervals = first..breaks[range].saturating_sub(1);
            let unit_base = match scheme {
                // ヒストグラムは入射エネルギーの間を補間しないので、方式によらない
                _ if law == Interpolation::Histogram => false,
                Scheme::Cartesian => false,
                Scheme::UnitBase => true,
                // 点が1つの範囲には区間がなく、補間しない
                Scheme::CorrespondingPoint if intervals.is_empty() => false,
                Scheme::CorrespondingPoint => {
                    if intervals.clone().all(|k| counts[k] == counts[k + 1]) {
                        return Err(Error::new(
                            line,
                            format!(
                                "corresponding-point interpolation (INT={code}) with equal grids \
                                 is not supported"
                            ),
                        ));
                    }
                    xlc_log::warn!(
                        code,
                        "corresponding-point interpolation between grids of different sizes; \
                         using unit-base interpolation as NJOY does"
                    );
                    true
                }
            };
            Ok(IncidentLaw { law, unit_base })
        })
        .collect()
}

/// 単位基底で補間する分布は、2点以上で、出射エネルギーの範囲の幅が正であること（範囲で正規化するため）。
pub(super) fn check_unit_base(
    laws: &[IncidentLaw],
    breaks: &[usize],
    spectra: &[Tab1],
    line: usize,
) -> Result<(), Error> {
    for (range, law) in laws.iter().enumerate() {
        if !law.unit_base {
            continue;
        }
        let first = if range == 0 { 0 } else { breaks[range - 1] - 1 };
        for spectrum in &spectra[first..breaks[range].min(spectra.len())] {
            let width = match spectrum.x[..] {
                [low, .., high] => high - low,
                _ => 0.0,
            };
            if width <= 0.0 || width.is_nan() {
                return Err(Error::new(
                    line,
                    "unit-base interpolation needs outgoing spectra of positive width",
                ));
            }
        }
    }
    Ok(())
}

/// 入射エネルギーの間の補間則が、補間する値の定義域にあること（`Interpolation::check_interval`）。補間する値は、入射
/// エネルギーと、両端の分布の密度（群ごとの確率と単位基底の密度は、負でない密度から負でない値になる）と、単位基底では
/// 出射エネルギーの範囲の両端。密度は最も小さい値で確かめる。
pub(super) fn check_incident_domain(
    laws: &[IncidentLaw],
    breaks: &[usize],
    incident: &[f64],
    spectra: &[Tab1],
    line: usize,
) -> Result<(), Error> {
    let lowest = |t: &Tab1| t.y.iter().copied().fold(f64::INFINITY, f64::min);
    let ends = |t: &Tab1| match t.x[..] {
        [low, .., high] => (low, high),
        [only] => (only, only),
        [] => (0.0, 0.0),
    };
    for k in 0..incident.len().saturating_sub(1) {
        let IncidentLaw { law, unit_base } = laws[breaks.partition_point(|&b| b < k + 2)];
        let (e0, e1) = (incident[k], incident[k + 1]);
        let (lo, hi) = (&spectra[k], &spectra[k + 1]);
        let ((a0, b0), (a1, b1)) = (ends(lo), ends(hi));
        let range_ends = if unit_base {
            vec![(a0, a1), (b0, b1)]
        } else {
            Vec::new()
        };
        std::iter::once((lowest(lo), lowest(hi)))
            .chain(range_ends)
            .try_for_each(|(v0, v1)| law.check_interval((e0, v0), (e1, v1)))
            .map_err(|message| {
                Error::new(
                    line,
                    format!("incident energies {} and {}: {message}", k + 1, k + 2),
                )
            })?;
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests;
