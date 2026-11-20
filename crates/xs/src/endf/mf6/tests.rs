use super::*;
use crate::endf::{IncidentLaw, Tape, rec};

fn r(f: [&str; 6]) -> String {
    rec(f, 825, 6, 16)
}

/// 中性子（LAW=1、2つの入射エネルギー）と、LAW=0 の別の生成粒子
fn section(lct: &str, neutron_law: &str, nd: &str) -> String {
    section_with_int(lct, neutron_law, nd, "22")
}

/// `section` の、入射エネルギーの補間則を `int` にしたもの（出射の点は 2 点と 3 点）
fn section_with_int(lct: &str, neutron_law: &str, nd: &str, int: &str) -> String {
    [
        r(["8016.0", "15.85751", "0", lct, "2", "0"]),
        // 別の生成粒子（LAW=0: 多重度だけ）
        r(["8015.0", "14.87", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.6e7", "1.0", "2.0e7", "1.0", "", ""]),
        // 中性子
        r(["1.0", "1.0", "0", neutron_law, "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.6e7", "2.0", "2.0e7", "2.0", "", ""]),
        r(["0.0", "0.0", "1", "2", "1", "2"]),
        r(["2", int, "", "", "", ""]),
        r(["0.0", "1.6e7", nd, "0", "4", "2"]),
        r(["0.0", "2.0e-6", "1.0e6", "0.0", "", ""]),
        r(["0.0", "2.0e7", "0", "1", "9", "3"]),
        // (E', b0, b1) の3つ組が3点
        r(["0.0", "0.0", "0.0", "1.0e6", "2.0e-6", "0.5"]),
        r(["2.0e6", "0.0", "0.0", "", "", ""]),
    ]
    .concat()
}

#[test]
fn reads_law1_neutrons_and_skips_other_products() {
    let tape = Tape::parse(&section("2", "1", "0")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert_eq!(mf6.neutrons.len(), 1);
    assert_eq!(mf6.awr, 15.85751);
    let n = &mf6.neutrons[0];
    assert_eq!(n.frame, Frame::CenterOfMass);
    assert_eq!(n.multiplicity.eval(1.8e7), Some(2.0));
    let NeutronSpectrum::Tabulated(d) = &n.spectrum else {
        panic!("LAW=1 is a table: {:?}", n.spectrum);
    };
    assert_eq!(d.incident, vec![1.6e7, 2.0e7]);
    // INT=22 は単位基底の線形補間
    assert_eq!(
        d.laws,
        vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true
        }]
    );
    assert_eq!(d.spectra[0].x, vec![0.0, 1.0e6]);
    assert_eq!(d.spectra[0].y, vec![2.0e-6, 0.0]);
    // NA=1 なので1点あたり3個の値。b0 だけを取る
    assert_eq!(d.spectra[1].x, vec![0.0, 1.0e6, 2.0e6]);
    assert_eq!(d.spectra[1].y, vec![0.0, 2.0e-6, 0.0]);
}

#[test]
fn two_body_neutron_law_is_marked() {
    // 中性子が LAW=2（2体反応の角度分布、入射エネルギー2点）で与えられた準位の非弾性散乱
    let text = [
        r(["50120.0", "118.872", "0", "2", "1", "0"]),
        r(["1.0", "1.0", "0", "2", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.2e6", "1.0", "2.0e7", "1.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        // E, LANG=0（Legendre 係数）, NW=1, NL=1
        r(["0.0", "1.2e6", "0", "0", "1", "1"]),
        r(["0.1", "", "", "", "", ""]),
        r(["0.0", "2.0e7", "0", "0", "1", "1"]),
        r(["0.3", "", "", "", "", ""]),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert!(mf6.two_body);
    assert!(mf6.neutrons.is_empty());
}

#[test]
fn isotropic_two_body_neutron_law_is_marked() {
    // 中性子が LAW=3（重心系で等方な2体反応。分布のレコードはない）の準位の非弾性散乱
    let text = [
        r(["24052.0", "51.5493", "0", "2", "1", "0"]),
        r(["1.0", "1.0", "0", "3", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.5e6", "1.0", "2.0e7", "1.0", "", ""]),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert!(mf6.two_body);
    assert!(mf6.neutrons.is_empty());
}

#[test]
fn incident_interpolation_schemes() {
    let law = |int: &str| {
        let tape = Tape::parse(&section_with_int("2", "1", "0", int)).unwrap();
        Mf6::parse(&tape.section(825, 6, 16).unwrap()).map(|mf6| match &mf6.neutrons[0].spectrum {
            NeutronSpectrum::Tabulated(d) => d.laws[0],
            other => panic!("LAW=1 is a table: {other:?}"),
        })
    };
    // ふつうの補間（ENDF のとおり、単位基底にはしない）
    assert_eq!(
        law("2").unwrap(),
        IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false
        }
    );
    // 対応点（INT=12）で隣の分布と点の数が違えば、NJOY と同じく単位基底に切り替える
    assert_eq!(
        law("12").unwrap(),
        IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true
        }
    );
    assert!(law("6").is_err());
}

/// 2つの入射エネルギーの出射分布がどちらも 2 点で、入射エネルギーの補間則が `int` の節
fn equal_grids(int: &str) -> String {
    [
        r(["8016.0", "15.85751", "0", "2", "1", "0"]),
        r(["1.0", "1.0", "0", "1", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["1.6e7", "2.0", "2.0e7", "2.0", "", ""]),
        r(["0.0", "0.0", "1", "2", "1", "2"]),
        r(["2", int, "", "", "", ""]),
        r(["0.0", "1.6e7", "0", "0", "4", "2"]),
        r(["0.0", "1.0e-6", "1.0e6", "1.0e-6", "", ""]),
        r(["0.0", "2.0e7", "0", "0", "4", "2"]),
        r(["0.0", "5.0e-7", "2.0e6", "5.0e-7", "", ""]),
    ]
    .concat()
}

#[test]
fn corresponding_points_with_equal_grids_are_not_supported() {
    let tape = Tape::parse(&equal_grids("12")).unwrap();
    assert!(Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_err());
    // ヒストグラムは入射エネルギーの間を補間しないので、対応点でも扱える
    let tape = Tape::parse(&equal_grids("11")).unwrap();
    assert!(Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_ok());
}

#[test]
fn other_neutron_law_is_error() {
    // LAW=5 は荷電粒子の弾性散乱で、中性子の分布ではない
    let tape = Tape::parse(&section("2", "5", "0")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 5);
}

#[test]
fn laboratory_frame() {
    let tape = Tape::parse(&section("1", "1", "0")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert_eq!(mf6.neutrons[0].frame, Frame::Laboratory);
}

#[test]
fn unsupported_frame_is_error() {
    // LCT=4 は生成粒子ごとに座標系が変わる（M0 では扱わない）
    let tape = Tape::parse(&section("4", "1", "0")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 1);
}

#[test]
fn discrete_lines_are_error() {
    let tape = Tape::parse(&section("2", "1", "1")).unwrap();
    assert!(Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_err());
}

/// H-2 の (n,2n) の形: 中性子と陽子が LAW=6（3 体、APSX = 2.99862）で、LCT は `lct`、粒子の数は `npsx`
fn phase_space(lct: &str, apsx: &str, npsx: &str) -> String {
    [
        r(["1002.0", "1.9968", "0", lct, "2", "0"]),
        r(["1.0", "1.0", "0", "6", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["3.339002e6", "2.0", "1.5e8", "2.0", "", ""]),
        r([apsx, "0.0", "0", "0", "0", npsx]),
        r(["1001.0", "0.99862", "0", "6", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["3.339002e6", "1.0", "1.5e8", "1.0", "", ""]),
        r([apsx, "0.0", "0", "0", "0", npsx]),
    ]
    .concat()
}

#[test]
fn law6_neutrons_are_phase_space_in_the_centre_of_mass() {
    // LAW=6 はいつも重心系（LCT=1 でも）
    let tape = Tape::parse(&phase_space("1", "2.99862", "3")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert_eq!(mf6.neutrons.len(), 1);
    assert!(mf6.has_phase_space());
    let n = &mf6.neutrons[0];
    assert_eq!(n.frame, Frame::CenterOfMass);
    assert_eq!(n.multiplicity.eval(1.0e7), Some(2.0));
    assert_eq!(
        n.spectrum,
        NeutronSpectrum::PhaseSpace(PhaseSpace {
            total_mass: 2.99862,
            particles: 3
        })
    );
}

#[test]
fn law1_neutrons_have_no_phase_space() {
    let tape = Tape::parse(&section("2", "1", "0")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert!(!mf6.has_phase_space());
}

#[test]
fn phase_space_outside_three_to_five_particles_is_error() {
    // ENDF-102 は n = 3〜5 の係数 C_n だけを与える
    for npsx in ["2", "6"] {
        let tape = Tape::parse(&phase_space("1", "2.99862", npsx)).unwrap();
        let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
        assert_eq!(err.line, 5, "NPSX={npsx}: {err}");
    }
}

#[test]
fn two_step_phase_space_is_error() {
    // APSX が生成粒子の質量の和（A + 1 ほど）より小さいのは 2 段の反応で、E_a は 1 段目の反跳のエネルギーになる
    // （ENDF-102 の 6.2.8）。式 6.28 の E_a では E_max が大きすぎる
    let tape = Tape::parse(&phase_space("1", "2.0", "3")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 5, "{err}");
    // 1 段の反応の APSX は A + 1 から Q/(m_n c²) ほど（0.002〜0.004）ずれるだけで、通る。境は中性子の質量の半分
    for (apsx, ok) in [("2.9928", true), ("2.5", true), ("2.49", false)] {
        let tape = Tape::parse(&phase_space("1", apsx, "3")).unwrap();
        assert_eq!(
            Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_ok(),
            ok,
            "APSX={apsx}"
        );
    }
}

#[test]
fn phase_space_lighter_than_the_neutron_is_error() {
    // A_p ≤ 1 では E_max の係数 (A_p − 1)/A_p が正にならない
    let tape = Tape::parse(&phase_space("1", "1.0", "3")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 5, "{err}");
}

#[test]
fn phase_space_max_energy_follows_endf_102() {
    let space = PhaseSpace {
        total_mass: 2.99862,
        particles: 3,
    };
    let (e, awr, q) = (1.0e7, 1.9968, -2.224e6);
    let expected = (2.99862 - 1.0) / 2.99862 * (awr / (awr + 1.0) * e + q);
    assert!((space.max_energy(e, awr, q) - expected).abs() < 1e-6 * expected);
}

#[test]
fn phase_space_density_integrates_to_one() {
    // u = sin²θ と置くと √u (1 − u)^m du = 2 sin²θ cos^(2m+1)θ dθ で、θ について滑らか
    for (particles, exponent) in [(3, 0.5), (4, 2.0), (5, 3.5)] {
        let space = PhaseSpace {
            total_mass: 5.0,
            particles,
        };
        assert_eq!(space.exponent(), exponent);
        let n = 100_000;
        let h = std::f64::consts::FRAC_PI_2 / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let t = (i as f64 + 0.5) * h;
                2.0 * t.sin().powi(2) * t.cos().powf(2.0 * exponent + 1.0) * h
            })
            .sum();
        let total = space.normalization() * integral;
        assert!((total - 1.0).abs() < 1e-9, "n = {particles}: {total}");
    }
}

/// Be-9 の (n,2n) の形: 中性子が LAW=7（実験室系の角度とエネルギー）で、入射エネルギーは 2 点。
/// 2e6 eV は余弦 −1、0、1 で lin-lin の E' の表、3e6 eV は余弦 −1、1 でヒストグラムの E' の表。
/// 余弦の補間則は `mu_int`
fn angle_energy(lct: &str, mu_int: &str) -> String {
    [
        r(["4009.0", "8.93478", "0", lct, "1", "0"]),
        r(["1.0", "1.0", "0", "7", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["2.0e6", "2.0", "3.0e6", "2.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        // E = 2e6: 余弦 3 点
        r(["0.0", "2.0e6", "0", "0", "1", "3"]),
        r(["3", mu_int, "", "", "", ""]),
        r(["0.0", "-1.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["0.0", "2.0e-7", "1.0e6", "2.0e-7", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["0.0", "4.0e-7", "1.0e6", "4.0e-7", "", ""]),
        r(["0.0", "1.0", "0", "0", "1", "3"]),
        r(["3", "2", "", "", "", ""]),
        r(["0.0", "0.0", "5.0e5", "2.5e-6", "8.0e5", "0.0"]),
        // E = 3e6: 余弦 2 点、E' はヒストグラム
        r(["0.0", "3.0e6", "0", "0", "1", "2"]),
        r(["2", mu_int, "", "", "", ""]),
        r(["0.0", "-1.0", "0", "0", "1", "3"]),
        r(["3", "1", "", "", "", ""]),
        r(["0.0", "1.0e-7", "1.0e6", "3.0e-7", "2.0e6", "0.0"]),
        r(["0.0", "1.0", "0", "0", "1", "3"]),
        r(["3", "1", "", "", "", ""]),
        r(["0.0", "6.0e-7", "5.0e5", "2.0e-7", "2.0e6", "0.0"]),
    ]
    .concat()
}

#[test]
fn law7_neutrons_integrate_the_cosines_into_a_laboratory_spectrum() {
    let tape = Tape::parse(&angle_energy("1", "2")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    assert!(!mf6.has_phase_space());
    let n = &mf6.neutrons[0];
    assert_eq!(n.frame, Frame::Laboratory);
    let NeutronSpectrum::Tabulated(d) = &n.spectrum else {
        panic!("LAW=7 becomes a table: {:?}", n.spectrum);
    };
    assert_eq!(d.incident, vec![2.0e6, 3.0e6]);
    assert_eq!(
        d.laws,
        vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false
        }]
    );
    // E = 2e6: 余弦の重み（台形）は 0.5、1、0.5。g(E') = 0.5 f(−1) + f(0) + 0.5 f(1)
    let close = |t: &Tab1, x: f64, expected: f64| {
        let value = t.eval(x).unwrap();
        assert!(
            (value - expected).abs() < 1e-18,
            "g({x}) = {value}, expected {expected}"
        );
    };
    let low = &d.spectra[0];
    close(low, 0.0, 5.0e-7);
    close(low, 2.5e5, 5.0e-7 + 0.5 * 1.25e-6);
    close(low, 5.0e5, 5.0e-7 + 0.5 * 2.5e-6);
    close(low, 6.5e5, 5.0e-7 + 0.5 * 1.25e-6);
    close(low, 9.0e5, 5.0e-7);
    close(low, 1.0e6, 5.0e-7);
    // E = 3e6: 重みは 1、1。ヒストグラムの段は E' の両方の表の区切りで変わる
    let high = &d.spectra[1];
    close(high, 0.0, 7.0e-7);
    close(high, 2.5e5, 7.0e-7);
    close(high, 7.5e5, 3.0e-7);
    close(high, 1.5e6, 5.0e-7);
    // 段の値は右の値（Tab1 の不連続点の規約）で、左の値は段の手前に残る
    close(high, 5.0e5, 3.0e-7);
    close(high, 1.0e6, 5.0e-7);
    close(high, 5.0e5 * (1.0 - 1e-12), 7.0e-7);
    close(high, 1.0e6 * (1.0 - 1e-12), 3.0e-7);
    let integral: f64 = high
        .x
        .windows(2)
        .zip(high.y.windows(2))
        .map(|(x, y)| 0.5 * (x[1] - x[0]) * (y[0] + y[1]))
        .sum();
    // 和（台形で厳密）: 0.5e6 × 7e-7 + 0.5e6 × 3e-7 + 1e6 × 5e-7 = 1（ENDF-102 の式 6.2 の規格化）
    assert!((integral - 1.0).abs() < 1e-12, "{integral}");
}

#[test]
fn law7_histogram_cosines_short_of_one_are_error() {
    // 余弦のヒストグラムでは、区間 [μ_i, μ_(i+1)) の値は左の表。E = 2e6 の重みは 1、1、0 で、
    // ∫∫ f = 2e-7 × 1e6 + 4e-7 × 1e6 = 0.6（台形なら 1）。規格化に足りないのでエラーになり、その値を示す
    let tape = Tape::parse(&angle_energy("1", "1")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 7, "{err}");
    assert!(err.message.contains("integrates to 0.600000"), "{err}");
}

#[test]
fn law7_outside_the_laboratory_frame_is_error() {
    // LAW=7 は実験室系の分布なので、重心系の LCT とは食い違う
    for lct in ["2", "3"] {
        let tape = Tape::parse(&angle_energy(lct, "2")).unwrap();
        assert!(
            Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_err(),
            "LCT={lct}"
        );
    }
}

#[test]
fn law7_cosine_interpolation_other_than_histogram_or_linear_is_error() {
    let tape = Tape::parse(&angle_energy("1", "3")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 7, "{err}");
}

/// 中性子が LAW=7 で、入射エネルギー 2e6 と 3e6 のどちらも余弦が `mu0` と 1 の 2 つ、E' の表は 0〜1e6 eV で一定の `value` の節。
/// ∫∫ f = (1 − mu0) × value × 1e6
fn cosine_pair(mu0: &str, value: &str) -> String {
    cosine_pair_with(mu0, value, value, "2", false)
}

/// `cosine_pair` の一般形: 余弦 `mu0` の表は `low`、余弦 1 の表は `high`、余弦の補間則は `mu_int`。`zero_at_threshold` なら、
/// 最初の入射エネルギー（しきい値）の表はどちらも 0
fn cosine_pair_with(
    mu0: &str,
    low: &str,
    high: &str,
    mu_int: &str,
    zero_at_threshold: bool,
) -> String {
    let energy = |e: &str, low: &str, high: &str| {
        let table = |mu: &str, value: &str| {
            [
                r(["0.0", mu, "0", "0", "1", "2"]),
                r(["2", "2", "", "", "", ""]),
                r(["0.0", value, "1.0e6", value, "", ""]),
            ]
            .concat()
        };
        [
            r(["0.0", e, "0", "0", "1", "2"]),
            r(["2", mu_int, "", "", "", ""]),
            table(mu0, low),
            table("1.0", high),
        ]
        .concat()
    };
    [
        r(["4009.0", "8.93478", "0", "1", "1", "0"]),
        r(["1.0", "1.0", "0", "7", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["2.0e6", "2.0", "3.0e6", "2.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        if zero_at_threshold {
            energy("2.0e6", "0.0", "0.0")
        } else {
            energy("2.0e6", low, high)
        },
        energy("3.0e6", low, high),
    ]
    .concat()
}

#[test]
fn law7_forward_cone_normalized_to_one_passes() {
    // ENDF-102 の 6.2.9 は、余弦が出射の粒子に開いた範囲（前方の円錐 μ_min〜1 でもよい）を覆うことを求める。
    // 0〜1 の円錐で f = 1e-6 なら ∫∫ f = 1 で、式 6.2 の規格化を満たす
    let tape = Tape::parse(&cosine_pair("0.0", "1.0e-6")).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    let NeutronSpectrum::Tabulated(d) = &mf6.neutrons[0].spectrum else {
        panic!("LAW=7 becomes a table");
    };
    assert!((d.spectra[0].eval(5.0e5).unwrap() - 1.0e-6).abs() < 1e-18);
}

#[test]
fn law7_histogram_cosines_of_a_normalized_distribution_take_the_left_table() {
    // 余弦のヒストグラムで 0〜1 の重みは左の表に 1、右の表に 0。左が 1e-6 なら ∫∫ f = 1 で通り、g は左の表になる
    // （lin-lin なら (1e-6 + 3e-6)/2 で、規格化も 2 になる）
    let tape = Tape::parse(&cosine_pair_with("0.0", "1.0e-6", "3.0e-6", "1", false)).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    let NeutronSpectrum::Tabulated(d) = &mf6.neutrons[0].spectrum else {
        panic!("LAW=7 becomes a table");
    };
    for spectrum in &d.spectra {
        assert!((spectrum.eval(5.0e5).unwrap() - 1.0e-6).abs() < 1e-18);
    }
}

#[test]
fn law7_zero_distribution_at_threshold_passes() {
    // しきい値の入射エネルギーでは断面積が 0 で、分布を 0 で与える評価がある。最初の入射エネルギーの 0 の分布は通す
    let tape = Tape::parse(&cosine_pair_with("0.0", "1.0e-6", "1.0e-6", "2", true)).unwrap();
    let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
    let NeutronSpectrum::Tabulated(d) = &mf6.neutrons[0].spectrum else {
        panic!("LAW=7 becomes a table");
    };
    assert_eq!(d.spectra[0].eval(5.0e5), Some(0.0));
}

#[test]
fn law7_distribution_short_of_one_is_error() {
    // 開いた範囲の一部を欠いた分布は和が 1 に足りず、中性子を黙って失うので、入射エネルギーと積分の値を示して止める
    let tape = Tape::parse(&cosine_pair("0.5", "1.0e-6")).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 7, "{err}");
    assert!(err.message.contains("integrates to 0.500000"), "{err}");
    assert!(err.message.contains("2e6"), "{err}");
}

#[test]
fn law7_cosines_outside_minus_one_to_one_are_error() {
    let text = angle_energy("1", "2").replacen("       -1.0", "       -1.5", 1);
    let tape = Tape::parse(&text).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 7, "{err}");
}

/// 中性子が LAW=7 で、入射エネルギー 2e6 と 3e6 のどちらも余弦が 1 つだけ（μ = 0）の節
fn single_cosine() -> String {
    let energy = |e: &str| {
        [
            r(["0.0", e, "0", "0", "1", "1"]),
            r(["1", "2", "", "", "", ""]),
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "5.0e-7", "1.0e6", "5.0e-7", "", ""]),
        ]
        .concat()
    };
    [
        r(["4009.0", "8.93478", "0", "1", "1", "0"]),
        r(["1.0", "1.0", "0", "7", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["2.0e6", "2.0", "3.0e6", "2.0", "", ""]),
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        energy("2.0e6"),
        energy("3.0e6"),
    ]
    .concat()
}

#[test]
fn law7_single_cosine_is_error() {
    // 余弦が 1 つでは積分の幅が 0 で、中性子を黙って失う
    let tape = Tape::parse(&single_cosine()).unwrap();
    let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
    assert_eq!(err.line, 7, "{err}");
}
