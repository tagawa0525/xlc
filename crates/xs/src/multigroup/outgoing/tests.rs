use super::*;
use crate::endf::IncidentLaw;

fn assert_close(actual: &[f64], expected: &[f64]) {
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() <= 1e-14, "{actual:?} vs {expected:?}");
    }
}

#[test]
fn evaporation_spectrum_is_cut_at_e_minus_u() {
    // E' ≤ E - U = 2 MeV。群ごとの割合は P(2, x) の差（x = E'/θ）
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let law = EnergyLaw::Evaporation {
        restriction: 1.0e6,
        theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 2.0e7], vec![1.0e5, 1.0e5])).unwrap(),
    };
    let grouped = GroupedLaw::new(&groups, &law);
    let mut out = vec![0.0; 3];
    grouped.accumulate(3.0e6, 2.0, &mut out);
    let p = |x: f64| gamma_p_2(x) / gamma_p_2(20.0);
    assert_close(
        &out,
        &[
            2.0 * (1.0 - p(10.0)),
            2.0 * (p(10.0) - p(1.0)),
            2.0 * p(1.0),
        ],
    );
    // E ≤ U では中性子を出さない
    let mut out = vec![0.0; 3];
    grouped.accumulate(1.0e6, 1.0, &mut out);
    assert_eq!(out, vec![0.0; 3]);
    // E - U が群境界を横切る点で、群ごとの割合の傾きが変わる
    let breakpoints = grouped.breakpoints();
    for e in [1.0e6, 1.1e6, 2.0e6] {
        assert!(breakpoints.contains(&e), "{e}: {breakpoints:?}");
    }
}

#[test]
fn maxwellian_interpolates_theta() {
    // θ(2 MeV) = 200 keV。E - U が群構造の上限を超える分は最上位の群に入る
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
    let law = EnergyLaw::Maxwellian {
        restriction: -3.0e7,
        theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 3.0e6], vec![1.0e5, 3.0e5])).unwrap(),
    };
    let grouped = GroupedLaw::new(&groups, &law);
    let mut out = vec![0.0; 3];
    grouped.accumulate(2.0e6, 1.0, &mut out);
    let p = |x: f64| gamma_p_3_2(x) / gamma_p_3_2(160.0);
    assert_close(&out, &[1.0 - p(5.0), p(5.0) - p(0.5), p(0.5)]);
    assert!(grouped.breakpoints().contains(&3.0e6));
}

#[test]
fn integrates_piecewise_linear_density_by_group() {
    let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
    let density = Tab1::linear(vec![1.0, 3.0], vec![0.5, 0.5]);
    assert_eq!(outgoing_group_integrals(&groups, &density), vec![0.5, 0.5]);
    // 群構造の外に出る分は端の群に入れる（NJOY と同じ。中性子を失わない）
    let density = Tab1::linear(vec![0.0, 3.0, 5.0], vec![0.2, 0.2, 0.2]);
    let r = outgoing_group_integrals(&groups, &density);
    assert!(
        (r[0] - 0.6).abs() < 1e-15 && (r[1] - 0.4).abs() < 1e-15,
        "{r:?}"
    );
}

#[test]
fn linear_density_split_inside_a_segment() {
    let groups = GroupStructure::new(vec![2.0, 1.0, 0.5]).unwrap();
    let density = Tab1::linear(vec![0.0, 2.0], vec![0.0, 1.0]);
    let r = outgoing_group_integrals(&groups, &density);
    // ∫ E/2 dE。0.5 eV 未満の分（0.0625）は最下位の群に入る
    assert!(
        (r[0] - 0.75).abs() < 1e-15 && (r[1] - 0.25).abs() < 1e-15,
        "{r:?}"
    );
}

#[test]
fn center_of_mass_energy_spreads_in_laboratory() {
    // 重心系で E'cm にほぼ単色の中性子は、実験室系で
    // [(√E'cm - √Ec)², (√E'cm + √Ec)²] に一様に分布する（Ec = E/(A+1)²）
    let (e, awr) = (1.6e7, 15.0);
    let ec = e / (awr + 1.0_f64).powi(2);
    assert_eq!(ec, 62500.0);
    let distributions = OutgoingDistributions {
        breaks: vec![1],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false,
        }],
        incident: vec![e],
        spectra: vec![Tab1 {
            breaks: vec![2],
            laws: vec![Interpolation::Histogram],
            x: vec![1.0e6 - 1.0, 1.0e6 + 1.0],
            y: vec![0.5, 0.5],
        }],
    };
    let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 1.0e5]).unwrap();
    let grouped = GroupedOutgoing::center_of_mass(&groups, &distributions, awr);
    let mut out = vec![0.0; 2];
    grouped.accumulate(e, 1.0, &mut out);
    // 実験室系の範囲は [562500, 1562500]
    assert!(
        (out[0] - 0.5625).abs() < 1e-4 && (out[1] - 0.4375).abs() < 1e-4,
        "{out:?}"
    );
}

/// 入射エネルギー 10 と 20 で、[0, 1] と [0, 3] の一様分布
fn widening(law: IncidentLaw) -> OutgoingDistributions {
    OutgoingDistributions {
        breaks: vec![2],
        laws: vec![law],
        incident: vec![10.0, 20.0],
        spectra: vec![
            Tab1::linear(vec![0.0, 1.0], vec![1.0, 1.0]),
            Tab1::linear(vec![0.0, 3.0], vec![1.0 / 3.0, 1.0 / 3.0]),
        ],
    }
}

#[test]
fn unit_base_keeps_the_shape_between_incident_energies() {
    // 単位基底では、中間の入射エネルギーで範囲が [0, 2] に伸びた一様分布になる
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let unit_base = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: true,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(unit_base));
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [0.0, 0.5, 0.5];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-6, "{out:?}");
    }
    // ふつうの線形補間なら、2つの分布の平均（[0, 1] に 2/3、[1, 3] に 1/3）
    let cartesian = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: false,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(cartesian));
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-6, "{out:?}");
    }
}

#[test]
fn unit_base_keeps_histogram_steps() {
    // 出射分布がヒストグラム（段差がある）でも、確率の合計と形を保つ。
    // [0, 0.5] に密度 2 と [0, 1.5] に密度 2/3 の中間は、[0, 1] に密度 1
    let step = |width: f64| Tab1 {
        breaks: vec![3],
        laws: vec![Interpolation::Histogram],
        x: vec![0.0, 0.5 * width, width],
        y: vec![2.0 / width, 0.0, 0.0],
    };
    let distributions = OutgoingDistributions {
        breaks: vec![2],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true,
        }],
        incident: vec![10.0, 20.0],
        spectra: vec![step(1.0), step(3.0)],
    };
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let grouped = GroupedOutgoing::new(&groups, &distributions);
    let mut out = vec![0.0; 3];
    grouped.accumulate(15.0, 1.0, &mut out);
    let expected = [0.0, 0.0, 1.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-9, "{out:?}");
    }
}

#[test]
fn unit_base_in_the_center_of_mass_matches_the_laboratory_for_a_heavy_target() {
    // 標的が十分重ければ、重心系と実験室系はほぼ同じ
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let law = IncidentLaw {
        law: Interpolation::LinLin,
        unit_base: true,
    };
    let lab = GroupedOutgoing::new(&groups, &widening(law));
    let cm = GroupedOutgoing::center_of_mass(&groups, &widening(law), 1.0e9);
    let (mut a, mut b) = (vec![0.0; 3], vec![0.0; 3]);
    lab.accumulate(15.0, 1.0, &mut a);
    cm.accumulate(15.0, 1.0, &mut b);
    for (x, y) in a.iter().zip(&b) {
        assert!((x - y).abs() < 1e-4, "{a:?} vs {b:?}");
    }
}

#[test]
fn incident_interpolation_follows_its_law() {
    // 線形・対数（LinLog）なら、入射エネルギーの対数の中点で2つの分布の平均
    let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
    let law = IncidentLaw {
        law: Interpolation::LinLog,
        unit_base: false,
    };
    let grouped = GroupedOutgoing::new(&groups, &widening(law));
    let mut out = vec![0.0; 3];
    grouped.accumulate((10.0_f64 * 20.0).sqrt(), 1.0, &mut out);
    let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
    for (o, e) in out.iter().zip(expected) {
        assert!((o - e).abs() < 1e-12, "{out:?}");
    }
}

#[test]
fn interpolates_between_incident_energies() {
    let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
    let distributions = OutgoingDistributions {
        breaks: vec![2],
        laws: vec![IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false,
        }],
        incident: vec![10.0, 20.0],
        spectra: vec![
            Tab1::linear(vec![2.0, 4.0], vec![0.5, 0.5]),
            Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]),
        ],
    };
    let grouped = GroupedOutgoing::new(&groups, &distributions);
    let mut out = vec![0.0; 2];
    grouped.accumulate(12.5, 2.0, &mut out);
    assert_eq!(out, vec![1.5, 0.5]);
    // 範囲外は端の分布
    let mut out = vec![0.0; 2];
    grouped.accumulate(30.0, 1.0, &mut out);
    assert_eq!(out, vec![0.0, 1.0]);
}
