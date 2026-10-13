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
    let d = &n.distributions;
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
        Mf6::parse(&tape.section(825, 6, 16).unwrap())
            .map(|mf6| mf6.neutrons[0].distributions.laws[0])
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
    let tape = Tape::parse(&section("2", "7", "0")).unwrap();
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
