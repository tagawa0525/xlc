use super::*;
use crate::endf::tape::tests::rec;
use crate::endf::{Interpolation, Tape};

fn mf2(lines: &[[&str; 6]]) -> Result<Mf2, Error> {
    let text: String = lines.iter().map(|f| rec(*f, 9237, 2, 151)).collect();
    let tape = Tape::parse(&text).unwrap();
    Mf2::parse(&tape.section(9237, 2, 151).unwrap())
}

const HEAD: [&str; 6] = ["9.223800+4", "2.360058+2", "0", "0", "1", "0"];

/// R-matrix limited（LRF=7）の範囲。`control` は IFG KRM NJS KRL の CONT
fn rml_lines(control: [&'static str; 6]) -> Vec<[&'static str; 6]> {
    vec![
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "1.9+5", "1", "7", "0", "0"],
        control,
        // 粒子対: MA MB ZA ZB IA IB / Q PNT SHF MT PA PB（光子と捕獲、中性子と標的）
        ["0.0", "0.0", "2", "0", "24", "4"],
        ["0.0", "5.743560+1", "0.0", "2.6+1", "1.0", "0.0"],
        ["0.0", "0.0", "0.0", "1.02+2", "0.0", "0.0"],
        ["1.0", "5.644629+1", "0.0", "2.6+1", "0.5", "-0.5"],
        ["0.0", "1.0", "0.0", "2.0", "0.0", "0.0"],
        // spin group: AJ PJ KBK KPS 6*NCH NCH と、チャネルごとに PPI L SCH BND APE APT
        ["1.0", "0.0", "0", "0", "12", "2"],
        ["1.0", "0.0", "0.0", "0.0", "0.0", "0.0"],
        ["2.0", "1.0", "1.0", "0.0", "0.63", "0.64"],
        // 共鳴: 0 0 0 NRS 6*NX NX と、共鳴ごとに ER とチャネルの幅
        ["0.0", "0.0", "0", "2", "12", "2"],
        ["-2.33+3", "1.0", "2.0", "0.0", "0.0", "0.0"],
        ["1.5+3", "0.5", "-3.0", "0.0", "0.0", "0.0"],
    ]
}

#[test]
fn reads_r_matrix_limited_range() {
    let data = mf2(&rml_lines(["0.0", "0.0", "0", "3", "1", "0"])).unwrap();
    let RangeKind::RMatrixLimited(rml) = &data.ranges[0].kind else {
        panic!("expected R-matrix limited");
    };
    assert_eq!(
        rml.pairs,
        vec![
            ParticlePair {
                ma: 0.0,
                mb: 57.4356,
                za: 0.0,
                zb: 26.0,
                ia: 1.0,
                ib: 0.0,
                q: 0.0,
                pnt: 0,
                shf: 0,
                mt: 102,
            },
            ParticlePair {
                ma: 1.0,
                mb: 56.44629,
                za: 0.0,
                zb: 26.0,
                ia: 0.5,
                ib: -0.5,
                q: 0.0,
                pnt: 1,
                shf: 0,
                mt: 2,
            },
        ]
    );
    assert_eq!(
        rml.spin_groups,
        vec![SpinGroup {
            aj: 1.0,
            pj: 0.0,
            channels: vec![
                RmlChannel {
                    pair: 0,
                    l: OrbitalL::S,
                    sch: 0.0,
                    bnd: 0.0,
                    ape: 0.0,
                    apt: 0.0,
                },
                RmlChannel {
                    pair: 1,
                    l: OrbitalL::P,
                    sch: 1.0,
                    bnd: 0.0,
                    ape: 0.63,
                    apt: 0.64,
                },
            ],
            resonances: vec![
                RmlResonance {
                    er: -2330.0,
                    widths: vec![1.0, 2.0],
                },
                RmlResonance {
                    er: 1500.0,
                    widths: vec![0.5, -3.0],
                },
            ],
        }]
    );
}

#[test]
fn rejects_r_matrix_limited_other_than_reich_moore() {
    let error = mf2(&rml_lines(["0.0", "0.0", "0", "2", "1", "0"]))
        .unwrap_err()
        .to_string();
    assert!(error.contains("KRM=2"), "{error}");
}

#[test]
fn rejects_r_matrix_limited_with_reduced_amplitudes() {
    let error = mf2(&rml_lines(["0.0", "0.0", "1", "3", "1", "0"]))
        .unwrap_err()
        .to_string();
    assert!(error.contains("IFG=1"), "{error}");
}

#[test]
fn reads_reich_moore_and_unresolved_ranges() {
    let data = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "2", "0"],
        // RRR: EL EH LRU LRF NRO NAPS
        ["1.0-5", "2.0+4", "1", "3", "0", "0"],
        ["0.0", "0.0", "0", "0", "1", "1"],
        // l=0: AWRI APL L 0 6*NRS NRS
        ["2.360060+2", "9.480000-1", "0", "0", "12", "2"],
        ["-4.0", "0.5", "2.0-2", "2.3-2", "0.0", "1.0-6"],
        ["6.674", "-0.5", "1.5-3", "2.3-2", "0.0", "0.0"],
        // URR: SPI AP LSSF 0 NLS 0
        ["2.0+4", "1.49+5", "2", "2", "0", "0"],
        ["0.0", "0.948", "1", "0", "1", "0"],
        // l=0: AWRI 0 L 0 NJS 0
        ["2.360060+2", "0.0", "0", "0", "1", "0"],
        // J: AJ 0 INT 0 6*NE+6 NE、続いて 0 0 AMUX AMUN AMUG AMUF と ES D GX GNO GG GF
        ["0.5", "0.0", "5", "0", "18", "2"],
        ["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"],
        ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
        ["1.49+5", "18.0", "0.0", "1.8-3", "2.3-2", "0.0"],
    ])
    .unwrap();
    assert_eq!(data.za, 92238.0);
    assert_eq!(data.ranges.len(), 2);
    assert_eq!((data.ranges[0].el, data.ranges[0].eh), (1.0e-5, 2.0e4));
    let RangeKind::ReichMoore(rm) = &data.ranges[0].kind else {
        panic!("expected Reich-Moore");
    };
    assert_eq!((rm.spi, rm.ap, rm.naps), (0.0, 0.0, 0));
    assert_eq!(rm.l_values.len(), 1);
    assert_eq!(rm.l_values[0].apl, 0.948);
    assert_eq!(
        rm.l_values[0].resonances[1],
        RmResonance {
            er: 6.674,
            aj: -0.5,
            gn: 1.5e-3,
            gg: 2.3e-2,
            gfa: 0.0,
            gfb: 0.0
        }
    );
    let RangeKind::Unresolved(urr) = &data.ranges[1].kind else {
        panic!("expected unresolved");
    };
    assert_eq!((urr.lssf, urr.spi, urr.ap, urr.naps), (1, 0.0, 0.948, 0));
    assert_eq!(urr.ap_energy, None);
    assert_eq!(urr.l_values.len(), 1);
    assert_eq!(
        (urr.l_values[0].awri, urr.l_values[0].l),
        (236.006, OrbitalL::S)
    );
    let j = &urr.l_values[0].j_values[0];
    assert_eq!((j.aj, j.int), (0.5, Interpolation::LogLog));
    assert_eq!(
        (j.amux, j.amun, j.amuf),
        (WidthDof::Fixed, WidthDof::One, WidthDof::Fixed)
    );
    assert_eq!(
        j.points[1],
        UrrPoint {
            e: 1.49e5,
            d: 18.0,
            gx: 0.0,
            gno: 1.8e-3,
            gg: 2.3e-2,
            gf: 0.0
        }
    );
}

#[test]
fn reads_energy_dependent_radius_of_unresolved_range() {
    let data = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["2.0+4", "1.49+5", "2", "2", "1", "0"],
        // AP(E) の TAB1: 0 0 0 0 NR NP、補間の区切り、(E, AP) の組
        ["0.0", "0.0", "0", "0", "1", "2"],
        ["2", "2", "", "", "", ""],
        ["2.0+4", "0.95", "1.49+5", "0.90", "", ""],
        ["0.0", "0.948", "1", "0", "1", "0"],
        ["2.360060+2", "0.0", "0", "0", "1", "0"],
        ["0.5", "0.0", "2", "0", "12", "1"],
        ["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"],
        ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
    ])
    .unwrap();
    let RangeKind::Unresolved(urr) = &data.ranges[0].kind else {
        panic!("expected unresolved");
    };
    let radius = urr.ap_energy.as_ref().unwrap();
    assert_eq!(
        (radius.x.clone(), radius.y.clone()),
        (vec![2.0e4, 1.49e5], vec![0.95, 0.90])
    );
    assert_eq!(urr.ap, 0.948);
}

/// 非分離共鳴 LRF=1 の範囲のあとに LRU=0 の範囲を置いて読み、LRF=1 の範囲と、
/// そのあとの範囲を正しく読めたこと（レコードの位置がずれないこと）を確かめる。
fn read_lrf1(lfw: &str, nro: &str, parameters: &[[&str; 6]]) -> Unresolved {
    let mut lines = vec![
        HEAD,
        ["9.423900+4", "1.0", "0", lfw, "2", "0"],
        ["5.7+3", "4.0+4", "2", "1", nro, "0"],
    ];
    lines.extend_from_slice(parameters);
    lines.push(["4.0+4", "2.0+7", "0", "0", "0", "0"]);
    lines.push(["0.0", "0.5", "0", "0", "0", "0"]);
    let data = mf2(&lines).unwrap();
    assert_eq!(data.ranges[1].kind, RangeKind::ScatteringRadius { ap: 0.5 });
    let RangeKind::Unresolved(urr) = &data.ranges[0].kind else {
        panic!("expected unresolved");
    };
    urr.clone()
}

/// LRF=1 の (l, J) の系列。競合反応の幅はなく、補間は NJOY と同じく線形
fn lrf1_sequence(amuf: WidthDof, points: Vec<UrrPoint>) -> UrrLValue {
    UrrLValue {
        awri: 237.9916,
        l: OrbitalL::S,
        j_values: vec![UrrJValue {
            aj: 0.5,
            int: Interpolation::LinLin,
            amux: WidthDof::Fixed,
            amun: WidthDof::One,
            amuf,
            points,
        }],
    }
}

fn lrf1_point(e: f64, gf: f64) -> UrrPoint {
    UrrPoint {
        e,
        d: 13.1,
        gx: 0.0,
        gno: 1.5e-3,
        gg: 3.1e-2,
        gf,
    }
}

#[test]
fn reads_unresolved_lrf1_without_fission_widths() {
    let urr = read_lrf1(
        "0",
        "0",
        &[
            // SPI AP LSSF 0 NLS 0、l ごとに AWRI 0 L 0 6*NJS NJS と D AJ AMUN GNO GG GF
            ["0.0", "0.888", "1", "0", "1", "0"],
            ["2.379916+2", "0.0", "0", "0", "6", "1"],
            ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
        ],
    );
    // エネルギーによらないパラメータは、範囲の両端で同じ値とする
    assert_eq!(
        urr,
        Unresolved {
            lssf: 1,
            spi: 0.0,
            ap: 0.888,
            ap_energy: None,
            naps: 0,
            l_values: vec![lrf1_sequence(
                WidthDof::Fixed,
                vec![lrf1_point(5.7e3, 0.0), lrf1_point(4.0e4, 0.0)],
            )],
        }
    );
}

#[test]
fn reads_unresolved_lrf1_with_fission_widths() {
    let urr = read_lrf1(
        "1",
        "0",
        &[
            // SPI AP LSSF 0 NE NLS と ES、l ごとに AWRI 0 L 0 NJS 0、
            // J ごとに 0 0 L MUF NE+6 0 と D AJ AMUN GNO GG 0 GF...
            ["0.0", "0.888", "1", "0", "2", "1"],
            ["5.7+3", "4.0+4", "", "", "", ""],
            ["2.379916+2", "0.0", "0", "0", "1", "0"],
            ["0.0", "0.0", "0", "1", "8", "0"],
            ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
            ["1.2-3", "1.4-3", "", "", "", ""],
        ],
    );
    // 核分裂幅のエネルギーの点に、ほかのパラメータを並べる
    assert_eq!(
        urr.l_values,
        vec![lrf1_sequence(
            WidthDof::One,
            vec![lrf1_point(5.7e3, 1.2e-3), lrf1_point(4.0e4, 1.4e-3)],
        )]
    );
}

/// LRF=1 の範囲 [5.7 keV, 40 keV] を読んだ結果のエラー
fn lrf1_error(lfw: &str, parameters: &[[&str; 6]]) -> String {
    let mut lines = vec![
        HEAD,
        ["9.423900+4", "1.0", "0", lfw, "1", "0"],
        ["5.7+3", "4.0+4", "2", "1", "0", "0"],
    ];
    lines.extend_from_slice(parameters);
    mf2(&lines).unwrap_err().to_string()
}

#[test]
fn rejects_invalid_lfw() {
    let error = lrf1_error(
        "2",
        &[
            ["0.0", "0.888", "1", "0", "1", "0"],
            ["2.379916+2", "0.0", "0", "0", "6", "1"],
            ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
        ],
    );
    assert!(error.contains("LFW=2"), "{error}");
}

#[test]
fn rejects_fission_width_energies_short_of_the_range() {
    // ES が範囲の上端 40 keV に届かない
    let error = lrf1_error(
        "1",
        &[
            ["0.0", "0.888", "1", "0", "2", "1"],
            ["5.7+3", "3.0+4", "", "", "", ""],
            ["2.379916+2", "0.0", "0", "0", "1", "0"],
            ["0.0", "0.0", "0", "1", "8", "0"],
            ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
            ["1.2-3", "1.4-3", "", "", "", ""],
        ],
    );
    assert!(error.contains("cover"), "{error}");
}

#[test]
fn reads_unresolved_lrf1_with_energy_dependent_radius() {
    let urr = read_lrf1(
        "0",
        "1",
        &[
            // AP(E) の TAB1 が先に来る
            ["0.0", "0.0", "0", "0", "1", "2"],
            ["2", "2", "", "", "", ""],
            ["5.7+3", "0.89", "4.0+4", "0.88", "", ""],
            ["0.0", "0.888", "1", "0", "1", "0"],
            ["2.379916+2", "0.0", "0", "0", "6", "1"],
            ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
        ],
    );
    let radius = urr.ap_energy.unwrap();
    assert_eq!((radius.x, radius.y), (vec![5.7e3, 4.0e4], vec![0.89, 0.88]));
}

#[test]
fn rounds_degrees_of_freedom_like_njoy() {
    // NJOY は自由度を最も近い整数に丸めて使う（Kr-83 の AMUX=1.0123 など）
    let data = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["2.0+4", "1.49+5", "2", "2", "0", "0"],
        ["0.0", "0.948", "1", "0", "1", "0"],
        ["2.360060+2", "0.0", "0", "0", "1", "0"],
        ["0.5", "0.0", "2", "0", "18", "2"],
        ["0.0", "0.0", "1.0123", "0.9", "0.0", "2.4"],
        // 同じエネルギーの繰り返しはパラメータの段差を表す（Th-229 など）
        ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
        ["2.0+4", "18.0", "0.0", "1.8-3", "2.3-2", "0.0"],
    ])
    .unwrap();
    let RangeKind::Unresolved(urr) = &data.ranges[0].kind else {
        panic!("expected unresolved");
    };
    let j = &urr.l_values[0].j_values[0];
    assert_eq!(
        (j.amux, j.amun, j.amuf),
        (WidthDof::One, WidthDof::One, WidthDof::Two)
    );
    assert_eq!(j.points.len(), 2);
}

#[test]
fn rejects_degrees_of_freedom_beyond_quadrature() {
    let result = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["2.0+4", "1.49+5", "2", "2", "0", "0"],
        ["0.0", "0.948", "1", "0", "1", "0"],
        ["2.360060+2", "0.0", "0", "0", "1", "0"],
        ["0.5", "0.0", "2", "0", "12", "1"],
        ["0.0", "0.0", "0.0", "5.0", "0.0", "0.0"],
        ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
    ]);
    assert!(result.is_err());
}

#[test]
fn reads_mlbw_range() {
    let data = mf2(&[
        ["4.009000+4", "8.913240+1", "0", "0", "1", "0"],
        ["4.009000+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "2.0+5", "1", "2", "0", "0"],
        ["0.0", "7.16-1", "0", "0", "1", "0"],
        // l=1: AWRI QX L LRX 6*NRS NRS、値は ER AJ GT GN GG GF
        ["8.913240+1", "0.0", "1", "0", "6", "1"],
        ["2.4e3", "1.5", "0.25", "0.2", "0.05", "0.0"],
    ])
    .unwrap();
    let RangeKind::Mlbw(bw) = &data.ranges[0].kind else {
        panic!("expected MLBW");
    };
    assert_eq!((bw.spi, bw.ap, bw.naps), (0.0, 0.716, 0));
    assert_eq!(bw.l_values[0].l, OrbitalL::P);
    assert_eq!(
        bw.l_values[0].resonances[0],
        BwResonance {
            er: 2.4e3,
            aj: 1.5,
            gn: 0.2,
            gg: 0.05,
            gf: 0.0
        }
    );
}

/// MLBW の範囲（1e-5〜50 eV、NRO=1）。`naps` は範囲の見出しの NAPS、`radius` は AP(E) の TAB1 の (E, AP) の 2 点
fn mlbw_with_radius(naps: &str, radius: [&str; 4]) -> Result<Mf2, Error> {
    mf2(&[
        ["6.014900+4", "1.476409+2", "0", "0", "1", "0"],
        ["6.014900+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "5.0+1", "1", "2", "1", naps],
        // AP(E) の TAB1: 0 0 0 0 NR NP、補間の区切り、(E, AP) の組
        ["0.0", "0.0", "0", "0", "1", "2"],
        ["2", "2", "", "", "", ""],
        [radius[0], radius[1], radius[2], radius[3], "", ""],
        ["3.5", "0.0", "0", "0", "1", "0"],
        ["1.476409+2", "0.0", "0", "0", "6", "1"],
        ["2.0", "3.0", "0.1", "0.05", "0.05", "0.0"],
    ])
}

#[test]
fn reads_energy_dependent_radius_of_mlbw_range() {
    let data = mlbw_with_radius("0", ["1.0-5", "0.70", "5.0+1", "0.60"]).unwrap();
    let RangeKind::Mlbw(bw) = &data.ranges[0].kind else {
        panic!("expected MLBW");
    };
    assert_eq!((bw.spi, bw.naps), (3.5, 0));
    assert_eq!(bw.l_values[0].resonances[0].er, 2.0);
    let radius = bw.ap_energy.as_ref().unwrap();
    assert_eq!(
        (radius.x.clone(), radius.y.clone()),
        (vec![1.0e-5, 50.0], vec![0.70, 0.60])
    );
}

#[test]
fn rejects_mlbw_radius_table_short_of_the_range() {
    // 範囲の外へ外挿しない
    for radius in [
        ["1.0-5", "0.70", "1.0+1", "0.60"],
        ["1.0", "0.70", "5.0+1", "0.60"],
    ] {
        let err = mlbw_with_radius("0", radius).unwrap_err();
        assert_eq!(err.line, 4, "{radius:?}");
    }
}

#[test]
fn rejects_mlbw_radius_table_with_a_step() {
    // 不連続（x の重複）では、範囲の中の値と、範囲の上端での内側の値が決まらない
    let err = mf2(&[
        ["6.014900+4", "1.476409+2", "0", "0", "1", "0"],
        ["6.014900+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "5.0+1", "1", "2", "1", "0"],
        ["0.0", "0.0", "0", "0", "1", "4"],
        ["4", "2", "", "", "", ""],
        ["1.0-5", "0.70", "2.5+1", "0.65", "2.5+1", "0.62"],
        ["5.0+1", "0.60", "", "", "", ""],
        ["3.5", "0.0", "0", "0", "1", "0"],
        ["1.476409+2", "0.0", "0", "0", "6", "1"],
        ["2.0", "3.0", "0.1", "0.05", "0.05", "0.0"],
    ])
    .unwrap_err();
    assert_eq!(err.line, 4);
}

#[test]
fn rejects_energy_dependent_mlbw_radius_unless_naps_is_zero() {
    for naps in ["1", "2"] {
        let err = mlbw_with_radius(naps, ["1.0-5", "0.70", "5.0+1", "0.60"]).unwrap_err();
        assert_eq!(err.line, 3, "NAPS={naps}");
    }
}

#[test]
fn competitive_width_is_error() {
    let err = mf2(&[
        ["4.009000+4", "8.913240+1", "0", "0", "1", "0"],
        ["4.009000+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "2.0+5", "1", "2", "0", "0"],
        ["0.0", "7.16-1", "0", "0", "1", "0"],
        ["8.913240+1", "-1.0+5", "1", "1", "6", "1"],
        ["2.4e3", "1.5", "0.25", "0.2", "0.05", "0.0"],
    ])
    .unwrap_err();
    assert_eq!(err.line, 5);
}

#[test]
fn reads_scattering_radius_only_range() {
    let data = mf2(&[
        ["1.001000+3", "9.991673-1", "0", "0", "1", "0"],
        ["1.001000+3", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "1.0+5", "0", "0", "0", "0"],
        ["0.5", "1.276553+0", "0", "0", "0", "0"],
    ])
    .unwrap();
    assert_eq!(
        data.ranges[0].kind,
        RangeKind::ScatteringRadius { ap: 1.276553 }
    );
}

#[test]
fn unsupported_formalism_is_error() {
    // LRF=1 (SLBW)
    let err = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "2.0+4", "1", "1", "0", "0"],
    ])
    .unwrap_err();
    assert_eq!(err.line, 3);
}

#[test]
fn energy_dependent_radius_in_reich_moore_is_error() {
    let err = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "2.0+4", "1", "3", "1", "0"],
    ])
    .unwrap_err();
    assert_eq!(err.line, 3);
}

#[test]
fn zero_resonance_energy_is_error() {
    // Er=0 では透過因子が 0 になり、中性子幅を換算できない
    let err = mf2(&[
        HEAD,
        ["9.223800+4", "1.0", "0", "0", "1", "0"],
        ["1.0-5", "2.0+4", "1", "3", "0", "0"],
        ["0.0", "0.0", "0", "0", "1", "1"],
        ["2.360060+2", "9.480000-1", "0", "0", "6", "1"],
        ["0.0", "0.5", "2.0-2", "2.3-2", "0.0", "0.0"],
    ])
    .unwrap_err();
    assert_eq!(err.line, 6);
}

#[test]
fn huge_range_count_is_error_not_abort() {
    let err = mf2(&[HEAD, ["9.223800+4", "1.0", "0", "0", "99999999999", "0"]]).unwrap_err();
    assert_eq!(err.line, 3);
}

#[test]
fn multiple_isotopes_are_error() {
    let err = mf2(&[["9.223800+4", "2.360058+2", "0", "0", "2", "0"]]).unwrap_err();
    assert_eq!(err.line, 1);
}
