use super::*;

/// 6 つの欄と制御欄の 1 行
fn rec(fields: [&str; 6], mat: i32, mf: u32, mt: u32) -> String {
    let data: String = fields.iter().map(|f| format!("{f:>11}")).collect();
    format!("{data}{mat:4}{mf:2}{mt:3}\n")
}

const BLANK: [&str; 6] = ["", "", "", "", "", ""];
const TAPE_ID: &str =
    " test tape                                                           1 0  0\n";

/// MF1/MT451 の最初の 2 行と、セクションの終わり
fn header(mat: i32, za: &str, excitation: &str, liso: &str) -> String {
    [
        rec([za, "2.399801+2", "-1", "0", "0", "0"], mat, 1, 451),
        rec([excitation, "1.0", "2", liso, "0", "6"], mat, 1, 451),
        rec(BLANK, mat, 1, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat()
}

fn am242m_decay() -> String {
    let mat = 3580;
    [
        TAPE_ID.to_string(),
        header(mat, "9.524200+4", "4.863000+4", "1"),
        rec(
            ["9.524200+4", "2.399801+2", "2", "1", "0", "4"],
            mat,
            8,
            457,
        ),
        rec(
            ["4.449622+9", "6.311520+7", "0", "0", "6", "0"],
            mat,
            8,
            457,
        ),
        rec(
            [
                "4.008322+4",
                "4.0+3",
                "5.373671+3",
                "3.6+2",
                "2.441003+4",
                "6.1+2",
            ],
            mat,
            8,
            457,
        ),
        rec(["5.0", "-1.0", "0", "0", "18", "3"], mat, 8, 457),
        rec(
            ["4.0", "0.0", "5.637120+6", "2.5+2", "4.590000-3", "1.2-4"],
            mat,
            8,
            457,
        ),
        rec(
            ["3.0", "0.0", "4.860000+4", "5.0+1", "9.954100-1", "1.2-4"],
            mat,
            8,
            457,
        ),
        rec(
            ["6.0", "0.0", "1.884000+8", "3.7+6", "1.60000-10", "6.0-11"],
            mat,
            8,
            457,
        ),
        // 続くスペクトルは読まない
        rec(["0.0", "0.0", "0", "0", "6", "4"], mat, 8, 457),
        rec(BLANK, mat, 8, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat()
}

#[test]
fn reads_the_material_of_a_metastable_state() {
    let tape = Tape::parse(&am242m_decay()).unwrap();
    let m = read_material(&tape).unwrap();
    assert_eq!((m.mat, m.z, m.a, m.state), (3580, 95, 242, 1));
    assert_eq!(m.excitation, 48_630.0);
}

#[test]
fn reads_decay_modes_ratios_and_energies() {
    let e = read_decay(&am242m_decay()).unwrap();
    assert_eq!((e.material.z, e.material.a, e.material.state), (95, 242, 1));
    let d = e.decay;
    assert!(!d.stable);
    assert_eq!(d.half_life, 4.449622e9);
    assert_eq!(d.energy, [4.008322e4, 5.373671e3, 2.441003e4]);
    assert_eq!(
        d.modes,
        vec![
            DecayMode {
                steps: vec![Step::Alpha],
                final_state: 0,
                ratio: 4.59e-3
            },
            DecayMode {
                steps: vec![Step::IsomericTransition],
                final_state: 0,
                ratio: 9.9541e-1
            },
            DecayMode {
                steps: vec![Step::SpontaneousFission],
                final_state: 0,
                ratio: 1.6e-10
            },
        ]
    );
}

fn decay_with_rtyp(rtyp: &str, rfs: &str) -> String {
    let mat = 100;
    [
        TAPE_ID.to_string(),
        header(mat, "3.508800+4", "0.0", "0"),
        rec(["3.508800+4", "8.7e1", "0", "0", "0", "0"], mat, 8, 457),
        rec(["1.0", "0.0", "0", "0", "6", "0"], mat, 8, 457),
        rec(["1.0", "0.0", "2.0", "0.0", "0.0", "0.0"], mat, 8, 457),
        rec(["0.0", "1.0", "0", "0", "6", "1"], mat, 8, 457),
        rec([rtyp, rfs, "1.0+7", "0.0", "1.0", "0.0"], mat, 8, 457),
        rec(BLANK, mat, 8, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat()
}

#[test]
fn reads_each_digit_of_rtyp_as_a_step() {
    let cases = [
        ("1.0", vec![Step::BetaMinus]),
        ("1.5", vec![Step::BetaMinus, Step::Neutron]),
        (
            "1.555",
            vec![Step::BetaMinus, Step::Neutron, Step::Neutron, Step::Neutron],
        ),
        ("2.7", vec![Step::EcBetaPlus, Step::Proton]),
        ("2.4", vec![Step::EcBetaPlus, Step::Alpha]),
        ("7.7", vec![Step::Proton, Step::Proton]),
        ("1.1", vec![Step::BetaMinus, Step::BetaMinus]),
        ("2.6", vec![Step::EcBetaPlus, Step::SpontaneousFission]),
    ];
    for (rtyp, steps) in cases {
        let d = read_decay(&decay_with_rtyp(rtyp, "1.0")).unwrap().decay;
        assert_eq!(d.modes[0].steps, steps, "{rtyp}");
        assert_eq!(d.modes[0].final_state, 1);
    }
}

#[test]
fn rejects_unknown_decay_types() {
    for rtyp in ["1.000000+1", "0.0", "8.0", "1.9"] {
        assert!(read_decay(&decay_with_rtyp(rtyp, "0.0")).is_err(), "{rtyp}");
    }
}

#[test]
fn rejects_a_fractional_daughter_state() {
    assert!(read_decay(&decay_with_rtyp("1.0", "0.5")).is_err());
}

#[test]
fn reads_a_stable_nuclide_without_modes() {
    let mat = 10;
    let text = [
        TAPE_ID.to_string(),
        header(mat, "2.004000+3", "0.0", "0"),
        rec(
            ["2.004000+3", "3.968219+0", "0", "0", "1", "0"],
            mat,
            8,
            457,
        ),
        rec(["0.0", "0.0", "0", "0", "6", "0"], mat, 8, 457),
        rec(["0.0", "0.0", "0.0", "0.0", "0.0", "0.0"], mat, 8, 457),
        rec(["0.0", "1.0", "0", "0", "6", "0"], mat, 8, 457),
        rec(["0.0", "0.0", "0.0", "0.0", "0.0", "0.0"], mat, 8, 457),
        rec(BLANK, mat, 8, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat();
    let d = read_decay(&text).unwrap().decay;
    assert!(d.stable);
    assert!(d.modes.is_empty());
}

fn mf3(mat: i32, mt: u32) -> String {
    mf3_q(mat, mt, "0.0")
}

/// QI を `qi` とした MF3 の節
fn mf3_q(mat: i32, mt: u32, qi: &str) -> String {
    [
        rec(["9.524100+4", "2.389860+2", "0", "0", "0", "0"], mat, 3, mt),
        rec(["2.0+8", qi, "0", "0", "1", "2"], mat, 3, mt),
        rec(["2", "2", "", "", "", ""], mat, 3, mt),
        rec(["1.0-5", "1.0", "2.0+7", "1.0", "", ""], mat, 3, mt),
        rec(BLANK, mat, 3, 0),
    ]
    .concat()
}

fn am241_neutron() -> String {
    let mat = 9543;
    [
        TAPE_ID.to_string(),
        header(mat, "9.524100+4", "0.0", "0"),
        mf3(mat, 1),
        mf3(mat, 16),
        mf3_q(mat, 18, "1.95+8"),
        mf3(mat, 102),
        rec(BLANK, mat, 0, 0),
        // (n,2n): NO = 1 の形（CONT）
        rec(["9.524100+4", "2.389860+2", "0", "0", "2", "1"], mat, 8, 16),
        rec(["9.524000+4", "0.0", "10", "0", "0", "0"], mat, 8, 16),
        rec(
            ["9.524000+4", "2.270000+3", "10", "1", "0", "0"],
            mat,
            8,
            16,
        ),
        rec(BLANK, mat, 8, 0),
        // 核分裂は読まない MT
        rec(["9.524100+4", "2.389860+2", "0", "0", "1", "1"], mat, 8, 18),
        rec(["0.0", "0.0", "6", "0", "0", "0"], mat, 8, 18),
        rec(BLANK, mat, 8, 0),
        // (n,γ): NO = 0 の形（LIST）。ZAP は 95242.01 と書かれている（ENDF/B-VIII.0 の Am-241 と同じ）
        rec(
            ["9.524100+4", "2.389860+2", "0", "0", "2", "0"],
            mat,
            8,
            102,
        ),
        rec(["9.524201+4", "0.0", "9", "0", "6", "0"], mat, 8, 102),
        rec(
            ["5.7672+4", "1.4", "9.624201+4", "0.37", "6.6+5", "2.01"],
            mat,
            8,
            102,
        ),
        rec(
            ["9.524201+4", "4.863000+4", "9", "2", "6", "0"],
            mat,
            8,
            102,
        ),
        rec(
            ["4.44+9", "3.1", "9.524201+4", "0.82", "4.8+4", "2.0"],
            mat,
            8,
            102,
        ),
        rec(BLANK, mat, 8, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat()
}

#[test]
fn reads_mf3_reactions_and_final_states() {
    let e = read_neutron(&am241_neutron(), &[16, 102, 103]).unwrap();
    assert_eq!((e.material.z, e.material.a, e.material.state), (95, 241, 0));
    assert_eq!(e.mts, vec![1, 16, 18, 102]);
    // 核分裂の Q は MT18 の QI（QM ではない）
    assert_eq!(e.fission_q, Some(1.95e8));
    assert_eq!(
        e.final_states,
        vec![
            (
                16,
                vec![
                    FinalStateEntry {
                        za: 95240,
                        excitation: 0.0,
                        mf: 10,
                        lfs: 0
                    },
                    FinalStateEntry {
                        za: 95240,
                        excitation: 2270.0,
                        mf: 10,
                        lfs: 1
                    },
                ]
            ),
            (
                102,
                vec![
                    FinalStateEntry {
                        za: 95242,
                        excitation: 0.0,
                        mf: 9,
                        lfs: 0
                    },
                    FinalStateEntry {
                        za: 95242,
                        excitation: 48_630.0,
                        mf: 9,
                        lfs: 2
                    },
                ]
            ),
        ]
    );
}

#[test]
fn rejects_a_zap_far_from_an_integer() {
    let text = am241_neutron().replace("9.524201+4", "9.524250+4");
    assert!(read_neutron(&text, &[102]).is_err());
}

fn u235_yields() -> String {
    let mat = 9228;
    [
        TAPE_ID.to_string(),
        header(mat, "9.223500+4", "0.0", "0"),
        rec(
            ["9.223500+4", "2.330250+2", "2", "0", "0", "0"],
            mat,
            8,
            454,
        ),
        rec(["2.530000-2", "0.0", "2", "0", "8", "2"], mat, 8, 454),
        rec(
            ["5.413500+4", "0.0", "2.0-3", "1.0-4", "5.313500+4", "0.0"],
            mat,
            8,
            454,
        ),
        rec(["6.0-2", "1.0-3", "", "", "", ""], mat, 8, 454),
        rec(["5.000000+5", "0.0", "0", "0", "4", "1"], mat, 8, 454),
        rec(["5.413500+4", "1.0", "1.0-3", "0.0", "", ""], mat, 8, 454),
        rec(BLANK, mat, 8, 0),
        rec(BLANK, mat, 0, 0),
    ]
    .concat()
}

#[test]
fn reads_independent_yields_per_energy() {
    let e = read_yields(&u235_yields()).unwrap();
    assert_eq!((e.material.z, e.material.a), (92, 235));
    assert_eq!(
        e.tables,
        vec![
            (
                0.0253,
                vec![
                    RawYield {
                        za: 54135,
                        state: 0,
                        value: 2e-3
                    },
                    RawYield {
                        za: 53135,
                        state: 0,
                        value: 6e-2
                    },
                ]
            ),
            (
                5e5,
                vec![RawYield {
                    za: 54135,
                    state: 1,
                    value: 1e-3
                }]
            ),
        ]
    );
}

#[test]
fn rejects_a_yield_count_that_differs_from_the_header() {
    let text = u235_yields().replace(
        &rec(
            ["9.223500+4", "2.330250+2", "2", "0", "0", "0"],
            9228,
            8,
            454,
        ),
        &rec(
            ["9.223500+4", "2.330250+2", "3", "0", "0", "0"],
            9228,
            8,
            454,
        ),
    );
    assert!(read_yields(&text).is_err());
}

#[test]
fn rejects_a_tape_with_several_materials() {
    let second = am242m_decay().replace("3580", "3581").replace(TAPE_ID, "");
    assert!(read_decay(&(am242m_decay() + &second)).is_err());
}
