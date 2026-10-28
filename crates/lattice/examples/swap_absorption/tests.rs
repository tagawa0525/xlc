use super::*;
use xlc_lattice::material::FissionTerm;

/// 3 群の物質（群 2 の散乱の行の和は全断面積を超える。(n,2n) のある群）
fn material(terms: &[[f64; 3]]) -> Material {
    Material {
        total: vec![0.5, 1.25, 2.0],
        scatter: vec![
            vec![0.25, 0.125, 0.0],
            vec![0.0, 1.0, 0.5],
            vec![0.0, 0.0, 1.5],
        ],
        anisotropic: vec![vec![
            vec![0.125, 0.0, 0.0],
            vec![0.0, 0.25, 0.0],
            vec![0.0, 0.0, 0.25],
        ]],
        fission: vec![0.0; 3],
        nu_fission: vec![0.0625 * terms.len() as f64; 3],
        fission_terms: terms
            .iter()
            .map(|chi| FissionTerm {
                chi: chi.to_vec(),
                nu_fission: vec![0.0625; 3],
            })
            .collect(),
    }
}

/// 中性子の過剰（(n,2n) の分。群 2 だけ）
const EXCESS: [f64; 3] = [0.0, 0.5, 0.0];

/// 物質の多重度を差し引かない吸収（正味の吸収 + 中性子の過剰）
fn gross(m: &Material) -> Vec<f64> {
    m.absorption()
        .iter()
        .zip(EXCESS)
        .map(|(a, e)| a + e)
        .collect()
}

#[test]
fn a_material_ratio_of_one_leaves_the_material_unchanged() {
    let m = material(&[]);
    // 中性子束を 2 の冪にして、比がちょうど 1 になるようにする
    let flux = [4.0, 2.0, 0.5];
    let absorption: Vec<f64> = gross(&m).iter().zip(flux).map(|(a, f)| a * f).collect();
    let ratio = material_ratio(&absorption, &flux, &gross(&m), &(0..3)).unwrap();
    assert_eq!(ratio, vec![1.0; 3]);
    assert_eq!(with_absorption(&m, &ratio, &EXCESS), m);
}

#[test]
fn a_material_ratio_moves_only_the_given_groups() {
    let m = material(&[]);
    let flux = [1.0; 3];
    let absorption: Vec<f64> = gross(&m).iter().map(|a| 2.0 * a).collect();
    let ratio = material_ratio(&absorption, &flux, &gross(&m), &(1..2)).unwrap();
    assert_eq!(ratio, vec![1.0, 2.0, 1.0]);
    let swapped = with_absorption(&m, &ratio, &EXCESS);
    // 群 2 の全断面積が、多重度を差し引かない吸収（1.25 − 1.5 + 0.5 = 0.25）だけ増える
    assert_eq!(swapped.total, vec![0.5, 1.5, 2.0]);
    assert_eq!(swapped.scatter, m.scatter);
    assert_eq!(swapped.anisotropic, m.anisotropic);
    assert_eq!(swapped.fission_terms, m.fission_terms);
}

#[test]
fn a_material_ratio_needs_absorption_and_flux_in_its_groups() {
    let m = material(&[]);
    let absorption = gross(&m);
    // 群 1 の多重度を差し引かない吸収を 0 にする
    let mut xlc = gross(&m);
    xlc[0] = 0.0;
    assert!(material_ratio(&absorption, &[1.0; 3], &xlc, &(0..1)).is_err());
    assert!(material_ratio(&absorption, &[1.0, 0.0, 1.0], &gross(&m), &(1..2)).is_err());
    // 範囲の外の群は比べない
    assert!(material_ratio(&absorption, &[1.0; 3], &xlc, &(1..3)).is_ok());
}

#[test]
fn the_same_chi_leaves_the_material_unchanged() {
    let mut m = material(&[[0.75, 0.25, 0.0], [0.5, 0.5, 0.0]]);
    let before = m.clone();
    let names = ["U-235".to_string(), "U-238".to_string()];
    let chi = BTreeMap::from([
        ("U-235".to_string(), vec![0.75, 0.25, 0.0]),
        ("U-238".to_string(), vec![0.5, 0.5, 0.0]),
    ]);
    replace_chi(&mut m, &names, &chi).unwrap();
    assert_eq!(m, before);
}

#[test]
fn chi_moves_only_the_terms_of_its_nuclides() {
    let mut m = material(&[[0.75, 0.25, 0.0], [0.5, 0.5, 0.0]]);
    let before = m.clone();
    let names = ["U-235".to_string(), "U-238".to_string()];
    let chi = BTreeMap::from([("U-238".to_string(), vec![0.625, 0.375, 0.0])]);
    replace_chi(&mut m, &names, &chi).unwrap();
    assert_eq!(m.fission_terms[0], before.fission_terms[0]);
    assert_eq!(m.fission_terms[1].chi, vec![0.625, 0.375, 0.0]);
    assert_eq!(
        m.fission_terms[1].nu_fission,
        before.fission_terms[1].nu_fission
    );
    assert_eq!(m.total, before.total);
}

#[test]
fn chi_of_a_nuclide_without_a_fission_term_is_an_error() {
    let mut m = material(&[[0.75, 0.25, 0.0]]);
    let names = ["U-235".to_string()];
    let chi = BTreeMap::from([("Pu-239".to_string(), vec![0.75, 0.25, 0.0])]);
    assert!(replace_chi(&mut m, &names, &chi).is_err());
    // 群の数の違う χ
    let chi = BTreeMap::from([("U-235".to_string(), vec![1.0, 0.0])]);
    assert!(replace_chi(&mut m, &names, &chi).is_err());
}
