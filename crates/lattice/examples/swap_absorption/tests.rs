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

/// 3 群の核種の実効断面積。群 1 に (n,2n)（0.5 b）と (n,α)（0.25 b）があり、非弾性散乱は無い
fn effective() -> xlc_format::library::Shielded {
    xlc_format::library::Shielded {
        total: vec![4.0, 6.0, 10.0],
        elastic: vec![3.0, 5.0, 8.0],
        fission: vec![0.125, 0.25, 1.0],
        nu_fission: vec![0.375, 0.625, 2.5],
        capture: vec![0.125, 0.75, 1.0],
    }
}

/// `effective` の核種の無限希釈の P0 の散乱行列（群 1 の行の和は弾性散乱 3.5 b と (n,2n) の 2 × 0.5 b）と弾性散乱
fn rows() -> (Vec<Vec<f64>>, Vec<f64>) {
    (
        vec![
            vec![2.5, 2.0, 0.0],
            vec![0.0, 4.5, 1.0],
            vec![0.0, 0.0, 9.0],
        ],
        vec![3.5, 5.5, 9.0],
    )
}

const NUCLIDE_EXCESS: [f64; 3] = [0.5, 0.0, 0.0];

#[test]
fn the_other_absorption_is_the_disappearance_besides_capture_and_fission() {
    let (scatter, elastic) = rows();
    let other = xlc_other(&effective(), &scatter, &elastic, &NUCLIDE_EXCESS);
    assert_eq!(other, vec![0.25, 0.0, 0.0]);
}

/// 中性子束 2、数密度 0.5 の OpenMC の反応率（`rates` は反応ごとの群の値）
fn reference(rates: &[(&str, [f64; 3])]) -> Reference {
    Reference {
        bounds: vec![1e7, 1e5, 1.0, 1e-5],
        flux: vec![2.0; 3],
        densities: BTreeMap::from([("U-238".to_string(), 0.5)]),
        rates: BTreeMap::from([(
            "U-238".to_string(),
            rates
                .iter()
                .map(|(s, v)| (s.to_string(), Score { mean: v.to_vec() }))
                .collect(),
        )]),
    }
}

#[test]
fn the_openmc_other_absorption_is_absorption_less_capture_and_fission() {
    let r = reference(&[
        ("absorption", [1.0, 2.0, 4.0]),
        ("(n,gamma)", [0.5, 1.5, 3.0]),
        ("fission", [0.25, 0.5, 1.0]),
    ]);
    // 反応率 / 中性子束 / 数密度 (barn)。範囲の外の群は 0
    assert_eq!(
        openmc_other(&r, "U-238", &(0..2), 3).unwrap(),
        vec![0.25, 0.0, 0.0]
    );
}

#[test]
fn the_openmc_excess_counts_the_extra_neutrons() {
    let r = reference(&[
        ("(n,2n)", [0.5, 0.25, 0.0]),
        ("(n,3n)", [0.125, 0.0, 0.0]),
        ("(n,4n)", [0.0625, 0.0, 0.0]),
    ]);
    // σ(n,2n) + 2σ(n,3n) + 3σ(n,4n)
    assert_eq!(
        openmc_excess(&r, "U-238", &(0..1), 3).unwrap(),
        vec![0.5 + 0.25 + 0.1875, 0.0, 0.0]
    );
}

#[test]
fn an_openmc_value_without_its_scores_or_flux_is_an_error() {
    let r = reference(&[("absorption", [1.0, 2.0, 4.0])]);
    assert!(openmc_other(&r, "U-238", &(0..1), 3).is_err());
    assert!(openmc_excess(&r, "U-238", &(0..1), 3).is_err());
    let mut r = reference(&[
        ("absorption", [1.0, 2.0, 4.0]),
        ("(n,gamma)", [0.5, 1.5, 3.0]),
        ("fission", [0.25, 0.5, 1.0]),
    ]);
    r.flux[1] = 0.0;
    assert!(openmc_other(&r, "U-238", &(0..2), 3).is_err());
    assert!(openmc_other(&r, "U-235", &(0..1), 3).is_err());
}

/// 1 つの燃料の核種（核分裂の項 1 つ）の物質。SPH 因子 1、数密度 0.5
fn fuel_material() -> Material {
    let s = effective();
    let density = 0.5;
    Material {
        total: s.total.iter().map(|v| density * v).collect(),
        scatter: rows()
            .0
            .iter()
            .map(|r| r.iter().map(|v| density * v).collect())
            .collect(),
        anisotropic: vec![],
        fission: s.fission.iter().map(|v| density * v).collect(),
        nu_fission: s.nu_fission.iter().map(|v| density * v).collect(),
        fission_terms: vec![FissionTerm {
            chi: vec![1.0, 0.0, 0.0],
            nu_fission: s.nu_fission.iter().map(|v| density * v).collect(),
        }],
    }
}

#[test]
fn a_ratio_of_one_and_no_shift_leave_the_fuel_unchanged() {
    let mut m = fuel_material();
    swap_nuclide(
        &mut m,
        Some(0),
        &[1.0; 3],
        0.5,
        &effective(),
        &Ratio::unchanged(3),
    );
    assert_eq!(m, fuel_material());
}

#[test]
fn a_shift_moves_only_the_total_of_its_groups() {
    let mut m = fuel_material();
    let mut r = Ratio::unchanged(3);
    r.shift = vec![0.25, 0.0, 0.0];
    let mu = [2.0, 1.0, 1.0];
    swap_nuclide(&mut m, Some(0), &mu, 0.5, &effective(), &r);
    let before = fuel_material();
    // μ N Δσ = 2 × 0.5 × 0.25
    assert_eq!(
        m.total,
        vec![before.total[0] + 0.25, before.total[1], before.total[2]]
    );
    assert_eq!(m.scatter, before.scatter);
    assert_eq!(m.fission, before.fission);
    assert_eq!(m.nu_fission, before.nu_fission);
    assert_eq!(m.fission_terms, before.fission_terms);
}

#[test]
fn a_shift_is_the_openmc_value_less_the_xlc_value_in_its_groups() {
    assert_eq!(
        shift(&[0.5, 0.25, 1.0], &[0.25, 0.25, 0.0], &(0..2)),
        vec![0.25, 0.0, 0.0]
    );
}

#[test]
fn other_and_excess_swaps_name_fuel_nuclides_and_groups() {
    let fuel = ["O-16".to_string(), "U-238".to_string()];
    let r = reference(&[]);
    assert_eq!(
        parse("*!1-2,U-238^1-1", &fuel, &r, 3).unwrap(),
        vec![
            Swap::Other("U-238".to_string(), 0..2),
            Swap::Excess("U-238".to_string(), 0..1),
        ]
    );
    assert!(parse("Pu-239!1-1", &fuel, &r, 3).is_err());
    assert!(parse("U-238^0-1", &fuel, &r, 3).is_err());
}
