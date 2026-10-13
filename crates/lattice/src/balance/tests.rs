use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::material::FissionTerm;
use crate::moc::{Options, solve};
use crate::tracking::Tracking;

fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
    Material {
        total: total.to_vec(),
        scatter: scatter.iter().map(|r| r.to_vec()).collect(),
        anisotropic: Vec::new(),
        fission: nu_fission.iter().map(|v| v / 2.5).collect(),
        nu_fission: nu_fission.to_vec(),
        fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
            vec![FissionTerm {
                chi: vec![1.0, 0.0],
                nu_fission: nu_fission.to_vec(),
            }]
        } else {
            vec![]
        },
    }
}

#[test]
fn the_balance_of_a_reflective_cell_closes() {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    // (n,2n) のように、散乱の行の和が全断面積を超える群のある減速材
    let moderator = two_group([0.6, 2.0], [[0.52, 0.1], [0.0, 1.98]], [0.0, 0.0]);
    let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
    let materials = [&fuel, &fuel, &moderator, &moderator];
    let names: Vec<String> = ["fuel", "fuel", "water", "water"]
        .into_iter()
        .map(String::from)
        .collect();
    let solution = solve(&tracking, &materials, &Options::default()).unwrap();
    // 物理的な中性子束は MOC の値の 2 倍とする（中性子束だけが変わり、反応率は変わらない）
    let physical: Vec<Vec<f64>> = solution
        .flux
        .iter()
        .map(|f| f.iter().map(|v| 2.0 * v).collect())
        .collect();
    let rates = by_material(&solution, &materials, &names, &physical).unwrap();
    // 生成の和は 1、正味の吸収の和は 1/k（外周が鏡面反射）
    let production: f64 = rates.values().flat_map(|r| &r.production).sum();
    let absorption: f64 = rates.values().flat_map(|r| &r.absorption).sum();
    assert!((production - 1.0).abs() < 1e-12, "{production}");
    assert!(
        (absorption * solution.k - 1.0).abs() < 1e-6,
        "{absorption} {}",
        solution.k
    );
    // 群から出る移行の和と正味の吸収の和は、全断面積の反応率
    let water = &rates["water"];
    let phi_moc: f64 = water.flux[0] / 2.0;
    let out: f64 = water.transfer[0].iter().sum();
    assert!((out + water.absorption[0] - 0.6 * phi_moc).abs() < 1e-12);
    // 中性子束だけが物理的な値（MOC の 2 倍）になる
    let fuel_flux: f64 = (0..2)
        .map(|r| solution.flux[r][1] * solution.volumes[r])
        .sum::<f64>();
    assert!(
        (rates["fuel"].flux[1] - 2.0 * fuel_flux / production_of(&solution, &materials)).abs()
            < 1e-12
    );
}

#[test]
fn mismatched_inputs_are_errors() {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let materials = [&fuel, &fuel];
    let names: Vec<String> = vec!["fuel".into(); 2];
    let solution = solve(&tracking, &materials, &Options::default()).unwrap();
    assert!(matches!(
        by_material(&solution, &materials[..1], &names, &solution.flux),
        Err(BalanceError::Regions { .. })
    ));
    let short: Vec<Vec<f64>> = solution.flux.iter().map(|f| f[..1].to_vec()).collect();
    assert!(matches!(
        by_material(&solution, &materials, &names, &short),
        Err(BalanceError::Groups(0))
    ));
}

/// 規格化の前の生成の和
fn production_of(solution: &Solution, materials: &[&Material]) -> f64 {
    materials
        .iter()
        .enumerate()
        .map(|(r, m)| {
            m.nu_fission
                .iter()
                .zip(&solution.flux[r])
                .map(|(n, f)| n * f)
                .sum::<f64>()
                * solution.volumes[r]
        })
        .sum()
}
