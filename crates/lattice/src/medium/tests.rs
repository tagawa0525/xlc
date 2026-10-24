use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::material::FissionTerm;
use crate::moc::{Options, solve};
use crate::tracking::{Boundary, Tracking};
use xlc_format::medium::Medium;
use xlc_format::table::Axis;

/// 2群の物質（P1 と P2 の散乱は P0 の 0.3 倍と 0.1 倍）。核分裂は2つの核種の項に分ける
fn material(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
    let moment = |c: f64| -> Vec<Vec<f64>> {
        scatter
            .iter()
            .map(|r| r.iter().map(|s| c * s).collect())
            .collect()
    };
    let fissile = nu_fission.iter().any(|&v| v > 0.0);
    Material {
        total: total.to_vec(),
        scatter: moment(1.0),
        anisotropic: vec![moment(0.3), moment(0.1)],
        fission: nu_fission.iter().map(|v| v / 2.5).collect(),
        nu_fission: nu_fission.to_vec(),
        fission_terms: if fissile {
            vec![
                FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.iter().map(|v| 0.75 * v).collect(),
                },
                FissionTerm {
                    chi: vec![0.6, 0.4],
                    nu_fission: nu_fission.iter().map(|v| 0.25 * v).collect(),
                },
            ]
        } else {
            vec![]
        },
    }
}

/// 燃料、真空の隙間、水のピンセル（外周は鏡面反射）の解と、領域の物質
fn pin_cell() -> (Solution, Vec<Material>) {
    let fuel = material([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = material([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let materials = vec![fuel, Material::void(2), water.clone(), water];
    let cell = PinCell::new(1.26, vec![0.4096, 0.418, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
    let refs: Vec<&Material> = materials.iter().collect();
    let options = Options {
        scattering_order: 2,
        ..Options::default()
    };
    (solve(&tracking, &refs, &options).unwrap(), materials)
}

#[test]
fn medium_reproduces_the_lattice_k_and_spectrum() {
    // 外周が鏡面反射の格子の媒質は、反応率と生成率を保存するので、その無限媒質の解は格子の k と平均の細かい群の
    // 中性子束に一致する（ADR-0128）
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let medium = homogenize_medium(&solution, &refs, &Boundaries::REFLECTIVE).unwrap();
    let m = Material::from_medium(&medium).unwrap();
    let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
    let options = Options {
        scattering_order: 2,
        ..Options::default()
    };
    let homogeneous = solve(&tracking, &[&m, &m], &options).unwrap();
    // 固有値の反復の許容差（k 1e-7）の程度
    assert!(
        (homogeneous.k - solution.k).abs() < 1e-6 * solution.k,
        "{} vs {}",
        homogeneous.k,
        solution.k
    );
    let spectrum = |s: &Solution| {
        let mut sum = [0.0; 2];
        for (phi, v) in s.flux.iter().zip(&s.volumes) {
            for g in 0..2 {
                sum[g] += v * phi[g];
            }
        }
        sum[1] / sum[0]
    };
    let (a, b) = (spectrum(&homogeneous), spectrum(&solution));
    assert!((a - b).abs() < 1e-5 * b, "{a} vs {b}");
}

#[test]
fn medium_weights_every_moment_and_fission_by_the_flux() {
    // どの次数の散乱もスカラー中性子束と体積で重みづけし、χ は核種の項の χ を生成率で重みづけした1つ
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let medium = homogenize_medium(&solution, &refs, &Boundaries::REFLECTIVE).unwrap();
    assert_eq!(medium.scatter.len(), 3);
    let weight = |g: usize| -> Vec<f64> {
        solution
            .flux
            .iter()
            .zip(&solution.volumes)
            .map(|(phi, v)| v * phi[g])
            .collect()
    };
    for g in 0..2 {
        let w = weight(g);
        let total: f64 = w.iter().sum();
        let average = |value: &dyn Fn(&Material) -> f64| {
            refs.iter().zip(&w).map(|(m, w)| w * value(m)).sum::<f64>() / total
        };
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-14 * b.abs(), "{a} vs {b}");
        close(medium.total[g], average(&|m| m.total[g]));
        close(medium.nu_fission[g], average(&|m| m.nu_fission[g]));
        close(medium.fission[g], average(&|m| m.fission[g]));
        for to in 0..2 {
            close(medium.scatter[0][g][to], average(&|m| m.scatter[g][to]));
            for l in 1..3 {
                close(
                    medium.scatter[l][g][to],
                    average(&|m| m.anisotropic.get(l - 1).map_or(0.0, |s| s[g][to])),
                );
            }
        }
    }
    // 燃料だけが核分裂し、項の χ は (1, 0) と (0.6, 0.4) で生成率は 3 : 1
    let chi_fast = 0.75 + 0.25 * 0.6;
    assert!((medium.chi[0] - chi_fast).abs() < 1e-14, "{:?}", medium.chi);
    assert!((medium.chi.iter().sum::<f64>() - 1.0).abs() < 1e-14);
}

#[test]
fn medium_needs_a_reflective_lattice() {
    // 外周に真空の辺があれば、媒質の無限媒質の解は格子の解に一致しない
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let open = Boundaries {
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let err = homogenize_medium(&solution, &refs, &open)
        .unwrap_err()
        .to_string();
    assert!(err.contains("reflective"), "{err}");
}

#[test]
fn materials_come_back_from_the_medium() {
    // 媒質の物質は、P0 を `scatter`、P1 以上を `anisotropic` に、核分裂を χ の1つの項に持つ
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let medium = homogenize_medium(&solution, &refs, &Boundaries::REFLECTIVE).unwrap();
    let m = Material::from_medium(&medium).unwrap();
    assert_eq!(m.scatter, medium.scatter[0]);
    assert_eq!(m.anisotropic, medium.scatter[1..].to_vec());
    assert_eq!(m.fission_terms.len(), 1);
    assert_eq!(m.fission_terms[0].chi, medium.chi);
    assert_eq!(m.fission_terms[0].nu_fission, medium.nu_fission);
    // 核分裂のない媒質は項を持たない
    let water = MediumState {
        fission: vec![0.0; 2],
        nu_fission: vec![0.0; 2],
        chi: vec![0.0; 2],
        ..medium
    };
    assert!(
        Material::from_medium(&water)
            .unwrap()
            .fission_terms
            .is_empty()
    );
}

#[test]
fn media_must_match_the_library_and_the_branches() {
    // 媒質の群境界はライブラリと同じ、軸（名前、点、密度、B-10 の割合）は平板の分岐の軸と同じ（ADR-0128）
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let state = homogenize_medium(&solution, &refs, &Boundaries::REFLECTIVE).unwrap();
    let axes = vec![
        Axis::boron(vec![1300.0], 0.199),
        Axis::moderator(vec![600.0], vec![0.743]),
        Axis::fuel(vec![600.0, 900.0]),
    ];
    let bounds = vec![1.0e7, 0.625, 1.0e-5];
    let medium = Medium {
        bounds: bounds.clone(),
        axes: axes.clone(),
        states: vec![state.clone(), state.clone()],
    };
    assert!(check_medium(&medium, &bounds, &axes).is_ok());
    // 平板の軸の公称の値（ADR-0200）は比べない（媒質は持たない）
    let nominal: Vec<Axis> = axes
        .iter()
        .map(|a| {
            let value = a.points()[0];
            a.clone().with_nominal(value)
        })
        .collect();
    assert!(check_medium(&medium, &bounds, &nominal).is_ok());
    let err = check_medium(&medium, &[2.0e7, 0.625, 1.0e-5], &axes).unwrap_err();
    assert!(err.contains("bounds"), "{err}");
    let mut other = axes.clone();
    other[1] = Axis::moderator(vec![600.0], vec![0.7]);
    for slab in [vec![], other] {
        let err = check_medium(&medium, &bounds, &slab).unwrap_err();
        assert!(err.contains("axes"), "{err}");
    }
    let single = Medium {
        bounds: bounds.clone(),
        axes: vec![],
        states: vec![state],
    };
    assert!(check_medium(&single, &bounds, &[]).is_ok());
    assert!(check_medium(&single, &bounds, &axes).is_err());
}

#[test]
fn media_must_have_the_scattering_order_of_the_other_materials() {
    let (solution, materials) = pin_cell();
    let refs: Vec<&Material> = materials.iter().collect();
    let state = homogenize_medium(&solution, &refs, &Boundaries::REFLECTIVE).unwrap();
    let fuel = Material::from_medium(&state).unwrap();
    let mut named = vec![
        ("fuel".to_string(), fuel.clone()),
        ("void".to_string(), Material::void(2)),
        ("water".to_string(), materials[2].clone()),
    ];
    assert!(check_orders(&named, &["fuel".to_string()]).is_ok());
    named[2].1.anisotropic.pop();
    let err = check_orders(&named, &["fuel".to_string()]).unwrap_err();
    assert!(err.contains("fuel") && err.contains("water"), "{err}");
}
