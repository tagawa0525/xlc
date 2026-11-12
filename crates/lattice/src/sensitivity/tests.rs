use super::*;

/// 1 群の、全断面積と自群の散乱 `scatter` だけを持つ物質（数密度 1 の核種の群定数の代わり）
fn one_group(total: f64, scatter: f64) -> Material {
    Material {
        total: vec![total],
        scatter: vec![vec![scatter]],
        anisotropic: Vec::new(),
        fission: vec![0.0],
        nu_fission: vec![0.0],
        fission_terms: Vec::new(),
    }
}

#[test]
fn sum_check_allows_rounding_of_cancelling_terms() {
    // 散乱の要素で核種の寄与が打ち消し合うと、和は項より桁が小さく、足す順の丸めの差が和に比べて大きい。許容は項の大きさで決める
    let list = vec![
        ("A".to_string(), 1.0, one_group(1.0, 0.1)),
        ("B".to_string(), 1.0, one_group(1.0, 0.2)),
        ("C".to_string(), 1.0, one_group(1.0, -0.3)),
    ];
    // 解いた物質は別の順に足したもの: 0.1 + (0.2 − 0.3) と (0.1 + 0.2) − 0.3 は相対で 1 ほど違う
    let solved = one_group(3.0, 0.1 + (0.2 - 0.3));
    assert_ne!(solved.scatter[0][0], (0.1 + 0.2) - 0.3);
    check_sum("fuel", &solved, &list).unwrap();
    // 項の大きさを超える食い違いは止める
    let wrong = one_group(3.0, 1e-6);
    assert!(check_sum("fuel", &wrong, &list).is_err());
}

#[test]
fn derivatives_tell_listed_nuclides_from_absent_ones() {
    // 数密度 0 で組成にある核種は、変化がなく感度 0（向きが定まらない）。どの物質の組成にもない核種は、名前の誤りの候補なので区別する
    let derivatives = Derivatives {
        materials: vec![
            Vec::new(),
            vec![
                ("U-235".to_string(), 7e-4, one_group(1.0, 0.1)),
                ("Xe-135".to_string(), 0.0, one_group(1.0, 0.1)),
            ],
        ],
    };
    assert_eq!(derivatives.changes("U-235").len(), 1);
    assert!(derivatives.changes("Xe-135").is_empty());
    assert!(derivatives.lists("Xe-135"));
    assert!(!derivatives.lists("Xe-136"));
}

/// 3 群の物質（細かい群 0、1 が少数群 1、群 2 が少数群 2）
fn three_group(total: [f64; 3], nu_fission: [f64; 3]) -> Material {
    use crate::material::FissionTerm;
    Material {
        total: total.to_vec(),
        scatter: vec![
            vec![0.2, 0.15, 0.05],
            vec![0.0, 0.3, 0.2],
            vec![0.0, 0.01, 0.9],
        ],
        anisotropic: Vec::new(),
        fission: nu_fission.iter().map(|v| v / 2.5).collect(),
        nu_fission: nu_fission.to_vec(),
        fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
            vec![FissionTerm {
                chi: vec![0.9, 0.1, 0.0],
                nu_fission: nu_fission.to_vec(),
            }]
        } else {
            Vec::new()
        },
    }
}

#[test]
fn response_sources_are_orthogonal_to_the_flux_and_weight_the_few_groups() {
    use crate::geometry::{Lattice, PinCell};
    let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let (f, m) = (
        three_group([0.5, 0.7, 1.2], [0.01, 0.02, 0.3]),
        three_group([0.4, 0.6, 1.0], [0.0; 3]),
    );
    let materials = [&f, &m];
    let solution = crate::moc::solve(&tracking, &materials, &Options::default()).unwrap();
    let group_of = [0, 0, 1];
    let responses = [
        Response::Scatter(0, 1),
        Response::NuFission(1),
        Response::Total(0),
    ];
    let sources = response_sources(&responses, &group_of, &materials, &solution).unwrap();
    let inner = |h: &[Vec<f64>]| -> f64 {
        h.iter()
            .zip(&solution.flux)
            .zip(&solution.volumes)
            .map(|((h, f), v)| v * h.iter().zip(f).map(|(a, b)| a * b).sum::<f64>())
            .sum()
    };
    for s in &sources {
        let size: f64 = s
            .iter()
            .zip(&solution.flux)
            .zip(&solution.volumes)
            .map(|((h, f), v)| v * h.iter().zip(f).map(|(a, b)| a.abs() * b).sum::<f64>())
            .sum();
        assert!(inner(s).abs() < 1e-14 * size);
    }
    // Σs,1→2 の分子は、少数群 1 の細かい群から群 2 への散乱: 群 0 は 0.05、群 1 は 0.2。少数群 2 の群には値がない
    let numerator = |r: usize, g: usize| {
        let h1 = [[0.05, 0.2, 0.0], [0.05, 0.2, 0.0]];
        h1[r][g]
    };
    let a: f64 = (0..2)
        .map(|r| {
            solution.volumes[r]
                * (0..3)
                    .map(|g| numerator(r, g) * solution.flux[r][g])
                    .sum::<f64>()
        })
        .sum();
    let b: f64 = (0..2)
        .map(|r| solution.volumes[r] * (solution.flux[r][0] + solution.flux[r][1]))
        .sum();
    for (r, s) in sources[0].iter().enumerate() {
        for (g, value) in s.iter().enumerate() {
            let expected = numerator(r, g) / a - if g < 2 { 1.0 / b } else { 0.0 };
            assert!((value - expected).abs() < 1e-12 * expected.abs().max(1.0 / b));
        }
    }
    // νΣf の少数群 2 は、燃料の領域だけが分子を持つ
    assert!(sources[1][1][2] < 0.0 && sources[1][0][2] > 0.0);
    // 核分裂のない格子の νΣf の応答は定まらない
    let water = [&m, &m];
    let solution_water = crate::moc::solve(&tracking, &[&f, &m], &Options::default()).unwrap();
    assert!(
        response_sources(
            &[Response::NuFission(0)],
            &group_of,
            &water,
            &solution_water
        )
        .is_err()
    );
}

#[test]
fn generalized_cmfd_takes_the_cmfd_options_of_the_solved_state() {
    let cmfd = CmfdOptions {
        flux_tolerance: 1e-9,
        ..CmfdOptions::default()
    };
    assert_eq!(generalized_cmfd(&Acceleration::Cmfd(cmfd)).unwrap(), cmfd);
    // 状態を CMFD なしで解いたなら、粗メッシュの判定を入力が決めていない
    assert!(generalized_cmfd(&Acceleration::None).is_err());
}

#[test]
fn table_responses_list_the_few_group_constants_of_the_table() {
    let names: Vec<String> = table_responses(2).iter().map(ToString::to_string).collect();
    assert_eq!(
        names,
        [
            "total-1",
            "total-2",
            "nu-fission-1",
            "nu-fission-2",
            "fission-1",
            "fission-2",
            "scatter-1-1",
            "scatter-1-2",
            "scatter-2-1",
            "scatter-2-2",
        ]
    );
}

/// 3 群の物質に P1 の散乱行列（P0 の `factor` 倍）を足したもの（拡散係数の B1 が P1 を要る）
fn with_p1(mut m: Material, factor: f64) -> Material {
    m.anisotropic = vec![
        m.scatter
            .iter()
            .map(|row| row.iter().map(|s| factor * s).collect())
            .collect(),
    ];
    m
}

/// 燃料と水の 3 群のピンセル（鏡面反射）の、トラッキングと物質
fn three_group_pin_cell() -> (Tracking, Material, Material) {
    use crate::geometry::{Lattice, PinCell};
    let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    // 水の熱群の散乱を燃料と違う値にする（すべての領域で同じなら、熱群の 1 つの細かい群だけの応答 Σs,2→H は中性子束の形によらず、
    // 源が 0 になる）
    let mut water = three_group([0.4, 0.6, 1.0], [0.0; 3]);
    water.scatter[2][1] = 0.03;
    water.scatter[2][2] = 0.85;
    (
        tracking,
        with_p1(three_group([0.5, 0.7, 1.2], [0.01, 0.02, 0.3]), 0.2),
        with_p1(water, 0.4),
    )
}

#[test]
fn response_values_are_the_table_few_group_constants() {
    use xlc_format::group::GroupStructure;
    let (tracking, f, m) = three_group_pin_cell();
    let materials = [&f, &m];
    let solution = crate::moc::solve(&tracking, &materials, &Options::default()).unwrap();
    let fine = GroupStructure::new(vec![2e7, 1e5, 0.625, 1e-5]).unwrap();
    let bounds = [2e7, 0.625, 1e-5];
    let group_of = coarse_groups(&fine, &bounds).unwrap();
    let responses = table_responses(2);
    assert_eq!(responses.len(), 10);
    let values = response_values(&responses, &group_of, &materials, &solution);
    let few = crate::fewgroup::homogenize(
        &solution,
        &materials,
        &fine,
        &bounds,
        &tracking.boundaries(),
    )
    .unwrap();
    for (response, (a, b)) in responses.iter().zip(&values) {
        let table = match *response {
            Response::Total(g) => few.total[g],
            Response::NuFission(g) => few.nu_fission[g],
            Response::Fission(g) => few.fission[g],
            Response::Scatter(g, h) => few.scatter[g][h],
        };
        assert!(
            (a / b - table).abs() <= 1e-14 * table.abs(),
            "{response}: {} vs {table}",
            a / b
        );
    }
}

#[test]
fn few_group_first_order_changes_match_central_differences() {
    // 燃料の全反応（全断面積、散乱、核分裂）を変える。δR は分子の断面積が変わる直接の項と、中性子束を通じた項の和
    let (tracking, f, m) = three_group_pin_cell();
    let materials = [&f, &m];
    let options = Options {
        k_tolerance: 1e-12,
        flux_tolerance: 1e-11,
        acceleration: Acceleration::Cmfd(CmfdOptions {
            k_tolerance: 1e-13,
            flux_tolerance: 1e-12,
            ..CmfdOptions::default()
        }),
        ..Options::default()
    };
    let (solution, forward) = crate::moc::solve_symmetric_from(
        &tracking,
        &materials,
        &options,
        &crate::symmetry::Detected::none(&tracking),
        None,
    )
    .unwrap();
    let (_, adjoint) =
        crate::moc::solve_adjoint_from(&tracking, &materials, &options, None, None).unwrap();
    let group_of = [0, 0, 1];
    let responses = table_responses(2);
    let settings = GeneralizedOptions {
        method: GeneralizedMethod::Gmres { restart: 100 },
        tolerance: 1e-12,
        max_iterations: options.max_iterations,
        batch: false,
    };
    let solve = FewGroupSolve::new(
        &tracking,
        &materials,
        &options,
        None,
        (&solution, &forward, &adjoint),
        &responses,
        &group_of,
        &settings,
    )
    .unwrap();
    let d = with_p1(three_group([0.1, 0.2, 0.3], [0.002, 0.004, 0.05]), 0.2);
    let first = solve
        .first_order(&responses, &group_of, &solution, &[(0, &d)])
        .unwrap();
    let at = |h: f64| -> Vec<f64> {
        let perturbed = f.added(&d, h).unwrap();
        let m2 = [&perturbed, &m];
        let s = crate::moc::solve(&tracking, &m2, &options).unwrap();
        response_values(&responses, &group_of, &m2, &s)
            .iter()
            .map(|(a, b)| a / b)
            .collect()
    };
    let central = |h: f64| -> Vec<f64> {
        at(h)
            .iter()
            .zip(at(-h))
            .map(|(p, q)| (p - q) / (2.0 * h))
            .collect()
    };
    let (coarse, fine) = (central(0.02), central(0.01));
    for (i, response) in responses.iter().enumerate() {
        let direct = fine[i] + (fine[i] - coarse[i]) / 3.0;
        assert!(
            (first[i] - direct).abs() < 1e-6 * direct.abs(),
            "{response}: {} vs {direct}",
            first[i]
        );
    }
}

#[test]
fn scaled_inputs_change_only_positive_densities_of_the_nuclide() {
    let materials = vec![
        MaterialInput {
            name: "fuel".into(),
            temperature: 900.0,
            composition: vec![("U-235".into(), 7e-4), ("Pu-240".into(), 0.0)],
        },
        MaterialInput {
            name: "clad".into(),
            temperature: 600.0,
            composition: vec![("Zr-90".into(), 2e-2)],
        },
    ];
    let scaled = scaled_inputs(&materials, "U-235", 0.1);
    assert_eq!(scaled[0].composition[0], ("U-235".into(), 7e-4 * 1.1));
    assert_eq!(scaled[0].composition[1], ("Pu-240".into(), 0.0));
    assert_eq!(scaled[1], materials[1]);
    assert_eq!(scaled_inputs(&materials, "Pu-240", 0.1), materials);
}

#[test]
fn zero_table_constants_change_only_through_their_numerators() {
    // 3 つの少数群（細かい群ごと）では、熱群から高速群への散乱 Σs,3→1 がどの物質でも 0。中性子束の形によらず 0 なので一般化
    // 随伴を解かず、変化は分子の断面積の直接の項 ⟨δh₁, φ⟩ / ⟨h₂, φ⟩ だけ
    let (tracking, f, m) = three_group_pin_cell();
    let materials = [&f, &m];
    let options = Options {
        k_tolerance: 1e-12,
        flux_tolerance: 1e-11,
        acceleration: Acceleration::Cmfd(CmfdOptions {
            k_tolerance: 1e-13,
            flux_tolerance: 1e-12,
            ..CmfdOptions::default()
        }),
        ..Options::default()
    };
    let (solution, forward) = crate::moc::solve_symmetric_from(
        &tracking,
        &materials,
        &options,
        &crate::symmetry::Detected::none(&tracking),
        None,
    )
    .unwrap();
    let (_, adjoint) =
        crate::moc::solve_adjoint_from(&tracking, &materials, &options, None, None).unwrap();
    let group_of = [0, 1, 2];
    let responses = table_responses(3);
    let settings = GeneralizedOptions {
        method: GeneralizedMethod::Gmres { restart: 100 },
        tolerance: 1e-10,
        max_iterations: options.max_iterations,
        batch: false,
    };
    let solve = FewGroupSolve::new(
        &tracking,
        &materials,
        &options,
        None,
        (&solution, &forward, &adjoint),
        &responses,
        &group_of,
        &settings,
    )
    .unwrap();
    let mut d = three_group([0.1, 0.2, 0.3], [0.0; 3]);
    d.scatter[2][0] = 0.004;
    let first = solve
        .first_order(&responses, &group_of, &solution, &[(0, &d)])
        .unwrap();
    let at = responses
        .iter()
        .position(|r| *r == Response::Scatter(2, 0))
        .unwrap();
    let b: f64 = solution
        .volumes
        .iter()
        .zip(&solution.flux)
        .map(|(v, phi)| v * phi[2])
        .sum();
    let expected = solution.volumes[0] * 0.004 * solution.flux[0][2] / b;
    assert!(
        (first[at] - expected).abs() < 1e-14 * expected,
        "{} vs {expected}",
        first[at]
    );
}
