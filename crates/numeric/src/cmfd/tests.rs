use super::*;

fn tolerances() -> Tolerances {
    Tolerances {
        k: 1e-12,
        flux: 1e-10,
        max_iterations: 10_000,
    }
}

/// 2 群、χ = (1, 0) の 1 つのセルの定数（除去 r1、吸収 a2、散乱 s12、νΣf の f1 と f2）
const ONE_CELL: (f64, f64, f64, f64, f64) = (0.027, 0.08, 0.017, 0.005, 0.12);

fn spectrum(f1: f64, f2: f64, volume: f64) -> Vec<Spectrum> {
    vec![Spectrum {
        chi: vec![1.0, 0.0],
        chi_sum: 1.0,
        nu_fission: vec![volume * f1, volume * f2],
    }]
}

#[test]
fn one_cell_gives_the_infinite_medium_k() {
    let (r1, a2, s12, f1, f2) = ONE_CELL;
    let coarse = Coarse::new(
        2,
        vec![0],
        vec![1.0, 0.3],
        vec![spectrum(f1, f2, 1.0)],
        vec![vec![r1, 0.0, -s12, a2]],
        vec![Vec::new()],
    )
    .unwrap();
    let solution = coarse.solve(1.0, &tolerances()).unwrap();
    // k∞ = νΣf1 / Σr1 + νΣf2 Σ12 / (Σr1 Σa2)
    let exact = f1 / r1 + f2 * s12 / (r1 * a2);
    assert!(
        (solution.k - exact).abs() < 1e-10,
        "{} vs {exact}",
        solution.k
    );
    // 比は、解の中性子束（φ2/φ1 = Σ12/Σa2）と組み立てた中性子束の、生成率 1 に規格化したものの比
    let shape = [1.0, s12 / a2];
    let given = [1.0, 0.3];
    let norm = |p: [f64; 2]| f1 * p[0] + f2 * p[1];
    for g in 0..2 {
        let expected = (shape[g] / norm(shape)) / (given[g] / norm(given));
        assert!(
            (solution.ratio[g] - expected).abs() < 1e-9 * expected,
            "group {g}: {} vs {expected}",
            solution.ratio[g]
        );
    }
}

/// 1 群の 1 次元の平板（幅 h のセル n 個、両端は中性子束 0 に近い強い漏れ）の有限差分の行列
fn slab(n: usize, leak_at_ends: f64) -> Coarse {
    let (d, h, absorption, nu_fission) = (1.2, 2.0, 0.02, 0.025);
    let coupling = d / h;
    let matrices = (0..n)
        .map(|c| {
            let ends = usize::from(c == 0) + usize::from(c + 1 == n);
            let inner = 2 - ends;
            vec![absorption * h + coupling * inner as f64 + leak_at_ends * ends as f64]
        })
        .collect();
    let neighbors = (0..n)
        .map(|c| {
            let mut list = Vec::new();
            if c > 0 {
                list.push((c - 1, vec![-coupling]));
            }
            if c + 1 < n {
                list.push((c + 1, vec![-coupling]));
            }
            list
        })
        .collect();
    let fission = (0..n)
        .map(|_| {
            vec![Spectrum {
                chi: vec![1.0],
                chi_sum: 1.0,
                nu_fission: vec![nu_fission * h],
            }]
        })
        .collect();
    Coarse::new(
        1,
        (0..n).collect(),
        vec![1.0; n],
        fission,
        matrices,
        neighbors,
    )
    .unwrap()
}

#[test]
fn slab_matches_the_dense_eigenvalue() {
    let n = 12;
    let coarse = slab(n, 0.4);
    let solution = coarse.solve(1.0, &tolerances()).unwrap();
    // 同じ行列の固有値を、密な逆反復（固有値のずっと上のシフトなしで十分収束させる）で求める
    let (d, h, absorption, nu_fission, leak) = (1.2, 2.0, 0.02, 0.025, 0.4);
    let c = d / h;
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        let ends = usize::from(i == 0) + usize::from(i + 1 == n);
        a[i * n + i] = absorption * h + c * (2 - ends) as f64 + leak * ends as f64;
        if i > 0 {
            a[i * n + i - 1] = -c;
        }
        if i + 1 < n {
            a[i * n + i + 1] = -c;
        }
    }
    let lu = Lu::new(a, n).unwrap();
    let mut phi = vec![1.0; n];
    let mut k = 1.0;
    for _ in 0..20_000 {
        let mut next: Vec<f64> = phi.iter().map(|p| nu_fission * h * p).collect();
        lu.solve(&mut next);
        let rate: f64 = next.iter().map(|p| nu_fission * h * p).sum();
        let old: f64 = phi.iter().map(|p| nu_fission * h * p).sum();
        k = rate / old;
        phi = next.iter().map(|p| p / rate).collect();
    }
    assert!((solution.k - k).abs() < 1e-10, "{} vs {k}", solution.k);
    // 対称な問題の解は対称
    for c in 0..n / 2 {
        let (a, b) = (solution.ratio[c], solution.ratio[n - 1 - c]);
        assert!(
            (a - b).abs() < 1e-8 * a,
            "cells {c} and {}: {a} vs {b}",
            n - 1 - c
        );
    }
}

#[test]
fn starting_from_the_solution_stops_at_once() {
    let coarse = slab(8, 0.4);
    let solution = coarse.solve(1.0, &tolerances()).unwrap();
    // 解の中性子束を組み立てた中性子束にすると、比は 1 で、反復は変化の見積もりに要る数回で止まる
    let flux: Vec<f64> = solution.ratio.clone();
    let again = Coarse { flux, ..coarse };
    let repeated = again.solve(solution.k, &tolerances()).unwrap();
    assert!(
        repeated.iterations <= 3,
        "{} iterations",
        repeated.iterations
    );
    assert!((repeated.k - solution.k).abs() < 1e-12);
    for r in &repeated.ratio {
        assert!((r - 1.0).abs() < 1e-9, "{r}");
    }
}

#[test]
fn iteration_limit_is_an_error() {
    let coarse = slab(40, 0.4);
    let limited = Tolerances {
        max_iterations: 3,
        ..tolerances()
    };
    assert!(matches!(
        coarse.solve(1.0, &limited),
        Err(CoarseError::NotConverged(3))
    ));
}

#[test]
fn singular_block_names_its_cell() {
    let result = Coarse::new(
        1,
        vec![7],
        vec![1.0],
        vec![Vec::new()],
        vec![vec![0.0]],
        vec![Vec::new()],
    );
    assert!(matches!(result, Err(CoarseError::Singular(7))));
}

#[test]
fn no_fission_is_an_error() {
    let coarse = Coarse::new(
        1,
        vec![0],
        vec![1.0],
        vec![Vec::new()],
        vec![vec![1.0]],
        vec![Vec::new()],
    )
    .unwrap();
    assert!(matches!(
        coarse.solve(1.0, &tolerances()),
        Err(CoarseError::NoFission)
    ));
}

#[test]
fn coupling_that_drives_a_flux_negative_is_an_error() {
    // 隣との結合が正（補正で符号が変わりうる）で大きく、2 つ目のセルの中性子束が負になる
    let fission = vec![spectrum(0.005, 0.12, 1.0), spectrum(0.005, 0.12, 1.0)];
    let (r1, a2, s12, _, _) = ONE_CELL;
    let coarse = Coarse::new(
        2,
        vec![3, 4],
        vec![1.0; 4],
        fission,
        vec![vec![r1, 0.0, -s12, a2]; 2],
        vec![vec![(1, vec![0.0, 0.0])], vec![(0, vec![1.0, 1.0])]],
    )
    .unwrap();
    assert!(matches!(
        coarse.solve(1.0, &tolerances()),
        Err(CoarseError::NonPositiveSolution { cell: 4, .. })
    ));
}
