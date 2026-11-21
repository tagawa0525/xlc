use super::*;
use xlc_numeric::lu::Lu;

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
        BlockSolve::Lu,
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

/// `Coarse::new` と `Coarse::colored` の引数（群の数、番号、組み立てた中性子束、核分裂の項、対角のブロック、隣との結合）
type Parts = (
    usize,
    Vec<usize>,
    Vec<f64>,
    Vec<Vec<Spectrum>>,
    Vec<Vec<f64>>,
    Vec<Vec<(usize, Vec<f64>)>>,
);

/// 1 群の 1 次元の平板（幅 h のセル n 個、両端は中性子束 0 に近い強い漏れ）の有限差分の行列
fn slab(n: usize, leak_at_ends: f64) -> Coarse {
    let (groups, labels, flux, fission, matrices, neighbors) = slab_parts(n, leak_at_ends);
    Coarse::new(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
    )
    .unwrap()
}

/// `slab` の引数
fn slab_parts(n: usize, leak_at_ends: f64) -> Parts {
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
    (
        1,
        (0..n).collect(),
        vec![1.0; n],
        fission,
        matrices,
        neighbors,
    )
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
        BlockSolve::Lu,
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
        BlockSolve::Lu,
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
        BlockSolve::Lu,
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

/// 2 群の nx × ny の平面（1 つのセルの定数 `ONE_CELL` に隣との結合と外周の漏れを足したもの）。組み立てた中性子束は位置で変え、
/// 解き始めの値が平らでない
fn plane_parts(nx: usize, ny: usize) -> Parts {
    let (r1, a2, s12, f1, f2) = ONE_CELL;
    let (c1, c2) = (0.004, 0.002);
    let cells = nx * ny;
    let mut matrices = Vec::new();
    let mut neighbors = Vec::new();
    for c in 0..cells {
        let (i, j) = (c % nx, c / nx);
        let mut list = Vec::new();
        let mut outer = 0;
        for (inside, other) in [
            (i > 0, c.wrapping_sub(1)),
            (i + 1 < nx, c + 1),
            (j > 0, c.wrapping_sub(nx)),
            (j + 1 < ny, c + nx),
        ] {
            if inside {
                list.push((other, vec![-c1, -c2]));
            } else {
                outer += 1;
            }
        }
        // 隣との結合の和と、外周の面からの漏れ（結合の 2 倍）を対角に足す
        let leak = (list.len() + 2 * outer) as f64;
        matrices.push(vec![r1 + c1 * leak, 0.0, -s12, a2 + c2 * leak]);
        neighbors.push(list);
    }
    let fission = (0..cells).map(|_| spectrum(f1, f2, 1.0)).collect();
    let flux = (0..cells)
        .flat_map(|c| [1.0 + 0.1 * (c % 7) as f64, 0.3 + 0.02 * (c % 5) as f64])
        .collect();
    (
        2,
        (100..100 + cells).collect(),
        flux,
        fission,
        matrices,
        neighbors,
    )
}

/// 平面の位置の市松の色
fn checkerboard(nx: usize, ny: usize) -> Vec<usize> {
    (0..nx * ny).map(|c| (c % nx + c / nx) % 2).collect()
}

#[test]
fn colored_solves_agree_with_the_cell_order() {
    // ADR-0234: 同じ色のセルを互いに読まずに解く順（市松の 2 色）は、固有値問題の解を変えない。どちらの順も同じ判定で止めるので、
    // 差は許容の 2 倍の中。比は与えたセルの順に返す
    for (nx, ny) in [(12, 1), (7, 5), (6, 6)] {
        let (groups, labels, flux, fission, matrices, neighbors) = plane_parts(nx, ny);
        let plain = Coarse::new(
            BlockSolve::Lu,
            groups,
            labels.clone(),
            flux.clone(),
            fission.clone(),
            matrices.clone(),
            neighbors.clone(),
        )
        .unwrap()
        .solve(1.0, &tolerances())
        .unwrap();
        let colored = Coarse::colored(
            BlockSolve::Lu,
            groups,
            labels,
            flux,
            fission,
            matrices,
            neighbors,
            &checkerboard(nx, ny),
        )
        .unwrap()
        .solve(1.0, &tolerances())
        .unwrap();
        let t = tolerances();
        assert!(
            (colored.k - plain.k).abs() < 2.0 * t.k,
            "{nx} x {ny}: k {} vs {}",
            colored.k,
            plain.k
        );
        assert_eq!(colored.ratio.len(), plain.ratio.len());
        for (cell, (a, b)) in colored.ratio.iter().zip(&plain.ratio).enumerate() {
            assert!(
                (a - b).abs() < 2.0 * t.flux * b,
                "{nx} x {ny}, cell {cell}: {a} vs {b}"
            );
        }
    }
}

#[test]
fn colored_solves_do_not_depend_on_the_thread_count() {
    // ADR-0234: 同じ色のセルは互いに読まず、和は決まった順に取るので、スレッドの数によらず解はビットで同じ
    let (nx, ny) = (40, 30);
    let solve_on = |threads: usize| {
        let (groups, labels, flux, fission, matrices, neighbors) = plane_parts(nx, ny);
        let coarse = Coarse::colored(
            BlockSolve::Lu,
            groups,
            labels,
            flux,
            fission,
            matrices,
            neighbors,
            &checkerboard(nx, ny),
        )
        .unwrap();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| coarse.solve(1.0, &tolerances()).unwrap())
    };
    let one = solve_on(1);
    for threads in [2, 3, 8] {
        assert_eq!(solve_on(threads), one, "{threads} threads");
    }
}

#[test]
fn colors_must_not_couple_cells_of_the_same_color() {
    // 同じ色の隣を読むと、並列に解いた結果が読む順で変わる。黙って順に解かず、組み立てた側の番号を付けたエラーにする
    let (groups, labels, flux, fission, matrices, neighbors) = slab_parts(4, 0.4);
    let result = Coarse::colored(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
        &[0, 1, 1, 0],
    );
    assert!(
        matches!(result, Err(CoarseError::SameColor { cell: 1, other: 2 })),
        "{result:?}"
    );
    let (groups, labels, flux, fission, matrices, neighbors) = slab_parts(4, 0.4);
    let result = Coarse::colored(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
        &[0, 1, 0],
    );
    assert!(
        matches!(result, Err(CoarseError::Colors { .. })),
        "{result:?}"
    );
}

#[test]
fn changes_at_rounding_level_count_as_converged() {
    let mut r = Remaining::default();
    r.update(1.0e-15);
    assert_eq!(r.update(2.0e-15), 0.0);
    let mut r = Remaining::default();
    assert_eq!(r.update(1.0), f64::INFINITY);
    // 縮小率 0.5 なら、最新の変化を含めて残りは 2 倍
    assert!((r.update(0.5) - 1.0).abs() < 1e-15);
}

/// 3 つのセルの 1 群の平板（`slab_parts` の 3 セル）を中央で鏡映した 2 つの代表（端のセルが 2 つの組、中央が 1 つ）の引数。
/// 組み立てた中性子束は端と中央で違え、全体の平板の [端, 中央, 端] と同じにする
fn mirrored_slab_parts() -> Parts {
    let (groups, _, _, fission, matrices, neighbors) = slab_parts(3, 0.4);
    let coupling = neighbors[1][0].1.clone();
    (
        groups,
        vec![0, 1],
        vec![1.0, 1.7],
        fission[..2].to_vec(),
        matrices[..2].to_vec(),
        vec![
            vec![(1, coupling.clone())],
            vec![(0, coupling.clone()), (0, coupling)],
        ],
    )
}

#[test]
fn weights_count_each_cell_as_many_times_as_it_stands_for() {
    // 解の比は生成率 1 に規格化した中性子束の比なので、代表の生成率を組のセルの数だけ数えると、全体を解いた平板の比と同じになる
    let (groups, labels, _, fission, matrices, neighbors) = slab_parts(3, 0.4);
    let full = Coarse::new(
        BlockSolve::Lu,
        groups,
        labels,
        vec![1.0, 1.7, 1.0],
        fission,
        matrices,
        neighbors,
    )
    .unwrap()
    .solve(1.0, &tolerances())
    .unwrap();
    let (groups, labels, flux, fission, matrices, neighbors) = mirrored_slab_parts();
    let plain = Coarse::new(
        BlockSolve::Lu,
        groups,
        labels.clone(),
        flux.clone(),
        fission.clone(),
        matrices.clone(),
        neighbors.clone(),
    )
    .unwrap();
    // 色ごとに解く形では、重みは与えたセルの順に受け取る
    let colored = Coarse::colored(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
        &[1, 0],
    )
    .unwrap();
    for coarse in [plain, colored] {
        let solution = coarse
            .weighted(vec![2.0, 1.0])
            .unwrap()
            .solve(1.0, &tolerances())
            .unwrap();
        assert!(
            (solution.k - full.k).abs() < 1e-11,
            "{} vs {}",
            solution.k,
            full.k
        );
        for (a, b) in solution.ratio.iter().zip(&full.ratio[..2]) {
            assert!((a - b).abs() < 1e-9 * b, "{a} vs {b}");
        }
    }
}

#[test]
fn unit_weights_leave_the_solution_bit_for_bit() {
    // core は重みを付けない。重み 1 の積は値を変えないので、重みを 1 にした問題の解はビットで同じ
    for colors in [None, Some(checkerboard(7, 5))] {
        let (groups, labels, flux, fission, matrices, neighbors) = plane_parts(7, 5);
        let cells = labels.len();
        let build = || {
            match &colors {
                None => Coarse::new(
                    BlockSolve::Lu,
                    groups,
                    labels.clone(),
                    flux.clone(),
                    fission.clone(),
                    matrices.clone(),
                    neighbors.clone(),
                ),
                Some(colors) => Coarse::colored(
                    BlockSolve::Lu,
                    groups,
                    labels.clone(),
                    flux.clone(),
                    fission.clone(),
                    matrices.clone(),
                    neighbors.clone(),
                    colors,
                ),
            }
            .unwrap()
        };
        let plain = build().solve(1.0, &tolerances()).unwrap();
        let weighted = build()
            .weighted(vec![1.0; cells])
            .unwrap()
            .solve(1.0, &tolerances())
            .unwrap();
        assert_eq!(weighted, plain);
    }
}

#[test]
fn weights_must_be_positive_and_one_per_cell() {
    for weights in [vec![2.0], vec![2.0, 0.0], vec![2.0, f64::NAN]] {
        let (groups, labels, flux, fission, matrices, neighbors) = mirrored_slab_parts();
        let result = Coarse::new(
            BlockSolve::Lu,
            groups,
            labels,
            flux,
            fission,
            matrices,
            neighbors,
        )
        .unwrap()
        .weighted(weights);
        assert!(
            matches!(
                result,
                Err(CoarseError::Weights { cells: 2, .. } | CoarseError::Weight { cell: 1, .. })
            ),
            "{result:?}"
        );
    }
}

#[test]
fn a_sweep_solves_each_block_by_the_lu_substitution() {
    // 対角のブロックは、xlc_numeric の LU 分解の前進・後退代入で解く（逆行列を掛ける形は、条件数の悪いブロックで誤差が大きい。
    // ADR-0275）。1 回は、先に解いたセルの新しい値と、まだ解いていないセルの前の値を読む
    let (groups, labels, flux, fission, matrices, neighbors) = plane_parts(3, 2);
    let ng = groups;
    let previous: Vec<f64> = (0..6 * ng).map(|j| 0.7 + 0.13 * (j % 5) as f64).collect();
    let source: Vec<f64> = (0..6 * ng).map(|j| 0.01 * (1 + j % 3) as f64).collect();
    let inverse_k = 1.0 / 1.07;
    let mut expected = previous.clone();
    for c in 0..6 {
        let mut b = source[c * ng..][..ng].to_vec();
        for s in &fission[c] {
            let rate = s
                .nu_fission
                .iter()
                .zip(&previous[c * ng..][..ng])
                .map(|(a, b)| a * b)
                .sum::<f64>()
                * inverse_k;
            for (r, x) in b.iter_mut().zip(&s.chi) {
                *r += x * rate;
            }
        }
        for (other, coefficients) in &neighbors[c] {
            for g in 0..ng {
                b[g] -= coefficients[g] * expected[other * ng + g];
            }
        }
        Lu::new(matrices[c].clone(), ng).unwrap().solve(&mut b);
        expected[c * ng..][..ng].copy_from_slice(&b);
    }
    let coarse = Coarse::new(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
    )
    .unwrap();
    let mut x = previous;
    coarse.sweep(&mut x, Some(&source), inverse_k);
    assert_eq!(x, expected);
}

#[test]
fn blocks_must_have_a_square_of_the_groups() {
    // 成分の足りないブロックを平らな配列に並べると、後ろのセルのブロックが黙ってずれる。組み立てた側の番号を付けたエラーにする
    for (groups, matrices) in [
        (2, vec![vec![1.0, 0.0, 0.0, 1.0], vec![1.0, 0.0, 1.0]]),
        (0, vec![Vec::new(), Vec::new()]),
    ] {
        let result = Coarse::new(
            BlockSolve::Lu,
            groups,
            vec![3, 4],
            vec![1.0; 2 * groups],
            vec![Vec::new(); 2],
            matrices,
            vec![Vec::new(); 2],
        );
        let message = result.map(|_| ()).unwrap_err().to_string();
        assert!(message.contains("entries"), "{message}");
    }
}

#[test]
fn a_weight_error_names_the_first_bad_weight_only() {
    // 全炉心の大きさの重みを並べると、エラーの文が長くなりすぎる
    let (groups, labels, flux, fission, matrices, neighbors) = mirrored_slab_parts();
    let message = Coarse::new(
        BlockSolve::Lu,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
    )
    .unwrap()
    .weighted(vec![2.0, -1.0])
    .map(|_| ())
    .unwrap_err()
    .to_string();
    assert!(message.contains("weight 1 is -1,"), "{message}");
}

#[test]
fn an_inverse_block_solve_multiplies_by_the_inverse_from_the_lu_factors() {
    // `BlockSolve::Inverse` は、LU 分解で単位ベクトルを解いて作った逆行列を掛ける（core の形。ADR-0275）
    let (groups, labels, flux, fission, matrices, neighbors) = plane_parts(3, 2);
    let ng = groups;
    let previous: Vec<f64> = (0..6 * ng).map(|j| 0.7 + 0.13 * (j % 5) as f64).collect();
    let source: Vec<f64> = (0..6 * ng).map(|j| 0.01 * (1 + j % 3) as f64).collect();
    let inverse_k = 1.0 / 1.07;
    let mut expected = previous.clone();
    for c in 0..6 {
        let lu = Lu::new(matrices[c].clone(), ng).unwrap();
        let mut inverse = vec![0.0; ng * ng];
        for j in 0..ng {
            let mut column = vec![0.0; ng];
            column[j] = 1.0;
            lu.solve(&mut column);
            for (i, v) in column.iter().enumerate() {
                inverse[i * ng + j] = *v;
            }
        }
        let mut b = source[c * ng..][..ng].to_vec();
        for s in &fission[c] {
            let rate = s.rate(&previous[c * ng..][..ng]) * inverse_k;
            for (r, x) in b.iter_mut().zip(&s.chi) {
                *r += x * rate;
            }
        }
        for (other, coefficients) in &neighbors[c] {
            for g in 0..ng {
                b[g] -= coefficients[g] * expected[other * ng + g];
            }
        }
        for i in 0..ng {
            expected[c * ng + i] = inverse[i * ng..][..ng]
                .iter()
                .zip(&b)
                .map(|(a, x)| a * x)
                .sum();
        }
    }
    let coarse = Coarse::new(
        BlockSolve::Inverse,
        groups,
        labels,
        flux,
        fission,
        matrices,
        neighbors,
    )
    .unwrap();
    let mut x = previous;
    coarse.sweep(&mut x, Some(&source), inverse_k);
    assert_eq!(x, expected);
}
