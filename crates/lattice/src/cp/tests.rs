use super::*;
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

#[test]
fn bickley_function_has_known_values() {
    // Ki3(0) = π/4、Ki3 は単調に減る
    assert!(close(ki3(0.0), PI / 4.0, 1e-12));
    // Ki3(x) = ∫_0^{π/2} cos²θ exp(−x/cosθ) dθ（SciPy の適応求積で求めた値）
    assert!(close(ki3(1.0), 0.237_845_082_192_855, 1e-8), "{}", ki3(1.0));
    assert!(close(ki3(2.0), 0.076_963_590_311_658, 1e-8), "{}", ki3(2.0));
    assert!(ki3(2.0) < ki3(1.0) && ki3(20.0) > 0.0);
}

/// 燃料・ギャップ・被覆管・水のセル
const RADII: [f64; 4] = [0.4096, 0.418, 0.475, 0.711];

#[test]
fn probabilities_are_conserved_and_reciprocal() {
    let totals = [0.6, 0.0, 0.3, 1.2];
    let p = Cylinder::new(RADII.to_vec())
        .unwrap()
        .probabilities(&totals)
        .unwrap();
    let volumes = Cylinder::new(RADII.to_vec()).unwrap().volumes();
    for i in [0, 2, 3] {
        // 白色境界では、どこかで必ず衝突する
        let sum: f64 = (0..4).map(|j| p[i][j]).sum();
        assert!(close(sum, 1.0, 1e-8), "row {i}: {sum}");
        for j in [0, 2, 3] {
            let (a, b) = (
                volumes[i] * totals[i] * p[i][j],
                volumes[j] * totals[j] * p[j][i],
            );
            assert!(close(a, b, 1e-8), "{i},{j}: {a} vs {b}");
        }
        // 真空では衝突しない
        assert_eq!(p[i][1], 0.0);
    }
}

#[test]
fn homogeneous_cell_is_an_infinite_medium() {
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let flux = cell.fixed_source(&[0.8; 4], &[2.4; 4]).unwrap();
    for f in &flux {
        assert!(close(*f, 3.0, 1e-8), "{f}");
    }
}

#[test]
fn thin_rim_conserves_neutrons() {
    // 広いリングの外の薄いリングは、弦の長さが y について急に変わるので求積の誤差が大きい。
    // その誤差によらず、どこかで必ず衝突し、均質なら無限媒質の中性子束になる
    let cell = Cylinder::new(vec![0.4, 0.4001, 0.5]).unwrap();
    let p = cell.probabilities(&[0.5; 3]).unwrap();
    for (i, row) in p.iter().enumerate() {
        let sum: f64 = row.iter().sum();
        assert!(close(sum, 1.0, 1e-12), "row {i}: {sum}");
    }
    for f in cell.fixed_source(&[0.5; 3], &[1.0; 3]).unwrap() {
        assert!(close(f, 2.0, 1e-12), "{f}");
    }
}

#[test]
fn void_ring_carries_the_surrounding_flux() {
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let flux = cell
        .fixed_source(&[0.8, 0.0, 0.8, 0.8], &[2.4, 0.0, 2.4, 2.4])
        .unwrap();
    // 真空の中では角度中性子束が弦に沿って変わらないので、真空のリングの中性子束も周りと同じになる
    for (i, f) in flux.iter().enumerate() {
        assert!(close(*f, 3.0, 1e-8), "{i}: {f}");
    }
}

#[test]
fn coarse_trapezoid_stays_within_its_bound() {
    // 刻み 0.25 を使う範囲（B < 4）の端の近くでも、刻み 0.125 との相対差が見積もり（3e-14）以下
    for b in [0.0, 1.0, 3.0, 3.99] {
        for (ta, tb) in [(1e-4, 1e-4), (1e-3, 0.5), (1e-6, 3.0)] {
            let f = |c: f64| (-b * c).exp() * phi1(ta * c) * phi1(tb * c) / c;
            let (coarse, fine) = (trapezoid(coarse_stride(b), f), trapezoid(1, f));
            assert!(
                (coarse - fine).abs() <= 3e-14 * fine,
                "B {b} τ ({ta}, {tb}): {coarse} vs {fine}"
            );
        }
    }
}

#[test]
fn quadrature_converges() {
    // 黒に近い燃料を細かいリングに分け、水の線源から入る中性子束をみる。2つ目は、広い燃料の外の薄いリム
    let cells = [
        (
            vec![0.1, 0.2, 0.3, 0.4096, 0.418, 0.475, 0.711],
            vec![30.0, 30.0, 30.0, 30.0, 0.0, 0.3, 1.2],
        ),
        (vec![0.4, 0.4001, 0.5], vec![20.0, 40.0, 1.2]),
    ];
    for (radii, totals) in cells {
        let mut sources = vec![0.0; radii.len()];
        sources[radii.len() - 1] = 1.0;
        let flux = |points| {
            Cylinder::with_points(radii.clone(), points)
                .unwrap()
                .fixed_source(&totals, &sources)
                .unwrap()
        };
        let (coarse, fine) = (flux(POINTS), flux(4 * POINTS));
        for (i, (c, f)) in coarse.iter().zip(&fine).enumerate() {
            assert!((c - f).abs() <= 1e-9 * f, "{radii:?} {i}: {c} vs {f}");
        }
    }
}

#[test]
fn bickley_function_vanishes_far_away() {
    // exp(−x) が 0 に丸められる距離でも止まる
    assert_eq!(ki3(750.0), 0.0);
    assert!(ki3(60.0) > 0.0 && ki3(60.0) < 1e-26);
}

#[test]
fn thin_gas_approaches_void() {
    // ギャップの Σ を 0 に近づけると、真空のときの中性子束に連続につながる
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let flux =
        |gap: f64, sources: &[f64]| cell.fixed_source(&[0.6, gap, 0.3, 1.2], sources).unwrap();
    for sources in [[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]] {
        let void = flux(0.0, &sources);
        for gap in [1e-6, 1e-10] {
            for (i, (f, v)) in flux(gap, &sources).iter().zip(&void).enumerate() {
                assert!((f - v).abs() <= 1e-5 * v, "{gap} {i}: {f} vs {v}");
            }
        }
    }
}

#[test]
fn rejects_bad_radii() {
    assert!(Cylinder::new(vec![]).is_err());
    assert!(Cylinder::new(vec![0.5, 0.4]).is_err());
    assert!(Cylinder::new(vec![0.0, 0.4]).is_err());
}
