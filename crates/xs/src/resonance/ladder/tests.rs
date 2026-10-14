use super::*;

/// n 個の標本の平均と分散
fn moments(n: usize, mut draw: impl FnMut() -> f64) -> (f64, f64) {
    let values: Vec<f64> = (0..n).map(|_| draw()).collect();
    let mean = values.iter().sum::<f64>() / n as f64;
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    (mean, variance)
}

#[test]
fn uniform_draws_are_flat_in_the_open_unit_interval() {
    let mut rng = Rng::new(1);
    let (mean, variance) = moments(200_000, || {
        let u = rng.uniform();
        assert!(u > 0.0 && u < 1.0);
        u
    });
    assert!((mean - 0.5).abs() < 3e-3, "{mean}");
    assert!((variance - 1.0 / 12.0).abs() < 1e-3, "{variance}");
}

#[test]
fn widths_follow_chi_square_with_mean_one() {
    // 自由度 ν の χ²/ν: 平均 1、分散 2/ν。自由度なし（Fixed）はゆらがない
    let mut rng = Rng::new(2);
    for (dof, nu) in [
        (WidthDof::One, 1.0),
        (WidthDof::Two, 2.0),
        (WidthDof::Three, 3.0),
        (WidthDof::Four, 4.0),
    ] {
        let (mean, variance) = moments(400_000, || chi_square(dof, &mut rng));
        assert!((mean - 1.0).abs() < 1e-2, "{dof:?}: {mean}");
        assert!(
            (variance - 2.0 / nu).abs() < 3e-2 * 2.0 / nu,
            "{dof:?}: {variance}"
        );
    }
    assert_eq!(chi_square(WidthDof::Fixed, &mut rng), 1.0);
}

#[test]
fn wigner_spacings_have_the_mean_spacing() {
    // Wigner の分布 p(s) = (π s / 2D²) exp(−π s²/(4D²)): 平均 D、2 次の積率 4D²/π
    let mut rng = Rng::new(3);
    let d = 3.0;
    let (mean, variance) = moments(400_000, || wigner_spacing(d, &mut rng));
    assert!((mean - d).abs() < 1e-2 * d, "{mean}");
    let second = variance + mean * mean;
    assert!(
        (second - 4.0 * d * d / std::f64::consts::PI).abs() < 1e-2 * second,
        "{second}"
    );
}

#[test]
fn the_same_seed_gives_the_same_draws() {
    let (mut a, mut b) = (Rng::new(42), Rng::new(42));
    for _ in 0..100 {
        assert_eq!(a.uniform(), b.uniform());
    }
    assert_ne!(Rng::new(42).uniform(), Rng::new(43).uniform());
}
