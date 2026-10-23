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

/// 系列 `q` のラダー `ladder` のエネルギー `energies` での [弾性散乱, 核分裂, 捕獲] を、エネルギーと共鳴の対ごとに 1 つずつ
/// 求める（`add_ladder` の式をそのまま点ごとに書いたもの）
fn pointwise(q: &LadderSequence, ladder: &[Resonance], energies: &[f64]) -> Vec<[f64; 3]> {
    use crate::special::faddeeva;
    let a = Complex64::new(1.0 - q.cos2, q.sin2p);
    let pole = |r: &Resonance| Complex64::new(r.energy, -0.5 * r.gamma);
    let root_pi = PI.sqrt();
    let elastic_scale = 0.5 * q.scale * root_pi / q.delta;
    let mut parts = vec![[0.0; 3]; energies.len()];
    for r in ladder {
        let p = pole(r);
        let c: Complex64 = ladder.iter().map(|s| s.gn / (pole(s).conj() - p)).sum();
        let ar = r.gn * (c + Complex64::i() * a);
        let y = 0.5 * r.gamma / q.delta;
        let c = q.scale * r.gn / r.gamma * root_pi * y / r.gamma;
        let (fission_re, capture_re) = (c * r.gf, c * q.gg);
        let iar = Complex64::i() * ar * elastic_scale;
        for (part, &e) in parts.iter_mut().zip(energies) {
            let x = (e - r.energy) / q.delta;
            let w = if x.abs() > ASYMPTOTIC || y > ASYMPTOTIC {
                Complex64::i() / (root_pi * Complex64::new(x, y))
            } else {
                faddeeva(Complex64::new(x, y))
            };
            part[0] += (iar * w).re;
            part[1] += fission_re * w.re;
            part[2] += capture_re * w.re;
        }
    }
    parts
}

#[test]
fn ladder_sums_give_the_bits_of_the_pointwise_formula() {
    // 計算の並べ方（ベクトル化）を変えても、点ごとの演算と足す順序が同じなので、ビットまで同じ（ADR-0188）。Δ = 0.2 eV では
    // 窓の中の点に Faddeeva 関数と漸近形の両方が、Δ = 1e-4 eV では y > 100 で漸近形だけが使われる
    for delta in [0.2, 1.0e-4] {
        let q = LadderSequence {
            scale: 1.0e3,
            gn: 0.05,
            gf: 0.1,
            gg: 0.03,
            gx: 0.0,
            dof_n: WidthDof::One,
            dof_f: WidthDof::Two,
            dof_x: WidthDof::Fixed,
            d: 1.0,
            delta,
            cos2: 0.6_f64.cos(),
            sin2p: 0.6_f64.sin(),
        };
        let mut rng = Rng::new(5);
        let half = 25.0;
        let ladder = sample_ladder(&q, half, 1.0, &mut rng);
        let energies: Vec<f64> = (0..250).map(|i| -half + 0.2 * (i as f64 + 0.37)).collect();
        let mut parts = Parts::zeroed(energies.len());
        add_ladder(&q, &ladder, &energies, &mut parts);
        let expected = pointwise(&q, &ladder, &energies);
        // 窓の端で Faddeeva 関数と漸近形が切り替わる共鳴（Δ = 0.2 eV）と、y > 100 で漸近形だけの共鳴（Δ = 1e-4 eV）があること
        let switching = ladder
            .iter()
            .filter(|r| {
                let x = energies.iter().map(|&e| ((e - r.energy) / delta).abs());
                let near = x.clone().filter(|&x| x <= ASYMPTOTIC).count();
                0 < near && near < energies.len()
            })
            .count();
        let wide = ladder
            .iter()
            .filter(|r| 0.5 * r.gamma / delta > ASYMPTOTIC)
            .count();
        if delta > 1e-3 {
            assert!(switching > 0 && wide == 0, "{switching}, {wide}");
        } else {
            assert_eq!(wide, ladder.len());
        }
        for (j, e) in expected.iter().enumerate() {
            let got = [parts.elastic[j], parts.fission[j], parts.capture[j]];
            assert_eq!(got.map(f64::to_bits), e.map(f64::to_bits), "sample {j}");
        }
    }
}
