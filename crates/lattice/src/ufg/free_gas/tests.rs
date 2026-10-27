use super::*;
use xlc_numeric::quadrature::GAUSS_LEGENDRE_8;

const KT: f64 = BOLTZMANN * 600.0;

/// 格子点 `start` からの中性子束 `flux`（全体の添字）の累積の積率
fn prefix_of(u: &[f64], flux: &[f64], start: usize) -> Prefix {
    let mut prefix = Prefix::new(start);
    for m in start..flux.len() {
        prefix.push(u, flux, m);
    }
    prefix
}

fn assert_rel(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs(),
        "actual {actual}, expected {expected}"
    );
}

/// [a, b] を `pieces` 等分した 8 点の Gauss-Legendre 則
fn integrate_fn(f: impl Fn(f64) -> f64, a: f64, b: f64, pieces: usize) -> f64 {
    let h = (b - a) / pieces as f64;
    (0..pieces)
        .map(|k| {
            let (mid, half) = (a + (k as f64 + 0.5) * h, 0.5 * h);
            GAUSS_LEGENDRE_8
                .iter()
                .map(|&(x, w)| w * half * f(mid + half * x))
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn kernel_is_the_integral_of_the_free_gas_law() {
    // S(E→E') = ((A+1)/A)² A/(4E) ∫ e^{−(α+β)²/(4α)}/√(4πα) dα（α は (√E ∓ √E')²/(A kT) の間、β = (E' − E)/kT）
    for (awr, e, e_out) in [
        (0.99917, 6.0, 5.0),
        (0.99917, 6.0, 6.05),
        (15.86, 6.0, 4.9),
        (15.86, 6.0, 6.08),
        (9.93, 100.0, 80.0),
        (236.0, 6.67, 6.6),
        (236.0, 6.67, 6.75),
    ] {
        let beta = (e_out - e) / KT;
        let (lo, hi) = (
            (e.sqrt() - f64::sqrt(e_out)).powi(2) / (awr * KT),
            (e.sqrt() + f64::sqrt(e_out)).powi(2) / (awr * KT),
        );
        // v = √α で積分する（被積分関数はガウス形）
        let integral = integrate_fn(
            |v| {
                let alpha = v * v;
                2.0 * v * (-(alpha + beta).powi(2) / (4.0 * alpha)).exp()
                    / (4.0 * PI * alpha).sqrt()
            },
            lo.sqrt(),
            hi.sqrt(),
            4000,
        );
        let expected = ((awr + 1.0) / awr).powi(2) * awr / (4.0 * e) * integral;
        assert_rel(kernel(e, e_out, awr, KT), expected, 1e-9);
    }
}

#[test]
fn hydrogen_kernel_has_the_wigner_wilkins_form() {
    // A = 1: 下方散乱は erf(√(E'/kT))/E（静止した標的の 1/E に近い）、上方散乱は e^{(E−E')/kT} erf(√(E/kT))/E
    let e = 6.0;
    for e_out in [3.0, 5.99, 6.0] {
        assert_rel(
            kernel(e, e_out, 1.0, KT),
            erf((e_out / KT).sqrt()) / e,
            1e-12,
        );
    }
    for e_out in [6.01, 6.05, 6.2] {
        assert_rel(
            kernel(e, e_out, 1.0, KT),
            ((e - e_out) / KT).exp() * erf((e / KT).sqrt()) / e,
            1e-10,
        );
    }
}

#[test]
fn kernel_integrates_to_the_effective_ratio() {
    for awr in [0.99917_f64, 15.86, 236.0] {
        for e in [2.0 * KT, 50.0 * KT, 1000.0 * KT] {
            let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
            let reach = e + 60.0 * KT + 30.0 * (e * KT / awr).sqrt();
            // v = √E' で積分する（E' → 0 の √E' の振る舞いを滑らかにする）
            let points = [0.0, alpha * e, e, reach].map(f64::sqrt);
            let total: f64 = points
                .windows(2)
                .map(|w| integrate_fn(|v| 2.0 * v * kernel(e, v * v, awr, KT), w[0], w[1], 2000))
                .sum();
            assert_rel(total, effective_ratio(e, awr, KT), 1e-8);
        }
    }
}

#[test]
fn kernel_approaches_the_target_at_rest_far_above_kt() {
    // T → 0: 規格化した核は [αE, E] で一様な 1/((1 − α)E)
    let kt = 1.0e-9;
    for awr in [0.99917_f64, 15.86] {
        let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
        let e = 10.0;
        let at_rest = 1.0 / ((1.0 - alpha) * e);
        for x in [0.1, 0.5, 0.9] {
            let e_out = e * (alpha + (1.0 - alpha) * x);
            assert_rel(
                kernel(e, e_out, awr, kt) / effective_ratio(e, awr, kt),
                at_rest,
                1e-6,
            );
        }
        assert!(kernel(e, 1.001 * e, awr, kt) < 1e-12 * at_rest);
        if alpha > 0.1 {
            assert!(kernel(e, 0.999 * alpha * e, awr, kt) < 1e-12 * at_rest);
        }
    }
}

#[test]
fn correction_vanishes_outside_its_cells() {
    // 区間の外の u' では D が静止した核の 1e−11 倍より小さい（TAIL と EDGE の選び方。1200 K でも範囲の下限で）
    let high: f64 = 2.0e4;
    let (mut bounds, mut cells) = (Vec::new(), Vec::new());
    for (awr, temperature) in [
        (0.99917, 600.0),
        (9.93, 600.0),
        (15.86, 600.0),
        (15.86, 1200.0),
        (0.99917, 1200.0),
        (89.13, 600.0),
        (236.0, 565.0),
        (236.0, 900.0),
    ] {
        let gas = Gas::new(awr, temperature, None);
        let grades = grades(4);
        for e in [4.0, 4.5, 50.0, 5000.0] {
            let u = (high / e).ln();
            gas.cells(high, u, e, &grades, 4, 9.0, &mut bounds, &mut cells);
            let width = gas.width.min(3.0);
            let samples = 20_000;
            for j in 0..=samples {
                let x = u - width - 0.5 + (width + 2.0) * j as f64 / samples as f64;
                if x < 0.0 || cells.iter().any(|&(a, b)| a <= x && x <= b) {
                    continue;
                }
                let inside = u - gas.width < x && x < u;
                let d = gas.correction(high, x, u, e, inside);
                assert!(d.abs() < 1e-11 * gas.gain, "A {awr}, E {e}: D({x}) = {d}");
            }
        }
    }
}

/// 試験の格子: 範囲 [0, 3]（レサジー）に、粗いところと細かいところのある格子
fn test_grid() -> Vec<f64> {
    let mut u: Vec<f64> = (0..=300).map(|k| 3.0 * k as f64 / 300.0).collect();
    // 1.2〜1.3 に細かい点（共鳴のくぼみの格子の代わり）
    u.extend((1..2000).map(|k| 1.2 + 0.1 * k as f64 / 2000.0));
    u.sort_by(f64::total_cmp);
    u.dedup();
    u
}

/// 構造のある中性子束: くぼみと振動
fn test_flux(u: f64) -> f64 {
    1.0 + 0.2 * (7.0 * u).sin() - 0.6 * (-((u - 1.25) / 0.01).powi(2)).exp()
}

/// 格子点 m の D の散乱源を、区間の機構（[`integrate`]）と、細かい直接の求積で求める（弾性散乱断面積は一定）
fn sources(gas: &Gas, cells: usize, u: &[f64], m: usize, high: f64) -> (f64, f64) {
    sources_with(gas, cells, u, m, high, &vec![3.8; u.len()])
}

/// [`sources`] の、格子点ごとの弾性散乱断面積 `sigma` のもの
fn sources_with(
    gas: &Gas,
    cells: usize,
    u: &[f64],
    m: usize,
    high: f64,
    sigma: &[f64],
) -> (f64, f64) {
    let flux: Vec<f64> = u.iter().map(|&x| test_flux(x)).collect();
    let prefix = prefix_of(u, &flux, 0);
    let history = History {
        flux: &flux[..m],
        prefix: &prefix,
        previous: Some((&flux, &prefix)),
    };
    let e = high * (-u[m]).exp();
    let (mut bounds, mut pieces) = (Vec::new(), Vec::new());
    let grades = grades(cells);
    let last = u[u.len() - 1];
    gas.cells(
        high,
        u[m],
        e,
        &grades,
        cells,
        last,
        &mut bounds,
        &mut pieces,
    );
    let (mut machine, mut direct) = (0.0, 0.0);
    for &(a, b) in &pieces {
        let inside = gas.inside(u[m], a, b);
        let f = |x: f64| {
            interpolate(u, sigma, locate(u, x), x) * gas.correction(high, x, u[m], e, inside)
        };
        let mid = 0.5 * (a + b);
        let q = Quadratic::through(a, b, f(a), f(mid), f(b));
        let (known, own) = integrate(&q, u, m, a, b, locate(u, a), locate(u, b), &history);
        machine += known + own * flux[m];
        direct += integrate_fn(|x| f(x) * interpolate(u, &flux, locate(u, x), x), a, b, 64);
    }
    (machine, direct)
}

#[test]
fn cells_integrate_the_correction_against_the_flux() {
    // 区間の 2 次式と累積の積率による D の散乱源は、細かい直接の求積に、区間を細かくするほど近づく
    let u = test_grid();
    let high = 20.0;
    for awr in [0.99917, 15.86] {
        let gas = Gas::new(awr, 600.0, None);
        for target in [0.6, 1.25, 1.28, 2.9] {
            let m = locate(&u, target);
            let errors: Vec<f64> = [2, 4, 8, 16]
                .iter()
                .map(|&cells| {
                    let (machine, direct) = sources(&gas, cells, &u, m, high);
                    (machine - direct).abs()
                })
                .collect();
            let (_, scale) = sources(&gas, 16, &u, m, high);
            assert!(
                errors[3] <= 1e-6 * scale.abs().max(1e-6),
                "A {awr}, u {target}: {errors:?} of {scale}"
            );
            assert!(errors[3] <= errors[0], "A {awr}, u {target}: {errors:?}");
        }
    }
}

#[test]
fn cells_integrate_the_correction_against_a_resonant_cross_section() {
    // 分離共鳴を持つ核種（U-238、900 K）の弾性散乱断面積が区間の中で共鳴の形に変わっても、区間の 2 次式による D の散乱源は、
    // 細かい直接の求積に、区間を細かくするほど近づく。共鳴は u = 1.25（E = 5.73 eV）で、幅は 900 K のドップラー幅ほど。
    // 2 次式は σs と D の積を近似するので、σs が一定のとき（`cells_integrate_the_correction_against_the_flux`）より遅く、
    // 区間の細かさを倍にするごとに誤差は 3 分の 1 より小さくなる
    let u = test_grid();
    let high = 20.0;
    let sigma: Vec<f64> = u
        .iter()
        .map(|&x| 9.0 + 2.0e3 / (1.0 + ((x - 1.25) / 0.01).powi(2)))
        .collect();
    let gas = Gas::new(236.0, 900.0, None);
    for target in [1.22, 1.25, 1.26, 1.29] {
        let m = locate(&u, target);
        let errors: Vec<f64> = [2, 4, 8, 16]
            .iter()
            .map(|&cells| {
                let (machine, direct) = sources_with(&gas, cells, &u, m, high, &sigma);
                (machine - direct).abs()
            })
            .collect();
        let (_, scale) = sources_with(&gas, 16, &u, m, high, &sigma);
        assert!(
            errors[3] <= 1e-5 * scale.abs(),
            "u {target}: {errors:?} of {scale}"
        );
        for w in errors.windows(2) {
            assert!(3.0 * w[1] <= w[0], "u {target}: {errors:?}");
        }
    }
}

#[test]
fn first_sweep_takes_the_flux_ahead_as_its_own() {
    // 前のスイープのないとき、m より後の中性子束は F_m とみなし、その分は F_m の係数に入る
    let u = test_grid();
    let flux = vec![1.0; u.len()];
    let prefix = prefix_of(&u, &flux, 0);
    let m = locate(&u, 1.0);
    let q = Quadratic::through(0.9, 1.1, 1.0, 2.0, 0.5);
    let first = History {
        flux: &flux[..m],
        prefix: &prefix,
        previous: None,
    };
    let later = History {
        flux: &flux[..m],
        prefix: &prefix,
        previous: Some((&flux, &prefix)),
    };
    let (qa, qb) = (locate(&u, 0.9), locate(&u, 1.1));
    let (k1, o1) = integrate(&q, &u, m, 0.9, 1.1, qa, qb, &first);
    let (k2, o2) = integrate(&q, &u, m, 0.9, 1.1, qa, qb, &later);
    // 一様な中性子束ではどちらも ∫ f
    let total = q.integral(0.9, 1.1);
    assert_rel(k1 + o1, total, 1e-12);
    assert_rel(k2 + o2, total, 1e-12);
    assert!(o1 > o2);
}

#[test]
fn idle_gases_have_no_cells() {
    // `idle` の出射では、上限より下の入射の区間が無い（区間を作らずに飛ばしても散乱源は変わらない）
    let high: f64 = 2.0e4;
    let (mut bounds, mut cells) = (Vec::new(), Vec::new());
    let grades = grades(4);
    let mut idle = 0;
    for (awr, temperature) in [
        (0.99917, 600.0),
        (15.86, 600.0),
        (89.13, 600.0),
        (236.0, 900.0),
    ] {
        let gas = Gas::new(awr, temperature, Some(400.0));
        for k in 0..400 {
            let e = 4.0 * 1.02_f64.powi(k);
            let u = (high / e).ln();
            gas.cells(high, u, e, &grades, 4, 9.0, &mut bounds, &mut cells);
            if gas.idle(high, u, e) {
                idle += 1;
                assert!(cells.is_empty(), "A {awr}, E {e}: {} cells", cells.len());
            }
        }
    }
    assert!(idle > 0);
}

#[test]
fn correction_shares_split_the_kernel_at_an_outgoing_energy() {
    // 入射 E' の D（エネルギーあたり (S − S_static)/N）を出射のエネルギー Ex の下と上に分けた積分は、核と静止した核をそれぞれ
    // 直接積分した値の差に合い、2 つの和は 1 − 1/N（∫ S dE = N、∫ S_static dE = 1）。区間の細かさは移行を数えるときの
    // [`SHARE_CELLS`]
    for (awr, temperature) in [
        (0.99917, 600.0),
        (15.86, 600.0),
        (15.86, 1200.0),
        (236.0, 900.0),
    ] {
        let gas = Gas::new(awr, temperature, None);
        let kt = gas.kt;
        let alpha = 1.0 - 1.0 / gas.gain;
        let grades = grades(SHARE_CELLS);
        for e in [4.0, 6.0, 30.0] {
            let n = effective_ratio(e, awr, kt);
            let edge = alpha * e;
            let mut outgoing = vec![
                edge - 0.01,
                edge + 0.005,
                0.5 * (edge + e),
                e - 0.3,
                e - 0.02,
                e + 0.01,
                e + 0.2,
                e + 5.0,
                2.0,
            ];
            outgoing.retain(|&x| x > 0.0);
            // 核の累積 ∫_0^Ex S dE を、αE' と E' と Ex の点で区切り、E' のまわり ± 50 kT を細かく積分する
            let mut breaks = vec![0.0, edge, e, (e - 50.0 * kt).max(0.0), e + 50.0 * kt];
            breaks.extend(&outgoing);
            breaks.retain(|&x| x >= 0.0);
            breaks.sort_by(f64::total_cmp);
            breaks.dedup();
            let mut cumulative = vec![0.0];
            for w in breaks.windows(2) {
                let (a, b) = (w[0], w[1]);
                let near = b > e - 50.0 * kt && a < e + 50.0 * kt;
                let pieces = if near { 20_000 } else { 2_000 };
                let piece = integrate_fn(|x| kernel(e, x, awr, kt), a, b, pieces);
                cumulative.push(cumulative.last().unwrap() + piece);
            }
            let total = *cumulative.last().unwrap();
            assert_rel(total, n, 1e-10);
            for &x in &outgoing {
                let below_kernel = cumulative[breaks.iter().position(|&b| b == x).unwrap()];
                let below_static = gas.gain / e * (x.clamp(edge, e) - edge);
                let expected = (below_kernel - below_static) / n;
                let [below, above] = gas.shares(e, x, &grades, SHARE_CELLS);
                assert!(
                    (below - expected).abs() < 1e-10,
                    "A {awr}, T {temperature}, E' {e}, Ex {x}: {below} vs {expected}"
                );
                assert!(
                    (below + above - (1.0 - 1.0 / n)).abs() < 1e-10,
                    "A {awr}, T {temperature}, E' {e}, Ex {x}: {below} + {above} vs {}",
                    1.0 - 1.0 / n
                );
            }
        }
    }
}

#[test]
fn locate_from_finds_the_point_of_locate_from_any_hint() {
    // 重なった格子点（不連続点の両側の極限）を含む昇順の格子で、起点によらず locate と同じ格子点
    let u = [0.0, 0.1, 0.2, 0.2, 0.35, 0.5, 0.5, 0.5, 0.9, 1.4];
    let xs = [-0.1, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 1.4, 2.0, f64::NAN];
    for &x in &xs {
        for hint in 0..=u.len() {
            assert_eq!(
                locate_from(&u, x, hint),
                locate(&u, x),
                "x = {x}, hint = {hint}"
            );
        }
    }
}

/// erfcx(x) の多倍長の値（x, hi, lo。hi + lo が値）。mpmath 1.4.1 の 60 桁で erfc(x)·exp(x²) を求め、hi に倍精度へ丸めた値、
/// lo に残りを丸めた値を置いた。x = 1e200 は mpmath の erfc が扱えないので、漸近展開の 2 項 (1 − 1/(2x²))/(√π x)（次の項は
/// 1e−800 の相対の大きさ）。x は区間の端と、2.5 ≤ x の一様な乱数の 10 万点で誤差が大きかった点を含む
const ERFCX_REFERENCE: &[(f64, f64, f64)] = &[
    (2.5, 0.2108063640611436, -5.6277259093102524e-18),
    (
        2.5124651723206233,
        0.20988347101234978,
        1.3445554062764545e-17,
    ),
    (
        2.5846515745793606,
        0.20468671679751543,
        -4.829544975244344e-18,
    ),
    (
        2.6050776195927656,
        0.20326047815980136,
        -2.2696291215348346e-18,
    ),
    (2.75, 0.1936620962790687, -1.2015846532739174e-17),
    (3.0, 0.17900115118138996, -5.4272175920200274e-18),
    (3.5, 0.1552936556088943, -1.355844542216092e-18),
    (
        3.9999999999999996,
        0.1369994576250614,
        -6.177841842557427e-18,
    ),
    (4.0, 0.13699945762506138, 7.196568139158719e-18),
    (4.000000000000001, 0.13699945762506135, 6.18981248696211e-18),
    (
        4.44892545588937,
        0.12382782538384467,
        6.6209316208332245e-18,
    ),
    (
        4.463089702624895,
        0.12345249057391196,
        1.6302355169082588e-18,
    ),
    (5.0, 0.11070463773306863, -1.832347493639739e-18),
    (6.0, 0.09277656780053835, 6.215364755528485e-18),
    (8.0, 0.06998516620088092, 3.2863406596468746e-18),
    (
        9.57370480820499,
        0.05861480779876944,
        2.9943406702944013e-18,
    ),
    (10.0, 0.05614099274382259, -1.6720611399896374e-18),
    (20.0, 0.02817434874105132, -4.829674968674034e-20),
    (
        36.82194217585145,
        0.015316458700546182,
        5.269510661225192e-20,
    ),
    (50.0, 0.011281536265323773, -8.629279370172108e-19),
    (100.0, 0.005641613782989433, -1.186182619499244e-19),
    (
        602.8467785386605,
        0.0009358743007719356,
        4.2957790157627714e-20,
    ),
    (1000.0, 0.0005641893014533876, 2.3118263082589717e-20),
    (10000.0, 5.641895807268084e-05, -1.484202149398048e-21),
    (1000000.0, 5.641895835474742e-07, 1.0958100640939896e-23),
    (100000000.0, 5.641895835477562e-09, 3.2112571487601016e-25),
    (
        1000000000000.0,
        5.641895835477563e-13,
        3.776500696651908e-30,
    ),
    (1e+100, 5.6418958354775625e-101, 2.76431524922385e-117),
    (1e200, 5.641895835477563e-201, 3.5981468851464104e-217),
];

#[test]
fn erfcx_is_within_two_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERFCX_REFERENCE {
        let ulp = f64::from_bits(hi.to_bits() + 1) - hi;
        let error = ((erfcx(x) - hi) - lo) / ulp;
        assert!(error.abs() <= 2.0, "x={x}: {error} ulp");
    }
}

/// erfcx(x) のビット。xs の `special::tests` の同じ名前のテストに同じ表を置き、2 つのクレートの erfcx が同じ値を返すことを確かめる
/// （2 つは同じ計算の複製。ADR-0203）
const ERFCX_BITS: &[(f64, u64)] = &[
    (0.0, 0x3ff0000000000000),
    (0.5, 0x3fe3b3bc3c98b0f1),
    (1.0, 0x3fdb5d8780f956ab),
    (2.0, 0x3fd058671b52c866),
    (2.4999999999999996, 0x3fcafbb3f3b72f80),
    (2.5, 0x3fcafbb3f3b7343b),
    (2.5124651723206233, 0x3fcadd7629fc0866),
    (2.5846515745793606, 0x3fca332ca14914a4),
    (2.6050776195927656, 0x3fca047079220171),
    (2.75, 0x3fc8c9eb68ff27d6),
    (3.0, 0x3fc6e9827d229d2d),
    (3.5, 0x3fc3e0a99a0ee914),
    (3.9999999999999996, 0x3fc18932bf08e155),
    (4.0, 0x3fc18932bf08e154),
    (4.000000000000001, 0x3fc18932bf08e153),
    (4.44892545588937, 0x3fbfb32e2c5bc0da),
    (4.463089702624895, 0x3fbf9a95199fec59),
    (5.0, 0x3fbc57239e943d1a),
    (6.0, 0x3fb7c0348489d721),
    (8.0, 0x3fb1ea8c4009b459),
    (9.57370480820499, 0x3fae02c2951d3468),
    (10.0, 0x3facbe831f997123),
    (20.0, 0x3f9cd9bc89b73548),
    (36.82194217585145, 0x3f8f5e3c49accfb1),
    (50.0, 0x3f871ac62a773371),
    (100.0, 0x3f771ba92b1d1777),
    (602.8467785386605, 0x3f4eaaaec1eaf0e6),
    (1000.0, 0x3f427cc34c9b7f2d),
    (10000.0, 0x3f0d946ca36c60e7),
    (1000000.0, 0x3ea2ee5a03c7620f),
    (100000000.0, 0x3e383b547fb79514),
    (1000000000000.0, 0x3d63d9c49b6a903e),
    (1e100, 0x2b1f974cd8b00f7d),
    (1e200, 0x165ba394ce53f796),
];

#[test]
fn erfcx_gives_the_same_bits_as_its_copy() {
    for &(x, bits) in ERFCX_BITS {
        assert_eq!(erfcx(x).to_bits(), bits, "x={x}");
    }
}
