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
        let gas = Gas::new(awr, temperature, f64::INFINITY);
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
        let gas = Gas::new(awr, 600.0, f64::INFINITY);
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
    let gas = Gas::new(236.0, 900.0, f64::INFINITY);
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
        let gas = Gas::new(awr, temperature, 400.0);
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
        let gas = Gas::new(awr, temperature, f64::INFINITY);
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

/// [`Thermal`] で、試験の格子の 3 回のスイープを模した散乱源。格子点ごとの [リング] の（既知の部分, F_m の係数）を、スイープの順に
/// 並べたものと、最初のスイープで持った区間の数を返す。`hold` は持つ区間の数の上限。中性子束は散乱源に依るようにする（持った区間を
/// 使うスイープでも、中性子束が前のスイープと違う）
fn swept_sources(hold: usize) -> (Vec<Vec<(f64, f64)>>, usize) {
    let (out, held, _) = swept_on(&test_grid(), 1.0, hold, None);
    (out, held)
}

/// [`swept_sources`] を、格子 `u`、弾性散乱断面積に掛ける `scale`、自由ガスの核の値を持つ表 `kernels` で解いたもの。3 つ目は
/// 核の評価（[`Gas::correction`]）の回数
fn swept_on(
    u: &[f64],
    scale: f64,
    hold: usize,
    kernels: Option<&KernelCache>,
) -> (Vec<Vec<(f64, f64)>>, usize, u64) {
    let (out, held, work) = swept_mixed(u, scale, hold, kernels, None);
    (out, held, work.evaluations)
}

/// [`swept_on`] に、質量比 30 の自由ガスの核種を、核種の並びの `extra` の位置に、2 つのリングの個数密度 `extra` で足したもの。
/// 3 つ目は自由ガスの核の区間の仕事
fn swept_mixed(
    u: &[f64],
    scale: f64,
    hold: usize,
    kernels: Option<&KernelCache>,
    extra: Option<(usize, [f64; 2])>,
) -> (Vec<Vec<(f64, f64)>>, usize, FreeGasWork) {
    swept_at(u, scale, hold, kernels, extra, [600.0, 900.0])
}

/// [`swept_mixed`] を、質量比 1 と 16 の核種（減速材）の温度と質量比 236 の核種（燃料）の温度 `temperatures` (K) で解いたもの
fn swept_at(
    u: &[f64],
    scale: f64,
    hold: usize,
    kernels: Option<&KernelCache>,
    extra: Option<(usize, [f64; 2])>,
    temperatures: [f64; 2],
) -> (Vec<Vec<(f64, f64)>>, usize, FreeGasWork) {
    let [moderator, fuel] = temperatures;
    let u = u.to_vec();
    let high = 20.0;
    let energies: Vec<f64> = u.iter().map(|&x| high * (-x).exp()).collect();
    let mut nuclides = vec![
        (
            Gas::new(0.99917, moderator, f64::INFINITY),
            20.0,
            vec![(1, 0.03)],
        ),
        (
            Gas::new(15.86, moderator, f64::INFINITY),
            3.8,
            vec![(0, 0.05), (1, 0.025)],
        ),
        (
            Gas::new(236.0, fuel, 400.0),
            9.0,
            vec![(0, 0.022), (1, 0.045)],
        ),
    ];
    if let Some((position, [first, second])) = extra {
        nuclides.insert(
            position,
            (
                Gas::new(30.0, 600.0, f64::INFINITY),
                5.0,
                vec![(0, first), (1, second)],
            ),
        );
    }
    let gases: Vec<Gas> = nuclides.iter().map(|(gas, ..)| *gas).collect();
    let sigma = nuclides
        .iter()
        .map(|(gas, value, _)| {
            let first = first_below(&energies, gas.top());
            (first, vec![scale * value; u.len() - first])
        })
        .collect();
    let rings = nuclides.into_iter().map(|(.., rings)| rings).collect();
    let n = 2;
    let tables = vec![None; gases.len()];
    let start = sweep_start(&gases, &tables, &energies).expect("free-gas nuclides");
    let mut thermal = Thermal::new(gases, sigma, rings, tables, 2, start, n, hold, kernels);
    let start = thermal.start();
    let mut flux: Vec<Vec<f64>> = vec![u[..start].iter().map(|&x| test_flux(x)).collect(); n];
    let mut out = Vec::new();
    let mut held = 0;
    for sweep in 0..3 {
        for f in &mut flux {
            f.truncate(start);
        }
        for m in start..u.len() {
            let (mut known, mut own) = (vec![0.0; n], vec![0.0; n]);
            thermal
                .source(m, &u, high, energies[m], &flux, &mut known, &mut own)
                .unwrap();
            for (i, f) in flux.iter_mut().enumerate() {
                f.push(test_flux(u[m]) * (1.0 + 0.01 * sweep as f64) + 1e-3 * (known[i] + own[i]));
            }
            thermal.record(m, &u, &flux);
            out.push(known.into_iter().zip(own).collect());
        }
        if sweep == 0 {
            held = thermal.held();
            // 持つ配列の区間、核種の区間の始まりと位置の表は、確保した分も上限を超えない
            let starts = (thermal.held_starts.capacity() * 4
                + thermal.held.capacity() * std::mem::size_of::<Option<Held>>())
            .div_ceil(HELD_CELL_BYTES);
            assert!(thermal.held_cells.capacity() + starts <= hold, "{hold}");
        }
        thermal.finish(&flux);
    }
    thermal.release();
    assert_eq!(thermal.held(), 0);
    (out, held, thermal.work().clone())
}

#[test]
fn held_cells_give_the_sources_of_rebuilt_cells_bit_for_bit() {
    // スイープの間で持った区間と 2 次式を積分だけに使っても、毎回作り直したときと散乱源がビットで同じ。上限が区間の一部だけを持つ
    // ときも（上限を超える格子点は作り直す）同じ
    let (rebuilt, none) = swept_sources(0);
    assert_eq!(none, 0);
    let (all, total) = swept_sources(usize::MAX);
    assert!(total > 0);
    let (part, some) = swept_sources(total / 2);
    assert!(0 < some && some <= total / 2, "{some} of {total}");
    let bits = |v: &[Vec<(f64, f64)>]| -> Vec<(u64, u64)> {
        v.iter()
            .flatten()
            .map(|&(k, o)| (k.to_bits(), o.to_bits()))
            .collect()
    };
    assert_eq!(bits(&all), bits(&rebuilt));
    assert_eq!(bits(&part), bits(&rebuilt));
}

#[test]
fn the_cells_a_solve_would_hold_do_not_depend_on_its_hold() {
    // 減速計算が上限なしで持つ区間の数（核種の区間の始まりと位置の表を 1 区間 48 B で数えたものを含む。表の上限を決める
    // ピンの 1 回の量）は、与えた上限によらず、すべてを持てたときに持った数と同じ
    let wanted = |hold: usize| swept_mixed(&test_grid(), 1.0, hold, None, None).2.wanted;
    let (_, total) = swept_sources(usize::MAX);
    let all = wanted(usize::MAX);
    assert!(all > total as u64, "{all} {total}");
    for hold in [0, total / 2, total] {
        assert_eq!(wanted(hold), all, "{hold}");
    }
    let (.., work) = swept_mixed(&test_grid(), 1.0, usize::MAX, None, None);
    assert_eq!(work.wanted, all);
}

#[test]
fn the_census_counts_the_table_bytes_of_its_points_by_nuclide() {
    // DEBUG の計数の格子点の大きさ（表に入れたときの大きさ）は、その格子点の核種ごとの値を表に入れた大きさの和
    let subscriber = tracing_subscriber::fmt()
        .with_max_level(tracing::Level::DEBUG)
        .with_writer(std::io::sink)
        .finish();
    tracing::subscriber::with_default(subscriber, || {
        let kernels = KernelCache::default();
        kernels.begin(|_| u64::MAX);
        let (.., work) = swept_mixed(&test_grid(), 1.0, usize::MAX, Some(&kernels), None);
        let counted: u64 = work.points.iter().map(|p| p.3).sum();
        assert!(counted > 0);
        assert_eq!(counted, kernels.bytes());
    });
}

#[test]
fn held_kernels_give_the_sources_of_fresh_kernels_bit_for_bit() {
    // 自由ガスの核の区間と D の値を、断面積の違う減速計算と、格子点の一部だけが同じ格子の減速計算の間で持っても、散乱源は持たない
    // ときとビットで同じ。同じ格子の 2 本目は核を評価せず、一部が同じ格子は残りの格子点だけを評価する
    let bits = |v: &[Vec<(f64, f64)>]| -> Vec<(u64, u64)> {
        v.iter()
            .flatten()
            .map(|&(k, o)| (k.to_bits(), o.to_bits()))
            .collect()
    };
    let grid = test_grid();
    // 最初と最後の点（範囲の上端と格子の最後の点）を残して、3 点に 1 点を次の点との中点に動かした格子
    let thinned: Vec<f64> = grid
        .iter()
        .enumerate()
        .map(|(k, &x)| {
            if k == 0 || k == grid.len() - 1 || k % 3 != 1 {
                x
            } else {
                0.5 * (x + grid[k + 1])
            }
        })
        .collect();
    for hold in [0, usize::MAX] {
        let kernels = KernelCache::default();
        kernels.begin(|_| u64::MAX);
        let (first, _, fresh) = swept_on(&grid, 1.0, hold, Some(&kernels));
        let (alone, _, alone_evaluations) = swept_on(&grid, 1.0, hold, None);
        assert_eq!(bits(&first), bits(&alone));
        // 持たない入力（hold 0）では、2 回目からのスイープも表の値を使う
        let sweeps = if hold == 0 { 3 } else { 1 };
        assert_eq!(fresh * sweeps, alone_evaluations, "hold {hold}");
        let (second, _, again) = swept_on(&grid, 1.7, hold, Some(&kernels));
        let (second_alone, ..) = swept_on(&grid, 1.7, hold, None);
        assert_eq!(bits(&second), bits(&second_alone));
        assert_eq!(again, 0, "hold {hold}");
        let (partial, _, some) = swept_on(&thinned, 1.3, hold, Some(&kernels));
        let (partial_alone, _, all) = swept_on(&thinned, 1.3, hold, None);
        assert_eq!(bits(&partial), bits(&partial_alone));
        assert!(0 < some && some < all, "{some} of {all}");
    }
}

#[test]
fn the_nuclides_whose_kernels_stay_reuse_their_values_across_temperature_branches() {
    // 表は核種ごとの核で引くので、燃料の温度の分岐（減速材の核は同じ）でも減速材の温度の分岐（燃料の核は同じ）でも、核の変わらない
    // 核種の値を使い、変わった核種だけを評価する。散乱源は表を持たないときとビットで同じ
    let bits = |v: &[Vec<(f64, f64)>]| -> Vec<(u64, u64)> {
        v.iter()
            .flatten()
            .map(|&(k, o)| (k.to_bits(), o.to_bits()))
            .collect()
    };
    let grid = test_grid();
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    swept_at(&grid, 1.0, usize::MAX, Some(&kernels), None, [600.0, 900.0]);
    for temperatures in [[600.0, 1200.0], [580.0, 900.0]] {
        let (shared, _, work) =
            swept_at(&grid, 1.1, usize::MAX, Some(&kernels), None, temperatures);
        let (fresh, _, alone) = swept_at(&grid, 1.1, usize::MAX, None, None, temperatures);
        assert_eq!(bits(&shared), bits(&fresh), "{temperatures:?}");
        assert!(
            0 < work.evaluations && work.evaluations < alone.evaluations,
            "{temperatures:?}: {} of {}",
            work.evaluations,
            alone.evaluations
        );
    }
    // 両方の分岐の核種の値を持ったので、減速材と燃料の両方を変えた状態も評価しない
    let (both, _, work) = swept_at(
        &grid,
        1.2,
        usize::MAX,
        Some(&kernels),
        None,
        [580.0, 1200.0],
    );
    let (fresh, ..) = swept_at(&grid, 1.2, usize::MAX, None, None, [580.0, 1200.0]);
    assert_eq!(bits(&both), bits(&fresh));
    assert_eq!(work.evaluations, 0);
}

#[test]
fn free_gas_nuclides_without_density_spare_their_work_and_keep_the_sources_bit_for_bit() {
    // 個数密度がちょうど 0 の核種の散乱源を飛ばしても、散乱源は計算したときとビットで同じ。どのリングでも 0 の核種は核を評価
    // せず、一部のリングで 0 の核種は、そのリングの積分を飛ばす
    let bits = |v: &[Vec<(f64, f64)>]| -> Vec<(u64, u64)> {
        v.iter()
            .flatten()
            .map(|&(k, o)| (k.to_bits(), o.to_bits()))
            .collect()
    };
    let integrations = |w: &FreeGasWork| w.integrations.iter().sum::<u64>();
    let grid = test_grid();
    // 並びの最後の核種は、区間の塊の分け方を変えないので、核種を外した計算とビットで同じで、仕事も同じ
    let (without, _, alone) = swept_mixed(&grid, 1.0, usize::MAX, None, None);
    let (last, _, work) = swept_mixed(&grid, 1.0, usize::MAX, None, Some((3, [0.0, 0.0])));
    assert_eq!(bits(&last), bits(&without));
    assert_eq!(work.evaluations, alone.evaluations);
    assert_eq!(integrations(&work), integrations(&alone));
    // 並びの途中の核種は、ほかの核種の区間の塊の分け方を変えず、無視できる個数密度（1e−300。散乱源の和の桁に届かない）で計算した
    // ものとビットで同じ
    for densities in [[0.0, 0.0], [0.0, 0.01]] {
        let negligible = densities.map(|d| if d == 0.0 { 1.0e-300 } else { d });
        let (skipped, _, work) = swept_mixed(&grid, 1.0, usize::MAX, None, Some((1, densities)));
        let (computed, _, full) = swept_mixed(&grid, 1.0, usize::MAX, None, Some((1, negligible)));
        assert_eq!(bits(&skipped), bits(&computed), "{densities:?}");
        assert!(integrations(&work) < integrations(&full), "{densities:?}");
        if densities == [0.0, 0.0] {
            // 核種の評価はすべて省く。残る差は、核種の区間がずらした塊の境で、ほかの核種の区間の端を評価し直す分だけ
            let own = full.evaluations - alone.evaluations;
            assert!(work.evaluations - alone.evaluations < own / 1000, "{own}");
        } else {
            assert_eq!(work.evaluations, full.evaluations);
        }
    }
    // 核を飛ばした減速計算と、その核種に個数密度のある減速計算が同じ表を引いても、どちらも表を持たないときとビットで同じ
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    for densities in [[0.0, 0.0], [0.01, 0.02], [0.0, 0.0]] {
        let (shared, ..) =
            swept_mixed(&grid, 1.0, usize::MAX, Some(&kernels), Some((1, densities)));
        let (fresh, ..) = swept_mixed(&grid, 1.0, usize::MAX, None, Some((1, densities)));
        assert_eq!(bits(&shared), bits(&fresh), "{densities:?}");
    }
}

/// 最初と最後の点（範囲の上端と格子の最後の点）を残して、ほかの格子点を次の点との中点に動かした格子（`grid` と重なる格子点は
/// 最後の点だけ）
fn moved(grid: &[f64]) -> Vec<f64> {
    grid.iter()
        .enumerate()
        .map(|(k, &x)| {
            if k == 0 || k == grid.len() - 1 {
                x
            } else {
                0.5 * (x + grid[k + 1])
            }
        })
        .collect()
}

/// 格子 `u` の格子点をすべて持った表の大きさ (B)
fn table_bytes(u: &[f64]) -> u64 {
    let probe = KernelCache::default();
    probe.begin(|_| u64::MAX);
    swept_on(u, 1.0, usize::MAX, Some(&probe));
    probe.bytes()
}

#[test]
fn a_full_kernel_cache_adds_no_points() {
    // 上限に届いた表は格子点を足さず、足さなかった格子点は減速計算ごとに評価し直す
    let grid = test_grid();
    let kernels = KernelCache::default();
    kernels.begin(|_| 0);
    let (_, _, all) = swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    assert_eq!(kernels.bytes(), 0);
    // どの格子点も表から引けず、作った格子点はどれも上限で入れなかったと数える
    let [reused, built, full] = kernels.counts();
    assert_eq!((reused, full), (0, built));
    assert!(built > 0);
    let (_, _, again) = swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    assert_eq!(again, all);
}

#[test]
fn a_kernel_table_keeps_the_points_of_its_burnup_point_through_its_shieldings() {
    // 同じ燃焼点の中では、上限に届かない限り、いくつの自己遮蔽（状態と摂動）で使わなくても格子点を捨てない
    let grid = test_grid();
    let other = moved(&grid);
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    for _ in 0..3 {
        kernels.begin(|_| u64::MAX);
        swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    }
    kernels.begin(|_| u64::MAX);
    let (_, _, again) = swept_on(&grid, 1.3, usize::MAX, Some(&kernels));
    assert_eq!(again, 0);
}

#[test]
fn a_new_burnup_point_drops_the_points_its_previous_point_did_not_use() {
    // 燃焼点の境（呼び出し側が渡す）では、直前の燃焼点で使った格子点（燃焼点の間の重なり）を残し、ほかを捨てる
    let grid = test_grid();
    let other = moved(&grid);
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    let first = kernels.bytes();
    kernels.next_point();
    kernels.begin(|_| u64::MAX);
    assert_eq!(kernels.bytes(), first);
    swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    kernels.next_point();
    kernels.begin(|_| u64::MAX);
    // 残るのは 2 つ目の燃焼点で使った格子点（動かした格子と、元の格子と重なる最後の点）だけ
    assert_eq!(kernels.bytes(), table_bytes(&other));
    let (_, _, again) = swept_on(&other, 1.3, usize::MAX, Some(&kernels));
    assert_eq!(again, 0);
}

#[test]
fn a_full_kernel_table_gives_way_only_to_the_points_of_a_later_burnup_point() {
    // 上限に届いた表は、今の燃焼点で使っていない格子点（前の燃焼点のもの）を捨てて今の燃焼点の格子点を入れ、今の燃焼点で使った
    // 格子点は、同じ燃焼点のほかの格子が上限に届いても捨てない
    let grid = test_grid();
    let other = moved(&grid);
    let cap = table_bytes(&grid);
    let kernels = KernelCache::default();
    kernels.begin(|_| cap);
    swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    assert_eq!(kernels.bytes(), cap);
    kernels.next_point();
    kernels.begin(|_| cap);
    swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    let (_, _, again) = swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    // 動かした格子の 2 本目は、ほとんどの格子点を表から引く（上限は元の格子の分で、動かした格子の格子点の数はほぼ同じ）
    let alone = swept_on(&other, 1.0, usize::MAX, None).2;
    assert!(4 * again < alone, "{again} of {alone}");
    // 同じ燃焼点の次の自己遮蔽（元の格子）は表に入れられず、動かした格子の格子点を締め出さない
    kernels.begin(|_| cap);
    swept_on(&grid, 1.3, usize::MAX, Some(&kernels));
    let (_, _, after) = swept_on(&other, 1.3, usize::MAX, Some(&kernels));
    assert_eq!(after, again);
}

#[test]
fn a_lower_cap_keeps_the_points_of_the_current_burnup_point() {
    // 上限が今の大きさより小さくなっても（使用可能メモリが表の分だけ減る）、今の燃焼点で使った格子点は捨てず、増やさない。前の
    // 燃焼点の格子点は、今の燃焼点の格子点を入れようとして上限に届いたときに捨てる
    let grid = test_grid();
    let other = moved(&grid);
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    let held = kernels.bytes();
    assert!(held > 0);
    assert_eq!(kernels.begin(|held| held - 1), 0);
    swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    assert_eq!(kernels.bytes(), held);
    let (_, _, again) = swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
    assert_eq!(again, 0);
    assert_eq!(kernels.begin(|held| held + 10), 10);
    kernels.next_point();
    assert_eq!(kernels.begin(|held| held - 1), 0);
    assert_eq!(kernels.bytes(), held);
    swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    assert!(kernels.bytes() < held);
    let (_, _, again) = swept_on(&other, 1.0, usize::MAX, Some(&kernels));
    let alone = swept_on(&other, 1.0, usize::MAX, None).2;
    assert!(again < alone, "{again} of {alone}");
}

#[test]
fn inserts_beside_the_drop_of_the_earlier_points_find_its_room() {
    // 上限に届いた表で、前の燃焼点の格子点を捨てている間に余地を確かめたほかの減速計算も、捨てた後の余地を使える
    let kernels = KernelCache::default();
    kernels.begin(|_| u64::MAX);
    let table = kernels.table(vec![1]);
    let point = || KernelPoint {
        d: vec![[0.0; 3]],
        evaluations: 0,
    };
    for k in 0..100_000 {
        kernels.insert(&table, [k, 0], point());
    }
    let full = kernels.bytes();
    kernels.next_point();
    kernels.begin(|_| full);
    let size = point().bytes();
    let threads = 32;
    let barrier = std::sync::Barrier::new(threads);
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                barrier.wait();
                assert!(kernels.has_room(size, 1));
            });
        }
    });
    assert_eq!(kernels.counts()[2], 0);
}

#[test]
fn the_kernel_table_grows_only_into_its_share_of_the_given_hold_memory() {
    // `ufg-hold-memory` の量は、前の自己遮蔽から残った表、表が増えうる量、ピンと塊に渡す量を合わせた全体（ADR-0248）
    let grid = test_grid();
    let other = moved(&grid);
    // 測ったピンの量（`demand`）で表に残りを渡すときも同じ
    for (whole, at_once, demand) in [
        (3, 1, None),
        (2, 4, None),
        (3, 2, Some(1)),
        (4, 3, Some(u64::MAX)),
    ] {
        let kernels = KernelCache::default();
        kernels.begin(|_| u64::MAX);
        swept_on(&grid, 1.0, usize::MAX, Some(&kernels));
        let whole = whole * kernels.bytes();
        let held = kernels.start();
        let total = super::super::hold_total(Some(whole), held).unwrap();
        let growth = kernels.limit(super::super::kernel_cap(total, held, at_once, demand));
        swept_on(&other, 1.0, usize::MAX, Some(&kernels));
        assert!(kernels.bytes() <= held + growth, "{at_once}");
        assert!(kernels.bytes() + (total - growth) <= whole, "{at_once}");
    }
}

#[test]
fn an_infinite_multiple_of_kt_leaves_the_free_gas_uncapped_at_any_temperature() {
    // 上限なし（`upper = inf`、`heavy = inf`）は、kT が 0 でもどの入射も上限より下にする（無限大 × 0 の NaN にしない）
    for temperature in [0.0, 600.0] {
        let gas = Gas::new(15.86, temperature, f64::INFINITY);
        assert_eq!(gas.cap(), f64::INFINITY, "{temperature} K");
    }
}
