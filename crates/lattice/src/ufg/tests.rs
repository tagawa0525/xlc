use super::*;
use xlc_format::group::GroupStructure;
use xlc_format::library::{NuclideData, Pointwise, Resolved, Shielded, TemperatureData};

const TEMPERATURE: f64 = 600.0;
const BOUNDS: [f64; 7] = [2.0e7, 2.0e4, 1.0e4, 1.0e3, 1.0e2, 4.0, 1.0e-5];

/// ポイントワイズの断面積 (energy, total, elastic, fission) を持つ核種。群定数は 1/E の重みの平均
struct Spec {
    name: &'static str,
    awr: f64,
    top: Option<f64>,
    points: Vec<[f64; 4]>,
}

fn constant(name: &'static str, awr: f64, total: f64, elastic: f64, fission: f64) -> Spec {
    let points = [4.0, 2.0e4].map(|e| [e, total, elastic, fission]).to_vec();
    Spec {
        name,
        awr,
        top: None,
        points,
    }
}

fn library(specs: Vec<Spec>) -> Library {
    let groups = GroupStructure::new(BOUNDS.to_vec()).unwrap();
    let ng = BOUNDS.len() - 1;
    let nuclides = specs
        .into_iter()
        .map(|spec| {
            let column = |k: usize| spec.points.iter().map(|p| p[k]).collect::<Vec<_>>();
            // 群定数はポイントワイズの断面積の両端の値を使う（範囲の外の減速源に使う）
            let first = spec.points[0];
            let group = |k: usize| vec![first[k]; ng];
            NuclideData {
                name: spec.name.into(),
                za: 0.0,
                awr: spec.awr,
                resolved: spec.top.map(|top| Resolved {
                    top,
                    potential: 10.0,
                }),
                thermal_scattering: false,
                temperatures: vec![TemperatureData {
                    reactions: Vec::new(),
                    temperature: TEMPERATURE,
                    chi: vec![0.0; ng],
                    scatter: vec![vec![vec![0.0; ng]; ng]],
                    shielded: vec![Shielded {
                        total: group(1),
                        elastic: group(2),
                        fission: group(3),
                        nu_fission: group(3).iter().map(|f| 2.4 * f).collect(),
                        capture: group(3),
                    }],
                    pointwise: Pointwise {
                        energy: column(0),
                        total: column(1),
                        elastic: column(2),
                        fission: column(3),
                        nu_fission: column(3).iter().map(|f| 2.4 * f).collect(),
                        capture: column(3),
                    },
                }],
            }
        })
        .collect();
    Library {
        groups,
        sigma0: vec![1.0e10],
        pointwise_range: Some([4.0, 2.0e4]),
        nuclides,
    }
}

fn material(composition: &[(&str, f64)]) -> MaterialInput {
    MaterialInput {
        name: "m".into(),
        temperature: TEMPERATURE,
        composition: composition
            .iter()
            .map(|&(n, d)| (n.to_string(), d))
            .collect(),
    }
}

/// 燃料・ギャップ・被覆管・水
const RADII: [f64; 4] = [0.4096, 0.418, 0.475, 0.711];

#[test]
fn range_follows_the_resolved_top_of_the_fuel() {
    let fuel = material(&[("U-235", 0.01)]);
    let clad = material(&[("Zr-90", 0.04)]);
    let water = material(&[("H-1", 0.05)]);
    let range = |top: f64| {
        let mut fissile = constant("U-235", 233.0, 10.0, 8.0, 1.0);
        fissile.top = Some(top);
        // 燃料の外の核種の分離共鳴は範囲を決めない
        let mut zirconium = constant("Zr-90", 89.0, 6.0, 6.0, 0.0);
        zirconium.top = Some(1.0e6);
        let lib = library(vec![
            fissile,
            zirconium,
            constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        ]);
        cell_range(&lib, &[Some(&fuel), None, Some(&clad), Some(&water)])
            .unwrap()
            .unwrap()
    };
    // 上端以下で最も高い群の境界
    assert_eq!(range(3.0e3), [4.0, 1.0e3]);
    // 20 keV で頭打ちにする
    assert_eq!(range(5.0e4), [4.0, 2.0e4]);
}

#[test]
fn no_range_without_resolved_fuel() {
    // 燃料がない、または燃料に分離共鳴を持つ核種がなければ、減速計算する範囲はない
    let lib = library(vec![
        constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        constant("U-235", 233.0, 10.0, 8.0, 1.0),
    ]);
    let water = material(&[("H-1", 0.05)]);
    let fuel = material(&[("U-235", 0.01)]);
    assert_eq!(cell_range(&lib, &[Some(&water)]).unwrap(), None);
    assert_eq!(
        cell_range(&lib, &[Some(&fuel), Some(&water)]).unwrap(),
        None
    );
}

#[test]
fn resolved_fuel_needs_pointwise_data() {
    // 燃料に分離共鳴があるのにポイントワイズの断面積のないライブラリはエラー。燃料のないセルは範囲がないだけでよい
    let mut fissile = constant("U-235", 233.0, 10.0, 8.0, 1.0);
    fissile.top = Some(3.0e3);
    let mut lib = library(vec![fissile, constant("H-1", 0.99917, 20.0, 20.0, 0.0)]);
    lib.pointwise_range = None;
    let water = material(&[("H-1", 0.05)]);
    let fuel = material(&[("U-235", 0.01)]);
    assert_eq!(cell_range(&lib, &[Some(&water)]).unwrap(), None);
    assert!(matches!(
        cell_range(&lib, &[Some(&fuel), Some(&water)]),
        Err(UfgError::NoPointwise)
    ));
}

#[test]
fn scatterers_keep_the_asymptotic_spectrum() {
    // 吸収がなく断面積が一定なら、どのリングでも中性子束はレサジーあたり一定（1/E）のまま。減速材は、窓が
    // 群の上端（中性子のない範囲）まで届かない H-2 にする（H-1 では範囲の上からの線源が e^{−7.6} だけ足りない）
    let lib = library(vec![
        constant("U-238", 236.0, 11.0, 11.0, 0.0),
        constant("Zr-90", 89.0, 6.5, 6.5, 0.0),
        constant("H-2", 1.9968, 3.4, 3.4, 0.0),
        constant("O-16", 15.86, 3.8, 3.8, 0.0),
    ]);
    let fuel = material(&[("U-238", 0.022), ("O-16", 0.045)]);
    let clad = material(&[("Zr-90", 0.043)]);
    let water = material(&[("H-2", 0.05), ("O-16", 0.025)]);
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let spectrum = slowing_down(
        &lib,
        &cell,
        &[Some(&fuel), None, Some(&clad), Some(&water)],
        [4.0, 1.0e4],
    )
    .unwrap();
    assert_eq!(spectrum.groups, 2..5);
    for (ring, flux) in spectrum.flux.iter().enumerate() {
        for (k, f) in flux.iter().enumerate() {
            let g = spectrum.groups.start + k;
            let per_lethargy = f / (BOUNDS[g] / BOUNDS[g + 1]).ln();
            assert!(
                (per_lethargy - 1.0).abs() < 1e-9,
                "{ring} {g}: {per_lethargy}"
            );
        }
    }
    // 群定数はポイントワイズの断面積の値のまま
    assert!((spectrum.nuclides[0][0].total[1] - 11.0).abs() < 1e-12);
    assert!(spectrum.nuclides[1].is_empty());
}

#[test]
fn discontinuities_keep_both_limits() {
    // 500 eV で 10 b から 20 b に跳ぶ微量の核種: 中性子束は 1/E のままなので、群平均は両側の値の
    // レサジーの重みの平均になる
    let step = Spec {
        name: "Pu-239",
        awr: 236.999,
        top: None,
        points: vec![
            [4.0, 10.0, 0.0, 0.0],
            [500.0, 10.0, 0.0, 0.0],
            [500.0, 20.0, 0.0, 0.0],
            [2.0e4, 20.0, 0.0, 0.0],
        ],
    };
    let lib = library(vec![
        constant("U-238", 236.0, 11.0, 11.0, 0.0),
        constant("H-2", 1.9968, 3.4, 3.4, 0.0),
        step,
    ]);
    let fuel = material(&[("U-238", 0.022), ("Pu-239", 1.0e-12)]);
    let water = material(&[("H-2", 0.05)]);
    let cell = Cylinder::new(vec![0.4, 0.7]).unwrap();
    let spectrum = slowing_down(&lib, &cell, &[Some(&fuel), Some(&water)], [4.0, 1.0e4]).unwrap();
    // 1 keV〜100 eV の群
    let k = 3 - spectrum.groups.start;
    let expected = (10.0 * 5.0f64.ln() + 20.0 * 2.0f64.ln()) / 10.0f64.ln();
    let total = spectrum.nuclides[0][1].total[k];
    assert!(
        (total - expected).abs() < 1e-9 * expected,
        "{total} vs {expected}"
    );
}

#[test]
fn hydrogen_slowing_down_matches_the_analytic_solution() {
    // H と吸収だけの核種の均質な媒質: 衝突密度 q は dq/du = −r q で減る。
    // r = 1 − Σs/((1 − α)Σt)、q(0) は範囲の上（中性子束 1/E、群の上端より上は 0）からの散乱
    let awr: f64 = 0.99917;
    let lib = library(vec![
        constant("H-1", awr, 20.0, 20.0, 0.0),
        constant("Pu-240", 238.0, 2.0, 0.0, 0.0),
    ]);
    let mix = material(&[("H-1", 0.05), ("Pu-240", 0.01)]);
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let spectrum = slowing_down(&lib, &cell, &[Some(&mix); 4], [4.0, 1.0e4]).unwrap();
    let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
    let (sigma_s, sigma_t) = (0.05 * 20.0, 0.05 * 20.0 + 0.01 * 2.0);
    let r = 1.0 - sigma_s / ((1.0 - alpha) * sigma_t);
    let q0 = sigma_s * (1.0 - 1.0e4 / 2.0e7) / (1.0 - alpha);
    for flux in &spectrum.flux {
        for (k, f) in flux.iter().enumerate() {
            let g = spectrum.groups.start + k;
            let (ua, ub) = ((1.0e4 / BOUNDS[g]).ln(), (1.0e4 / BOUNDS[g + 1]).ln());
            let expected = q0 / (sigma_t * r) * ((-r * ua).exp() - (-r * ub).exp());
            assert!(
                (f - expected).abs() < 1e-7 * expected,
                "{g}: {f} vs {expected}"
            );
        }
    }
}

#[test]
fn reports_the_lethargy_points_of_the_grid() {
    // 核種の点は範囲の端（4 eV、2e4 eV）にしかないので、格子は群の境界（1e4、1e3、1e2、4 eV）と、その間を最も重い核種
    // （Pu-240）の窓の幅の 1/WINDOW_POINTS 以下の間隔で埋めた点
    let lib = library(vec![
        constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        constant("Pu-240", 238.0, 2.0, 0.0, 0.0),
    ]);
    let mix = material(&[("H-1", 0.05), ("Pu-240", 0.01)]);
    let cell = Cylinder::new(RADII.to_vec()).unwrap();
    let spectrum = slowing_down(&lib, &cell, &[Some(&mix); 4], [4.0, 1.0e4]).unwrap();
    let step = -(237.0f64 / 239.0).powi(2).ln() / WINDOW_POINTS as f64;
    let between = |ratio: f64| (ratio.ln() / step).ceil() as usize - 1;
    let expected = 4 + 2 * between(10.0) + between(25.0);
    assert_eq!(spectrum.points, expected);
}

#[test]
fn resonance_is_shielded_more_inside_the_fuel() {
    // 200 eV の強い吸収の共鳴: 燃料の内側ほど中性子束がへこみ、実効断面積が小さい
    let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
    let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
    energies.splice(0..0, [4.0, 100.0]);
    energies.extend([1.0e3, 2.0e4]);
    let points = energies
        .iter()
        .map(|&e| [e, 11.0 + resonance(e), 11.0, 0.0])
        .collect();
    let lib = library(vec![
        Spec {
            name: "U-238",
            awr: 236.0,
            top: Some(2.0e4),
            points,
        },
        constant("U-235", 233.0, 10.0, 8.0, 1.0),
        constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        constant("O-16", 15.86, 3.8, 3.8, 0.0),
    ]);
    let fuel = material(&[("U-238", 0.022), ("U-235", 0.0005), ("O-16", 0.045)]);
    let water = material(&[("H-1", 0.05), ("O-16", 0.025)]);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let spectrum = slowing_down(&lib, &cell, &rings, [4.0, 1.0e4]).unwrap();
    // 100 eV〜1 keV の群
    let k = 3 - spectrum.groups.start;
    let total = |ring: usize| spectrum.nuclides[ring][0].total[k];
    // 無限希釈（1/E の重み）の平均
    let dilute = {
        let (mut rate, mut flux) = (0.0, 0.0);
        for w in energies
            .windows(2)
            .filter(|w| w[0] >= 100.0 && w[1] <= 1.0e3)
        {
            let du = (w[1] / w[0]).ln();
            rate += 0.5 * du * (resonance(w[0]) + resonance(w[1]));
            flux += du;
        }
        11.0 + rate / flux
    };
    assert!(
        total(0) < total(1) && total(1) < total(2),
        "{:?}",
        [total(0), total(1), total(2)]
    );
    assert!(total(2) < dilute, "{} vs {dilute}", total(2));
}

#[test]
fn rejects_a_range_off_the_group_bounds() {
    let lib = library(vec![constant("H-1", 0.99917, 20.0, 20.0, 0.0)]);
    let water = material(&[("H-1", 0.05)]);
    let cell = Cylinder::new(vec![0.5]).unwrap();
    for range in [[4.0, 5.0e3], [3.0, 1.0e4], [1.0e2, 1.0e2]] {
        assert!(
            slowing_down(&lib, &cell, &[Some(&water)], range).is_err(),
            "{range:?}"
        );
    }
    // リングの数が合わない
    assert!(slowing_down(&lib, &cell, &[Some(&water); 2], [4.0, 1.0e4]).is_err());
}

/// 200 eV の強い吸収の共鳴を持つ燃料、真空のギャップ、水のセル（[`resonance_is_shielded_more_inside_the_fuel`] と同じ）
fn resonant_cell(total_at_peak: f64) -> (Library, MaterialInput, MaterialInput) {
    let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
    let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
    energies.splice(0..0, [4.0, 100.0]);
    energies.extend([1.0e3, 2.0e4]);
    // 共鳴の頂点（180 eV から 200 番目の点）の全断面積だけを `total_at_peak` にする
    let peak = 2 + 200;
    let points = energies
        .iter()
        .enumerate()
        .map(|(i, &e)| {
            let total = if i == peak {
                total_at_peak
            } else {
                11.0 + resonance(e)
            };
            [e, total, 11.0, 0.0]
        })
        .collect();
    let lib = library(vec![
        Spec {
            name: "U-238",
            awr: 236.0,
            top: Some(2.0e4),
            points,
        },
        constant("U-235", 233.0, 10.0, 8.0, 1.0),
        constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        constant("O-16", 15.86, 3.8, 3.8, 0.0),
    ]);
    let fuel = material(&[("U-238", 0.022), ("U-235", 0.0005), ("O-16", 0.045)]);
    let water = material(&[("H-1", 0.05), ("O-16", 0.025)]);
    (lib, fuel, water)
}

#[test]
fn table_tracks_give_the_group_constants_of_the_exact_tracks() {
    // ADR-0156: 既定の表（10 倍あたり 16 点）で内挿した飛程の群定数は、格子点ごとに求めた飛程の群定数と近い
    let (lib, fuel, water) = resonant_cell(11.0 + 1.0e4);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let solve = |collision| {
        slowing_down_with(
            &lib,
            &cell,
            &rings,
            &no_background(&rings),
            [4.0, 1.0e4],
            collision,
        )
        .unwrap()
    };
    let (table, exact) = (solve(Collision::default()), solve(Collision::Exact));
    for (ring, (t, e)) in table.nuclides.iter().zip(&exact.nuclides).enumerate() {
        for (k, (t, e)) in t.iter().zip(e).enumerate() {
            let fields = [
                (&t.total, &e.total),
                (&t.elastic, &e.elastic),
                (&t.fission, &e.fission),
                (&t.nu_fission, &e.nu_fission),
            ];
            for (f, (a, b)) in fields.into_iter().enumerate() {
                for (g, (a, b)) in a.iter().zip(b).enumerate() {
                    let d = if *b == 0.0 {
                        a.abs()
                    } else {
                        (a / b - 1.0).abs()
                    };
                    assert!(
                        d < 1e-6,
                        "ring {ring} nuclide {k} field {f} group {g}: {a} vs {b}"
                    );
                }
            }
        }
    }
}

#[test]
fn non_finite_totals_stop_with_the_table_as_without() {
    // 格子点の全断面積が NaN なら、表を作るときも、格子点ごとに求めるときと同じく止まる
    let (lib, fuel, water) = resonant_cell(f64::NAN);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    for collision in [Collision::Exact, Collision::default()] {
        let result = slowing_down_with(
            &lib,
            &cell,
            &rings,
            &no_background(&rings),
            [4.0, 1.0e4],
            collision,
        );
        assert!(
            matches!(result, Err(UfgError::Cp(CpError::Total(0)))),
            "{collision:?}: {:?}",
            result.err()
        );
    }
}

#[test]
fn zero_totals_stop_the_table_with_their_ring() {
    // ln Σ の格子は Σ = 0 を覆えないので、表にできないリングを番号で示して止める
    // 燃料は U-238 だけにして、共鳴の頂点の 1 点で全断面積を 0 にする
    let (lib, _, water) = resonant_cell(0.0);
    let fuel = material(&[("U-238", 0.022)]);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let result = slowing_down_with(
        &lib,
        &cell,
        &rings,
        &no_background(&rings),
        [4.0, 1.0e4],
        Collision::default(),
    );
    assert!(
        matches!(result, Err(UfgError::ZeroTotal(0))),
        "{:?}",
        result.err()
    );
}

#[test]
fn tables_with_as_many_nodes_as_points_give_the_exact_tracks() {
    // ADR-0182: 格子点ごとに求めるより点の多い表は作らず、そのセルは格子点ごとに求める（10 倍あたり 10 万点で、燃料の全断面積の
    // 3 桁ほどの範囲に約 27 万点。格子点は U-238 の散乱の窓で決まり、その数より多い）
    let (lib, fuel, water) = resonant_cell(11.0 + 1.0e4);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let solve = |collision| {
        slowing_down_with(
            &lib,
            &cell,
            &rings,
            &no_background(&rings),
            [4.0, 1.0e4],
            collision,
        )
        .unwrap()
    };
    assert_eq!(
        solve(Collision::Table { density: 100_000 }),
        solve(Collision::Exact)
    );
}

/// [`resonant_cell`] の燃料の U-238 の数密度を、内側から `u238` の 3 つのリングの物質にしたセルを、既定の表と格子点ごとに解く
fn ring_cells(u238: [f64; 3]) -> (CellSpectrum, CellSpectrum) {
    let (lib, _, water) = resonant_cell(11.0 + 1.0e4);
    let fuel = u238.map(|n| material(&[("U-238", n), ("U-235", 0.0005), ("O-16", 0.045)]));
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [
        Some(&fuel[0]),
        Some(&fuel[1]),
        Some(&fuel[2]),
        None,
        Some(&water),
    ];
    let solve = |collision| {
        slowing_down_with(
            &lib,
            &cell,
            &rings,
            &no_background(&rings),
            [4.0, 1.0e4],
            collision,
        )
        .unwrap()
    };
    (solve(Collision::default()), solve(Collision::Exact))
}

#[test]
fn close_ring_materials_give_the_group_constants_of_the_exact_tracks() {
    // ADR-0182: 核種の並びが同じで数密度の少し違うリングの物質は、燃料の 1 つの軸と 1 次の展開で解け、群定数は格子点ごとに
    // 求めたものと近い（燃焼で外側の U-238 が減ったリング）。このセルで表だけの差（3 リングを同じ組成にしたとき）は 8e-7 で、
    // 展開の差（上限で衝突確率に 1e-5 の水準）がそれと同じ程度に加わる
    let (table, exact) = ring_cells([0.0221, 0.0220, 0.0218]);
    // 表と展開を通っている（格子点ごとに求めたものと丸めで一致しない）
    assert_ne!(table, exact);
    for (ring, (t, e)) in table.nuclides.iter().zip(&exact.nuclides).enumerate() {
        for (k, (t, e)) in t.iter().zip(e).enumerate() {
            let fields = [
                (&t.total, &e.total),
                (&t.elastic, &e.elastic),
                (&t.fission, &e.fission),
                (&t.nu_fission, &e.nu_fission),
            ];
            for (f, (a, b)) in fields.into_iter().enumerate() {
                for (g, (a, b)) in a.iter().zip(b).enumerate() {
                    let d = if *b == 0.0 {
                        a.abs()
                    } else {
                        (a / b - 1.0).abs()
                    };
                    assert!(
                        d < 3e-6,
                        "ring {ring} nuclide {k} field {f} group {g}: {a} vs {b}"
                    );
                }
            }
        }
    }
}

#[test]
fn distant_ring_materials_are_solved_point_by_point() {
    // ADR-0182: どの格子点でもリングの全断面積の差が展開の上限を超えるセルは、表を作っても衝突確率の数が減らないので、
    // 格子点ごとに求める
    let (table, exact) = ring_cells([0.044, 0.022, 0.011]);
    assert_eq!(table, exact);
}

#[test]
fn capture_is_shielded_with_the_same_flux_as_fission() {
    // テストのライブラリは捕獲に核分裂と同じポイントワイズの列を持たせる。共鳴を核分裂にも置けば、遮蔽した捕獲は核分裂に一致する
    let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
    let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
    energies.splice(0..0, [4.0, 100.0]);
    energies.extend([1.0e3, 2.0e4]);
    let points = energies
        .iter()
        .map(|&e| [e, 11.0 + resonance(e), 11.0, 0.1 * resonance(e)])
        .collect();
    let lib = library(vec![
        Spec {
            name: "U-238",
            awr: 236.0,
            top: Some(2.0e4),
            points,
        },
        constant("H-1", 0.99917, 20.0, 20.0, 0.0),
        constant("O-16", 15.86, 3.8, 3.8, 0.0),
    ]);
    let fuel = material(&[("U-238", 0.022), ("O-16", 0.045)]);
    let water = material(&[("H-1", 0.05), ("O-16", 0.025)]);
    let cell = Cylinder::new(vec![0.2, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), None, Some(&water)];
    let spectrum = slowing_down(&lib, &cell, &rings, [4.0, 1.0e4]).unwrap();
    for ring in [0, 1] {
        let u238 = &spectrum.nuclides[ring][0];
        assert_eq!(u238.capture, u238.fission);
        assert!(u238.capture.iter().any(|&c| c > 0.0));
    }
}

#[test]
fn points_beyond_the_expansion_do_not_widen_the_axes() {
    // 展開の上限を超える格子点は表を読まないので、軸の範囲を広げない（Gd の共鳴の頂点で燃料の軸が伸びて表が大きくならない）
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let layout = Layout::new(
        &cell,
        vec![Some(0), Some(1), Some(2), None, Some(3)],
        vec![0, 0, 0, 1],
    )
    .unwrap();
    let mut ranges = vec![[f64::INFINITY, f64::NEG_INFINITY]; 2];
    assert!(widen_ranges(&layout, &[0.5, 0.5, 0.5, 1.4], &mut ranges));
    assert!(!widen_ranges(
        &layout,
        &[10.0, 100.0, 1000.0, 1.5],
        &mut ranges
    ));
    // 燃料の軸は ln の重みの平均を exp で戻すので丸めの範囲で 0.5、水の軸はその Σ そのもの
    for x in ranges[0] {
        assert!((x - 0.5).abs() <= 1e-15, "{ranges:?}");
    }
    assert_eq!(ranges[1], [1.4, 1.4]);
}

/// [`resonant_cell`] の燃料のリングで、U-235（組成の 2 番目）だけを背景にして群定数 `u235` を与えた `background`
fn u235_background(
    rings: &[Option<&MaterialInput>],
    u235: &Shielded,
) -> Vec<Vec<Option<Shielded>>> {
    let mut background = no_background(rings);
    for (ring, m) in background.iter_mut().zip(rings) {
        if m.is_some_and(|m| m.composition[1].0 == "U-235") {
            ring[1] = Some(u235.clone());
        }
    }
    background
}

#[test]
fn a_background_of_the_constant_cross_sections_gives_the_pointwise_solution() {
    // ADR-0185: 断面積が一定の U-235 を、同じ値の群定数の背景にしても、中性子束とほかの核種の群定数はポイントワイズで解いたものと
    // 同じ（背景の項は U-235 だけなので、まとめた α も U-235 のもの）
    let (lib, fuel, water) = resonant_cell(11.0 + 1.0e4);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let dilute = lib.lookup("U-235", TEMPERATURE).unwrap().1.shielded[0].clone();
    let range = [4.0, 1.0e4];
    let pointwise = slowing_down_with(
        &lib,
        &cell,
        &rings,
        &no_background(&rings),
        range,
        Collision::Exact,
    )
    .unwrap();
    let background = u235_background(&rings, &dilute);
    let lumped =
        slowing_down_with(&lib, &cell, &rings, &background, range, Collision::Exact).unwrap();
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs();
    for (ring, (a, b)) in lumped.flux.iter().zip(&pointwise.flux).enumerate() {
        for (g, (a, b)) in a.iter().zip(b).enumerate() {
            assert!(close(*a, *b), "ring {ring} group {g}: {a} vs {b}");
        }
    }
    for ring in 0..3 {
        // U-238 は共鳴で遮蔽した値が、U-235 は与えた値がそのまま一致する
        for k in 0..3 {
            let (a, b) = (&lumped.nuclides[ring][k], &pointwise.nuclides[ring][k]);
            for (x, y) in a
                .total
                .iter()
                .zip(&b.total)
                .chain(a.fission.iter().zip(&b.fission))
            {
                assert!(close(*x, *y), "ring {ring} nuclide {k}: {x} vs {y}");
            }
        }
    }
}

#[test]
fn background_nuclides_keep_the_given_group_constants() {
    // ADR-0185: 背景の核種の群定数は、群の中で一定の断面積を中性子束で縮約したものなので、与えた値のまま
    let (lib, fuel, water) = resonant_cell(11.0 + 1.0e4);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let ng = BOUNDS.len() - 1;
    let by_group = |base: f64| (0..ng).map(|g| base + g as f64).collect::<Vec<_>>();
    let given = Shielded {
        total: by_group(10.0),
        elastic: by_group(8.0),
        fission: by_group(1.0),
        nu_fission: by_group(2.4),
        capture: by_group(1.0),
    };
    let background = u235_background(&rings, &given);
    let spectrum = slowing_down_with(
        &lib,
        &cell,
        &rings,
        &background,
        [4.0, 1.0e4],
        Collision::default(),
    )
    .unwrap();
    let groups = spectrum.groups.clone();
    for ring in 0..3 {
        let u235 = &spectrum.nuclides[ring][1];
        assert_eq!(u235.total, given.total[groups.clone()]);
        assert_eq!(u235.capture, given.capture[groups.clone()]);
        assert_eq!(u235.nu_fission, given.nu_fission[groups.clone()]);
    }
}

#[test]
fn the_background_follows_the_rings_and_their_nuclides() {
    // 背景は、リングごとに物質の核種の数だけ、ライブラリのすべての群の値を持つ
    let (lib, fuel, water) = resonant_cell(11.0 + 1.0e4);
    let cell = Cylinder::new(vec![0.2, 0.35, 0.4096, 0.418, 0.711]).unwrap();
    let rings = [Some(&fuel), Some(&fuel), Some(&fuel), None, Some(&water)];
    let solve = |background: &[Vec<Option<Shielded>>]| {
        slowing_down_with(
            &lib,
            &cell,
            &rings,
            background,
            [4.0, 1.0e4],
            Collision::default(),
        )
    };
    let mut short = no_background(&rings);
    short.pop();
    assert!(matches!(solve(&short), Err(UfgError::Background(4))));
    let mut extra = no_background(&rings);
    extra[1].push(None);
    assert!(matches!(solve(&extra), Err(UfgError::Background(1))));
    let mut coarse = no_background(&rings);
    let mut one = lib.lookup("U-235", TEMPERATURE).unwrap().1.shielded[0].clone();
    one.capture.pop();
    coarse[2][1] = Some(one);
    assert!(matches!(solve(&coarse), Err(UfgError::Background(2))));
}

#[test]
fn lumped_background_keeps_the_slowing_down_power() {
    // ADR-0185: まとめた項は、群ごとの巨視的な断面積の和と、核種の減速能の和 Σ ξ N σs を保つ α を持つ
    let lib = library(vec![
        constant("O-16", 15.86, 3.8, 3.8, 0.0),
        constant("U-238", 236.0, 11.0, 9.0, 0.0),
    ]);
    let ng = BOUNDS.len() - 1;
    let flat = |v: f64| vec![v; ng];
    let shielded = |t: f64, s: f64| Shielded {
        total: flat(t),
        elastic: flat(s),
        fission: flat(0.0),
        nu_fission: flat(0.0),
        capture: flat(t - s),
    };
    let (oxygen, uranium) = (shielded(3.8, 3.8), shielded(11.0, 9.0));
    let members = [("O-16", 0.045, &oxygen), ("U-238", 0.022, &uranium)];
    let data = lump(&lib, &members, &BOUNDS, 1..5).unwrap();
    let xi = |awr: f64| mean_lethargy_gain(scattering_alpha("x", awr).unwrap());
    let (so, su) = (0.045 * 3.8, 0.022 * 9.0);
    let expected = (so * xi(15.86) + su * xi(236.0)) / (so + su);
    let alpha = 1.0 - 1.0 / data.gain;
    assert!((mean_lethargy_gain(alpha) - expected).abs() < 1e-12);
    assert!((data.width + alpha.ln()).abs() < 1e-12);
    let mut cursor = 0;
    let row = evaluate(&data, &mut cursor, 500.0, false).unwrap();
    assert!((row[0] - (0.045 * 3.8 + 0.022 * 11.0)).abs() < 1e-15);
    assert!((row[1] - (so + su)).abs() < 1e-15);
}

#[test]
fn group_backgrounds_jump_at_the_group_bounds() {
    // 背景の断面積は群の境界で跳び、境界の格子点の高い側は上の群、低い側は下の群の値になる
    let lib = library(vec![constant("U-238", 236.0, 11.0, 9.0, 0.0)]);
    let ng = BOUNDS.len() - 1;
    let by_group = (0..ng).map(|g| g as f64).collect::<Vec<_>>();
    let s = Shielded {
        total: by_group.clone(),
        elastic: by_group.clone(),
        fission: vec![0.0; ng],
        nu_fission: vec![0.0; ng],
        capture: vec![0.0; ng],
    };
    let data = lump(&lib, &[("U-238", 1.0, &s)], &BOUNDS, 1..5).unwrap();
    let mut cursor = data.first_cursor();
    let at =
        |cursor: &mut usize, e: f64, below: bool| evaluate(&data, cursor, e, below).unwrap()[0];
    // 群 1 は 2e4〜1e4 eV、群 2 は 1e4〜1e3 eV
    assert_eq!(at(&mut cursor, 2.0e4, true), 1.0);
    assert_eq!(at(&mut cursor, 1.5e4, false), 1.0);
    assert_eq!(at(&mut cursor, 1.0e4, false), 1.0);
    assert_eq!(at(&mut cursor, 1.0e4, true), 2.0);
    assert_eq!(at(&mut cursor, 5.0, false), 4.0);
}
