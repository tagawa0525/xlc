use std::path::Path;

use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, product_points};

use super::*;
use crate::input::{Boron, BoronSearch, Case, Cycle, Inventory, Plant};

const PLANT: &str = r#"
pitch = 20.0
height = 20.0
map = ["*  *", "*  *"]
"#;

const INVENTORY: &str = r#"
[design.A]
bounds = [0.0, 20.0]
layers = ["a.msgpack"]

[design.B]
bounds = [0.0, 20.0]
layers = ["b.msgpack"]

[assembly.A1]
design = "A"
[assembly.A2]
design = "A"
[assembly.B1]
design = "B"
[assembly.B2]
design = "B"
"#;

const CYCLE: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["A1  B1", "B2  A2"]
"#;

/// +x の面だけ中性子束 0 で、ほかは反射（漏れのある 2 × 2 の集合体）
const CASE: &str = r#"
cycle = "cycle.toml"
axial-mesh = [0.0, 10.0, 20.0]
radial-nodes = 1
boundary = { x = ["reflective", "zero-flux"], y = ["reflective", "reflective"], z = ["reflective", "reflective"] }
"#;

const PPM: [f64; 2] = [1000.0, 1500.0];
const TM: [f64; 2] = [560.0, 570.0];
const TF: [f64; 2] = [560.0, 570.0];

/// 2群の定数。熱群の吸収はホウ素濃度で、減速は減速材の温度で、速群の吸収は燃料の温度で変わる。`fissile` は
/// 核分裂の強さ
fn few(fissile: f64, ppm: f64, tm: f64, tf: f64) -> FewGroup {
    let absorption = vec![0.01 + 1.0e-6 * (tf - 565.0), 0.06 + 4.0e-5 * (ppm - 1000.0)];
    let removal = 0.02 - 1.0e-5 * (tm - 565.0);
    let scatter = vec![vec![0.5, removal], vec![0.0, 1.0]];
    let nu_fission = vec![0.005 * fissile, 0.1 * fissile];
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.0,
        total: vec![absorption[0] + 0.5 + removal, absorption[1] + 1.0],
        fission: nu_fission.iter().map(|x| x / 2.43).collect(),
        nu_fission,
        absorption,
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter,
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

fn table(fissile: f64) -> Table {
    let mut states = Vec::new();
    for ppm in PPM {
        for tm in TM {
            for tf in TF {
                states.push(few(fissile, ppm, tm, tf));
            }
        }
    }
    let mut table = Table {
        points: Vec::new(),
        split: 1,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        axes: vec![
            Axis::boron(PPM.to_vec(), 0.199),
            Axis::moderator(TM.to_vec(), vec![0.75, 0.74]),
            Axis::fuel(TF.to_vec()),
        ],
        states,
    };
    table.points = product_points(&table.axes);
    table
}

/// 炉心の並びと、`layout.files` の順のテーブル
fn core() -> (Layout, Vec<Table>) {
    let base = Path::new("dir");
    let case = Case::parse(CASE, base).unwrap();
    let cycle = Cycle::parse(CYCLE, base).unwrap();
    let plant = Plant::parse(PLANT, base).unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = layout
        .files
        .iter()
        .map(|f| table(if f.ends_with("a.msgpack") { 1.0 } else { 0.9 }))
        .collect();
    (layout, tables)
}

fn given(name: &str, boron: f64, temperature: f64) -> CaseState {
    CaseState {
        name: name.into(),
        boron: Boron::Given(boron),
        b10: 0.199,
        moderator: temperature,
        fuel: temperature,
        burnup: None,
        banks: Default::default(),
    }
}

fn searching(bracket: [f64; 2], tolerance: f64) -> CaseState {
    CaseState {
        name: "critical".into(),
        boron: Boron::Search(BoronSearch {
            bracket,
            k: 1.0,
            tolerance,
            max_solves: 20,
        }),
        b10: 0.199,
        moderator: 565.0,
        fuel: 565.0,
        burnup: None,
        banks: Default::default(),
    }
}

#[test]
fn boron_search_converges_to_the_target_k() {
    let (layout, tables) = core();
    let options = Options::default();
    let state = searching([1000.0, 1500.0], 1.0e-7);
    let solved = solve_state(&layout, &tables, &state, &options).unwrap();
    let searched = solved.searched.unwrap();
    assert!(
        (solved.solution.k - 1.0).abs() <= 1.0e-7,
        "{}",
        solved.solution.k
    );
    assert!(
        1000.0 < searched.ppm && searched.ppm < 1500.0,
        "{}",
        searched.ppm
    );
    // この定数の k は範囲の中で 1.13 から 0.88 まで曲がって下がる。傾き（約 −4.6e-4 /ppm）から、二分法では許容に
    // 入るまで範囲を 2e-4 ppm まで狭める（両端と約 21 回）。Illinois 法はその半分より少ない
    assert!(searched.solves <= 10, "{} solves", searched.solves);
    // 求めたホウ素濃度で解き直しても同じ k
    let direct = solve_at(&layout, &tables, Some(&state.at(searched.ppm)), &options).unwrap();
    assert_eq!(direct.solution.k, solved.solution.k);
}

#[test]
fn search_errors_are_reported() {
    let (layout, tables) = core();
    let options = Options::default();
    let error = |state: &CaseState| solve_state(&layout, &tables, state, &options).unwrap_err();
    // 範囲がテーブルの外（内挿のエラー）、範囲が解を挟まない、固有値の判定より細かい許容
    let e = error(&searching([900.0, 1500.0], 1.0e-7));
    assert!(
        matches!(e, CaseError::State(ref m) if m.contains("outside")),
        "{e}"
    );
    let e = error(&searching([1400.0, 1500.0], 1.0e-7));
    assert!(
        matches!(e, CaseError::Search(ref m) if m.contains("do not bracket")),
        "{e}"
    );
    let e = error(&searching([1000.0, 1500.0], options.k_tolerance));
    assert!(
        matches!(e, CaseError::Search(ref m) if m.contains("eigenvalue tolerance")),
        "{e}"
    );
    // 許容の検査は解く前にできる
    assert!(check_search(&searching([1000.0, 1500.0], options.k_tolerance), &options).is_err());
    assert!(check_search(&searching([1000.0, 1500.0], 1.0e-7), &options).is_ok());
}

#[test]
fn searched_states_solved_together_match_solved_alone() {
    let (layout, tables) = core();
    let options = Options::default();
    let states = [
        given("aro", 1100.0, 565.0),
        searching([1000.0, 1500.0], 1.0e-7),
    ];
    for (state, solved) in states
        .iter()
        .zip(solve_states(&layout, &tables, &states, &options))
    {
        let solved = solved.unwrap();
        let alone = solve_state(&layout, &tables, state, &options).unwrap();
        assert_eq!(solved.solution.k, alone.solution.k, "{}", state.name);
        assert_eq!(solved.searched, alone.searched, "{}", state.name);
    }
}

#[test]
fn states_solved_together_match_states_solved_alone() {
    let (layout, tables) = core();
    let options = Options::default();
    let states = [
        given("aro-1100", 1100.0, 565.0),
        given("aro-1300", 1300.0, 565.0),
        given("cold", 1100.0, 560.0),
        given("hot", 1100.0, 570.0),
    ];
    let together = solve_states(&layout, &tables, &states, &options);
    assert_eq!(together.len(), states.len());
    let mut k = Vec::new();
    for (state, solved) in states.iter().zip(together) {
        let solved = solved.unwrap();
        let alone = solve_state(&layout, &tables, state, &options).unwrap();
        assert_eq!(solved.solution.k, alone.solution.k, "{}", state.name);
        k.push(solved.solution.k);
    }
    // 価値と係数は2つの状態の k の差（ρ = 1 − 1/k）: ホウ素を増やすと k が下がる
    let rho = |k: f64| 1.0 - 1.0 / k;
    let boron_worth = (rho(k[1]) - rho(k[0])) / 200.0;
    assert!(boron_worth < 0.0, "{boron_worth}");
    // 温度を上げると減速が減って k が下がる
    assert!(rho(k[3]) < rho(k[2]), "{} {}", k[2], k[3]);
}

#[test]
fn errors_are_reported_per_state() {
    let (layout, tables) = core();
    let options = Options::default();
    // テーブルの範囲の外の状態だけがエラーになり、ほかの状態は解ける
    let states = [
        given("inside", 1100.0, 565.0),
        given("outside", 900.0, 565.0),
    ];
    let solved = solve_states(&layout, &tables, &states, &options);
    assert!(solved[0].is_ok());
    assert!(matches!(solved[1], Err(CaseError::State(_))));
}

#[test]
fn states_give_the_same_k_on_one_and_many_threads() {
    let (layout, tables) = core();
    let options = Options::default();
    let states = [
        given("aro-1100", 1100.0, 565.0),
        given("aro-1300", 1300.0, 565.0),
        given("cold", 1100.0, 560.0),
        given("hot", 1100.0, 570.0),
    ];
    let k = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve_states(&layout, &tables, &states, &options))
            .into_iter()
            .map(|s| s.unwrap().solution.k)
            .collect::<Vec<_>>()
    };
    assert_eq!(k(1), k(4));
}

#[test]
fn tables_must_match_the_layout_files() {
    let (layout, mut tables) = core();
    let options = Options::default();
    let state = given("aro", 1100.0, 565.0);
    tables.push(tables[0].clone());
    assert!(solve_at(&layout, &tables, Some(&state.at(1100.0)), &options).is_err());
    tables.truncate(layout.files.len() - 1);
    assert!(solve_at(&layout, &tables, Some(&state.at(1100.0)), &options).is_err());
}

#[test]
fn tables_with_history_residuals_need_micro_depleted_nodes() {
    // 残り R は、履歴の変数を進めるミクロ燃焼のノードにだけ足せる（ADR-0215）。燃焼させない状態の解では黙って無視せず止める
    let (layout, mut tables) = core();
    let options = Options::default();
    let state = given("aro", 1100.0, 565.0);
    tables[0]
        .histories
        .push(xlc_format::table::HistoryCorrection {
            variable: xlc_format::table::MODERATOR_DENSITY.into(),
            nominal: 0.75,
            branch: 0.74,
            rods: Vec::new(),
        });
    let err = solve_at(&layout, &tables, Some(&state.at(1100.0)), &options)
        .err()
        .unwrap()
        .to_string();
    assert!(err.contains("history"), "{err}");
}

// 履歴の補正（ADR-0215）

use xlc_format::table::{
    HistoryCorrection, HistoryRod, MODERATOR_DENSITY, Micro, MicroDerivative, MicroNuclide,
    MicroReaction, ROD_OUT, RateResidual, Residual, add_linear_correction,
};

use crate::state::State;

/// 2 群の少数群定数（`s` で値を変える。10 進で短く書けない値）
fn history_few(s: f64) -> FewGroup {
    let x = |v: f64| v * (1.0 + s / 7.0) / 3.0;
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: x(3.5),
        total: vec![x(0.55), x(1.3)],
        absorption: vec![x(0.01), x(0.1)],
        nu_fission: vec![x(0.007), x(0.14)],
        fission: vec![x(0.0028), x(0.057)],
        chi: vec![1.0 - s / 1000.0, s / 1000.0],
        diffusion: vec![x(1.4), x(0.4)],
        scatter: vec![vec![x(0.52), x(0.017)], vec![x(0.0015), x(1.2)]],
        adf: vec![vec![x(1.0), x(1.01)]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

/// 2 群の微視的断面積（`s` で値を変える）。`derivative` なら数密度の微分も持つ
fn sigma(s: f64, derivative: bool) -> MicroNuclide {
    let x = |v: f64| v * (1.0 + s / 5.0) / 3.0;
    MicroNuclide {
        total: vec![x(11.0), x(130.0)],
        nu_fission: vec![x(1.1), x(170.0)],
        fission: vec![x(0.4), x(70.0)],
        scatter: vec![vec![x(9.0), x(0.1)], vec![x(0.01), x(8.0)]],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![vec![x(0.3), x(20.0)]],
        }],
        thermal_fission: Some(0.9),
        derivative: derivative.then(|| MicroDerivative {
            total: vec![x(9.0), x(110.0)],
            nu_fission: vec![x(1.0), x(-150.0)],
            fission: vec![x(0.35), x(-60.0)],
            scatter: vec![vec![x(8.5), x(0.12)], vec![x(0.02), x(7.0)]],
        }),
    }
}

/// 公称の燃焼の均質化した数密度（燃焼度 0、10、20）
const NOMINAL_DENSITIES: [[f64; 2]; 3] = [[7.0e-4, 0.0], [6.0e-4, 1.0e-5], [5.0e-4, 2.0e-5]];
/// 履歴の分岐の均質化した数密度（燃焼度 0、20）
const BRANCH_DENSITIES: [[f64; 2]; 2] = [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]];

/// 履歴の分岐の自分の条件の定数（燃焼度 0、20。燃焼度 0 は公称と同じ組成なので、公称のテーブルの条件の状態と同じ定数）
fn branch_constants(table: &Table) -> [FewGroup; 2] {
    // 制御棒の軸を持てば、その値 out（最初の値）
    let mut at = vec![0, 0, 1, 0];
    at.resize(table.axes.len(), 0);
    let p = table.point_index(&at).unwrap();
    [table.states[p].clone(), history_few(9.5)]
}

/// 履歴の分岐の自分の条件の、反応率の微視的断面積（燃焼度 0、20 の点 `k`。微分は持たない）。収率の境界より下の割合も公称と違う
fn branch_sigma(k: usize) -> [MicroNuclide; 2] {
    let mut out = [sigma(12.0 + k as f64, false), sigma(24.0 + k as f64, false)];
    out[0].thermal_fission = Some(0.85);
    out
}

/// 反応率の微視的断面積の差 `a − b`（lattice histories と同じ式。収率の境界より下は割合 × 群の和の σf の差）
fn rate_difference(a: &MicroNuclide, b: &MicroNuclide) -> RateResidual {
    let minus = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p - q).collect::<Vec<_>>();
    let below = |n: &MicroNuclide| n.thermal_fission.unwrap() * n.fission.iter().sum::<f64>();
    RateResidual {
        fission: minus(&a.fission, &b.fission),
        reactions: a
            .reactions
            .iter()
            .zip(&b.reactions)
            .map(|(x, y)| MicroReaction {
                mt: x.mt,
                channels: x
                    .channels
                    .iter()
                    .zip(&y.channels)
                    .map(|(p, q)| minus(p, q))
                    .collect(),
            })
            .collect(),
        thermal_fission: Some(below(a) - below(b)),
    }
}

/// 少数群定数の差 `a − b` の残りの項目
fn difference(a: &FewGroup, b: &FewGroup) -> Residual {
    let minus = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p - q).collect::<Vec<_>>();
    let rows = |x: &[Vec<f64>], y: &[Vec<f64>]| {
        x.iter()
            .zip(y)
            .map(|(p, q)| minus(p, q))
            .collect::<Vec<_>>()
    };
    Residual {
        total: minus(&a.total, &b.total),
        absorption: minus(&a.absorption, &b.absorption),
        nu_fission: minus(&a.nu_fission, &b.nu_fission),
        fission: minus(&a.fission, &b.fission),
        chi: minus(&a.chi, &b.chi),
        diffusion: minus(&a.diffusion, &b.diffusion),
        scatter: rows(&a.scatter, &b.scatter),
        adf: rows(&a.adf, &b.adf),
    }
}

/// 公称のテーブル: 燃焼度 0、10、20、ホウ素 1300 ppm、減速材 565 K と 0.743 g/cm³（公称）と 600 K と 0.661 g/cm³（履歴の分岐の
/// 条件）、燃料 565 K。公称でない状態は 0 と 20 で解く。U-235 と Pu-240（数密度の微分を持つ）を追う。`rods` が空でなければ制御棒の軸を
/// 持つ。減速材の密度の履歴の残り R（制御棒の値 out、燃焼度 0 と 20）は、lattice と同じく分岐の定数 − (条件の状態 + 線形の補正)
fn history_table(rods: &[&str]) -> Table {
    let mut axes = vec![
        Axis::Burnup {
            points: vec![0.0, 10.0, 20.0],
        },
        Axis::boron(vec![1300.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.661]).with_nominal(565.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ];
    if !rods.is_empty() {
        axes.push(Axis::Rod {
            values: rods.iter().map(|r| r.to_string()).collect(),
        });
    }
    let points: Vec<Vec<usize>> = product_points(&axes)
        .into_iter()
        .filter(|p| p[0] != 1 || p[2] == 0)
        .collect();
    let n = points.len();
    let mut table = Table {
        axes,
        split: 1,
        states: (0..n).map(|p| history_few(p as f64)).collect(),
        heavy_metal: Some(400.0 / 3.0),
        micro: Some(Micro {
            nuclides: vec!["U-235".into(), "Pu-240".into()],
            densities: NOMINAL_DENSITIES.iter().map(|d| d.to_vec()).collect(),
            states: (0..n)
                .map(|p| vec![sigma(p as f64, false), sigma(2.0 * p as f64, true)])
                .collect(),
        }),
        histories: Vec::new(),
        points,
    };
    let branch = branch_constants(&table);
    let micro = table.micro.as_ref().unwrap();
    let residuals = [0usize, 2]
        .iter()
        .zip(&branch)
        .zip(&BRANCH_DENSITIES)
        .map(|((&u, fb), nj)| {
            let mut at = vec![u, 0, 1, 0];
            if !rods.is_empty() {
                at.push(0);
            }
            let p = table.point_index(&at).unwrap();
            let delta: Vec<f64> = nj
                .iter()
                .zip(&NOMINAL_DENSITIES[u])
                .map(|(a, b)| a - b)
                .collect();
            let mut corrected = table.states[p].clone();
            add_linear_correction(&mut corrected, &micro.states[p], &delta);
            difference(fb, &corrected)
        })
        .collect();
    // 反応率の微視的断面積の差: 履歴の分岐 − 公称の条件の状態（lattice と同じ）
    let rates = [0usize, 2]
        .iter()
        .enumerate()
        .map(|(k, &u)| {
            let mut at = vec![u, 0, 1, 0];
            if !rods.is_empty() {
                at.push(0);
            }
            let p = table.point_index(&at).unwrap();
            branch_sigma(k)
                .iter()
                .zip(&micro.states[p])
                .map(|(b, n)| rate_difference(b, n))
                .collect()
        })
        .collect();
    table.histories = vec![HistoryCorrection {
        variable: MODERATOR_DENSITY.into(),
        nominal: 0.743,
        branch: 0.661,
        rods: vec![HistoryRod {
            rod: ROD_OUT.into(),
            burnup: vec![0.0, 20.0],
            residuals,
            micro: rates,
        }],
    }];
    table.validate().unwrap();
    table
}

fn at(moderator: f64, burnup: f64) -> State {
    State {
        boron: 1300.0,
        b10: 0.199,
        moderator,
        fuel: 565.0,
        burnup: Some(burnup),
    }
}

/// 相対 1e-12 で同じ（丸めは項の大きさの範囲）
fn close(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| (x - y).abs() <= 1.0e-12 * x.abs().max(y.abs()).max(1.0e-6))
}

#[test]
fn history_branch_points_return_the_branch_constants() {
    // w = 1、条件 j、ノードの数密度が履歴の分岐の均質化した数密度なら、第 1 項の格子計算の定数をそのまま返す（ADR-0215 の確認）
    let table = history_table(&[]);
    let branch = branch_constants(&table);
    for (k, &b) in [0.0, 20.0].iter().enumerate() {
        let f =
            node_constants(&table, &at(600.0, b), None, &BRANCH_DENSITIES[k], &[0.661]).unwrap();
        let fb = &branch[k];
        assert!(
            close(&f.total, &fb.total),
            "{b}: {:?} {:?}",
            f.total,
            fb.total
        );
        assert!(close(&f.absorption, &fb.absorption), "{b}");
        assert!(close(&f.nu_fission, &fb.nu_fission), "{b}");
        assert!(close(&f.fission, &fb.fission), "{b}");
        assert!(close(&f.chi, &fb.chi), "{b}");
        assert!(close(&f.diffusion, &fb.diffusion), "{b}");
        assert!(close(&f.scatter.concat(), &fb.scatter.concat()), "{b}");
        assert!(close(&f.adf.concat(), &fb.adf.concat()), "{b}");
    }
}

#[test]
fn nominal_histories_leave_the_constants_bit_identical() {
    // w = 0（公称の履歴）では、R を持たないテーブルと同じ定数（ビットで同じ）
    let with = history_table(&[]);
    let mut without = with.clone();
    without.histories.clear();
    let n = [5.5e-4, 1.6e-5];
    for state in [
        at(565.0, 10.0),
        at(565.0, 20.0),
        at(600.0, 7.0),
        at(580.0, 15.0),
    ] {
        assert_eq!(
            node_constants(&with, &state, None, &n, &[0.743]).unwrap(),
            node_constants(&without, &state, None, &n, &[0.743]).unwrap(),
            "{state:?}"
        );
    }
}

#[test]
fn history_residuals_scale_with_the_weight_of_the_variable() {
    // w = (H − H⁰) / (Hᵇ − H⁰) を掛けた R を足す。H = 0.702 は w = 0.5
    let table = history_table(&[]);
    let n = [5.5e-4, 1.6e-5];
    let state = at(600.0, 20.0);
    let base = node_constants(&table, &state, None, &n, &[0.743]).unwrap();
    let half = node_constants(&table, &state, None, &n, &[0.702]).unwrap();
    let r = &table.histories[0].rods[0].residuals[1];
    let w = (0.702 - 0.743) / (0.661 - 0.743);
    let added = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a - b).collect::<Vec<_>>();
    let scaled = |x: &[f64]| x.iter().map(|x| w * x).collect::<Vec<_>>();
    assert!(close(&added(&half.total, &base.total), &scaled(&r.total)));
    assert!(close(
        &added(&half.nu_fission, &base.nu_fission),
        &scaled(&r.nu_fission)
    ));
    assert!(close(
        &added(&half.diffusion, &base.diffusion),
        &scaled(&r.diffusion)
    ));
    assert!(close(
        &added(&half.scatter.concat(), &base.scatter.concat()),
        &scaled(&r.scatter.concat())
    ));
}

#[test]
fn history_residuals_that_the_table_lacks_are_errors() {
    // 制御棒を入れた値の R を持たないテーブルで、制御棒の値 aic のノードが公称でない履歴を持つ
    let table = history_table(&[ROD_OUT, "aic"]);
    let n = [5.5e-4, 1.6e-5];
    let state = at(600.0, 20.0);
    let err = node_constants(&table, &state, Some("aic"), &n, &[0.70]).unwrap_err();
    assert!(err.contains("aic"), "{err}");
    // 公称の履歴なら R は要らない
    assert!(node_constants(&table, &state, Some("aic"), &n, &[0.743]).is_ok());
    assert!(node_constants(&table, &state, None, &n, &[0.70]).is_ok());
    // ノードの履歴の変数の数が core の変数と違う
    let err = node_constants(&table, &state, None, &n, &[]).unwrap_err();
    assert!(err.contains("history"), "{err}");
    // core の持たない履歴の変数の R（検査の後に書き換えたもの）
    let mut unknown = history_table(&[]);
    unknown.histories[0].variable = "rods".into();
    let err = node_constants(&unknown, &state, None, &n, &[0.70]).unwrap_err();
    assert!(err.contains("rods"), "{err}");
}

#[test]
fn history_branch_points_deplete_with_the_branch_rate_cross_sections() {
    // w = 1、条件 j では、反応率の微視的断面積は履歴の分岐の自分の条件の値（ADR-0215 の案 S2）
    let table = history_table(&[]);
    for (k, &b) in [0.0, 20.0].iter().enumerate() {
        let m = node_micro(&table, &at(600.0, b), None, &[0.661]).unwrap();
        for (got, want) in m.nuclides.iter().zip(&branch_sigma(k)) {
            assert!(close(&got.fission, &want.fission), "{b}");
            assert!(
                close(
                    &got.reactions[0].channels[0],
                    &want.reactions[0].channels[0]
                ),
                "{b}"
            );
            assert!(
                close(
                    &[got.thermal_fission.unwrap()],
                    &[want.thermal_fission.unwrap()]
                ),
                "{b}: {:?} {:?}",
                got.thermal_fission,
                want.thermal_fission
            );
        }
        // 基準の数密度は公称の履歴のまま
        assert_eq!(m.densities, NOMINAL_DENSITIES[2 * k].to_vec());
    }
}

#[test]
fn nominal_histories_leave_the_rate_cross_sections_bit_identical() {
    // w = 0 では、反応率の微視的断面積も履歴の補正のないテーブルとビットで同じ
    let with = history_table(&[]);
    let mut without = with.clone();
    without.histories.clear();
    for state in [
        at(565.0, 10.0),
        at(565.0, 20.0),
        at(600.0, 7.0),
        at(580.0, 15.0),
    ] {
        assert_eq!(
            node_micro(&with, &state, None, &[0.743]).unwrap(),
            interpolate_micro(&without, &state, None).unwrap(),
            "{state:?}"
        );
    }
}

#[test]
fn rate_residuals_scale_with_the_weight_and_leave_the_linear_correction() {
    // w = 0.5 で σ + w Δσ。収率の境界より下の割合は、境界より下の断面積 + w Δ を補正した σf の和で割ったもの
    let table = history_table(&[]);
    let state = at(600.0, 20.0);
    let base = interpolate_micro(&table, &state, None).unwrap();
    let half = node_micro(&table, &state, None, &[0.702]).unwrap();
    let w = (0.702 - 0.743) / (0.661 - 0.743);
    let rates = &table.histories[0].rods[0].micro[1];
    for ((h, b), r) in half.nuclides.iter().zip(&base.nuclides).zip(rates) {
        let expected: Vec<f64> = b
            .fission
            .iter()
            .zip(&r.fission)
            .map(|(x, d)| x + w * d)
            .collect();
        assert!(close(&h.fission, &expected));
        let capture: Vec<f64> = b.reactions[0].channels[0]
            .iter()
            .zip(&r.reactions[0].channels[0])
            .map(|(x, d)| x + w * d)
            .collect();
        assert!(close(&h.reactions[0].channels[0], &capture));
        let below = b.thermal_fission.unwrap() * b.fission.iter().sum::<f64>()
            + w * r.thermal_fission.unwrap();
        let fraction = below / expected.iter().sum::<f64>();
        assert!(close(&[h.thermal_fission.unwrap()], &[fraction]));
        // 線形の補正の係数（全断面積、νσf、散乱、微分）は補正しない
        assert_eq!(h.total, b.total);
        assert_eq!(h.nu_fission, b.nu_fission);
        assert_eq!(h.scatter, b.scatter);
        assert_eq!(h.derivative, b.derivative);
    }
    // ノードの少数群定数は反応率の差によらない（線形の補正の係数と R は変えない）
    let mut zero = table.clone();
    for row in &mut zero.histories[0].rods[0].micro {
        for r in row.iter_mut() {
            r.fission = vec![0.0; 2];
            r.reactions[0].channels[0] = vec![0.0; 2];
            r.thermal_fission = Some(0.0);
        }
    }
    let n = [5.5e-4, 1.6e-5];
    assert_eq!(
        node_constants(&table, &state, None, &n, &[0.702]).unwrap(),
        node_constants(&zero, &state, None, &n, &[0.702]).unwrap()
    );
}

#[test]
fn rate_residuals_that_give_negative_cross_sections_or_fractions_are_errors() {
    let state = at(600.0, 20.0);
    let mut table = history_table(&[]);
    table.histories[0].rods[0].micro[1][0].fission[1] = -1.0e3;
    let err = node_micro(&table, &state, None, &[0.661]).unwrap_err();
    assert!(err.contains("negative") && err.contains("U-235"), "{err}");
    let mut table = history_table(&[]);
    table.histories[0].rods[0].micro[1][1].reactions[0].channels[0][0] = -1.0e3;
    let err = node_micro(&table, &state, None, &[0.661]).unwrap_err();
    assert!(err.contains("negative") && err.contains("Pu-240"), "{err}");
    let mut table = history_table(&[]);
    table.histories[0].rods[0].micro[1][0].thermal_fission = Some(1.0e3);
    let err = node_micro(&table, &state, None, &[0.661]).unwrap_err();
    assert!(err.contains("fraction"), "{err}");
    // 制御棒の値の残りを持たないのに w が 0 でない
    let table = history_table(&[ROD_OUT, "aic"]);
    let err = node_micro(&table, &state, Some("aic"), &[0.70]).unwrap_err();
    assert!(err.contains("aic"), "{err}");
}
