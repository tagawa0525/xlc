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
