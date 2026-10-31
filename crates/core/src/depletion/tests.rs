use std::path::Path;

use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, Table, product_points};

use super::*;
use crate::case::{solve_at, solve_burnups};
use crate::input::{BoronSearch, Case, CaseState, Cycle, Inventory, Layout, MidpointBoron, Plant};
use crate::nem::Options;

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
[assembly.A3]
design = "A"
[assembly.A4]
design = "A"
[assembly.B1]
design = "B"
[assembly.B2]
design = "B"
"#;

/// 4 体とも設計 A（一様な炉心）
const UNIFORM: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["A1  A2", "A3  A4"]
"#;

/// A と B を互い違いに
const MIXED: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["A1  B1", "B2  A2"]
"#;

const CASE: &str = r#"
cycle = "cycle.toml"
axial-mesh = [0.0, 10.0, 20.0]
radial-nodes = 1
boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], z = ["reflective", "reflective"] }

[state.hfp]
boron = 1200
b10 = 0.199
moderator = 565
fuel = 565

[depletion]
specific-power = 40.0
steps = [{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]
history = "history.msgpack"
"#;

const BURNUP: [f64; 3] = [0.0, 10.0, 20.0];
const PPM: [f64; 2] = [1000.0, 1500.0];
const TEMPERATURE: [f64; 2] = [560.0, 570.0];

/// 熱群の吸収のホウ素濃度の係数 (1/cm/ppm)
const WORTH: f64 = 4.0e-5;

/// 2 群の定数。核分裂は燃焼度で減り（`fissile` × (1 − 0.02 B)）、熱群の吸収はホウ素濃度で変わる（`worth` × (ppm − 1000)）
fn few(fissile: f64, burnup: f64, ppm: f64, worth: f64) -> FewGroup {
    let absorption = vec![0.01, 0.06 + worth * (ppm - 1000.0)];
    let removal = 0.02;
    let strength = fissile * (1.0 - 0.02 * burnup);
    let nu_fission = vec![0.005 * strength, 0.1 * strength];
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.0,
        total: vec![absorption[0] + 0.5 + removal, absorption[1] + 1.0],
        fission: nu_fission.iter().map(|x| x / 2.43).collect(),
        nu_fission,
        absorption,
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.5, removal], vec![0.0, 1.0]],
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

/// 燃焼度の軸を持つテーブル（燃焼度、ホウ素、減速材、燃料の順）
fn table(fissile: f64, heavy_metal: Option<f64>) -> Table {
    table_with_worth(fissile, heavy_metal, WORTH)
}

/// `table` の、熱群の吸収のホウ素濃度の係数を `worth` にしたもの
fn table_with_worth(fissile: f64, heavy_metal: Option<f64>, worth: f64) -> Table {
    let mut states = Vec::new();
    for burnup in BURNUP {
        for ppm in PPM {
            for _ in TEMPERATURE {
                for _ in TEMPERATURE {
                    states.push(few(fissile, burnup, ppm, worth));
                }
            }
        }
    }
    let mut table = Table {
        points: Vec::new(),
        axes: vec![
            Axis::Burnup {
                points: BURNUP.to_vec(),
            },
            Axis::boron(PPM.to_vec(), 0.199),
            Axis::moderator(TEMPERATURE.to_vec(), vec![0.75, 0.74]),
            Axis::fuel(TEMPERATURE.to_vec()),
        ],
        split: 1,
        states,
        heavy_metal,
        micro: None,
        histories: Vec::new(),
    };
    table.points = product_points(&table.axes);
    table
}

/// 炉心の並び、`layout.files` の順のテーブル（A は核分裂 1.0、B は 0.8。重金属は A が 400、B が 380 g/cm）、ケース
fn core(case: &str, cycle: &str, heavy_metal: Option<f64>) -> (Layout, Vec<Table>, Case) {
    let base = Path::new("dir");
    let case = Case::parse(case, base).unwrap();
    let cycle = Cycle::parse(cycle, base).unwrap();
    let plant = Plant::parse(PLANT, base).unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = layout
        .files
        .iter()
        .map(|f| {
            if f.ends_with("a.msgpack") {
                table(1.0, heavy_metal)
            } else {
                table(0.8, heavy_metal.map(|m| 0.95 * m))
            }
        })
        .collect();
    (layout, tables, case)
}

/// 固有値の反復を詰める（出力の割合の反復の誤差が、比べる燃焼度の差より小さくなるように）
fn tight() -> Options {
    Options {
        k_tolerance: 1.0e-13,
        flux_tolerance: 1.0e-12,
        ..Options::default()
    }
}

fn run(case: &str, cycle: &str) -> Result<CoreHistory, CaseError> {
    let (layout, tables, case) = core(case, cycle, Some(400.0));
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |_| {},
    )
}

#[test]
fn advance_shares_the_core_increment_by_power_per_heavy_metal() {
    // 質量 1、3、0 g、出力の割合 1/2、1/2、0。炉心の平均の増分 1 MWd/kg は 4 g に 4 g MWd/kg のエネルギー
    let b = advance(&[5.0, 6.0, 0.0], &[1.0, 3.0, 0.0], 1.0, &[2.0, 2.0, 0.0]).unwrap();
    assert_eq!(b[0], 5.0 + 0.5 * 4.0 / 1.0);
    assert!((b[1] - (6.0 + 0.5 * 4.0 / 3.0)).abs() <= 1e-15, "{}", b[1]);
    assert_eq!(b[2], 0.0);
}

#[test]
fn advance_rejects_fission_without_heavy_metal_and_cores_without_fission() {
    assert!(advance(&[0.0, 0.0], &[1.0, 0.0], 1.0, &[1.0, 1.0]).is_err());
    assert!(advance(&[0.0, 0.0], &[1.0, 1.0], 1.0, &[0.0, 0.0]).is_err());
}

#[test]
fn midpoint_steps_converge_at_second_order() {
    // 2 ノード（質量 1 g ずつ）。核分裂率がノードの燃焼度で exp(−a B) と減る（出力分布が燃焼で動く）
    let rates = |b: &[f64]| -> Result<Vec<f64>, CaseError> {
        Ok(vec![(-0.05 * b[0]).exp(), (-0.02 * b[1]).exp()])
    };
    let burn = |steps: usize| {
        let mut b = vec![0.0, 0.0];
        let delta = 20.0 / steps as f64;
        for _ in 0..steps {
            let start = rates(&b).unwrap();
            b = midpoint(&b, &[1.0, 1.0], delta, &start, rates).unwrap();
        }
        b
    };
    let reference = burn(4096);
    let error = |steps: usize| (burn(steps)[0] - reference[0]).abs();
    let ratio = error(8) / error(16);
    assert!((3.5..4.5).contains(&ratio), "{ratio}");
    // 炉心の平均の燃焼度はステップによらず 20 MWd/kg（エネルギーの保存）
    let b = burn(3);
    assert!((0.5 * (b[0] + b[1]) - 20.0).abs() <= 1e-12, "{b:?}");
}

#[test]
fn uniform_cores_burn_at_the_specific_power() {
    let history = run(CASE, UNIFORM).unwrap();
    assert_eq!(history.assemblies, ["A1", "A2", "A3", "A4"]);
    assert_eq!(history.axial_mesh, [0.0, 10.0, 20.0]);
    assert_eq!(history.specific_power, 40.0);
    // 解いた状態（比べるツールがテーブルの同じ状態を引く）
    assert_eq!(history.states.len(), 1);
    assert_eq!(
        history.states["hfp"],
        CoreState {
            boron: Some(1200.0),
            search: None,
            b10: 0.199,
            moderator: 565.0,
            fuel: 565.0,
            banks: Default::default(),
        }
    );
    // 集合体ごとの 2 層。質量は 400 g/cm × 10 cm
    assert_eq!(history.heavy_metal, vec![vec![4000.0; 2]; 4]);
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let efpd = [0.0, 50.0, 125.0, 250.0];
    assert_eq!(history.points.len(), efpd.len());
    for (point, e) in history.points.iter().zip(efpd) {
        // 燃焼度 = 出力密度 (W/g) × EFPD / 1000。ノードの値は出力の割合の反復の誤差（判定 1e-12）の分だけずれうる
        let want = 40.0 * e / 1000.0;
        assert_eq!(point.efpd, e);
        // 与えたホウ素濃度で 1 回ずつ解く
        assert_eq!((point.boron, point.solves), (1200.0, 1));
        assert!(
            (point.burnup - want).abs() <= 1e-12 * want.max(1.0),
            "{point:?}"
        );
        let close = |b: f64| (b - want).abs() <= 1e-10 * want.max(1.0);
        assert!(
            point.assembly_burnups.iter().all(|&b| close(b)),
            "{point:?}"
        );
        assert!(
            point.burnups.iter().flatten().all(|&b| close(b)),
            "{point:?}"
        );
        // 出力の割合は 4 体 × 2 層で一様
        for p in point.powers.iter().flatten() {
            assert!((p - 0.125).abs() <= 1e-9, "{point:?}");
        }
        // k は、その燃焼度の一様な状態で解いた k
        let mut state = case.states[0].at(1200.0);
        state.burnup = Some(want);
        let uniform = solve_at(&layout, &tables, Some(&state), &tight()).unwrap();
        assert!(
            (point.k - uniform.solution.k).abs() <= 1e-11,
            "{} vs {}",
            point.k,
            uniform.solution.k
        );
    }
}

#[test]
fn nodes_burn_by_their_own_power() {
    let history = run(CASE, MIXED).unwrap();
    assert_eq!(history.assemblies, ["A1", "A2", "B1", "B2"]);
    let mass: f64 = history.heavy_metal.iter().flatten().sum();
    for pair in history.points.windows(2) {
        let (before, after) = (&pair[0], &pair[1]);
        // 増えたエネルギー Σ M ΔB は、出力密度 × 炉心の質量 × EFPD
        let energy: f64 = history
            .heavy_metal
            .iter()
            .flatten()
            .zip(
                after
                    .burnups
                    .iter()
                    .flatten()
                    .zip(before.burnups.iter().flatten()),
            )
            .map(|(m, (b1, b0))| m * (b1 - b0))
            .sum();
        let want = 40.0 * (after.efpd - before.efpd) / 1000.0 * mass;
        assert!((energy - want).abs() <= 1e-12 * want, "{energy} vs {want}");
        assert!((after.burnup - 40.0 * after.efpd / 1000.0).abs() <= 1e-12 * after.burnup);
    }
    // 核分裂の強い A は B より出力の割合が高く、早く燃える
    let last = history.points.last().unwrap();
    assert!(
        last.assembly_burnups[0] > last.assembly_burnups[2],
        "{last:?}"
    );
    assert!(last.powers[0][0] > last.powers[2][0], "{last:?}");
    let total: f64 = last.powers.iter().flatten().sum();
    assert!((total - 1.0).abs() <= 1e-12, "{total}");
}

#[test]
fn fissile_tables_without_heavy_metal_cannot_deplete() {
    // B のテーブルだけが質量を持たない。解く前に、そのファイルを示して止まる
    let (layout, mut tables, case) = core(CASE, MIXED, Some(400.0));
    for (file, table) in layout.files.iter().zip(&mut tables) {
        if file.ends_with("b.msgpack") {
            table.heavy_metal = None;
        }
    }
    let mut points = 0;
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &Options::default(),
        |_| points += 1,
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("b.msgpack") && error.contains("heavy_metal"),
        "{error}"
    );
    assert_eq!(points, 0);
}

#[test]
fn solving_by_node_burnups_rejects_a_state_burnup() {
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let mut state = case.states[0].at(1200.0);
    state.burnup = Some(5.0);
    let burnups = vec![0.0; layout.node_material.len()];
    assert!(crate::case::solve_burnups(&layout, &tables, &state, &burnups, &tight()).is_err());
    state.burnup = None;
    assert!(crate::case::solve_burnups(&layout, &tables, &state, &burnups, &tight()).is_ok());
}

#[test]
fn reflector_sites_take_no_burnup() {
    // 集合体 1 体の両側に反射体。反射体のテーブルは燃焼度の軸も核分裂も持たない
    let base = Path::new("dir");
    let plant = Plant::parse(
        "pitch = 20.0\nheight = 20.0\nmap = [\"R  *  R\"]\n[reflector.R]\nbounds = [0.0, 20.0]\nlayers = [\"r.msgpack\"]\n",
        base,
    )
    .unwrap();
    let cycle = Cycle::parse(
        "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\nloading = [\".  A1  .\"]\n",
        base,
    )
    .unwrap();
    let case = Case::parse(
        &CASE.replace(
            "[\"reflective\", \"reflective\"], y",
            "[\"zero-flux\", \"zero-flux\"], y",
        ),
        base,
    )
    .unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables: Vec<Table> = layout
        .files
        .iter()
        .map(|f| {
            if f.ends_with("r.msgpack") {
                let mut t = table(0.0, None);
                // 燃焼度の軸を除いた、最初の燃焼度の点の状態
                t.axes.remove(0);
                t.points = product_points(&t.axes);
                t.states.truncate(t.states.len() / BURNUP.len());
                t
            } else {
                table(1.0, Some(400.0))
            }
        })
        .collect();
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    // 履歴は集合体だけを持つ
    assert_eq!(history.assemblies, ["A1"]);
    let last = history.points.last().unwrap();
    assert!((last.burnup - 40.0 * 250.0 / 1000.0).abs() <= 1e-12 * last.burnup);
    let total: f64 = last.powers.iter().flatten().sum();
    assert!((total - 1.0).abs() <= 1e-12, "{total}");
}

#[test]
fn burnups_past_the_table_are_errors() {
    // 40 W/g × 600 EFPD = 24 MWd/kg はテーブルの 20 MWd/kg を越える（外挿しない。ADR-0187）
    let case = CASE.replace("{ efpd = 250.0 }", "{ efpd = 600.0 }");
    assert!(run(&case, UNIFORM).is_err());
}

#[test]
fn solve_errors_name_the_point_or_midpoint() {
    // 中点（14.5 MWd/kg）はテーブルの中で、終わりの点（24 MWd/kg）で内挿のエラーになる。探索のエラーでなくても、どの解かを示す
    let case = CASE.replace("{ efpd = 250.0 }", "{ efpd = 600.0 }");
    let error = run(&case, UNIFORM).unwrap_err().to_string();
    assert!(error.contains("the point at 600 EFPD"), "{error}");
    let error = run(&searching(&case, 1.0e-9), UNIFORM)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("the midpoint of the step to 600 EFPD"),
        "{error}"
    );
    // 表示が原因を含むので、エラーの連鎖（main が `{:#}` で出す）で原因を繰り返さない
    let error = run(&case, UNIFORM).unwrap_err();
    let mut chain = vec![error.to_string()];
    let mut cause = std::error::Error::source(&error);
    while let Some(e) = cause {
        chain.push(e.to_string());
        cause = e.source();
    }
    let shown = chain.join(": ");
    let inner = chain[0]
        .strip_prefix("the point at 600 EFPD: ")
        .unwrap_or_default();
    assert!(
        !inner.is_empty() && shown.matches(inner).count() == 1,
        "{shown}"
    );
}

#[test]
fn histories_round_trip_through_their_file() {
    let history = run(CASE, MIXED).unwrap();
    // 省いた中点のホウ素は、ADR-0191 の探索
    assert_eq!(history.midpoint_boron, MidpointBoron::Search);
    let bytes = history.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), history);
    let start = CoreHistory {
        midpoint_boron: MidpointBoron::Start,
        ..history.clone()
    };
    let bytes = start.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), start);
    let value = xlc_format::msgpack::read_value(&bytes).unwrap();
    let format = value
        .as_map()
        .unwrap()
        .iter()
        .find(|(k, _)| k.as_str() == Some("format"))
        .and_then(|(_, v)| v.as_str());
    assert_eq!(format, Some("xlc-core-history"));
}

#[test]
fn core_average_burnups_do_not_drift_over_many_steps() {
    // 集合体 1 体・ノード 1 つ。EFPD は 2A の履歴の 0.1〜20 MWd/kg の点（38.41 W/g）。ステップの増分を足し合わせると、丸めで
    // 20.000000000000004 とテーブルの最後の点 20 を越える。終わりの累積の EFPD から求めると越えない
    let efpd = [
        2.603488674824265,
        6.508721687060661,
        13.017443374121326,
        19.526165061181988,
        26.034886748242652,
        39.052330122363976,
        52.069773496485304,
        65.08721687060662,
        78.10466024472794,
        91.12210361884927,
        104.13954699297058,
        130.17443374121325,
        156.20932048945588,
        182.24420723769853,
        208.27909398594116,
        234.31398073418382,
        260.34886748242644,
        286.38375423066907,
        312.41864097891175,
        338.4535277271544,
        364.488414475397,
        390.5233012236397,
        416.5581879718823,
        468.62796146836763,
        520.6977349648529,
    ];
    let steps: Vec<String> = efpd.iter().map(|e| format!("{{ efpd = {e:?} }}")).collect();
    let case = CASE
        .replace("axial-mesh = [0.0, 10.0, 20.0]", "axial-mesh = [0.0, 20.0]")
        .replace("specific-power = 40.0", "specific-power = 38.41")
        .replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            &format!("[{}]", steps.join(", ")),
        );
    let base = Path::new("dir");
    let case = Case::parse(&case, base).unwrap();
    let cycle = Cycle::parse(
        "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\nloading = [\"A1\"]\n",
        base,
    )
    .unwrap();
    let plant = Plant::parse("pitch = 20.0\nheight = 20.0\nmap = [\"*\"]\n", base).unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = vec![table(1.0, Some(400.0))];
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &Options::default(),
        |_| {},
    )
    .unwrap();
    for (point, e) in history.points.iter().skip(1).zip(efpd) {
        let want = 38.41 * e * 1.0e-3;
        assert!(
            (point.burnup - want).abs() <= 4.0 * f64::EPSILON * want,
            "{} vs {want}",
            point.burnup
        );
    }
}

#[test]
fn points_are_reported_as_they_are_solved() {
    // 最後のステップ（24 MWd/kg）はテーブルを越えて止まるが、それまでの点は解いた順に渡してある
    let case = CASE.replace("{ efpd = 250.0 }", "{ efpd = 600.0 }");
    let (layout, tables, case) = core(&case, UNIFORM, Some(400.0));
    let mut seen = Vec::new();
    let result = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |p| seen.push(p.efpd),
    );
    assert!(result.is_err());
    assert_eq!(seen, [0.0, 50.0, 125.0]);
    // 解き終えたケースの点は、返す履歴の点と同じ
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let mut points = Vec::new();
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |p| points.push(p.clone()),
    )
    .unwrap();
    assert_eq!(points, history.points);
}

/// 状態のホウ素を探索にしたケース（目標の k = 1、範囲はテーブルのホウ素の軸の両端）
fn searching(case: &str, tolerance: f64) -> String {
    case.replace(
        "boron = 1200\n",
        &format!(
            "search = {{ boron = [1000.0, 1500.0], k = 1.0, tolerance = {tolerance:e}, max-solves = 30 }}\n"
        ),
    )
}

/// 集合体 1 体（設計 A）、ノード 1 つの炉心
fn single(case: &str) -> (Layout, Vec<Table>, Case) {
    let base = Path::new("dir");
    let case = case.replace("axial-mesh = [0.0, 10.0, 20.0]", "axial-mesh = [0.0, 20.0]");
    let case = Case::parse(&case, base).unwrap();
    let cycle = Cycle::parse(
        "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\nloading = [\"A1\"]\n",
        base,
    )
    .unwrap();
    let plant = Plant::parse("pitch = 20.0\nheight = 20.0\nmap = [\"*\"]\n", base).unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    (layout, vec![table(1.0, Some(400.0))], case)
}

#[test]
fn searched_points_are_critical_and_record_their_boron() {
    // 40 W/g で 2、5 MWd/kg。k はホウ素濃度でも燃焼度でも下がるので、求めるホウ素は燃焼で下がる
    let case = searching(CASE, 1.0e-9).replace(
        "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
        "[{ efpd = 50.0 }, { efpd = 125.0 }]",
    );
    let history = run(&case, MIXED).unwrap();
    // 履歴の状態は探索の設定を持ち、与えた濃度を持たない
    let state = &history.states["hfp"];
    assert_eq!(state.boron, None);
    assert_eq!(
        state.search,
        Some(BoronSearch {
            bracket: [1000.0, 1500.0],
            k: 1.0,
            tolerance: 1.0e-9,
            max_solves: 30,
        })
    );
    assert_eq!(history.points.len(), 3);
    for point in &history.points {
        assert!((point.k - 1.0).abs() <= 1.0e-9, "{point:?}");
        assert!((1000.0..=1500.0).contains(&point.boron), "{point:?}");
        // 両端と、間の少なくとも 1 点
        assert!((3..=30).contains(&point.solves), "{point:?}");
    }
    assert!(
        history.points.windows(2).all(|w| w[1].boron < w[0].boron),
        "{:?}",
        history.points
    );
    // 燃焼度の進め方は ADR-0190 のまま（エネルギーの保存）
    let last = history.points.last().unwrap();
    assert!((last.burnup - 40.0 * 125.0 / 1000.0).abs() <= 1e-12 * last.burnup);
}

#[test]
fn searched_points_match_the_state_search_at_the_same_uniform_burnup() {
    // 集合体 1 体では燃焼度が一様なので、点の探索は、その燃焼度を状態に与えた探索（ADR-0145）と同じ計算になる
    let case = searching(CASE, 1.0e-9).replace(
        "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
        "[{ efpd = 50.0 }, { efpd = 125.0 }]",
    );
    let (layout, tables, case) = single(&case);
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    for point in &history.points {
        // 解いたのはノードの燃焼度（質量の重みの平均の `burnup` とビットで同じとは限らない）
        let state = CaseState {
            burnup: Some(point.burnups[0][0]),
            ..case.states[0].clone()
        };
        let solved = crate::case::solve_state(&layout, &tables, &state, &tight()).unwrap();
        let searched = solved.searched.unwrap();
        assert_eq!(
            (point.boron, point.solves, point.k),
            (searched.ppm, searched.solves, solved.solution.k)
        );
    }
}

#[test]
fn searched_midpoints_keep_second_order() {
    // B のホウ素の係数を A の 4 倍にして、出力の割合がホウ素濃度で動く炉心。中点でも臨界にすれば、ステップを半分にした
    // 燃焼度の差の比は 2 次の 4 に近い（中点を初めの点のホウ素で解くと、ホウ素を通した部分が 1 次になる）
    let burn = |steps: usize| {
        let efpd: Vec<String> = (1..=steps)
            .map(|i| format!("{{ efpd = {:?} }}", 150.0 * i as f64 / steps as f64))
            .collect();
        let case = searching(CASE, 1.0e-12).replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            &format!("[{}]", efpd.join(", ")),
        );
        let (layout, mut tables, case) = core(&case, MIXED, Some(400.0));
        for (file, table) in layout.files.iter().zip(&mut tables) {
            if file.ends_with("b.msgpack") {
                *table = table_with_worth(0.8, Some(380.0), 4.0 * WORTH);
            }
        }
        let history = deplete(
            &layout,
            &tables,
            &case.states,
            case.depletion.as_ref().unwrap(),
            None,
            None,
            &tight(),
            |_| {},
        )
        .unwrap();
        history.points.last().unwrap().assembly_burnups[0]
    };
    let (b2, b4, b8) = (burn(2), burn(4), burn(8));
    let ratio = (b2 - b4) / (b4 - b8);
    assert!((3.5..4.5).contains(&ratio), "{ratio} ({b2}, {b4}, {b8})");
}

#[test]
fn searches_that_do_not_bracket_stop_at_the_step() {
    // 10 MWd/kg のステップ（125 → 250 EFPD）の中点 7.5 MWd/kg で、1000 ppm でも k < 1。ステップと中点を示して止まり、それまでの
    // 点は渡してある
    let case = searching(CASE, 1.0e-9);
    let (layout, tables, case) = core(&case, MIXED, Some(400.0));
    let mut seen = Vec::new();
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |p| seen.push(p.efpd),
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("the midpoint of the step to 250 EFPD") && error.contains("bracket"),
        "{error}"
    );
    assert_eq!(seen, [0.0, 50.0, 125.0]);
}

#[test]
fn search_tolerances_within_the_eigenvalue_tolerance_are_errors_before_solving() {
    // 固有値の反復の判定（既定の 1e-8）以下の許容は区別できない（ADR-0145）。解く前に止まる
    let case = searching(CASE, 1.0e-9);
    let (layout, tables, case) = core(&case, MIXED, Some(400.0));
    let mut points = 0;
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &Options::default(),
        |_| points += 1,
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("tolerance"), "{error}");
    assert_eq!(points, 0);
}

#[test]
fn history_states_hold_either_the_boron_or_the_search() {
    let given = run(CASE, MIXED).unwrap();
    let searched = run(
        &searching(CASE, 1.0e-9).replace(", { efpd = 250.0 }", ""),
        MIXED,
    )
    .unwrap();
    for history in [&given, &searched] {
        let bytes = history.to_msgpack().unwrap();
        assert_eq!(&CoreHistory::from_msgpack(&bytes).unwrap(), history);
    }
    // 両方を持つ状態、どちらもない状態は読めない
    let with_state = |history: &CoreHistory, state: CoreState| CoreHistory {
        states: [("hfp".to_string(), state)].into(),
        ..history.clone()
    };
    let both = with_state(
        &searched,
        CoreState {
            boron: Some(1200.0),
            ..searched.states["hfp"].clone()
        },
    );
    assert!(CoreHistory::from_msgpack(&both.to_msgpack().unwrap()).is_err());
    let neither = with_state(
        &given,
        CoreState {
            boron: None,
            search: None,
            ..given.states["hfp"].clone()
        },
    );
    assert!(CoreHistory::from_msgpack(&neither.to_msgpack().unwrap()).is_err());
    // 探索の設定は入力と同じ検査を通る（範囲の下端 < 上端、正の目標と許容、両端と間の 1 点を解ける回数）
    let bad = |search: BoronSearch| {
        with_state(
            &searched,
            CoreState {
                search: Some(search),
                ..searched.states["hfp"].clone()
            },
        )
    };
    let good = searched.states["hfp"].search.unwrap();
    for search in [
        BoronSearch {
            bracket: [1500.0, 1000.0],
            ..good
        },
        BoronSearch { k: 0.0, ..good },
        BoronSearch {
            tolerance: -1.0,
            ..good
        },
        BoronSearch {
            max_solves: 2,
            ..good
        },
    ] {
        let bytes = bad(search).to_msgpack().unwrap();
        assert!(CoreHistory::from_msgpack(&bytes).is_err(), "{search:?}");
    }
}

#[test]
fn midpoints_at_the_starting_boron_fall_to_first_order() {
    // ADR-0191 の案 S1（中点を、ステップの初めの点で探したホウ素で解く）。`searched_midpoints_keep_second_order` と同じ炉心で、
    // ホウ素を通した中点の出力の割合の誤差がステップの長さに比例し、ステップを半分にした燃焼度の差の比が 1 次の 2 に近い
    let burn = |steps: usize| {
        let efpd: Vec<String> = (1..=steps)
            .map(|i| format!("{{ efpd = {:?} }}", 150.0 * i as f64 / steps as f64))
            .collect();
        let case = searching(CASE, 1.0e-12).replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            &format!("[{}]", efpd.join(", ")),
        );
        let (layout, mut tables, case) = core(&case, MIXED, Some(400.0));
        for (file, table) in layout.files.iter().zip(&mut tables) {
            if file.ends_with("b.msgpack") {
                *table = table_with_worth(0.8, Some(380.0), 4.0 * WORTH);
            }
        }
        let state = &case.states[0];
        let Boron::Search(search) = state.boron else {
            panic!("the state searches")
        };
        let options = tight();
        let critical = |burnups: &[f64]| {
            let (solved, found) = crate::case::search_boron(&search, |ppm| {
                solve_burnups(&layout, &tables, &state.at(ppm), burnups, &options)
            })
            .unwrap();
            (fission_rates(&solved.problem, &solved.solution), found.ppm)
        };
        let [nx, ny, nz] = layout.mesh.counts();
        let masses: Vec<f64> = (0..nx * ny * nz)
            .map(|site| {
                let [_, _, z] = layout.mesh.position(site);
                layout
                    .site_file(site)
                    .and_then(|f| tables[f].heavy_metal)
                    .map_or(0.0, |m| m * layout.mesh.dz[z])
            })
            .collect();
        let mut burnups = vec![0.0; masses.len()];
        let (mut rates, mut ppm) = critical(&burnups);
        for i in 1..=steps {
            let target = 40.0 * 150.0 * i as f64 / steps as f64 * 1.0e-3;
            let increment = target - weighted(&burnups, &masses);
            burnups = midpoint(&burnups, &masses, increment, &rates, |b| {
                let solved = solve_burnups(&layout, &tables, &state.at(ppm), b, &options)?;
                Ok(fission_rates(&solved.problem, &solved.solution))
            })
            .unwrap();
            (rates, ppm) = critical(&burnups);
        }
        // 入力で案 S1 を選んだ燃焼は、上の組み立てと同じ燃焼度になる（ステップの目標の燃焼度の掛け算の順だけが違う）
        let mut start = case.depletion.clone().unwrap();
        start.midpoint_boron = MidpointBoron::Start;
        let selected = deplete(
            &layout,
            &tables,
            &case.states,
            &start,
            None,
            None,
            &options,
            |_| {},
        )
        .unwrap();
        assert_eq!(selected.midpoint_boron, MidpointBoron::Start);
        let selected_a1 = &selected.points.last().unwrap().burnups[0];
        for (z, b) in selected_a1.iter().enumerate() {
            let built = burnups[layout.mesh.site([0, ny - 1, z])];
            assert!(
                (b - built).abs() <= 1.0e-12 * built,
                "node {z}: {b} vs {built}"
            );
        }
        assert_eq!(selected_a1.len(), nz);
        // 集合体 A1 の燃焼度（履歴の `assembly_burnups[0]` と同じ量）
        let history = deplete(
            &layout,
            &tables,
            &case.states,
            case.depletion.as_ref().unwrap(),
            None,
            None,
            &options,
            |_| {},
        )
        .unwrap();
        let a1 = |b: &[f64]| {
            let sites: Vec<usize> = (0..b.len())
                .filter(|&site| {
                    let [x, y, _] = layout.mesh.position(site);
                    // 図の上の行が +y。A1 は左上
                    x == 0 && y == ny - 1
                })
                .collect();
            sites.iter().map(|&s| b[s] * masses[s]).sum::<f64>()
                / sites.iter().map(|&s| masses[s]).sum::<f64>()
        };
        (
            a1(&burnups),
            history.points.last().unwrap().assembly_burnups[0],
        )
    };
    let ((s2, c2), (s4, c4), (s8, c8)) = (burn(2), burn(4), burn(8));
    let ratio = (s2 - s4) / (s4 - s8);
    assert!((1.5..2.5).contains(&ratio), "{ratio} ({s2}, {s4}, {s8})");
    // 中点でも探す形（採用）との差は、8 ステップでも燃焼度の 0.1% を超える
    assert!((s8 - c8).abs() > 1.0e-3 * c8, "{s8} vs {c8} ({c2}, {c4})");
}

// ミクロ燃焼（ADR-0216）

use xlc_depletion::chain::{Channel, Fission, Nuclide, Product, Reaction, Yield};
use xlc_format::table::{Micro, MicroNuclide, MicroReaction};

use crate::case::Compositions;
use crate::input::{MicroInput, RateResidualWeight, Xenon};
use crate::micro::MicroChain;

/// U-235 の捕獲 / 核分裂の比（群ごとに同じなので、反応率の比は中性子束の形によらない）
const CAPTURE_TO_FISSION: f64 = 0.2;
/// U-235 の核分裂の Q 値 (eV)
const Q: f64 = 193.0e6;
/// 1 eV (J)
const EV: f64 = 1.602_176_634e-19;

fn chain_nuclide(name: &str, z: u32, a: u32) -> Nuclide {
    Nuclide {
        name: name.into(),
        z,
        a,
        state: 0,
        excitation: 0.0,
        decay: None,
        reactions: Vec::new(),
        fission: None,
    }
}

fn capture_to(nuclide: usize) -> Reaction {
    Reaction {
        mt: 102,
        channels: vec![Channel {
            final_state: None,
            products: vec![Product { nuclide, count: 1 }],
        }],
    }
}

/// Ba-140、U-235（捕獲で U-236、核分裂で Ba-140）、U-236（捕獲で U-237）、U-237 のチェーン
fn small_chain() -> Chain {
    let mut nuclides = vec![
        chain_nuclide("Ba-140", 56, 140),
        chain_nuclide("U-235", 92, 235),
        chain_nuclide("U-236", 92, 236),
        chain_nuclide("U-237", 92, 237),
    ];
    nuclides[1].reactions = vec![capture_to(2)];
    nuclides[2].reactions = vec![capture_to(3)];
    let yields = vec![Yield {
        nuclide: 0,
        value: 0.06,
    }];
    nuclides[1].fission = Some(Fission {
        q: Q,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5.0e5],
        yields: vec![yields.clone(), yields],
    });
    let chain = Chain { nuclides };
    chain.check().unwrap();
    chain
}

/// 2 群の微視的断面積。捕獲は核分裂の `capture` 倍
fn micro_nuclide(fission: [f64; 2], capture: [f64; 2], thermal: Option<f64>) -> MicroNuclide {
    MicroNuclide {
        total: vec![
            fission[0] + capture[0] + 4.0,
            fission[1] + capture[1] + 10.0,
        ],
        nu_fission: fission.iter().map(|f| 2.43 * f).collect(),
        fission: fission.to_vec(),
        scatter: vec![vec![3.9, 0.1], vec![0.0, 10.0]],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![capture.to_vec()],
        }],
        thermal_fission: thermal,
        derivative: None,
    }
}

/// 燃焼の初めの U-235 の数密度 (1/(barn cm))
const U235: f64 = 7.0e-4;

/// テーブルに、U-235 と U-236 の状態によらない微視的断面積と、燃焼度の点（0、10、20 MWd/kg）ごとの数密度を足す
fn with_micro(mut table: Table, nuclides: &[&str]) -> Table {
    let f = [1.0, 50.0];
    let u235 = micro_nuclide(f, f.map(|x| CAPTURE_TO_FISSION * x), Some(0.95));
    let u236 = micro_nuclide([0.0, 0.0], [0.0, 0.0], None);
    table.micro = Some(Micro {
        nuclides: nuclides.iter().map(|s| s.to_string()).collect(),
        densities: vec![vec![U235, 0.0], vec![6.2e-4, 1.4e-5], vec![5.5e-4, 2.6e-5]],
        states: vec![vec![u235, u236]; table.points.len()],
    });
    table.validate().unwrap();
    table
}

fn micro_case(case: &str) -> Case {
    let mut case = Case::parse(case, Path::new("dir")).unwrap();
    case.depletion.as_mut().unwrap().micro = Some(MicroInput {
        chain: "chain.msgpack".into(),
        xenon: Xenon::Tracked,
        history_weight: 1.0,
        rate_residual_weight: RateResidualWeight::History,
    });
    case
}

fn micro_run(cycle: &str, nuclides: [&[&str]; 2]) -> Result<CoreHistory, CaseError> {
    let (layout, tables, _) = core(CASE, cycle, Some(400.0));
    let case = micro_case(CASE);
    let tables: Vec<Table> = layout
        .files
        .iter()
        .zip(tables)
        .map(|(f, t)| {
            let n = if f.ends_with("a.msgpack") {
                nuclides[0]
            } else {
                nuclides[1]
            };
            with_micro(t, n)
        })
        .collect();
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
}

const TRACKED: [&str; 2] = ["U-235", "U-236"];

#[test]
fn micro_depletion_follows_the_cecm_steps_of_the_power_normalized_rates() {
    // 一様な炉心の各ノードの出力密度は 40 W/g × 400 g/cm / (20 cm × 20 cm)。核分裂性の追う核種は U-235 だけなので、1 原子あたりの
    // 核分裂率は 出力密度 / (N₂₃₅ Q) で、消滅の率は (1 + 捕獲 / 核分裂) 倍
    let history = micro_run(UNIFORM, [&TRACKED, &TRACKED]).unwrap();
    assert_eq!(history.nuclides, Some(vec!["U-235".into(), "U-236".into()]));
    let power = 40.0 * 400.0 / (20.0 * 20.0);
    let loss = |n: f64| (1.0 + CAPTURE_TO_FISSION) * power * 1.0e-24 / (n * Q * EV);
    let (mut u235, mut u236) = (U235, 0.0);
    let mut efpd = 0.0;
    for point in &history.points {
        // ステップ（ADR-0173 の CE/CM）: 初めの組成の率で中点の組成へ、中点の組成の率で初めの組成から終わりへ
        let h = (point.efpd - efpd) * 86400.0;
        let half = u235 * (-loss(u235) * 0.5 * h).exp();
        let end = u235 * (-loss(half) * h).exp();
        u236 += CAPTURE_TO_FISSION / (1.0 + CAPTURE_TO_FISSION) * (u235 - end);
        u235 = end;
        efpd = point.efpd;
        let densities = point.densities.as_ref().unwrap();
        assert_eq!(densities.len(), 4);
        for node in densities.iter().flatten() {
            let close = |a: f64, b: f64| (a - b).abs() <= 1e-10 * b.abs().max(1e-30);
            assert!(close(node[0], u235), "{} {node:?} vs {u235}", point.efpd);
            assert!(close(node[1], u236), "{} {node:?} vs {u236}", point.efpd);
        }
    }
}

#[test]
fn micro_depletion_starts_from_the_table_and_keeps_its_constants_there() {
    // 燃焼度 0 の点は、テーブルの数密度から始めるので補正が 0 で、ミクロ燃焼しない燃焼と同じ k（ADR-0216 の確認 (a)）
    let micro = micro_run(UNIFORM, [&TRACKED, &TRACKED]).unwrap();
    let macroscopic = run(CASE, UNIFORM).unwrap();
    assert_eq!(micro.points[0].k, macroscopic.points[0].k);
    assert_eq!(
        micro.points[0].densities.as_ref().unwrap()[0][0],
        [U235, 0.0]
    );
    // 燃焼すると、数密度がテーブルの公称の履歴（ここでは任意の値）から外れた分だけ k が変わる
    assert_ne!(micro.points[1].k, macroscopic.points[1].k);
    // ファイルを通して同じ履歴に戻る
    let bytes = micro.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), micro);
}

#[test]
fn micro_depletion_needs_micro_tables_of_the_same_nuclides_and_its_chain() {
    // 設計 B のテーブルが別の核種を追う
    assert!(micro_run(MIXED, [&TRACKED, &["U-235", "U-237"]]).is_err());
    // ミクロ燃焼の微視的断面積のないテーブル
    let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let case = micro_case(CASE);
    let depletion = case.depletion.as_ref().unwrap();
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        depletion,
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap_err();
    assert!(error.to_string().contains("micro"), "{error}");
    // ミクロ燃焼のケースにチェーンがない。ミクロ燃焼しないケースにチェーンを渡す
    let micro_tables: Vec<Table> = tables
        .into_iter()
        .map(|t| with_micro(t, &TRACKED))
        .collect();
    assert!(
        deplete(
            &layout,
            &micro_tables,
            &case.states,
            depletion,
            None,
            None,
            &tight(),
            |_| {}
        )
        .is_err()
    );
    let (_, _, plain) = core(CASE, UNIFORM, Some(400.0));
    assert!(
        deplete(
            &layout,
            &micro_tables,
            &plain.states,
            plain.depletion.as_ref().unwrap(),
            Some(&small_chain()),
            None,
            &tight(),
            |_| {}
        )
        .is_err()
    );
}

#[test]
fn micro_corrector_rates_take_the_cross_sections_at_the_midpoint_burnup() {
    // U-235 の捕獲 / 核分裂の比を燃焼度で 0.2 + 0.01 u（燃焼度の点 0、10、20 の間は線形に内挿される）。予測子は初めの燃焼度の比、修正子は
    // 中点の燃焼度の比で消滅する
    let ratio = |u: f64| CAPTURE_TO_FISSION + 0.01 * u;
    let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let case = micro_case(CASE);
    let tables: Vec<Table> = tables
        .into_iter()
        .map(|t| {
            let mut t = with_micro(t, &TRACKED);
            let points = t.points.clone();
            for (p, state) in points.iter().zip(&mut t.micro.as_mut().unwrap().states) {
                let f = [1.0, 50.0];
                state[0] = micro_nuclide(f, f.map(|x| ratio(BURNUP[p[0]]) * x), Some(0.95));
            }
            t.validate().unwrap();
            t
        })
        .collect();
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    let power = 40.0 * 400.0 / (20.0 * 20.0);
    let loss = |n: f64, u: f64| (1.0 + ratio(u)) * power * 1.0e-24 / (n * Q * EV);
    let mut u235 = U235;
    let (mut efpd, mut burnup) = (0.0, 0.0);
    for point in &history.points[1..] {
        let h = (point.efpd - efpd) * 86400.0;
        let middle = burnup + 0.5 * 40.0 * (point.efpd - efpd) / 1000.0;
        let half = u235 * (-loss(u235, burnup) * 0.5 * h).exp();
        u235 *= (-loss(half, middle) * h).exp();
        (efpd, burnup) = (point.efpd, point.burnup);
        for node in point.densities.as_ref().unwrap().iter().flatten() {
            assert!(
                (node[0] - u235).abs() <= 1e-10 * u235,
                "{} {node:?} vs {u235}",
                point.efpd
            );
        }
    }
}

/// 一様な炉心のテーブル（微視的断面積つき）と、その追う核種のチェーン
fn micro_layout() -> (Layout, Vec<Table>, Case, MicroChain) {
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let tables: Vec<Table> = tables
        .into_iter()
        .map(|t| with_micro(t, &TRACKED))
        .collect();
    let micro = tables[0].micro.as_ref().unwrap();
    let chain = MicroChain::new(&small_chain(), &micro.nuclides, &micro.states[0]).unwrap();
    (layout, tables, case, chain)
}

#[test]
fn corrected_constants_are_checked_and_compositions_cover_every_site() {
    let (layout, tables, case, chain) = micro_layout();
    let sites = layout.node_material.len();
    let burnups = vec![0.0; sites];
    let state = case.states[0].at(1200.0);
    // テーブルは履歴の残りを持たないので、ノードの履歴の変数の値は定数を変えない
    let histories = vec![Some(vec![0.745]); sites];
    let solve = |densities: &[Option<Vec<f64>>]| {
        crate::case::solve_compositions(
            &layout,
            &tables,
            &state,
            &burnups,
            Some(&Compositions {
                chain: &chain,
                densities,
                histories: &histories,
            }),
            &tight(),
        )
    };
    // テーブルの数密度なら補正は 0 で解ける
    let nominal = vec![Some(chain.initial(&[U235, 0.0])); sites];
    assert!(solve(&nominal).is_ok());
    // 補正で吸収が負になる組成は、解く前にエラー
    let negative = vec![Some(chain.initial(&[-1.0, 0.0])); sites];
    let error = solve(&negative).unwrap_err().to_string();
    assert!(error.contains("micro"), "{error}");
    // 位置の数と違う組成
    assert!(solve(&nominal[1..]).is_err());
}

#[test]
fn histories_hold_densities_exactly_when_they_name_nuclides() {
    let history = micro_run(UNIFORM, [&TRACKED, &TRACKED]).unwrap();
    let mut missing = history.clone();
    missing.points[1].densities = None;
    assert!(CoreHistory::from_msgpack(&missing.to_msgpack().unwrap()).is_err());
    let mut short = history.clone();
    short.points[1].densities.as_mut().unwrap()[0][0].pop();
    assert!(CoreHistory::from_msgpack(&short.to_msgpack().unwrap()).is_err());
    // 集合体と軸方向のノードの数、有限の値
    let mut rows = history.clone();
    rows.points[1].densities.as_mut().unwrap().pop();
    assert!(CoreHistory::from_msgpack(&rows.to_msgpack().unwrap()).is_err());
    let mut layers = history.clone();
    layers.points[1].densities.as_mut().unwrap()[0].pop();
    assert!(CoreHistory::from_msgpack(&layers.to_msgpack().unwrap()).is_err());
    let mut nan = history.clone();
    nan.points[1].densities.as_mut().unwrap()[0][0][0] = f64::NAN;
    assert!(CoreHistory::from_msgpack(&nan.to_msgpack().unwrap()).is_err());
    let mut unnamed = history.clone();
    unnamed.nuclides = None;
    assert!(CoreHistory::from_msgpack(&unnamed.to_msgpack().unwrap()).is_err());
    // Xe の扱い（#0200）は、ミクロ燃焼の履歴だけが持つ
    assert_eq!(
        CoreHistory::from_msgpack(&history.to_msgpack().unwrap())
            .unwrap()
            .xenon,
        Some(Xenon::Tracked)
    );
    let mut no_xenon = history.clone();
    no_xenon.xenon = None;
    assert!(CoreHistory::from_msgpack(&no_xenon.to_msgpack().unwrap()).is_err());
    let plain = run(CASE, UNIFORM).unwrap();
    assert_eq!(plain.xenon, None);
    let mut extra = plain;
    extra.xenon = Some(Xenon::Equilibrium);
    assert!(CoreHistory::from_msgpack(&extra.to_msgpack().unwrap()).is_err());
}

#[test]
fn fissile_tables_off_the_depleted_nodes_are_errors() {
    // 燃焼度で引かない使い方（反射体の位置など）の核分裂のあるテーブルは、補正も燃焼もしないので止める
    let (mut layout, tables, case, _) = micro_layout();
    let m = layout.node_material[0].unwrap();
    let t = layout.materials[m].table;
    layout.tables[t].node = None;
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        micro_case(CASE).depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("not a depleted node"), "{error}");
}

// 部分挿入のノードのミクロ燃焼（ADR-0217）

/// 集合体 1 体の上から制御棒（バンク D、AIC だけ）を入れる炉心の形。先端は 0.001 × ステップ (cm)
const ROD_PLANT: &str = r#"
pitch = 20.0
height = 20.0
map = ["*"]

[rods]
step = 0.001
bottom = 0.0
steps = 20000
top = 20.0
map = ["D"]

[rods.bank.D]
type = "rcca"

[rods.type.rcca]
segments = [{ value = "aic", length = 20.0 }]
"#;

const ROD_INVENTORY: &str = r#"
[design.A]
bounds = [0.0, 20.0]
layers = ["a.msgpack"]

[assembly.A1]
design = "A"
"#;

const ROD_CYCLE: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["A1"]
"#;

/// 軸方向に 2 つのノード（0〜10、10〜20 cm）。上のノードに先端が来る
const ROD_CASE: &str = r#"
cycle = "cycle.toml"
axial-mesh = [0.0, 10.0, 20.0]
radial-nodes = 1
boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], z = ["reflective", "reflective"] }

[state.hfp]
boron = 1200
b10 = 0.199
moderator = 565
fuel = 565
banks = { D = 20000 }

[depletion]
specific-power = 40.0
steps = [{ efpd = 50.0 }, { efpd = 125.0 }]
history = "history.msgpack"
"#;

/// `with_micro` のテーブルに制御棒の軸（`out`、`aic`）を足す。`aic` は熱群の吸収を足し、U-235 の微視的断面積の熱群を小さくする
fn with_rod_axis(table: Table) -> Table {
    let micro = table.micro.clone().unwrap();
    let mut states = Vec::new();
    let mut micro_states = Vec::new();
    for (state, nuclides) in table.states.iter().zip(&micro.states) {
        states.push(state.clone());
        micro_states.push(nuclides.clone());
        let mut rodded = state.clone();
        rodded.absorption[1] += 0.03;
        rodded.total[1] += 0.03;
        states.push(rodded);
        let mut rodded = nuclides.clone();
        let u235 = &mut rodded[0];
        // 熱群の核分裂は 0.7 倍、捕獲は 1.5 倍（捕獲 / 核分裂の比が変わり、U-236 の生まれ方が微視的断面積の混ぜ方で変わる）
        let (fission, capture) = (u235.fission[1], u235.reactions[0].channels[0][1]);
        u235.fission[1] = 0.7 * fission;
        u235.nu_fission[1] *= 0.7;
        u235.reactions[0].channels[0][1] = 1.5 * capture;
        u235.total[1] += u235.fission[1] - fission + u235.reactions[0].channels[0][1] - capture;
        u235.thermal_fission = Some(0.9);
        micro_states.push(rodded);
    }
    let mut out = Table {
        points: Vec::new(),
        axes: table.axes.clone(),
        split: table.split,
        states,
        heavy_metal: table.heavy_metal,
        micro: Some(Micro {
            states: micro_states,
            ..micro
        }),
        histories: Vec::new(),
    };
    out.axes.push(Axis::Rod {
        values: vec!["out".into(), "aic".into()],
    });
    out.points = product_points(&out.axes);
    out.validate().unwrap();
    out
}

/// バンク D をステップ `step` に置いて、集合体 1 体の 2 つのノードをミクロ燃焼させる
fn rodded_micro_run(step: u32) -> CoreHistory {
    rodded_micro_run_with(ROD_CASE, step)
}

/// `rodded_micro_run` をケース `case` で
fn rodded_micro_run_with(case: &str, step: u32) -> CoreHistory {
    let base = Path::new("dir");
    let mut case = micro_case(case);
    let cycle = Cycle::parse(ROD_CYCLE, base).unwrap();
    let plant = Plant::parse(ROD_PLANT, base).unwrap();
    let inventory = Inventory::parse(ROD_INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = vec![with_rod_axis(with_micro(table(1.0, Some(400.0)), &TRACKED))];
    case.states[0].banks.insert("D".into(), step);
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap()
}

/// 上のノードが、燃焼の終わりまでに消えた U-235 の原子あたりに作った U-236 の原子（捕獲の割合）。ノードの出力によらず、
/// 微視的断面積（とノードの中性子束のスペクトル）で決まる
fn capture_share(history: &CoreHistory) -> f64 {
    let n = &history.points.last().unwrap().densities.as_ref().unwrap()[0][1];
    n[1] / (U235 - n[0])
}

/// 最後の点の k (pcm)
fn last_k(history: &CoreHistory) -> f64 {
    1.0e5 * history.points.last().unwrap().k
}

#[test]
fn partial_node_micro_depletion_approaches_the_uniform_nodes_at_the_insertion_limits() {
    // 上のノードの挿入の割合が 0 と 1 に近づくと、部分挿入のノードのミクロ燃焼は、そのノードを棒のないノードと棒の入った
    // ノードにした燃焼に、割合の 1 次で近づく（区間の重みが 0 に近づく。ADR-0127 の「境界がノードの面に近づくと」）
    let out = rodded_micro_run(20000);
    let rodded = rodded_micro_run(10000);
    // 棒の効きと、捕獲の割合の違いは区別できる大きさ
    assert!(last_k(&out) - last_k(&rodded) > 100.0);
    assert!(capture_share(&rodded) > 1.2 * capture_share(&out));
    for (limit, near, nearer) in [(&out, 19990, 19999), (&rodded, 10010, 10001)] {
        let (near, nearer) = (rodded_micro_run(near), rodded_micro_run(nearer));
        let share = |h: &CoreHistory| capture_share(h) / capture_share(limit) - 1.0;
        let k = |h: &CoreHistory| last_k(h) - last_k(limit);
        let shares = format!("{} at 1e-3, {} at 1e-4", share(&near), share(&nearer));
        let ks = format!("{} pcm at 1e-3, {} pcm at 1e-4", k(&near), k(&nearer));
        // 割合を 1/10 にすると差もおよそ 1/10 になる（k は、区間の定数の z の面の不連続因子が 1 からずれる分（ADR-0127 の
        // `identical_intervals_reproduce_the_uniform_node`）を含むので、1/5 より小さいことだけを求める）
        assert!(share(&nearer).abs() < 1.0e-4, "{shares}");
        assert!(
            (5.0..20.0).contains(&(share(&near) / share(&nearer))),
            "{shares}"
        );
        assert!(k(&nearer).abs() < k(&near).abs() / 5.0, "{ks}");
    }
}

#[test]
fn partial_node_micro_rates_weigh_the_rodded_segment_by_its_flux() {
    // 棒の入った区間は中性子束が低いので、中性子束の重み（既定）は体積の重み（比べるためだけの案）より棒の入った区間の
    // 微視的断面積を小さく数え、上のノードの捕獲の割合は棒のないノードの値に近い
    let out = capture_share(&rodded_micro_run(20000));
    let flux = capture_share(&rodded_micro_run(15000));
    let volume = capture_share(&rodded_micro_run_with(
        &format!("{ROD_CASE}\n[partial-rods]\nweighting = \"volume\"\n"),
        15000,
    ));
    assert!(out < flux && flux < volume, "{out} {flux} {volume}");
}

#[test]
fn partial_node_micro_depletion_is_monotonic_in_the_insertion() {
    // 先端を下げる（挿入の割合を増やす）ほど、最後の点の k は下がる
    let ks: Vec<f64> = [20000, 17500, 15000, 12500, 10000]
        .into_iter()
        .map(|s| rodded_micro_run(s).points.last().unwrap().k)
        .collect();
    assert!(ks.windows(2).all(|w| w[1] < w[0]), "{ks:?}");
}

// サイクルのステップ（ADR-0218）: 出力の割合、崩壊だけのステップ、ステップごとの状態

use xlc_depletion::chain::{Decay, DecayEnergy, DecayMode};

/// I-135 と Xe-135 の半減期 (s)
const I135_HALF_LIFE: f64 = 6.57 * 3600.0;
const XE135_HALF_LIFE: f64 = 9.14 * 3600.0;
/// U-235 の核分裂あたりの I-135 と Xe-135 の独立収率
const YIELD_I135: f64 = 0.06;
const YIELD_XE135: f64 = 0.003;
/// Xe-135 の捕獲 / U-235 の核分裂の微視的断面積の比（群ごとに同じなので、Xe-135 の 1 原子あたりの捕獲率は U-235 の核分裂率の
/// この倍で、中性子束の形によらない）
const XE_TO_FISSION: f64 = 1.0e4;

/// 半減期 `half_life` (s) の β- 崩壊で `product` になる
fn beta_to(half_life: f64, product: usize) -> Decay {
    Decay {
        half_life,
        energy: DecayEnergy {
            light: 0.0,
            electromagnetic: 0.0,
            heavy: 0.0,
        },
        modes: vec![DecayMode {
            mode: "beta-".into(),
            ratio: 1.0,
            products: vec![Product {
                nuclide: product,
                count: 1,
            }],
        }],
    }
}

/// `small_chain` の核種に、I-135（β- で Xe-135）、Xe-135（捕獲で Xe-136、β- で Cs-135）、Xe-136、Cs-135 を足し、U-235 の核分裂の
/// 収率に I-135 と Xe-135 を足したチェーン（核種は (Z, A) の順）
fn xenon_chain() -> Chain {
    let mut nuclides = vec![
        chain_nuclide("I-135", 53, 135),
        chain_nuclide("Xe-135", 54, 135),
        chain_nuclide("Xe-136", 54, 136),
        chain_nuclide("Cs-135", 55, 135),
        chain_nuclide("Ba-140", 56, 140),
        chain_nuclide("U-235", 92, 235),
        chain_nuclide("U-236", 92, 236),
        chain_nuclide("U-237", 92, 237),
    ];
    nuclides[0].decay = Some(beta_to(I135_HALF_LIFE, 1));
    nuclides[1].decay = Some(beta_to(XE135_HALF_LIFE, 3));
    nuclides[1].reactions = vec![capture_to(2)];
    nuclides[5].reactions = vec![capture_to(6)];
    nuclides[6].reactions = vec![capture_to(7)];
    let yields = vec![
        Yield {
            nuclide: 0,
            value: YIELD_I135,
        },
        Yield {
            nuclide: 1,
            value: YIELD_XE135,
        },
        Yield {
            nuclide: 4,
            value: 0.06,
        },
    ];
    nuclides[5].fission = Some(Fission {
        q: Q,
        yield_source: "U-235".into(),
        energies: vec![0.0253, 5.0e5],
        yields: vec![yields.clone(), yields],
    });
    let chain = Chain { nuclides };
    chain.check().unwrap();
    chain
}

/// U-235 の核分裂の微視的断面積（`with_micro` と同じ）
const U235_FISSION: [f64; 2] = [1.0, 50.0];

/// テーブルに、U-235、U-236（`with_micro` と同じ）、I-135（反応なし）、Xe-135（捕獲は U-235 の核分裂の `XE_TO_FISSION` 倍）の
/// 状態によらない微視的断面積と、燃焼度の点ごとの数密度（I-135 と Xe-135 は 0）を足す
fn with_xenon(table: Table) -> Table {
    let mut table = with_micro(table, &TRACKED);
    let micro = table.micro.as_mut().unwrap();
    micro.nuclides.extend(["I-135".into(), "Xe-135".into()]);
    for densities in &mut micro.densities {
        densities.extend([0.0, 0.0]);
    }
    let iodine = MicroNuclide {
        total: vec![4.0, 10.0],
        nu_fission: vec![0.0; 2],
        fission: vec![0.0; 2],
        scatter: vec![vec![3.9, 0.1], vec![0.0, 10.0]],
        reactions: Vec::new(),
        thermal_fission: None,
        derivative: None,
    };
    let xenon = micro_nuclide([0.0, 0.0], U235_FISSION.map(|f| XE_TO_FISSION * f), None);
    for nuclides in &mut micro.states {
        nuclides.extend([iodine.clone(), xenon.clone()]);
    }
    table.validate().unwrap();
    table
}

/// 追う核種 U-235、U-236、I-135、Xe-135 の数密度 `n` を、U-235 の 1 原子あたりの核分裂率 `r` (1/s) で時間 `t` (s) 進めた解析解
/// （反応率が一定の 1 つのステップ。`r` = 0 は崩壊だけ）
fn xenon_bateman(n: [f64; 4], r: f64, t: f64) -> [f64; 4] {
    let [u235, u236, i0, x0] = n;
    let (li, lx) = (
        std::f64::consts::LN_2 / I135_HALF_LIFE,
        std::f64::consts::LN_2 / XE135_HALF_LIFE,
    );
    // U-235 の消滅の率、Xe-135 の消滅の率
    let a = (1.0 + CAPTURE_TO_FISSION) * r;
    let b = lx + XE_TO_FISSION * r;
    let u = u235 * (-a * t).exp();
    let iodine = i0 * (-li * t).exp()
        + YIELD_I135 * r * u235 * ((-a * t).exp() - (-li * t).exp()) / (li - a);
    // Xe-135 の源は e^(−a t) と e^(−λ_I t) の項の和: s_a e^(−a t) + s_i e^(−λ_I t)
    let s_a = YIELD_XE135 * r * u235 + li * YIELD_I135 * r * u235 / (li - a);
    let s_i = li * i0 - li * YIELD_I135 * r * u235 / (li - a);
    let xenon = x0 * (-b * t).exp()
        + s_a * ((-a * t).exp() - (-b * t).exp()) / (b - a)
        + s_i * ((-li * t).exp() - (-b * t).exp()) / (b - li);
    [
        u,
        u236 + CAPTURE_TO_FISSION / (1.0 + CAPTURE_TO_FISSION) * (u235 - u),
        iodine,
        xenon,
    ]
}

/// 一様な炉心（`UNIFORM`）を、Xe のチェーンとステップ `steps` でミクロ燃焼させる
fn xenon_run(steps: &str) -> Result<CoreHistory, CaseError> {
    xenon_run_in(UNIFORM, "tracked", steps)
}

/// `xenon_run` と同じく、装荷 `cycle` と Xe の扱い `xenon`（`micro` の `xenon`）で
fn xenon_run_in(cycle: &str, xenon: &str, steps: &str) -> Result<CoreHistory, CaseError> {
    let case = CASE
        .replace(
            "[depletion]",
            &format!("[depletion]\nmicro = {{ chain = \"chain.msgpack\", xenon = \"{xenon}\" }}"),
        )
        .replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            steps,
        );
    let (layout, tables, case) = core(&case, cycle, Some(400.0));
    let tables: Vec<Table> = tables.into_iter().map(with_xenon).collect();
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&xenon_chain()),
        None,
        &tight(),
        |_| {},
    )
}

/// 一様な炉心のノードの定格の出力密度 (W/cm³): 40 W/g × 400 g/cm / (20 cm × 20 cm)
const UNIFORM_POWER: f64 = 40.0 * 400.0 / (20.0 * 20.0);

/// ノードの組成 `n` の CE/CM の 1 ステップ（出力の割合 `power`、時間 `t` (s)）の解析の参照: 初めの組成で規格化した率で中点の
/// 組成へ、中点の組成で規格化した率で初めの組成から終わりへ（ADR-0216 の案 S1。核分裂するのは U-235 だけなので、1 原子あたりの
/// 核分裂率は 出力密度 / (N₂₃₅ Q)）
fn xenon_cecm(n: [f64; 4], power: f64, t: f64) -> [f64; 4] {
    let rate = |n: [f64; 4]| power * UNIFORM_POWER * 1.0e-24 / (n[0] * Q * EV);
    let half = xenon_bateman(n, rate(n), 0.5 * t);
    xenon_bateman(n, rate(half), t)
}

/// 点のノードの数密度がすべて `expected` に、`scale` に対して相対 `tolerance` で近い
fn assert_densities(point: &CorePoint, expected: [f64; 4], tolerance: f64) {
    for node in point.densities.as_ref().unwrap().iter().flatten() {
        for (i, (&got, &want)) in node.iter().zip(&expected).enumerate() {
            // Xe と I は、U-235 の数密度の大きさの丸め（CRAM の誤差）を受ける。それぞれの大きさか、その 1e-6 倍の U-235 の大きさ
            let scale = want.abs().max(1.0e-6 * expected[0]);
            assert!(
                (got - want).abs() <= tolerance * scale,
                "{} days, nuclide {i}: {got} vs {want}",
                point.days
            );
        }
    }
}

#[test]
fn power_fractions_stretch_the_step_time_and_keep_the_burnup() {
    // 出力の割合 0.5 のステップは、同じ EFPD で燃焼度の増分が同じで、時間が 2 倍、ノードの出力密度が半分（ADR-0218）
    let half = xenon_run("[{ efpd = 50.0, power = 0.5 }, { efpd = 125.0 }]").unwrap();
    let full = xenon_run("[{ efpd = 50.0 }, { efpd = 125.0 }]").unwrap();
    let days: Vec<f64> = half.points.iter().map(|p| p.days).collect();
    assert_eq!(days, [0.0, 100.0, 175.0]);
    let powers: Vec<Option<f64>> = half.points.iter().map(|p| p.power).collect();
    assert_eq!(powers, [None, Some(0.5), Some(1.0)]);
    let mut n = [U235, 0.0, 0.0, 0.0];
    let mut efpd = 0.0;
    for (h, f) in half.points.iter().zip(&full.points).skip(1) {
        assert_eq!(h.efpd, f.efpd);
        // 燃焼度は EFPD だけで決まる
        let burnup = 40.0 * h.efpd * 1.0e-3;
        assert!(
            (h.burnup - burnup).abs() <= 1.0e-12 * burnup,
            "{}",
            h.burnup
        );
        assert!((h.burnup - f.burnup).abs() <= 1.0e-12 * burnup);
        let power = h.power.unwrap();
        n = xenon_cecm(n, power, (h.efpd - efpd) / power * 86400.0);
        efpd = h.efpd;
        assert_densities(h, n, 1.0e-9);
    }
    // U-235 は燃焼度で決まるので同じで、出力の水準の低いステップの終わりの Xe-135 は少ない
    let node = |h: &CoreHistory, k: usize| h.points[k].densities.as_ref().unwrap()[0][0].clone();
    assert!((node(&half, 1)[0] / node(&full, 1)[0] - 1.0).abs() < 1.0e-12);
    assert!(node(&half, 1)[3] < 0.9 * node(&full, 1)[3]);
}

#[test]
fn decay_steps_keep_the_burnup_and_decay_iodine_and_xenon() {
    // 崩壊だけのステップ（ADR-0218）は燃焼度と EFPD を進めず、I-135 と Xe-135 を崩壊で進める（反応率 0 の解析解）
    let history = xenon_run(
        "[{ efpd = 50.0 }, { decay-days = 0.25 }, { decay-days = 18.0 }, { efpd = 60.0 }]",
    )
    .unwrap();
    let efpd: Vec<f64> = history.points.iter().map(|p| p.efpd).collect();
    assert_eq!(efpd, [0.0, 50.0, 50.0, 50.0, 60.0]);
    let days: Vec<f64> = history.points.iter().map(|p| p.days).collect();
    assert_eq!(days, [0.0, 50.0, 50.25, 68.25, 78.25]);
    let powers: Vec<Option<f64>> = history.points.iter().map(|p| p.power).collect();
    assert_eq!(powers, [None, Some(1.0), Some(0.0), Some(0.0), Some(1.0)]);
    let burnups: Vec<f64> = history.points.iter().map(|p| p.burnup).collect();
    assert_eq!(burnups[1], burnups[2]);
    assert_eq!(burnups[2], burnups[3]);
    let mut n = xenon_cecm([U235, 0.0, 0.0, 0.0], 1.0, 50.0 * 86400.0);
    assert_densities(&history.points[1], n, 1.0e-9);
    // 停止の 6 時間後は I-135 の崩壊で Xe-135 が増え、18 日の後は I-135 も Xe-135 も崩壊で消える
    let node = |k: usize| history.points[k].densities.as_ref().unwrap()[0][0].clone();
    assert!(node(2)[3] > node(1)[3], "{:?} {:?}", node(1), node(2));
    for (k, d) in [(2, 0.25), (3, 18.0)] {
        n = xenon_bateman(n, 0.0, d * 86400.0);
        assert_densities(&history.points[k], n, 1.0e-10);
    }
    assert!(node(3)[3] < 1.0e-12 * node(1)[3]);
    // Xe-135 が抜けると k は上がる（Xe-135 の公称の数密度は 0 で、線形の補正が吸収を足す）
    assert!(history.points[3].k > history.points[1].k);
    // 出力のステップに戻ると、崩壊した組成から CE/CM で進む
    n = xenon_cecm(n, 1.0, 10.0 * 86400.0);
    assert_densities(&history.points[4], n, 1.0e-9);
}

#[test]
fn decay_steps_need_micro_depletion() {
    // 組成を持たない燃焼に、崩壊だけのステップで変わるものはない（入力の検査を通らない形を直接渡しても止まる）
    let case = CASE
        .replace(
            "[depletion]",
            "[depletion]\nmicro = { chain = \"chain.msgpack\", xenon = \"tracked\" }",
        )
        .replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            "[{ efpd = 50.0 }, { decay-days = 1.0 }]",
        );
    let (layout, tables, case) = core(&case, UNIFORM, Some(400.0));
    let mut depletion = case.depletion.clone().unwrap();
    depletion.micro = None;
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        &depletion,
        None,
        None,
        &tight(),
        |_| {},
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("decay step"), "{error}");
}

/// `ROD_CASE` の状態を、バンク D の位置の違う 2 つ（out は抜いた位置、in は上のノードの半分まで）にしたもの
fn two_rod_states(steps: &str) -> String {
    ROD_CASE
        .replace(
            "[state.hfp]\nboron = 1200\nb10 = 0.199\nmoderator = 565\nfuel = 565\nbanks = { D = 20000 }\n",
            "[state.out]\nboron = 1200\nb10 = 0.199\nmoderator = 565\nfuel = 565\nbanks = { D = 20000 }\n\n\
             [state.in]\nboron = 1200\nb10 = 0.199\nmoderator = 565\nfuel = 565\nbanks = { D = 15000 }\n",
        )
        .replace(
            "[depletion]",
            "[depletion]\nmicro = { chain = \"chain.msgpack\", xenon = \"tracked\" }",
        )
        .replace("[{ efpd = 50.0 }, { efpd = 125.0 }]", steps)
}

/// `rodded_micro_run_with` と同じ集合体 1 体を、ケース `case` の状態とステップで燃焼させる。並び、テーブル、ケースも返す
fn rod_states_run(case: &str) -> (CoreHistory, Layout, Vec<Table>, Case) {
    let base = Path::new("dir");
    let case = Case::parse(case, base).unwrap();
    let cycle = Cycle::parse(ROD_CYCLE, base).unwrap();
    let plant = Plant::parse(ROD_PLANT, base).unwrap();
    let inventory = Inventory::parse(ROD_INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = vec![with_rod_axis(with_micro(table(1.0, Some(400.0)), &TRACKED))];
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    (history, layout, tables, case)
}

#[test]
fn steps_solve_in_the_states_they_name() {
    // バンク D を途中で入れて抜く。点は、そこから始まるステップの状態で解く（最後の点は最後のステップの状態。ADR-0218 の案 P2）
    let steps = "[{ efpd = 50.0, state = \"out\" }, { efpd = 125.0, state = \"in\" }, { efpd = 200.0, state = \"out\" }]";
    let (history, layout, tables, case) = rod_states_run(&two_rod_states(steps));
    let states: Vec<&str> = history.points.iter().map(|p| p.state.as_str()).collect();
    assert_eq!(states, ["out", "in", "out", "out"]);
    assert_eq!(history.states.keys().collect::<Vec<_>>(), ["in", "out"]);
    assert_eq!(history.states["in"].banks["D"], 15000);
    // 点の k は、その状態の並びで、点のノードの燃焼度と組成を解いた値
    let micro = tables[0].micro.as_ref().unwrap();
    let chain = MicroChain::new(&small_chain(), &micro.nuclides, &micro.states[0]).unwrap();
    for point in &history.points {
        let state = case.states.iter().find(|s| s.name == point.state).unwrap();
        let banked = layout.with_banks(&state.banks).unwrap();
        let burnups: Vec<f64> = point.burnups[0].clone();
        let densities: Vec<Option<Vec<f64>>> = point.densities.as_ref().unwrap()[0]
            .iter()
            .map(|n| Some(chain.initial(n)))
            .collect();
        let histories: Vec<Option<Vec<f64>>> = point.histories.as_ref().unwrap()[0]
            .iter()
            .map(|h| Some(h.clone()))
            .collect();
        let solved = crate::case::solve_compositions(
            &banked,
            &tables,
            &state.at(1200.0),
            &burnups,
            Some(&Compositions {
                chain: &chain,
                densities: &densities,
                histories: &histories,
            }),
            &tight(),
        )
        .unwrap();
        assert_eq!(point.k, solved.solution.k, "{} EFPD", point.efpd);
    }
    // 棒を入れたステップは、すべて抜いた燃焼より上のノードの燃焼度を下げ、入れた状態で解いた点の k は低い
    let (out, ..) = rod_states_run(&two_rod_states(
        "[{ efpd = 50.0, state = \"out\" }, { efpd = 200.0 }, { efpd = 201.0, state = \"in\" }]",
    ));
    assert_eq!(history.points[0].k, out.points[0].k);
    assert!(1.0e5 * (out.points[1].k - history.points[1].k) > 100.0);
    let upper = |h: &CoreHistory| h.points.last().unwrap().burnups[0][1];
    assert!(upper(&history) < upper(&out));
    // 状態を省いたステップは、前のステップの状態
    let (omitted, ..) = rod_states_run(&two_rod_states(
        "[{ efpd = 50.0, state = \"out\" }, { efpd = 125.0, state = \"in\" }, { efpd = 200.0 }]",
    ));
    let (named, ..) = rod_states_run(&two_rod_states(
        "[{ efpd = 50.0, state = \"out\" }, { efpd = 125.0, state = \"in\" }, { efpd = 200.0, state = \"in\" }]",
    ));
    assert_eq!(omitted, named);
}

#[test]
fn history_points_name_their_state_and_step_power() {
    let history = xenon_run("[{ efpd = 50.0, power = 0.5 }, { decay-days = 2.0 }]").unwrap();
    let bytes = history.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), history);
    // 点の状態は履歴の状態にあること。燃焼点 0 だけが出力の割合を持たない
    let with_point = |k: usize, edit: &dyn Fn(&mut CorePoint)| {
        let mut h = history.clone();
        edit(&mut h.points[k]);
        CoreHistory::from_msgpack(&h.to_msgpack().unwrap())
    };
    assert!(with_point(1, &|p| p.state = "cold".into()).is_err());
    assert!(with_point(0, &|p| p.power = Some(1.0)).is_err());
    assert!(with_point(2, &|p| p.power = None).is_err());
    // 出力の割合は 0 以上の有限の値で、0 は EFPD の進まない崩壊だけのステップだけ。EFPD と日数は減らない
    assert!(with_point(1, &|p| p.power = Some(-1.0)).is_err());
    assert!(with_point(1, &|p| p.power = Some(f64::NAN)).is_err());
    assert!(with_point(1, &|p| p.power = Some(0.0)).is_err());
    assert!(with_point(2, &|p| p.power = Some(1.0)).is_err());
    assert!(with_point(2, &|p| p.days = 1.0).is_err());
    // 丸めで日数の増えない短いステップの履歴も、書いたものを読める
    assert!(with_point(2, &|p| p.days = 100.0).is_ok());
    assert!(with_point(2, &|p| p.efpd = 1.0).is_err());
}

#[test]
fn depletions_without_states_or_steps_are_errors() {
    // 入力の検査を通らない形を直接渡しても、添字で止まらずにエラーを返す
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let mut empty = case.depletion.clone().unwrap();
    empty.steps.clear();
    for (states, depletion) in [
        (&[][..], case.depletion.as_ref().unwrap()),
        (&[][..], &empty),
        (&case.states[..], &empty),
    ] {
        assert!(
            deplete(
                &layout,
                &tables,
                states,
                depletion,
                None,
                None,
                &tight(),
                |_| {}
            )
            .is_err()
        );
    }
}

/// 一様な炉心のノードの、U-235 の数密度 `u235` と出力の割合 `power` で釣り合う I-135 と Xe-135 の数密度（1 原子あたりの核分裂率は
/// 出力密度 / (N₂₃₅ Q)）
fn xenon_equilibrium(u235: f64, power: f64) -> [f64; 2] {
    let r = power * UNIFORM_POWER * 1.0e-24 / (u235 * Q * EV);
    let (li, lx) = (
        std::f64::consts::LN_2 / I135_HALF_LIFE,
        std::f64::consts::LN_2 / XE135_HALF_LIFE,
    );
    let iodine = YIELD_I135 * r * u235 / li;
    let xenon = (YIELD_XE135 + YIELD_I135) * r * u235 / (lx + XE_TO_FISSION * r);
    [iodine, xenon]
}

#[test]
fn equilibrium_xenon_balances_the_operating_points_and_decays_after_shutdowns() {
    // 平衡 Xe（#0200）: 崩壊だけのステップの前の出力のステップの終わりの点は I-135 と Xe-135 がそのステップの出力で釣り合い、燃焼点 0 と
    // 崩壊だけのステップの後の点は、組成を時間で進めた値のまま。ほかの核種はステップの CE/CM（Xe は核分裂しないので、U-235 は時間で追う燃焼と同じ）
    let steps = "[{ efpd = 50.0, power = 0.5 }, { decay-days = 0.25 }, { efpd = 60.0 }]";
    let history = xenon_run_in(UNIFORM, "equilibrium", steps).unwrap();
    let tracked = xenon_run_in(UNIFORM, "tracked", steps).unwrap();
    assert_eq!(history.xenon, Some(Xenon::Equilibrium));
    assert_eq!(tracked.xenon, Some(Xenon::Tracked));
    let node = |h: &CoreHistory, k: usize| h.points[k].densities.as_ref().unwrap()[0][0].clone();
    assert_densities(&history.points[0], [U235, 0.0, 0.0, 0.0], 0.0);
    for k in 1..4 {
        assert!(
            (node(&history, k)[0] / node(&tracked, k)[0] - 1.0).abs() <= 1.0e-12,
            "point {k}"
        );
    }
    let at = |k: usize, power: f64| {
        let n = node(&history, k);
        let [iodine, xenon] = xenon_equilibrium(n[0], power);
        [n[0], n[1], iodine, xenon]
    };
    let balanced = at(1, 0.5);
    assert_densities(&history.points[1], balanced, 1.0e-9);
    // 停止の 6 時間後は、平衡の組成から崩壊で進めた値
    assert_densities(
        &history.points[2],
        xenon_bateman(balanced, 0.0, 0.25 * 86400.0),
        1.0e-10,
    );
    assert_densities(&history.points[3], at(3, 1.0), 1.0e-9);
    // 点の固有値計算の回数は、平衡の固定点のすべての解の和（一様な炉心の反応率は Xe によらないので、平衡にした組成を 1 度解き直して
    // 止まる）
    let solves: Vec<usize> = history.points.iter().map(|p| p.solves).collect();
    assert_eq!(solves, [1, 2, 1, 2]);
    // 出力のステップが続く点は、次のステップの出力で釣り合う（点は次のステップの条件で解く。ADR-0218 の案 P2。ステップの初めの
    // 点、中点、終わりの点が同じ出力の平衡になる）
    let changed = xenon_run_in(
        UNIFORM,
        "equilibrium",
        "[{ efpd = 50.0, power = 0.5 }, { efpd = 60.0, power = 0.8 }]",
    )
    .unwrap();
    let n = node(&changed, 1);
    let [iodine, xenon] = xenon_equilibrium(n[0], 0.8);
    assert_densities(&changed.points[1], [n[0], n[1], iodine, xenon], 1.0e-9);
    let n = node(&changed, 2);
    let [iodine, xenon] = xenon_equilibrium(n[0], 0.8);
    assert_densities(&changed.points[2], [n[0], n[1], iodine, xenon], 1.0e-9);
}

#[test]
fn equilibrium_xenon_steps_keep_second_order() {
    // A と B の炉心は出力の割合が燃焼と Xe で動く。中点と点の解で Xe を平衡にしても、ステップを半分にした燃焼度の差の比は 2 次の 4 に近い
    let burn = |steps: usize| {
        let efpd: Vec<String> = (1..=steps)
            .map(|i| format!("{{ efpd = {:?} }}", 150.0 * i as f64 / steps as f64))
            .collect();
        let history =
            xenon_run_in(MIXED, "equilibrium", &format!("[{}]", efpd.join(", "))).unwrap();
        history.points.last().unwrap().assembly_burnups[0]
    };
    let (b2, b4, b8) = (burn(2), burn(4), burn(8));
    let ratio = (b2 - b4) / (b4 - b8);
    assert!((3.5..4.5).contains(&ratio), "{ratio} ({b2}, {b4}, {b8})");
}

/// `two_rod_states` のケースを、Xe のチェーンと平衡 Xe で燃焼させる（`rod_states_run` と同じ集合体 1 体）
fn xenon_rod_states_run(steps: &str) -> CoreHistory {
    let case = two_rod_states(steps).replace("xenon = \"tracked\"", "xenon = \"equilibrium\"");
    let base = Path::new("dir");
    let case = Case::parse(&case, base).unwrap();
    let cycle = Cycle::parse(ROD_CYCLE, base).unwrap();
    let plant = Plant::parse(ROD_PLANT, base).unwrap();
    let inventory = Inventory::parse(ROD_INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables = vec![with_rod_axis(with_xenon(table(1.0, Some(400.0))))];
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&xenon_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap()
}

#[test]
fn equilibrium_before_a_decay_step_balances_in_the_operating_state() {
    // 崩壊だけのステップの前の点は、崩壊だけのステップの状態で解く（ADR-0218 の案 P2）が、Xe は停止の前に運転していた状態の中性子束で
    // 釣り合う（停止の状態の中性子束で運転したことはない）
    let changed = xenon_rod_states_run(
        "[{ efpd = 50.0, state = \"out\" }, { decay-days = 1.0, state = \"in\" }, { efpd = 60.0 }]",
    );
    let kept = xenon_rod_states_run(
        "[{ efpd = 50.0, state = \"out\" }, { decay-days = 1.0 }, { efpd = 60.0, state = \"in\" }]",
    );
    assert_eq!(changed.points[1].state, "in");
    assert_eq!(kept.points[1].state, "out");
    let node = |h: &CoreHistory| h.points[1].densities.clone().unwrap();
    for (a, b) in node(&changed)
        .iter()
        .flatten()
        .zip(node(&kept).iter().flatten())
    {
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() <= 1.0e-9 * y.abs(), "{x} vs {y}");
        }
    }
    // 棒を抜いた状態では上下のノードの Xe-135 が同じ（棒の入った状態の中性子束では上のノードが少ない）
    let rows = node(&changed);
    assert!((rows[0][0][3] / rows[0][1][3] - 1.0).abs() <= 1.0e-9);
}

#[test]
fn equilibrium_xenon_without_xenon_and_iodine_stops_before_solving() {
    // 平衡にする Xe-135 と I-135 を解かないチェーンと追う核種は、最初の点を解く前に止める
    let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let mut case = micro_case(CASE);
    case.depletion
        .as_mut()
        .unwrap()
        .micro
        .as_mut()
        .unwrap()
        .xenon = Xenon::Equilibrium;
    let tables: Vec<Table> = tables
        .into_iter()
        .map(|t| with_micro(t, &TRACKED))
        .collect();
    let mut points = 0;
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| points += 1,
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("Xe-135"), "{error}");
    assert_eq!(points, 0);
}

// restart と装荷替え（ADR-0220）

/// テーブルに、U-235、U-236（`with_micro` と同じ）と Xe-135 を追う微視的断面積を足す。I-135 は追わないので、Xe-135 へ崩壊で流れ込む
/// 補助の核種になる
fn with_auxiliary_iodine(table: Table) -> Table {
    let mut table = with_micro(table, &TRACKED);
    let micro = table.micro.as_mut().unwrap();
    micro.nuclides.push("Xe-135".into());
    for densities in &mut micro.densities {
        densities.push(0.0);
    }
    let xenon = micro_nuclide([0.0, 0.0], U235_FISSION.map(|f| XE_TO_FISSION * f), None);
    for nuclides in &mut micro.states {
        nuclides.push(xenon.clone());
    }
    table.validate().unwrap();
    table
}

/// `cycle` の装荷に restart の行 `restart` を足したもの
fn restarted_cycle(cycle: &str, restart: &str) -> String {
    cycle.replace(
        "inventory = \"assemblies.toml\"\n",
        &format!("inventory = \"assemblies.toml\"\n{restart}\n"),
    )
}

/// 装荷 `cycle` を、Xe の扱い `xenon`（`None` はミクロ燃焼しない）とステップ `steps` で、前のサイクルの履歴 `restart` から燃焼させる。
/// ミクロ燃焼は `xenon_chain` と、I-135 を補助の核種にするテーブル
fn cycle_run(
    cycle: &str,
    xenon: Option<&str>,
    steps: &str,
    restart: Option<&CoreHistory>,
) -> Result<CoreHistory, CaseError> {
    let micro = xenon.map_or_else(String::new, |x| {
        format!("micro = {{ chain = \"chain.msgpack\", xenon = \"{x}\" }}\n")
    });
    let case = CASE
        .replace("[depletion]\n", &format!("[depletion]\n{micro}"))
        .replace(
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            steps,
        );
    let (layout, tables, case) = core(&case, cycle, Some(400.0));
    let (tables, chain): (Vec<Table>, _) = match xenon {
        Some(_) => (
            tables.into_iter().map(with_auxiliary_iodine).collect(),
            Some(xenon_chain()),
        ),
        None => (tables, None),
    };
    deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        chain.as_ref(),
        restart,
        &tight(),
        |_| {},
    )
}

/// 相対 `tolerance` で近い（`scale` は 0 に近い値の比べ方の大きさ）
fn near(a: f64, b: f64, tolerance: f64, scale: f64) -> bool {
    (a - b).abs() <= tolerance * b.abs().max(scale)
}

#[test]
fn restarts_continue_the_depletion_of_one_case() {
    // サイクルの途中で履歴に書いて restart した燃焼は、切らずに燃焼させた 1 つのケースの続きと、k、ノードの燃焼度、追う核種と補助の
    // 核種の数密度が丸めの範囲で一致する（ADR-0220）。A と B の炉心は出力の割合が燃焼で動く
    for xenon in ["tracked", "equilibrium"] {
        let full = cycle_run(
            MIXED,
            Some(xenon),
            "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
            None,
        )
        .unwrap();
        let first = cycle_run(
            MIXED,
            Some(xenon),
            "[{ efpd = 50.0 }, { efpd = 125.0 }]",
            None,
        )
        .unwrap();
        // 履歴はファイルを通して読む
        let first = CoreHistory::from_msgpack(&first.to_msgpack().unwrap()).unwrap();
        assert_eq!(first.auxiliary, Some(vec!["I-135".to_string()]));
        let second = cycle_run(
            &restarted_cycle(MIXED, "restart = { history = \"first.msgpack\" }"),
            Some(xenon),
            "[{ efpd = 125.0 }]",
            Some(&first),
        )
        .unwrap();
        // EFPD と日数は 0 から数え直し、燃焼度はノードごとに引き継ぐ
        let efpd: Vec<f64> = second.points.iter().map(|p| p.efpd).collect();
        assert_eq!(efpd, [0.0, 125.0], "{xenon}");
        assert_eq!(second.assemblies, full.assemblies);
        assert_eq!(second.heavy_metal, full.heavy_metal);
        for (got, want) in second.points.iter().zip(&full.points[2..]) {
            assert!(
                near(got.k, want.k, 1.0e-12, 1.0),
                "{xenon}: {} vs {}",
                got.k,
                want.k
            );
            assert!(near(got.burnup, want.burnup, 1.0e-12, 1.0), "{xenon}");
            for (g, w) in got
                .burnups
                .iter()
                .flatten()
                .zip(want.burnups.iter().flatten())
            {
                assert!(near(*g, *w, 1.0e-12, 1.0), "{xenon}: {g} vs {w}");
            }
            let (g, w) = (
                got.densities.as_ref().unwrap(),
                want.densities.as_ref().unwrap(),
            );
            for (g, w) in g
                .iter()
                .flatten()
                .flatten()
                .zip(w.iter().flatten().flatten())
            {
                assert!(near(*g, *w, 1.0e-9, 1.0e-20), "{xenon}: {g} vs {w}");
            }
        }
        // 最後の点の補助の核種（I-135）も続きと一致し、運転中の I-135 は 0 でない
        let (g, w) = (
            second.auxiliary_densities.as_ref().unwrap(),
            full.auxiliary_densities.as_ref().unwrap(),
        );
        for (g, w) in g
            .iter()
            .flatten()
            .flatten()
            .zip(w.iter().flatten().flatten())
        {
            assert!(
                *w > 0.0 && near(*g, *w, 1.0e-9, 1.0e-20),
                "{xenon}: {g} vs {w}"
            );
        }
    }
}

#[test]
fn restarts_without_the_auxiliary_nuclides_would_not_continue() {
    // 補助の核種を 0 に戻すと（今の初めの組成と同じ）、停止の後の Xe-135 が I-135 の崩壊の分だけ違う。restart が補助の核種を
    // 引き継ぐことの対照
    let steps = "[{ efpd = 50.0 }, { decay-days = 0.25 }]";
    let full = cycle_run(MIXED, Some("tracked"), "[{ efpd = 50.0 }]", None).unwrap();
    let decayed = cycle_run(MIXED, Some("tracked"), steps, None).unwrap();
    let mut zeroed = full.clone();
    for v in zeroed
        .auxiliary_densities
        .as_mut()
        .unwrap()
        .iter_mut()
        .flatten()
        .flatten()
    {
        *v = 0.0;
    }
    let cycle = restarted_cycle(MIXED, "restart = { history = \"first.msgpack\" }");
    let decay = "[{ decay-days = 0.25 }, { efpd = 1.0 }]";
    let kept = cycle_run(&cycle, Some("tracked"), decay, Some(&full)).unwrap();
    let lost = cycle_run(&cycle, Some("tracked"), decay, Some(&zeroed)).unwrap();
    let xenon = |h: &CoreHistory, k: usize| h.points[k].densities.as_ref().unwrap()[0][0][2];
    assert!(near(xenon(&kept, 1), xenon(&decayed, 2), 1.0e-9, 1.0e-20));
    assert!(xenon(&lost, 1) < 0.9 * xenon(&kept, 1));
}

/// 装荷替えの後の燃焼点 0 の、集合体 `id` のノードの燃焼度と数密度（追う核種）
fn assembly_state(history: &CoreHistory, k: usize, id: &str) -> (Vec<f64>, Vec<Vec<f64>>) {
    let a = history.assemblies.iter().position(|a| a == id).unwrap();
    let point = &history.points[k];
    (
        point.burnups[a].clone(),
        point.densities.as_ref().unwrap()[a].clone(),
    )
}

#[test]
fn reloads_carry_the_states_by_assembly_id() {
    // A1 と A2 を入れ替えて A2 を回し、B2 を取り出して新しい A3 を装荷する。引き継ぐ集合体は、履歴の最後の点の集合体の状態をその
    // まま持ち（回転はノードの状態を回さない）、新しい集合体は燃焼度 0 とテーブルの燃焼度 0 の組成（ADR-0220）
    let first = cycle_run(
        MIXED,
        Some("tracked"),
        "[{ efpd = 50.0 }, { efpd = 125.0 }]",
        None,
    )
    .unwrap();
    let reloaded = "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\n\
                    restart = { history = \"first.msgpack\", fresh = [\"A3\"] }\n\
                    loading = [\"A2@90  B1\", \"A1  A3\"]\n";
    let second = cycle_run(reloaded, Some("tracked"), "[{ efpd = 50.0 }]", Some(&first)).unwrap();
    assert_eq!(second.assemblies, ["A1", "A2", "A3", "B1"]);
    let last = first.points.len() - 1;
    for id in ["A1", "A2", "B1"] {
        assert_eq!(
            assembly_state(&second, 0, id),
            assembly_state(&first, last, id),
            "{id}"
        );
    }
    let (burnups, densities) = assembly_state(&second, 0, "A3");
    assert_eq!(burnups, [0.0, 0.0]);
    assert_eq!(densities, [[U235, 0.0, 0.0], [U235, 0.0, 0.0]]);
    // 炉心の平均の燃焼度は装荷した集合体の質量の重みで、ステップの増分は定格の出力密度 × EFPD
    let mass: f64 = second.heavy_metal.iter().flatten().sum();
    let start: f64 = second.points[0]
        .burnups
        .iter()
        .flatten()
        .zip(second.heavy_metal.iter().flatten())
        .map(|(b, m)| b * m)
        .sum::<f64>()
        / mass;
    assert!(near(second.points[0].burnup, start, 1.0e-15, 1.0));
    assert!(near(
        second.points[1].burnup,
        start + 40.0 * 50.0e-3,
        1.0e-12,
        1.0
    ));
    // 新しい集合体の補助の核種は 0、引き継いだ集合体は履歴の値から進む
    let a3 = second.assemblies.iter().position(|a| a == "A3").unwrap();
    assert!(second.auxiliary_densities.as_ref().unwrap()[a3][0][0] > 0.0);
}

#[test]
fn restart_mismatches_are_errors() {
    // 黙って合わせない（ADR-0220 の「確かめてエラーにするもの」）
    let steps = "[{ efpd = 50.0 }]";
    let first = cycle_run(MIXED, Some("tracked"), "[{ efpd = 50.0 }]", None).unwrap();
    let restart = |fresh: &str| {
        restarted_cycle(
            MIXED,
            &format!("restart = {{ history = \"first.msgpack\", fresh = [{fresh}] }}"),
        )
    };
    let error = |cycle: &str, xenon: Option<&str>, history: Option<&CoreHistory>| {
        cycle_run(cycle, xenon, steps, history)
            .unwrap_err()
            .to_string()
    };
    let same = restart("");
    // サイクルと燃焼の組み合わせ: restart のサイクルに履歴がない、restart のないサイクルに履歴がある
    assert!(error(&same, Some("tracked"), None).contains("restart"));
    assert!(error(MIXED, Some("tracked"), Some(&first)).contains("restart"));
    // 履歴にも fresh にもない ID、両方にある ID
    let reloaded = "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\n\
                    restart = { history = \"first.msgpack\" }\nloading = [\"A1  B1\", \"B2  A3\"]\n";
    assert!(error(reloaded, Some("tracked"), Some(&first)).contains("A3"));
    assert!(error(&restart("\"A1\""), Some("tracked"), Some(&first)).contains("A1"));
    // 軸方向のメッシュ、重金属の質量
    let mut h = first.clone();
    h.axial_mesh[1] = 10.000000000000002;
    assert!(error(&same, Some("tracked"), Some(&h)).contains("axial"));
    let mut h = first.clone();
    h.heavy_metal[0][0] *= 1.0 + f64::EPSILON;
    assert!(error(&same, Some("tracked"), Some(&h)).contains("heavy"));
    // ミクロ燃焼の有無、追う核種、補助の核種
    let plain = cycle_run(MIXED, None, "[{ efpd = 50.0 }]", None).unwrap();
    assert!(error(&same, Some("tracked"), Some(&plain)).contains("micro"));
    assert!(error(&same, None, Some(&first)).contains("micro"));
    let mut h = first.clone();
    h.nuclides = Some(vec!["U-235".into(), "U-236".into(), "Xe-136".into()]);
    assert!(error(&same, Some("tracked"), Some(&h)).contains("nuclides"));
    let mut h = first.clone();
    h.auxiliary = Some(vec!["Te-135".into()]);
    assert!(error(&same, Some("tracked"), Some(&h)).contains("auxiliary"));
    // 数密度を持つノードがミクロ燃焼するノードと違う、点がない
    let mut h = first.clone();
    h.points.last_mut().unwrap().densities.as_mut().unwrap()[0][0].clear();
    assert!(error(&same, Some("tracked"), Some(&h)).contains("densities"));
    let mut h = first.clone();
    h.points.clear();
    assert!(error(&same, Some("tracked"), Some(&h)).contains("point"));
    // 燃焼しない restart の燃焼は、重金属の質量と燃焼度だけを引き継ぐ
    let again = cycle_run(&same, None, "[{ efpd = 50.0 }]", Some(&plain)).unwrap();
    assert_eq!(again.points[0].burnups, plain.points[1].burnups);
}

#[test]
fn micro_histories_hold_the_auxiliary_densities_of_their_last_point() {
    // ミクロ燃焼の履歴は補助の核種の名前と最後の点の数密度を必ず持ち、ミクロ燃焼しない履歴は持たない（ADR-0220）
    let micro = cycle_run(MIXED, Some("tracked"), "[{ efpd = 50.0 }]", None).unwrap();
    let plain = cycle_run(MIXED, None, "[{ efpd = 50.0 }]", None).unwrap();
    assert_eq!(plain.auxiliary, None);
    assert_eq!(plain.auxiliary_densities, None);
    for history in [&micro, &plain] {
        let bytes = history.to_msgpack().unwrap();
        assert_eq!(&CoreHistory::from_msgpack(&bytes).unwrap(), history);
    }
    let rejected = |h: &CoreHistory| CoreHistory::from_msgpack(&h.to_msgpack().unwrap()).is_err();
    let mut h = micro.clone();
    h.auxiliary = None;
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.auxiliary_densities = None;
    assert!(rejected(&h));
    let mut h = plain.clone();
    h.auxiliary = Some(Vec::new());
    h.auxiliary_densities = Some(Vec::new());
    assert!(rejected(&h));
    // 形: [集合体][軸方向のノード][補助の核種]。核分裂のないノードだけが空で、有限の値
    let mut h = micro.clone();
    h.auxiliary_densities.as_mut().unwrap()[0].pop();
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.auxiliary_densities.as_mut().unwrap()[0][0].clear();
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.auxiliary_densities.as_mut().unwrap()[0][0][0] = f64::NAN;
    assert!(rejected(&h));
    // 点のない履歴は、補助の核種の形でなく点がないことで止める
    let mut h = micro.clone();
    h.points.clear();
    let error = CoreHistory::from_msgpack(&h.to_msgpack().unwrap())
        .unwrap_err()
        .to_string();
    assert!(error.contains("no point"), "{error}");
}

// 履歴の変数（ADR-0215）

#[test]
fn history_variables_average_the_step_values_by_burnup() {
    // 一定の値なら、どの燃焼度でもその値（ビットで同じ）
    let mut h = advance_history(&[0.743], 0.0, 0.1, &[0.743], 1.0);
    assert_eq!(h, [0.743]);
    for (start, end) in [(0.1, 0.25), (0.25, 13.7), (13.7, 60.0)] {
        h = advance_history(&h, start, end, &[0.743], 1.0);
        assert_eq!(h, [0.743]);
    }
    // 燃焼度 0 からは、ステップの値（燃焼度 0 の値は瞬時の値）
    assert_eq!(advance_history(&[0.743], 0.0, 2.0, &[0.661], 1.0), [0.661]);
    // 途中で値を変えると、燃焼度で重みを付けた平均 (1/B) ∫ x dB'
    let h = advance_history(&[0.743], 0.0, 20.0, &[0.743], 1.0);
    let h = advance_history(&h, 20.0, 30.0, &[0.661], 1.0);
    let h = advance_history(&h, 30.0, 60.0, &[0.661], 1.0);
    let expected = (0.743 * 20.0 + 0.661 * 40.0) / 60.0;
    assert!((h[0] - expected).abs() <= 1.0e-15, "{h:?} vs {expected}");
    // 燃焼度が進まなければ（崩壊だけのステップ）、そのまま
    assert_eq!(advance_history(&[0.7], 12.0, 12.0, &[0.661], 1.0), [0.7]);
}

#[test]
fn history_steps_weight_the_step_values_by_the_history_weight() {
    // PARCS の形 H(B₀ + ΔB) = (H(B₀) B₀ + ω x ΔB) / (B₀ + ω ΔB)（ADR-0215 の検討した案）
    for omega in [0.5, 1.25, 2.0, 2.5] {
        let h = advance_history(&[0.743], 20.0, 30.0, &[0.661], omega);
        let expected = (0.743 * 20.0 + omega * 0.661 * 10.0) / (20.0 + omega * 10.0);
        assert!(
            (h[0] - expected).abs() <= 1.0e-15,
            "{omega}: {h:?} vs {expected}"
        );
        // 燃焼度 0 からはステップの値、燃焼度が進まなければそのまま、一定の値ならその値（ビットで同じ）
        assert_eq!(
            advance_history(&[0.743], 0.0, 2.0, &[0.661], omega),
            [0.661]
        );
        assert_eq!(advance_history(&[0.7], 12.0, 12.0, &[0.661], omega), [0.7]);
        assert_eq!(
            advance_history(&[0.743], 13.7, 60.0, &[0.743], omega),
            [0.743]
        );
    }
    // ω が大きいほど、最近のステップの値に寄る
    let h = |omega| advance_history(&[0.743], 20.0, 30.0, &[0.661], omega)[0];
    assert!(h(2.0) < h(1.25) && h(1.25) < h(1.0));
    // ω = 1 は燃焼度で重みを付けた平均の形 H + (x − H) ΔB / B とビットで同じ
    for (hv, start, end, x) in [
        (0.743, 20.0, 22.0, 0.661),
        (0.7183, 0.1, 0.25, 0.70212),
        (0.69, 13.7, 60.0, 0.743),
    ] {
        assert_eq!(
            advance_history(&[hv], start, end, &[x], 1.0),
            [hv + (x - hv) * ((end - start) / end)]
        );
    }
}

#[test]
fn micro_depletion_records_the_moderator_density_history_of_the_nodes() {
    // ノードの減速材の密度は、テーブルの格子点の（温度、密度）を状態の温度で内挿した値で、一定の状態では一定
    let history = micro_run(UNIFORM, [&TRACKED, &TRACKED]).unwrap();
    assert_eq!(
        history.history_variables,
        Some(vec!["moderator-density".to_string()])
    );
    let (_, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let x = crate::state::moderator_density(&tables[0], 565.0).unwrap();
    for point in &history.points {
        for node in point.histories.as_ref().unwrap().iter().flatten() {
            assert_eq!(node, &[x], "{} EFPD", point.efpd);
        }
    }
    let bytes = history.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), history);
    // ミクロ燃焼しない燃焼は履歴の変数を持たない
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let plain = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        None,
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    assert_eq!(plain.history_variables, None);
    assert!(plain.points.iter().all(|p| p.histories.is_none()));
}

#[test]
fn node_histories_follow_the_states_of_the_steps() {
    // 560 K（0.75 g/cm³）で燃焼させ、途中から 570 K（0.74 g/cm³）にすると、ノードの値は燃焼度で重みを付けた平均
    let case = micro_case(
        &CASE
            .replace(
                "[state.hfp]\nboron = 1200\nb10 = 0.199\nmoderator = 565\nfuel = 565\n",
                "[state.cold]\nboron = 1200\nb10 = 0.199\nmoderator = 560\nfuel = 565\n\n\
                 [state.hot]\nboron = 1200\nb10 = 0.199\nmoderator = 570\nfuel = 565\n",
            )
            .replace(
                "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
                "[{ efpd = 50.0, state = \"cold\" }, { efpd = 125.0, state = \"hot\" }, { efpd = 250.0 }]",
            ),
    );
    let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let tables: Vec<Table> = tables
        .into_iter()
        .map(|t| with_micro(t, &TRACKED))
        .collect();
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    let values = |k: usize| -> Vec<(f64, f64)> {
        let p = &history.points[k];
        p.burnups
            .iter()
            .flatten()
            .zip(p.histories.as_ref().unwrap().iter().flatten())
            .map(|(&b, h)| (b, h[0]))
            .collect()
    };
    // 燃焼点 0 は最初のステップの状態の瞬時の値、最初のステップの後もその値
    assert!(values(0).iter().all(|&(_, h)| h == 0.75));
    assert!(values(1).iter().all(|&(_, h)| h == 0.75));
    for ((b1, _), (b3, h3)) in values(1).into_iter().zip(values(3)) {
        let expected = (0.75 * b1 + 0.74 * (b3 - b1)) / b3;
        assert!((h3 - expected).abs() <= 1.0e-15, "{h3} vs {expected}");
    }
}

#[test]
fn node_histories_step_with_the_history_weight_of_the_micro_input() {
    // 560 K（0.75 g/cm³）から 570 K（0.74 g/cm³）に変えた燃焼で、ノードの値は燃焼点ごとに ω の式で進む
    let mut case = micro_case(
        &CASE
            .replace(
                "[state.hfp]\nboron = 1200\nb10 = 0.199\nmoderator = 565\nfuel = 565\n",
                "[state.cold]\nboron = 1200\nb10 = 0.199\nmoderator = 560\nfuel = 565\n\n\
                 [state.hot]\nboron = 1200\nb10 = 0.199\nmoderator = 570\nfuel = 565\n",
            )
            .replace(
                "[{ efpd = 50.0 }, { efpd = 125.0 }, { efpd = 250.0 }]",
                "[{ efpd = 50.0, state = \"cold\" }, { efpd = 125.0, state = \"hot\" }, { efpd = 250.0 }]",
            ),
    );
    let omega = 2.0;
    case.depletion
        .as_mut()
        .unwrap()
        .micro
        .as_mut()
        .unwrap()
        .history_weight = omega;
    let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
    let tables: Vec<Table> = tables
        .into_iter()
        .map(|t| with_micro(t, &TRACKED))
        .collect();
    let history = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap();
    let values = |k: usize| -> Vec<(f64, f64)> {
        let p = &history.points[k];
        p.burnups
            .iter()
            .flatten()
            .zip(p.histories.as_ref().unwrap().iter().flatten())
            .map(|(&b, h)| (b, h[0]))
            .collect()
    };
    assert!(values(1).iter().all(|&(_, h)| h == 0.75));
    for k in 2..=3 {
        for ((b0, h0), (b1, h1)) in values(k - 1).into_iter().zip(values(k)) {
            let expected = (h0 * b0 + omega * 0.74 * (b1 - b0)) / (b0 + omega * (b1 - b0));
            assert!(
                (h1 - expected).abs() <= 1.0e-15,
                "point {k}: {h1} vs {expected}"
            );
        }
    }
}

#[test]
fn restarts_carry_the_node_histories_by_assembly() {
    // 引き継ぐ集合体は前のサイクルの最後の点の履歴の変数を、新しい集合体は燃焼点 0 の瞬時の値を持つ（ADR-0215、ADR-0036）
    let mut first = cycle_run(
        MIXED,
        Some("tracked"),
        "[{ efpd = 50.0 }, { efpd = 125.0 }]",
        None,
    )
    .unwrap();
    // 前のサイクルが別の状態で燃焼したことにする
    let a1 = first.assemblies.iter().position(|a| a == "A1").unwrap();
    let last = first.points.len() - 1;
    for node in &mut first.points[last].histories.as_mut().unwrap()[a1] {
        node[0] = 0.7425;
    }
    let first = CoreHistory::from_msgpack(&first.to_msgpack().unwrap()).unwrap();
    let reloaded = "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\n\
                    restart = { history = \"first.msgpack\", fresh = [\"A3\"] }\n\
                    loading = [\"A2@90  B1\", \"A1  A3\"]\n";
    let second = cycle_run(reloaded, Some("tracked"), "[{ efpd = 50.0 }]", Some(&first)).unwrap();
    let node_values = |h: &CoreHistory, k: usize, id: &str| -> Vec<Vec<f64>> {
        let a = h.assemblies.iter().position(|a| a == id).unwrap();
        h.points[k].histories.as_ref().unwrap()[a].clone()
    };
    for id in ["A1", "A2", "B1"] {
        assert_eq!(
            node_values(&second, 0, id),
            node_values(&first, last, id),
            "{id}"
        );
    }
    assert_eq!(node_values(&second, 0, "A1"), [[0.7425], [0.7425]]);
    let fresh = node_values(&second, 0, "A3");
    assert_eq!(fresh, node_values(&second, 0, "A2"));
    // 履歴の変数が core の変数と違う restart はエラー
    let mut other = first.clone();
    other.history_variables = Some(vec!["rods".into()]);
    let cycle = restarted_cycle(MIXED, "restart = { history = \"first.msgpack\" }");
    let error = cycle_run(&cycle, Some("tracked"), "[{ efpd = 50.0 }]", Some(&other))
        .unwrap_err()
        .to_string();
    assert!(error.contains("history variables"), "{error}");
}

#[test]
fn micro_histories_hold_node_histories_of_the_densities_shape() {
    // 履歴の変数の値は、追う核種の数密度と同じ [集合体][軸方向のノード] の形で、数密度を持つノードだけが変数の数の有限の値を持つ
    let micro = cycle_run(MIXED, Some("tracked"), "[{ efpd = 50.0 }]", None).unwrap();
    let rejected = |h: &CoreHistory| CoreHistory::from_msgpack(&h.to_msgpack().unwrap()).is_err();
    let mut h = micro.clone();
    h.history_variables = None;
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.points[1].histories = None;
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.points[0].histories.as_mut().unwrap()[0].pop();
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.points[1].histories.as_mut().unwrap()[0][0].clear();
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.points[1].histories.as_mut().unwrap()[0][0].push(0.7);
    assert!(rejected(&h));
    let mut h = micro.clone();
    h.points[1].histories.as_mut().unwrap()[0][0][0] = f64::INFINITY;
    assert!(rejected(&h));
}

#[test]
fn micro_tables_hold_the_same_history_residuals() {
    // R を持つテーブルと持たないテーブルを混ぜると、持たないテーブルのノードの履歴を補正しない（黙って落とさない。ADR-0215）
    let (layout, tables, _) = core(CASE, MIXED, Some(400.0));
    let case = micro_case(CASE);
    let tables: Vec<Table> = layout
        .files
        .iter()
        .zip(tables)
        .map(|(f, t)| {
            let mut t = with_micro(t, &TRACKED);
            if f.ends_with("a.msgpack") {
                t.histories = vec![xlc_format::table::HistoryCorrection {
                    variable: xlc_format::table::MODERATOR_DENSITY.into(),
                    nominal: 0.75,
                    branch: 0.74,
                    rods: Vec::new(),
                }];
            }
            t
        })
        .collect();
    let error = deplete(
        &layout,
        &tables,
        &case.states,
        case.depletion.as_ref().unwrap(),
        Some(&small_chain()),
        None,
        &tight(),
        |_| {},
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("history residuals"), "{error}");
}

#[test]
fn node_histories_that_do_not_match_the_step_values_are_errors() {
    // 値を持つ位置や変数の数が食い違えば、履歴を黙って落とさない
    let h = vec![Some(vec![0.743]), None];
    let start = [1.0, 0.0];
    let end = [2.0, 0.0];
    let advanced = advance_histories(&h, &start, &end, &[Some(vec![0.661]), None], 1.0).unwrap();
    assert_eq!(advanced[1], None);
    assert!(advance_histories(&h, &start, &end, &[None, None], 1.0).is_err());
    assert!(
        advance_histories(&h, &start, &end, &[Some(vec![0.661]), Some(vec![0.7])], 1.0).is_err()
    );
    assert!(advance_histories(&h, &start, &end, &[Some(vec![0.661, 0.5]), None], 1.0).is_err());
    assert!(advance_histories(&h, &start, &end, &[Some(vec![0.661])], 1.0).is_err());
}

/// `with_micro` のテーブルに、減速材の密度の履歴の補正（H⁰ = 0.75 の 560 K、Hᵇ = 0.74 の 570 K。燃焼度 0、10、20 の点）を持たせる。
/// U と L は 0（履歴の分岐の数密度は公称のもの）、反応率の微視的断面積の差は U-235 の熱群の捕獲だけ `capture`。軸に公称の値（ホウ素 1000 ppm、560 K、560 K）を置く
fn with_histories(mut table: Table, capture: f64) -> Table {
    table.axes[1] = table.axes[1].clone().with_nominal(1000.0);
    table.axes[2] = table.axes[2].clone().with_nominal(560.0);
    table.axes[3] = table.axes[3].clone().with_nominal(560.0);
    let zero = xlc_format::table::Residual {
        total: vec![0.0; 2],
        absorption: vec![0.0; 2],
        nu_fission: vec![0.0; 2],
        fission: vec![0.0; 2],
        chi: vec![0.0; 2],
        diffusion: vec![0.0; 2],
        scatter: vec![vec![0.0; 2]; 2],
        adf: vec![vec![0.0; 2]; 6],
    };
    let rate = |capture: [f64; 2], thermal: Option<f64>| xlc_format::table::RateResidual {
        fission: vec![0.0; 2],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![capture.to_vec()],
        }],
        thermal_fission: thermal,
    };
    table.histories = vec![xlc_format::table::HistoryCorrection {
        variable: xlc_format::table::MODERATOR_DENSITY.into(),
        nominal: 0.75,
        branch: 0.74,
        rods: vec![xlc_format::table::HistoryRod {
            rod: xlc_format::table::ROD_OUT.into(),
            burnup: BURNUP.to_vec(),
            untracked: vec![zero.clone(); 3],
            linearization: vec![zero; 3],
            densities: table.micro.as_ref().unwrap().densities.clone(),
            micro: vec![vec![rate([0.0, capture], Some(0.0)), rate([0.0, 0.0], None)]; 3],
        }],
    }];
    table.validate().unwrap();
    table
}

#[test]
fn micro_depletion_rates_take_the_weighted_rate_residuals_of_the_history() {
    // 565 K のノードは H = 0.745 g/cm³ で w = 0.5。反応率の U-235 の熱群の捕獲は σ + w Δσ で、補正のないテーブルの捕獲を
    // 0.5 Δσ 増やした燃焼と同じ数密度になる（R は 0 で、線形の補正の係数は補正しないので定数も同じ。ADR-0215 の案 S2）
    let capture = 4.0;
    let deplete_with = |case: &str, change: &dyn Fn(Table) -> Table| {
        let (layout, tables, _) = core(case, UNIFORM, Some(400.0));
        let case = micro_case(case);
        let tables: Vec<Table> = tables
            .into_iter()
            .map(|t| change(with_micro(t, &TRACKED)))
            .collect();
        deplete(
            &layout,
            &tables,
            &case.states,
            case.depletion.as_ref().unwrap(),
            Some(&small_chain()),
            None,
            &tight(),
            |_| {},
        )
        .unwrap()
    };
    let corrected = deplete_with(CASE, &|t| with_histories(t, capture));
    let shifted = deplete_with(CASE, &|mut t| {
        for state in &mut t.micro.as_mut().unwrap().states {
            state[0].reactions[0].channels[0][1] += 0.5 * capture;
        }
        t
    });
    let last = |h: &CoreHistory| h.points.last().unwrap().densities.clone().unwrap();
    for (a, b) in last(&corrected)
        .iter()
        .flatten()
        .zip(last(&shifted).iter().flatten())
    {
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() <= 1e-10 * y.abs(), "{a:?} vs {b:?}");
        }
    }
    // 補正しない燃焼とは違う（U-236 が増える）
    let plain = deplete_with(CASE, &|t| t);
    assert!(last(&corrected)[0][0][1] > 1.001 * last(&plain)[0][0][1]);
    // 公称の履歴（560 K。w = 0）の燃焼は、補正のないテーブルとビットで同じ
    let nominal = CASE.replace("moderator = 565", "moderator = 560");
    let with = deplete_with(&nominal, &|t| with_histories(t, capture));
    let reference = deplete_with(&nominal, &|t| {
        let mut t = with_histories(t, capture);
        t.histories[0].rods[0]
            .micro
            .iter_mut()
            .flatten()
            .for_each(|r| {
                r.reactions[0].channels[0] = vec![0.0; 2];
                r.thermal_fission = r.thermal_fission.map(|_| 0.0);
            });
        t
    });
    assert_eq!(with.points, reference.points);
}

#[test]
fn micro_depletion_rates_take_the_rate_residual_weight_of_the_micro_input() {
    // `rate-residual-weight = "production-index"` の燃焼は、Δσ を P / P_j で足す（ADR-0215 の案 R6）。履歴の分岐の追う核種の数密度が
    // 燃焼度 0 より後も公称と同じテーブルでは、P_j が 0 で比を作れないので止める（Δσ は追わない核種の違いも含み、0 とは限らない）。
    // 数密度が違えば補正しない燃焼と違う
    let capture = 4.0;
    let deplete_with = |weight: RateResidualWeight, change: &dyn Fn(Table) -> Table| {
        let (layout, tables, _) = core(CASE, UNIFORM, Some(400.0));
        let mut case = micro_case(CASE);
        case.depletion
            .as_mut()
            .unwrap()
            .micro
            .as_mut()
            .unwrap()
            .rate_residual_weight = weight;
        let tables: Vec<Table> = tables
            .into_iter()
            .map(|t| change(with_micro(t, &TRACKED)))
            .collect();
        deplete(
            &layout,
            &tables,
            &case.states,
            case.depletion.as_ref().unwrap(),
            Some(&small_chain()),
            None,
            &tight(),
            |_| {},
        )
    };
    let zeroed = |mut t: Table| {
        t.histories[0].rods[0]
            .micro
            .iter_mut()
            .flatten()
            .for_each(|r| {
                r.reactions[0].channels[0] = vec![0.0; 2];
                r.thermal_fission = r.thermal_fission.map(|_| 0.0);
            });
        t
    };
    let index = RateResidualWeight::ProductionIndex;
    let err = deplete_with(index, &|t| with_histories(t, capture)).unwrap_err();
    assert!(err.to_string().contains("index"), "{err}");
    // 同じテーブルでも w の重みは比を使わない
    assert!(deplete_with(RateResidualWeight::History, &|t| with_histories(t, capture)).is_ok());
    // 履歴の分岐の U-235 が公称の半分（燃焼度 0 は公称と同じ。比 P / P_j が Δσ で断面積を負にしない大きさ）
    let shifted = |t: Table| {
        let mut t = with_histories(t, capture);
        for row in t.histories[0].rods[0].densities.iter_mut().skip(1) {
            row[0] *= 0.5;
        }
        t.validate().unwrap();
        t
    };
    let by_index = deplete_with(index, &|t| shifted(t)).unwrap();
    let plain = deplete_with(index, &|t| zeroed(shifted(t))).unwrap();
    assert_ne!(by_index.points, plain.points);
}
