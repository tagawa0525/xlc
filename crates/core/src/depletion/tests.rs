use std::path::Path;

use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, Table};

use super::*;
use crate::case::solve_at;
use crate::input::{Case, Cycle, Inventory, Layout, Plant};
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

/// 2 群の定数。核分裂は燃焼度で減り（`fissile` × (1 − 0.02 B)）、熱群の吸収はホウ素濃度で変わる
fn few(fissile: f64, burnup: f64, ppm: f64) -> FewGroup {
    let absorption = vec![0.01, 0.06 + 4.0e-5 * (ppm - 1000.0)];
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
    let mut states = Vec::new();
    for burnup in BURNUP {
        for ppm in PPM {
            for _ in TEMPERATURE {
                for _ in TEMPERATURE {
                    states.push(few(fissile, burnup, ppm));
                }
            }
        }
    }
    Table {
        axes: vec![
            Axis::Burnup {
                points: BURNUP.to_vec(),
            },
            Axis::Boron {
                ppm: PPM.to_vec(),
                b10: 0.199,
            },
            Axis::Moderator {
                temperature: TEMPERATURE.to_vec(),
                density: vec![0.75, 0.74],
            },
            Axis::Fuel {
                temperature: TEMPERATURE.to_vec(),
            },
        ],
        split: 1,
        states,
        heavy_metal,
    }
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
        &case.states[0],
        case.depletion.as_ref().unwrap(),
        &tight(),
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
    // 集合体ごとの 2 層。質量は 400 g/cm × 10 cm
    assert_eq!(history.heavy_metal, vec![vec![4000.0; 2]; 4]);
    let (layout, tables, case) = core(CASE, UNIFORM, Some(400.0));
    let efpd = [0.0, 50.0, 125.0, 250.0];
    assert_eq!(history.points.len(), efpd.len());
    for (point, e) in history.points.iter().zip(efpd) {
        // 燃焼度 = 出力密度 (W/g) × EFPD / 1000。ノードの値は出力の割合の反復の誤差（判定 1e-12）の分だけずれうる
        let want = 40.0 * e / 1000.0;
        assert_eq!(point.efpd, e);
        assert!(
            (point.burnup - want).abs() <= 1e-12 * want.max(1.0),
            "{point:?}"
        );
        let close = |b: f64| (b - want).abs() <= 1e-10 * want.max(1.0);
        assert!(close(point.burnup), "{} vs {want}", point.burnup);
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
    let (layout, tables, case) = core(CASE, UNIFORM, None);
    let result = deplete(
        &layout,
        &tables,
        &case.states[0],
        case.depletion.as_ref().unwrap(),
        &Options::default(),
    );
    assert!(result.is_err());
}

#[test]
fn burnups_past_the_table_are_errors() {
    // 40 W/g × 600 EFPD = 24 MWd/kg はテーブルの 20 MWd/kg を越える（外挿しない。ADR-0187）
    let case = CASE.replace("{ efpd = 250.0 }", "{ efpd = 600.0 }");
    assert!(run(&case, UNIFORM).is_err());
}

#[test]
fn histories_round_trip_through_their_file() {
    let history = run(CASE, MIXED).unwrap();
    let bytes = history.to_msgpack().unwrap();
    assert_eq!(CoreHistory::from_msgpack(&bytes).unwrap(), history);
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
        &case.states[0],
        case.depletion.as_ref().unwrap(),
        &Options::default(),
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
