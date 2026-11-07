use std::path::Path;

use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Table, product_points};

use super::*;
use crate::input::{Case, Cycle, Inventory, Layout, Plant};

const PLANT: &str = r#"
pitch = 20.0
height = 30.0
map = ["R  *  *"]

[reflector.R]
bounds = [0.0, 30.0]
layers = ["r.msgpack"]
"#;

/// 反射体の位置のない炉心（`radial-nodes` で分けられる）
const SPLIT_PLANT: &str = r#"
pitch = 20.0
height = 30.0
map = ["*  *"]
"#;

const SPLIT_CYCLE: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["A1  A2"]
"#;

/// 下の層は核分裂のない層（ノズル）、上の 2 層は燃料
const INVENTORY: &str = r#"
[design.A]
bounds = [0.0, 10.0, 30.0]
layers = ["nozzle.msgpack", "fuel.msgpack"]

[assembly.A1]
design = "A"
[assembly.A2]
design = "A"
"#;

const CYCLE: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = [".  A1  A2"]
"#;

const CASE: &str = r#"
cycle = "cycle.toml"
axial-mesh = [0.0, 10.0, 20.0, 30.0]
radial-nodes = 1
boundary = { x = ["vacuum", "reflective"], y = ["reflective", "reflective"], z = ["vacuum", "vacuum"] }

[thermal-hydraulics]
power = 0.2
flow = 2.0
pressure = 15.5
fuel-heat-fraction = 0.974
effective-temperature = "rowlands"
rod = { pellet-radius = 0.4096, clad-inner-radius = 0.418, clad-outer-radius = 0.475, pellet-density = 10.257, gap-conductance = 1.0e4, rods = 120 }
channel = { flow-area = 2.86, wetted-perimeter = 10.4 }

[state.hfp]
boron = 1100
b10 = 0.199
inlet = 565
"#;

fn few(fissile: f64) -> FewGroup {
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.0,
        total: vec![0.53, 1.06],
        absorption: vec![0.01, 0.06],
        nu_fission: vec![0.005 * fissile, 0.1 * fissile],
        fission: vec![0.002 * fissile, 0.04 * fissile],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.5, 0.02], vec![0.0, 1.0]],
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

fn table(fissile: f64) -> Table {
    let mut table = Table {
        points: Vec::new(),
        split: 1,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        axes: Vec::new(),
        states: vec![few(fissile)],
    };
    table.points = product_points(&table.axes);
    table
}

fn channels() -> (Layout, Channels) {
    layout_channels(CASE, CYCLE, PLANT)
}

/// 反射体のない 2 体を `radial-nodes` = 2 で分けた並び
fn split_channels() -> (Layout, Channels) {
    let case = CASE.replace("radial-nodes = 1", "radial-nodes = 2");
    layout_channels(&case, SPLIT_CYCLE, SPLIT_PLANT)
}

fn layout_channels(case: &str, cycle: &str, plant: &str) -> (Layout, Channels) {
    let base = Path::new("dir");
    let case = Case::parse(case, base).unwrap();
    let cycle = Cycle::parse(cycle, base).unwrap();
    let plant = Plant::parse(plant, base).unwrap();
    let inventory = Inventory::parse(INVENTORY, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables: Vec<Table> = layout
        .files
        .iter()
        .map(|f| {
            table(if f.ends_with("fuel.msgpack") {
                1.0
            } else {
                0.0
            })
        })
        .collect();
    let channels = Channels::new(&layout, &tables).unwrap();
    (layout, channels)
}

/// 2 体の集合体の位置が流路になり、`radial-nodes` = 2 の 4 つのノードが同じ層の状態を持つ
#[test]
fn assembly_positions_become_channels() {
    let (layout, channels) = split_channels();
    assert_eq!(channels.count(), 2);
    let th = layout.thermal().unwrap();
    // 位置ごとの出力: 燃料のノードだけ。下の燃料の層は 1 ノード 5 kW、上は 10 kW（流路あたり 4 × 5 + 4 × 10 = 60 kW）
    let mut powers = vec![0.0; layout.node_material.len()];
    for (site, power) in powers.iter_mut().enumerate() {
        let [_, _, z] = layout.mesh.position(site);
        if channels.channel_of(site).is_some() && channels.fissile(site) {
            *power = if z == 1 { 5.0e3 } else { 10.0e3 };
        }
    }
    let solved = channels.solve(th, 565.0, &powers, None).unwrap();
    let h_in = water::region_one(15.5, 565.0).unwrap().enthalpy;
    let m = 1.0; // 2 kg/s を 2 体に
    for channel in &solved.channels {
        assert!((channel.rise - 60.0e3 / m / 1e3).abs() < 1e-9);
        // ノードの中点のエンタルピー: 下のノズル h_in、下の燃料 h_in + 20/2、上の燃料 h_in + 20 + 40/2（kJ/kg）
        for &site in &channel.sites {
            let [_, _, z] = layout.mesh.position(site);
            let h = h_in + [0.0, 10.0, 40.0][z];
            let t = water::temperature(15.5, h).unwrap();
            let node = solved.nodes[site].unwrap();
            assert!(
                (node.moderator - t).abs() < 1e-9,
                "{site}: {} vs {t}",
                node.moderator
            );
            // 核分裂のないノードの燃料の温度は減速材の温度
            if z == 0 {
                assert_eq!(node.fuel, node.moderator);
            } else {
                assert!(node.fuel > node.moderator);
            }
        }
    }
}

/// 反射体の位置は流路にならず、入口の温度
#[test]
fn reflector_positions_take_the_inlet_temperature() {
    let (layout, channels) = channels();
    assert_eq!(channels.count(), 2);
    let th = layout.thermal().unwrap();
    let powers = vec![0.0; layout.node_material.len()];
    let solved = channels.solve(th, 565.0, &powers, None).unwrap();
    let reflector = (0..powers.len())
        .find(|&s| layout.node_material[s].is_some() && channels.channel_of(s).is_none())
        .unwrap();
    assert_eq!(
        solved.nodes[reflector],
        Some(NodeTemperatures {
            moderator: 565.0,
            fuel: 565.0
        })
    );
}

/// 下に出力のない流路のノードと、出力のない流路の出口は、入口の温度そのもの（エンタルピーからの逆の解きの丸めで入口の温度の
/// 格子点の下に出ると、テーブルの外になる。Problem 6 の出力を下げた解で 564.9999999999999 K になった）
#[test]
fn channel_nodes_without_power_below_take_the_inlet_temperature() {
    let (layout, channels) = channels();
    let th = layout.thermal().unwrap();
    let powers = vec![0.0; layout.node_material.len()];
    for inlet in [555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 590.0, 600.0] {
        let solved = channels.solve(th, inlet, &powers, None).unwrap();
        for site in (0..powers.len()).filter(|&s| channels.channel_of(s).is_some()) {
            assert_eq!(
                solved.nodes[site],
                Some(NodeTemperatures {
                    moderator: inlet,
                    fuel: inlet
                }),
                "site {site} at the inlet {inlet} K"
            );
        }
        for c in &solved.channels {
            assert_eq!((c.rise, c.outlet), (0.0, inlet), "inlet {inlet} K");
        }
    }
}

/// 線出力はノードの出力 × 燃料の発熱の割合 / (棒の本数 / radial-nodes² × 高さ)
#[test]
fn rods_take_the_fuel_share_of_the_node_power() {
    let (layout, channels) = split_channels();
    let th = layout.thermal().unwrap();
    let site = (0..layout.node_material.len())
        .find(|&s| channels.channel_of(s).is_some() && channels.fissile(s))
        .unwrap();
    let mut powers = vec![0.0; layout.node_material.len()];
    powers[site] = 8.0e3;
    let solved = channels.solve(th, 565.0, &powers, None).unwrap();
    let node = solved.nodes[site].unwrap();
    let q = 8.0e3 * 0.974 / (120.0 / 4.0 * 0.1);
    let area = th.flow_area * 1e-4;
    let rod = th
        .rod
        .temperatures(
            q,
            node.moderator,
            15.5,
            1.0 / area,
            4.0 * th.flow_area / th.wetted_perimeter * 1e-2,
            0.0,
            th.effective,
        )
        .unwrap();
    assert_eq!(node.fuel, rod.effective);
    assert!(solved.fuel_max >= node.fuel);
}

/// 燃料棒の熱伝導度は、ノードごとの燃焼度で引く（燃焼のケース。ADR-0229 の燃料棒の 4.）
#[test]
fn rods_take_the_burnup_of_their_node() {
    let (layout, channels) = split_channels();
    let th = layout.thermal().unwrap();
    let site = (0..layout.node_material.len())
        .find(|&s| channels.channel_of(s).is_some() && channels.fissile(s))
        .unwrap();
    let mut powers = vec![0.0; layout.node_material.len()];
    powers[site] = 8.0e3;
    let mut burnups = vec![0.0; layout.node_material.len()];
    burnups[site] = 30.0;
    let solved = channels.solve(th, 565.0, &powers, Some(&burnups)).unwrap();
    let fresh = channels.solve(th, 565.0, &powers, None).unwrap();
    let node = solved.nodes[site].unwrap();
    let q = 8.0e3 * 0.974 / (120.0 / 4.0 * 0.1);
    let rod = |burnup: f64| {
        th.rod
            .temperatures(
                q,
                node.moderator,
                15.5,
                1.0 / (th.flow_area * 1e-4),
                4.0 * th.flow_area / th.wetted_perimeter * 1e-2,
                burnup,
                th.effective,
            )
            .unwrap()
            .effective
    };
    assert_eq!(node.fuel, rod(30.0));
    // 燃焼した UO2 は熱伝導度が下がり、燃料が熱い
    assert!(node.fuel > fresh.nodes[site].unwrap().fuel + 1.0);
    // 燃焼度の数が位置の数と違えばエラー
    assert!(
        channels
            .solve(th, 565.0, &powers, Some(&burnups[1..]))
            .is_err()
    );
}

#[test]
fn negative_or_misplaced_powers_are_errors() {
    let (layout, channels) = channels();
    let th = layout.thermal().unwrap();
    let mut powers = vec![0.0; layout.node_material.len()];
    assert!(channels.solve(th, 565.0, &powers[1..], None).is_err());
    let site = (0..powers.len())
        .find(|&s| channels.channel_of(s).is_some() && channels.fissile(s))
        .unwrap();
    powers[site] = -1.0;
    assert!(channels.solve(th, 565.0, &powers, None).is_err());
    // 核分裂のないノードの出力
    powers[site] = 0.0;
    let nozzle = (0..powers.len())
        .find(|&s| channels.channel_of(s).is_some() && !channels.fissile(s))
        .unwrap();
    powers[nozzle] = 1.0;
    assert!(channels.solve(th, 565.0, &powers, None).is_err());
}

/// 核分裂のある層のない集合体の位置（水の集合体など）は流路にならず、流量は燃料の集合体だけに分ける（ADR-0229）
#[test]
fn positions_without_fuel_take_no_channel_and_no_flow() {
    let inventory = format!(
        "{INVENTORY}\n[design.W]\nbounds = [0.0, 30.0]\nlayers = [\"water.msgpack\"]\n\n[assembly.W1]\ndesign = \"W\"\n"
    );
    let plant = SPLIT_PLANT.replace("map = [\"*  *\"]", "map = [\"*  *  *\"]");
    let cycle = SPLIT_CYCLE.replace("loading = [\"A1  A2\"]", "loading = [\"A1  A2  W1\"]");
    let base = Path::new("dir");
    let case = Case::parse(CASE, base).unwrap();
    let cycle = Cycle::parse(&cycle, base).unwrap();
    let plant = Plant::parse(&plant, base).unwrap();
    let inventory = Inventory::parse(&inventory, base).unwrap();
    let layout = Layout::new(&case, &cycle, &plant, &inventory).unwrap();
    let tables: Vec<Table> = layout
        .files
        .iter()
        .map(|f| {
            table(if f.ends_with("fuel.msgpack") {
                1.0
            } else {
                0.0
            })
        })
        .collect();
    let channels = Channels::new(&layout, &tables).unwrap();
    assert_eq!(channels.count(), 2);
    let th = layout.thermal().unwrap();
    let mut powers = vec![0.0; layout.node_material.len()];
    for (site, power) in powers.iter_mut().enumerate() {
        if channels.channel_of(site).is_some() && channels.fissile(site) {
            *power = 10.0e3;
        }
    }
    let solved = channels.solve(th, 565.0, &powers, None).unwrap();
    // 2 kg/s を燃料の 2 体に: 流路あたり 1 kg/s、2 層の燃料の 20 kW で 20 kJ/kg
    for channel in &solved.channels {
        assert!((channel.rise - 20.0).abs() < 1e-9, "{}", channel.rise);
    }
    // 水の集合体の位置は入口の温度
    let water = (0..powers.len())
        .find(|&s| layout.node_material[s].is_some() && channels.channel_of(s).is_none())
        .unwrap();
    assert_eq!(solved.nodes[water].unwrap().moderator, 565.0);
}
