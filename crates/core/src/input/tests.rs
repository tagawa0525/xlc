use super::*;
use crate::nem::{Boundary, Solution};

const PLANT: &str = r#"
pitch = 21.5
height = 10.0
map = [
  "R  *  R",
  "R  *  R",
]

[reflector.R]
bounds = [0.0, 10.0]
layers = ["r.msgpack"]
"#;

const INVENTORY: &str = r#"
[design.A]
bounds = [0.0, 4.0, 10.0]
layers = ["a-bottom.msgpack", "a-top.msgpack"]

[assembly.J01]
design = "A"

[assembly.J02]
design = "A"

[assembly.J03]
design = "A"
"#;

const CYCLE: &str = r#"
plant = "plant.toml"
inventory = "assemblies.toml"
loading = [
  ".  J01      .",
  ".  J02@90   .",
]
"#;

const CASE: &str = r#"
cycle = "cycle-1.toml"
axial-mesh = [0.0, 4.0, 7.0, 10.0]
radial-nodes = 1
boundary = { x = ["vacuum", "vacuum"], y = ["reflective", "reflective"], z = ["reflective", "reflective"] }
"#;

fn edited(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replace(from, to)
}

fn base() -> &'static Path {
    Path::new("dir")
}

fn plant() -> Plant {
    Plant::parse(PLANT, base()).unwrap()
}

fn inventory() -> Inventory {
    Inventory::parse(INVENTORY, base()).unwrap()
}

fn cycle() -> Cycle {
    Cycle::parse(CYCLE, base()).unwrap()
}

/// 分けていないテーブルの、ファイルごとの部分ノードの定数（1つずつ）
fn unsplit(constants: &[FewGroup]) -> Vec<Vec<FewGroup>> {
    constants.iter().map(|c| vec![c.clone()]).collect()
}

fn case() -> Case {
    Case::parse(CASE, base()).unwrap()
}

fn layout_of(case: &Case, cycle: &Cycle, plant: &Plant) -> Result<Layout, InputError> {
    Layout::new(case, cycle, plant, &inventory())
}

#[test]
fn parses_plant() {
    let p = plant();
    assert_eq!(p.pitch, 21.5);
    assert_eq!(p.height, 10.0);
    let r = Site::Reflector {
        name: "R".to_string(),
        turns: 0,
    };
    assert_eq!(
        p.map,
        vec![
            vec![r.clone(), Site::Assembly, r.clone()],
            vec![r.clone(), Site::Assembly, r]
        ]
    );
    let column = &p.reflectors["R"];
    assert_eq!(column.bounds, vec![0.0, 10.0]);
    // パスはファイルのあるディレクトリから辿る
    assert_eq!(column.layers, vec![PathBuf::from("dir/r.msgpack")]);
}

#[test]
fn invalid_plants_are_errors() {
    let bad = |from: &str, to: &str| Plant::parse(&edited(PLANT, from, to), base()).is_err();
    assert!(bad("pitch = 21.5", "pitch = 0.0"));
    assert!(bad("height = 10.0", "height = -1.0"));
    // 行の長さの違い、空の図、知らない名前
    assert!(bad("\"R  *  R\",\n]", "\"R  *\",\n]"));
    assert!(bad(
        "map = [\n  \"R  *  R\",\n  \"R  *  R\",\n]",
        "map = []"
    ));
    assert!(bad("\"R  *  R\",\n]", "\"R  *  W\",\n]"));
    // 層の境界: 0 から始まり、増え、最後が炉心の高さ。層は境界より1つ少ない
    assert!(bad("bounds = [0.0, 10.0]", "bounds = [1.0, 10.0]"));
    assert!(bad("bounds = [0.0, 10.0]", "bounds = [0.0, 9.0]"));
    assert!(bad(
        "bounds = [0.0, 10.0]",
        "bounds = [0.0, 5.0, 5.0, 10.0]"
    ));
    assert!(bad(
        "layers = [\"r.msgpack\"]",
        "layers = [\"r.msgpack\", \"s.msgpack\"]"
    ));
    // 反射体の名前は図の記号と区別できること
    assert!(bad("[reflector.R]", "[reflector.\"R@1\"]"));
    assert!(bad("[reflector.R]", "[reflector.\"*\"]"));
    // 知らないキー
    assert!(bad("height = 10.0", "height = 10.0\nsize = 1"));
}

#[test]
fn reflectors_can_be_rotated() {
    // 反射体の向き（燃料に面する辺）は回転で合わせる
    let text = edited(PLANT, "\"R  *  R\",\n]", "\"R  *  R@180\",\n]");
    let p = Plant::parse(&text, base()).unwrap();
    assert_eq!(
        p.map[1][2],
        Site::Reflector {
            name: "R".to_string(),
            turns: 2
        }
    );
    assert!(
        Plant::parse(
            &edited(PLANT, "\"R  *  R\",\n]", "\"R  *  R@45\",\n]"),
            base()
        )
        .is_err()
    );
    assert!(
        Plant::parse(
            &edited(PLANT, "\"R  *  R\",\n]", "\"R  *@90  R\",\n]"),
            base()
        )
        .is_err()
    );
    let layout = layout_of(&case(), &cycle(), &p).unwrap();
    let m = layout.materials[layout.node_material[2].unwrap()];
    assert_eq!(
        (layout.files[m.file].clone(), m.turns),
        (PathBuf::from("dir/r.msgpack"), 2)
    );
}

/// 外形の面の条件を書いたケース（ADR-0132）
fn outlined_case() -> Case {
    let text = edited(
        CASE,
        "z = [\"reflective\", \"reflective\"] }",
        "z = [\"reflective\", \"reflective\"], outline = \"vacuum\" }",
    );
    Case::parse(&text, base()).unwrap()
}

#[test]
fn dashes_mark_positions_without_nodes() {
    // `-` はノードを置かない位置（ADR-0132）。サイクルの装荷図ではその位置を `.` にし、ケースに外形の面の条件を書く
    let p = Plant::parse(&edited(PLANT, "\"R  *  R\",\n]", "\"-  *  R\",\n]"), base()).unwrap();
    let layout = layout_of(&outlined_case(), &cycle(), &p).unwrap();
    // 図の下の行（y = 0）の左の列（x = 0）。並びは x が最も速く、次に y、z
    for z in 0..3 {
        assert_eq!(layout.node_material[6 * z], None, "z = {z}");
        assert!(layout.node_material[1 + 6 * z].is_some());
        assert!(layout.node_material[3 + 6 * z].is_some());
    }
    assert_eq!(layout.ids[1][0], None);
    // ノードのない位置は、反射体の位置と分けて示す（上の行が +y）
    assert_eq!(
        layout.empty,
        vec![vec![false, false, false], vec![true, false, false]]
    );
    let constants: Vec<FewGroup> = layout.reflector.iter().map(|&r| sided(!r)).collect();
    let problem = layout.problem(&unsplit(&constants)).unwrap();
    assert_eq!(problem.outline, Some(Boundary::Vacuum));
    assert_eq!(problem.node_material, layout.node_material);
}

#[test]
fn invalid_outlines_are_errors() {
    let dashed =
        Plant::parse(&edited(PLANT, "\"R  *  R\",\n]", "\"-  *  R\",\n]"), base()).unwrap();
    // `-` があるのに外形の条件がない、`-` がないのに外形の条件がある
    let err = layout_of(&case(), &cycle(), &dashed)
        .unwrap_err()
        .to_string();
    assert!(err.contains("outline"), "{err}");
    let err = layout_of(&outlined_case(), &cycle(), &plant())
        .unwrap_err()
        .to_string();
    assert!(err.contains("outline"), "{err}");
    // `-` に集合体を置く
    let loaded = Cycle::parse(
        &edited(CYCLE, "\".  J02@90   .\"", "\"J03  J02@90   .\""),
        base(),
    )
    .unwrap();
    let err = layout_of(&outlined_case(), &loaded, &dashed)
        .unwrap_err()
        .to_string();
    assert!(err.contains("J03"), "{err}");
    // ノードに囲まれた `-`（図の端まで `-` だけをたどって出られない）は、炉心の中の真空の穴になる
    let enclosed = edited(
        PLANT,
        "map = [\n  \"R  *  R\",\n  \"R  *  R\",\n]",
        "map = [\n  \"R  *  R\",\n  \"R  -  R\",\n  \"R  *  R\",\n]",
    );
    let err = Plant::parse(&enclosed, base()).unwrap_err().to_string();
    assert!(err.contains("enclosed"), "{err}");
    // 端まで `-` でつながれば受け付ける
    let open = edited(
        PLANT,
        "map = [\n  \"R  *  R\",\n  \"R  *  R\",\n]",
        "map = [\n  \"-  -  R\",\n  \"R  -  R\",\n  \"R  *  R\",\n]",
    );
    assert!(Plant::parse(&open, base()).is_ok());
    // `-` で炉心のほかのノードから切り離されたノードは、別の問題になる
    let cut = edited(
        PLANT,
        "map = [\n  \"R  *  R\",\n  \"R  *  R\",\n]",
        "map = [\n  \"R  -  R\",\n  \"-  -  *\",\n]",
    );
    let err = Plant::parse(&cut, base()).unwrap_err().to_string();
    assert!(err.contains("cut off"), "{err}");
}

#[test]
fn parses_inventory() {
    let inv = inventory();
    assert_eq!(inv.designs["A"].bounds, vec![0.0, 4.0, 10.0]);
    assert_eq!(
        inv.designs["A"].layers,
        vec![
            PathBuf::from("dir/a-bottom.msgpack"),
            PathBuf::from("dir/a-top.msgpack")
        ]
    );
    assert_eq!(inv.assemblies["J02"], "A");
}

#[test]
fn invalid_inventories_are_errors() {
    let bad =
        |from: &str, to: &str| Inventory::parse(&edited(INVENTORY, from, to), base()).is_err();
    // 台帳にない設計
    assert!(bad(
        "[assembly.J03]\ndesign = \"A\"",
        "[assembly.J03]\ndesign = \"B\""
    ));
    // ID は空白、`@`、`.` を含まない
    assert!(bad("[assembly.J03]", "[assembly.\"J 03\"]"));
    assert!(bad("[assembly.J03]", "[assembly.\"J03@1\"]"));
    assert!(bad("[assembly.J03]", "[assembly.\".\"]"));
    // ID の重複は TOML の誤り
    assert!(bad("[assembly.J03]", "[assembly.J02]"));
    assert!(bad("bounds = [0.0, 4.0, 10.0]", "bounds = [0.0, 4.0]"));
    assert!(bad(
        "design = \"A\"\n\n[assembly.J02]",
        "design = \"A\"\nbatch = 1\n\n[assembly.J02]"
    ));
}

#[test]
fn parses_cycle_with_rotations() {
    let c = cycle();
    assert_eq!(c.plant, PathBuf::from("dir/plant.toml"));
    assert_eq!(c.inventory, PathBuf::from("dir/assemblies.toml"));
    let placed = |id: &str, turns| {
        Some(Placed {
            id: id.to_string(),
            turns,
        })
    };
    assert_eq!(
        c.loading,
        vec![
            vec![None, placed("J01", 0), None],
            vec![None, placed("J02", 1), None]
        ]
    );
    for (text, turns) in [("@0", 0), ("@180", 2), ("@270", 3)] {
        let c = Cycle::parse(&edited(CYCLE, "J02@90", &format!("J02{text}")), base()).unwrap();
        assert_eq!(c.loading[1][1].as_ref().map(|p| p.turns), Some(turns));
    }
}

#[test]
fn invalid_cycles_are_errors() {
    let bad = |from: &str, to: &str| Cycle::parse(&edited(CYCLE, from, to), base()).is_err();
    assert!(bad("J02@90", "J02@45"));
    assert!(bad("J02@90", "J02@"));
    assert!(bad("J02@90", "J02@-90"));
    assert!(bad("\".  J01      .\"", "\".  J01\""));
    assert!(bad("inventory = \"assemblies.toml\"\n", ""));
}

#[test]
fn parses_case() {
    let c = case();
    assert_eq!(c.cycle, PathBuf::from("dir/cycle-1.toml"));
    assert_eq!(c.axial_mesh, vec![0.0, 4.0, 7.0, 10.0]);
    assert_eq!(c.radial_nodes, 1);
    assert_eq!(c.boundaries[0], Boundary::Vacuum);
    assert_eq!(c.boundaries[2], Boundary::Reflective);
}

#[test]
fn invalid_cases_are_errors() {
    let bad = |from: &str, to: &str| Case::parse(&edited(CASE, from, to), base()).is_err();
    assert!(bad("[0.0, 4.0, 7.0, 10.0]", "[1.0, 4.0, 7.0, 10.0]"));
    assert!(bad("[0.0, 4.0, 7.0, 10.0]", "[0.0, 7.0, 4.0, 10.0]"));
    assert!(bad("[0.0, 4.0, 7.0, 10.0]", "[0.0]"));
    assert!(bad("radial-nodes = 1", "radial-nodes = 0"));
    assert!(bad("\"vacuum\", \"vacuum\"", "\"vacuum\""));
}

const STATES: &str = r#"
[state.cold]
boron = 1170
b10 = 0.199
moderator = 560
fuel = 560

[state.aro]
boron = 1291
b10 = 0.199
moderator = 565
fuel = 565
"#;

fn with_states() -> String {
    format!("{CASE}{STATES}")
}

fn named(name: &str, boron: f64, temperature: f64) -> CaseState {
    CaseState {
        name: name.into(),
        state: State {
            boron,
            b10: 0.199,
            moderator: temperature,
            fuel: temperature,
        },
    }
}

#[test]
fn parses_named_states_in_name_order() {
    let c = Case::parse(&with_states(), base()).unwrap();
    assert_eq!(
        c.states,
        vec![named("aro", 1291.0, 565.0), named("cold", 1170.0, 560.0)]
    );
    // 状態を書かなければ、軸のないテーブルで1つ解く
    assert!(case().states.is_empty());
}

#[test]
fn invalid_states_are_errors() {
    let bad =
        |from: &str, to: &str| Case::parse(&edited(&with_states(), from, to), base()).is_err();
    // 4つとも必須で、知らない項目はエラー
    assert!(bad("boron = 1170\n", ""));
    assert!(bad("b10 = 0.199\nmoderator = 560", "moderator = 560"));
    assert!(bad("fuel = 560", "fuel = 560\nrods = 1"));
    // 名前は図の記号と区別できること
    assert!(bad("[state.cold]", "[state.\"a b\"]"));
    assert!(bad("[state.cold]", "[state.\"*\"]"));
    // 名前のない1つの [state] と、中身のない [state]
    let single = "[state]\nboron = 1170\nb10 = 0.199\nmoderator = 565\nfuel = 565\n";
    let error = Case::parse(&format!("{CASE}{single}"), base()).unwrap_err();
    assert!(error.to_string().contains("[state.<name>]"), "{error}");
    assert!(Case::parse(&format!("{CASE}[state]\n"), base()).is_err());
}

#[test]
fn pin_powers_need_at_most_one_state() {
    let pins = "pin-powers = \"pins.msgpack\"\n";
    assert!(Case::parse(&format!("{CASE}{pins}{STATES}"), base()).is_err());
    let one = &STATES[..STATES.find("[state.aro]").unwrap()];
    let c = Case::parse(&format!("{CASE}{pins}{one}"), base()).unwrap();
    assert_eq!(c.states, vec![named("cold", 1170.0, 560.0)]);
    assert_eq!(c.pin_powers, Some(PathBuf::from("dir/pins.msgpack")));
}

#[test]
fn split_nodes_wait_for_reflector_constants() {
    // 反射体の定数の分割（ADR-0108 の次の段階）まで、反射体のある炉心は分けられない
    let c = Case::parse(
        &edited(CASE, "radial-nodes = 1", "radial-nodes = 2"),
        base(),
    )
    .unwrap();
    assert_eq!(c.radial_nodes, 2);
    let err = layout_of(&c, &cycle(), &plant()).unwrap_err().to_string();
    assert!(
        err.contains("reflector") && err.contains("ADR-0108"),
        "{err}"
    );
}

#[test]
fn split_layout_places_sub_nodes() {
    // 燃料だけの炉心を 2 × 2 に分ける（ADR-0108）。部分ノードの番号はテーブルの並び（-y の行から、各行は -x から）。
    // 90° 回した集合体は、元の部分ノード 0（-x -y）が +x -y に来る
    let plant = Plant::parse("pitch = 21.5\nheight = 10.0\nmap = [\"*  *\"]\n", base()).unwrap();
    let cycle = Cycle::parse(
        "plant = \"plant.toml\"\ninventory = \"assemblies.toml\"\nloading = [\"J01  J02@90\"]\n",
        base(),
    )
    .unwrap();
    let case = Case::parse(
        &edited(CASE, "radial-nodes = 1", "radial-nodes = 2"),
        base(),
    )
    .unwrap();
    let layout = layout_of(&case, &cycle, &plant).unwrap();
    assert_eq!(layout.mesh.dx, vec![10.75; 4]);
    assert_eq!(layout.mesh.dy, vec![10.75; 2]);
    let at = |x: usize, y: usize| {
        let m = layout.materials[layout.node_material[x + 4 * y].unwrap()];
        (m.turns, m.sub)
    };
    assert_eq!(
        [at(0, 0), at(1, 0), at(0, 1), at(1, 1)],
        [(0, 0), (0, 1), (0, 2), (0, 3)]
    );
    assert_eq!(
        [at(2, 0), at(3, 0), at(2, 1), at(3, 1)],
        [(1, 2), (1, 0), (1, 3), (1, 1)]
    );
    // 位置ごとの ID は集合体の単位
    assert_eq!(
        layout.ids,
        vec![vec![Some("J01".to_string()), Some("J02".to_string())]]
    );
    // 部分ノードの定数を回して置く。ファイルの定数の数が 4 でなければエラー
    let subs: Vec<FewGroup> = (0..4)
        .map(|q| FewGroup {
            k_inf: 1.0 + q as f64,
            ..sided(true)
        })
        .collect();
    let constants = vec![subs.clone(); layout.files.len()];
    let problem = layout.problem(&constants).unwrap();
    for (material, m) in problem.materials.iter().zip(&layout.materials) {
        assert_eq!(*material, rotated(&subs[m.sub], m.turns).unwrap());
    }
    assert!(
        layout
            .problem(&unsplit(&subs[..layout.files.len()]))
            .is_err()
    );
}

#[test]
fn layout_places_layers_and_rotations() {
    let layout = layout_of(&case(), &cycle(), &plant()).unwrap();
    assert_eq!(layout.mesh.dx, vec![21.5; 3]);
    assert_eq!(layout.mesh.dy, vec![21.5; 2]);
    assert_eq!(layout.mesh.dz, vec![4.0, 3.0, 3.0]);
    // ノードの並びは x が最も速く、次に y（下から）、z。図の上の行が +y
    let at = |x: usize, y: usize, z: usize| {
        let m = layout.materials[layout.node_material[x + 3 * (y + 2 * z)].unwrap()];
        (layout.files[m.file].clone(), m.turns)
    };
    let file = |name: &str| PathBuf::from(format!("dir/{name}"));
    // 下の行（y = 0）の中央は J02@90、上の行は J01
    assert_eq!(at(1, 0, 0), (file("a-bottom.msgpack"), 1));
    assert_eq!(at(1, 1, 0), (file("a-bottom.msgpack"), 0));
    assert_eq!(at(1, 1, 1), (file("a-top.msgpack"), 0));
    assert_eq!(at(1, 0, 2), (file("a-top.msgpack"), 1));
    assert_eq!(at(0, 0, 0), (file("r.msgpack"), 0));
    assert_eq!(at(2, 1, 2), (file("r.msgpack"), 0));
    // 同じファイルと回転は1つの物質にまとめる
    assert_eq!(layout.materials.len(), 5);
    assert_eq!(layout.files.len(), 3);
    assert!(
        layout.reflector[layout
            .files
            .iter()
            .position(|f| *f == file("r.msgpack"))
            .unwrap()]
    );
    // 位置ごとの ID（出力の図に使う）
    assert_eq!(
        layout.ids,
        vec![
            vec![None, Some("J01".to_string()), None],
            vec![None, Some("J02".to_string()), None]
        ]
    );
}

#[test]
fn invalid_layouts_are_errors() {
    let with_cycle = |from: &str, to: &str| {
        let c = Cycle::parse(&edited(CYCLE, from, to), base()).unwrap();
        layout_of(&case(), &c, &plant())
    };
    // 同じ ID が2か所、台帳にない ID、`*` の位置に `.`、`*` でない位置に ID
    assert!(with_cycle("J02@90", "J01").is_err());
    assert!(with_cycle("J02@90", "J09").is_err());
    assert!(with_cycle("J02@90", ".").is_err());
    assert!(with_cycle("\".  J01      .\"", "\"J03  J01  .\"").is_err());
    // 行と列の数が炉心の形と違う
    assert!(with_cycle("  \".  J02@90   .\",\n", "").is_err());
    let wide = "loading = [\n  \".  J01  .  .\",\n  \".  J02  .  .\",\n]";
    let c = Cycle::parse(
        &edited(
            CYCLE,
            "loading = [\n  \".  J01      .\",\n  \".  J02@90   .\",\n]",
            wide,
        ),
        base(),
    )
    .unwrap();
    assert!(layout_of(&case(), &c, &plant()).is_err());
    let with_case = |from: &str, to: &str| {
        let c = Case::parse(&edited(CASE, from, to), base()).unwrap();
        layout_of(&c, &cycle(), &plant())
    };
    // 層の境界がメッシュに乗っていない、メッシュの最後が炉心の高さと違う
    assert!(with_case("[0.0, 4.0, 7.0, 10.0]", "[0.0, 5.0, 10.0]").is_err());
    assert!(with_case("[0.0, 4.0, 7.0, 10.0]", "[0.0, 4.0, 7.0, 12.0]").is_err());
}

/// 辺ごとの量が区別できる2群の少数群定数
fn sided(fissile: bool) -> FewGroup {
    let nu = if fissile { 0.1 } else { 0.0 };
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        // 少数群定数の検査は k∞ > 0 を求める（反射体の値も正にしておく）
        k_inf: if fissile { 1.0 } else { 0.01 },
        total: vec![0.5, 1.0],
        absorption: vec![0.01, 0.08],
        nu_fission: vec![nu, nu],
        fission: vec![nu / 2.4, nu / 2.4],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.47, 0.02], vec![0.0, 0.92]],
        // 辺（左、右、上、下、−z、+z）
        adf: (1..=6).map(|s| vec![s as f64, 10.0 * s as f64]).collect(),
        vacuum_ratio: vec![
            None,
            Some(vec![0.4, 0.3]),
            None,
            None,
            Some(vec![0.2, 0.1]),
            None,
        ],
        normalized_adf: vec![true, false, false, false, false, true],
        height: None,
        pins: None,
    }
}

#[test]
fn rotation_turns_sides_counterclockwise() {
    let m = sided(true);
    let r = rotated(&m, 1).unwrap();
    // 90° 回すと、右の辺が上に、上が左に、左が下に、下が右に来る。z の辺は変わらない（ADR-0128）
    let order = |v: &[Vec<f64>]| v.iter().map(|g| g[0]).collect::<Vec<_>>();
    assert_eq!(order(&r.adf), vec![3.0, 4.0, 2.0, 1.0, 5.0, 6.0]);
    assert_eq!(
        r.vacuum_ratio,
        vec![
            None,
            None,
            Some(vec![0.4, 0.3]),
            None,
            Some(vec![0.2, 0.1]),
            None
        ]
    );
    assert_eq!(
        r.normalized_adf,
        vec![false, false, false, true, false, true]
    );
    assert_eq!(rotated(&m, 0).unwrap(), m);
    assert_eq!(rotated(&rotated(&m, 1).unwrap(), 3).unwrap(), m);
    assert_eq!(rotated(&rotated(&m, 2).unwrap(), 2).unwrap(), m);
}

#[test]
fn rotation_turns_pins_counterclockwise() {
    // 90° 回すと (x, y) が (H - y, x) に移る（H はノードの高さ）。列 2 × 行 1 のセルは、列 1 × 行 2 になり、
    // 角（-x -y、+x -y、+x +y、-x +y）は1つずつ先に進む（ADR-0096）
    let m = FewGroup {
        pins: Some(Pins {
            dx: vec![1.0, 2.0],
            dy: vec![2.0],
            form: vec![vec![vec![1.0, 2.0]], vec![vec![3.0, 4.0]]],
            corners: vec![
                vec![1.0, 1.5],
                vec![2.0, 2.5],
                vec![3.0, 3.5],
                vec![4.0, 4.5],
            ],
            normalized_corners: vec![true, false, false, false],
            halves: vec![false, true, false, false],
        }),
        ..sided(true)
    };
    let r = rotated(&m, 1).unwrap().pins.unwrap();
    assert_eq!(r.dx, vec![2.0]);
    assert_eq!(r.dy, vec![1.0, 2.0]);
    assert_eq!(
        r.form,
        vec![vec![vec![1.0], vec![2.0]], vec![vec![3.0], vec![4.0]]]
    );
    assert_eq!(
        r.corners,
        vec![
            vec![4.0, 4.5],
            vec![1.0, 1.5],
            vec![2.0, 2.5],
            vec![3.0, 3.5]
        ]
    );
    // 燃料の CDF で割った角の印は角の量と同じく回す（ADR-0131）
    assert_eq!(r.normalized_corners, vec![false, true, false, false]);
    // 分け合ったセルの印は辺の量と同じく回す（右の辺の印は上の辺に来る）
    assert_eq!(r.halves, vec![false, false, true, false]);
    // 2 行 × 2 列で、行の向きも確かめる: 元の左上（行 1、列 0）は、回すと左下（行 0、列 0）
    let square = FewGroup {
        pins: Some(Pins {
            dx: vec![1.0, 2.0],
            dy: vec![1.5, 0.5],
            form: vec![vec![vec![1.0, 2.0], vec![3.0, 4.0]]; 2],
            corners: vec![vec![1.0, 1.0]; 4],
            normalized_corners: vec![false; 4],
            halves: vec![false; 4],
        }),
        ..sided(true)
    };
    let r = rotated(&square, 1).unwrap().pins.unwrap();
    assert_eq!(r.dx, vec![0.5, 1.5]);
    assert_eq!(r.form[0], vec![vec![3.0, 1.0], vec![4.0, 2.0]]);
    for turns in 1..4 {
        assert_eq!(
            rotated(&rotated(&square, turns).unwrap(), 4 - turns).unwrap(),
            square
        );
        assert_eq!(rotated(&rotated(&m, turns).unwrap(), 4 - turns).unwrap(), m);
    }
}

#[test]
fn problem_rotates_materials_and_rejects_fissile_reflectors() {
    let layout = layout_of(&case(), &cycle(), &plant()).unwrap();
    let constants: Vec<FewGroup> = layout
        .files
        .iter()
        .zip(&layout.reflector)
        .map(|(_, &r)| sided(!r))
        .collect();
    let problem = layout.problem(&unsplit(&constants)).unwrap();
    assert_eq!(problem.mesh, layout.mesh);
    assert_eq!(problem.node_material, layout.node_material);
    assert_eq!(problem.boundaries, case().boundaries);
    for (material, m) in problem.materials.iter().zip(&layout.materials) {
        assert_eq!(*material, rotated(&constants[m.file], m.turns).unwrap());
    }
    // 反射体の層が核分裂すればエラー（燃料は ID を付けて台帳に書く）
    let fissile: Vec<FewGroup> = layout.files.iter().map(|_| sided(true)).collect();
    assert!(layout.problem(&unsplit(&fissile)).is_err());
    // 定数の数がファイルの数と違う
    assert!(layout.problem(&unsplit(&constants[1..])).is_err());
}

#[test]
fn reflectors_and_assemblies_may_share_layer_files() {
    // 下部の反射体の層を、径方向の反射体と集合体の設計で共有する
    let text = edited(
        PLANT,
        "bounds = [0.0, 10.0]\nlayers = [\"r.msgpack\"]",
        "bounds = [0.0, 4.0, 10.0]\nlayers = [\"a-bottom.msgpack\", \"r.msgpack\"]",
    );
    let p = Plant::parse(&text, base()).unwrap();
    let layout = layout_of(&case(), &cycle(), &p).unwrap();
    let shared = layout
        .files
        .iter()
        .position(|f| f == Path::new("dir/a-bottom.msgpack"))
        .unwrap();
    // 反射体に使うファイルは、核分裂の検査の対象になる
    assert!(layout.reflector[shared]);
    let constants: Vec<FewGroup> = layout.reflector.iter().map(|&r| sided(!r)).collect();
    let problem = layout.problem(&unsplit(&constants)).unwrap();
    // 反射体の位置（図の左と右の列）に置いた物質だけが、f* の向きを確かめる反射体の物質（ADR-0095）。
    // 同じファイルでも、集合体の位置のものは別の物質
    for (n, m) in problem.node_material.iter().enumerate() {
        assert_eq!(
            problem.reflectors.contains(&m.unwrap()),
            n % 3 != 1,
            "node {n}"
        );
    }
    let fissile: Vec<FewGroup> = (0..layout.files.len())
        .map(|f| sided(f == shared))
        .collect();
    assert!(layout.problem(&unsplit(&fissile)).is_err());
}

#[test]
fn layer_paths_are_normalized() {
    // 同じファイルを別の綴りで書いても1つのファイルになる
    let text = edited(
        PLANT,
        "layers = [\"r.msgpack\"]",
        "layers = [\"./sub/../r.msgpack\"]",
    );
    let p = Plant::parse(&text, base()).unwrap();
    assert_eq!(
        p.reflectors["R"].layers,
        vec![PathBuf::from("dir/r.msgpack")]
    );
    let inv = Inventory::parse(
        &edited(INVENTORY, "\"a-top.msgpack\"", "\"../dir/a-top.msgpack\""),
        base(),
    )
    .unwrap();
    assert_eq!(
        inv.designs["A"].layers[1],
        PathBuf::from("dir/a-top.msgpack")
    );
}

#[test]
fn rotating_assembly_sites_in_the_plant_points_to_the_loading() {
    let err = Plant::parse(
        &edited(PLANT, "\"R  *  R\",\n]", "\"R  *@90  R\",\n]"),
        base(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("loading"), "{err}");
}

#[test]
fn assembly_fission_rates_follow_the_map() {
    let layout = layout_of(&case(), &cycle(), &plant()).unwrap();
    let constants: Vec<FewGroup> = layout.reflector.iter().map(|&r| sided(!r)).collect();
    let problem = layout.problem(&unsplit(&constants)).unwrap();
    // ノード n の高速群の中性子束を n + 1 にする（熱群は 0）
    let solution = Solution {
        k: 1.0,
        flux: (0..layout.node_material.len())
            .map(|n| vec![n as f64 + 1.0, 0.0])
            .collect(),
        leakage: vec![0.0, 0.0],
        currents: std::array::from_fn(|_| Vec::new()),
        iterations: 1,
    };
    let rates = layout.assembly_fission_rates(&problem, &solution);
    // 上の行の J01 は (x, y) = (1, 1) のノード 4、10、16、下の行の J02 は (1, 0) のノード 1、7、13。
    // 層の厚さは 4、3、3 cm なので、J01 は 5×4 + 11×3 + 17×3 = 104、J02 は 2×4 + 8×3 + 14×3 = 74 に比例する
    let mean = (104.0 + 74.0) / 2.0;
    let expected = [
        vec![None, Some(104.0 / mean), None],
        vec![None, Some(74.0 / mean), None],
    ];
    for (row, want) in rates.iter().zip(&expected) {
        for (got, want) in row.iter().zip(want) {
            match (got, want) {
                (Some(g), Some(w)) => assert!((g - w).abs() < 1e-12, "{g} vs {w}"),
                (None, None) => {}
                _ => panic!("{rates:?}"),
            }
        }
    }
}

#[test]
fn case_may_write_pin_powers() {
    // 再構成したピン出力のファイル（ADR-0096）。パスはケースのファイルのあるディレクトリから
    let text = edited(
        CASE,
        "radial-nodes = 1",
        "radial-nodes = 1\npin-powers = \"pins.msgpack\"",
    );
    let c = Case::parse(&text, base()).unwrap();
    assert_eq!(c.pin_powers, Some(PathBuf::from("dir/pins.msgpack")));
    assert_eq!(case().pin_powers, None);
}

#[test]
fn layers_may_put_marshak_on_their_vacuum_faces() {
    // 径方向の真空の比 r を持たない層（軸方向の反射体のスラブの定数）は、その層の真空の面を Marshak の条件
    // （r = 1/2）にすることを入力で明示する（ADR-0132）。r を持たない径方向の辺に r = 1/2 を置く
    let text = edited(
        PLANT,
        "layers = [\"r.msgpack\"]",
        "layers = [{ file = \"r.msgpack\", vacuum = \"marshak\" }]",
    );
    let p = Plant::parse(&text, base()).unwrap();
    assert_eq!(
        p.reflectors["R"].layers,
        vec![PathBuf::from("dir/r.msgpack")]
    );
    let layout = layout_of(&case(), &cycle(), &p).unwrap();
    let plain = |fissile: bool| {
        let mut m = sided(fissile);
        m.vacuum_ratio = vec![None; 6];
        m.normalized_adf = vec![false; 6];
        m
    };
    let constants: Vec<FewGroup> = layout.reflector.iter().map(|&r| plain(!r)).collect();
    let problem = layout.problem(&unsplit(&constants)).unwrap();
    for (n, m) in problem.node_material.iter().enumerate() {
        let m = &problem.materials[m.unwrap()];
        // 径方向の 4 辺だけ。z の辺は変えない
        let expected = if n % 3 == 1 {
            vec![None; 6]
        } else {
            vec![
                Some(vec![0.5, 0.5]),
                Some(vec![0.5, 0.5]),
                Some(vec![0.5, 0.5]),
                Some(vec![0.5, 0.5]),
                None,
                None,
            ]
        };
        assert_eq!(m.vacuum_ratio, expected, "node {n}");
    }
    // 定数が径方向の辺に r を持てば、Marshak に替えずにエラー
    let mut with_ratio = constants.clone();
    for (c, &r) in with_ratio.iter_mut().zip(&layout.reflector) {
        if r {
            c.vacuum_ratio[1] = Some(vec![0.4, 0.3]);
        }
    }
    let err = layout
        .problem(&unsplit(&with_ratio))
        .unwrap_err()
        .to_string();
    assert!(err.contains("marshak"), "{err}");
    // 書けるのは marshak だけ。知らない項目はエラー
    for bad in [
        "layers = [{ file = \"r.msgpack\", vacuum = \"vacuum\" }]",
        "layers = [{ file = \"r.msgpack\", vacuum = \"marshak\", size = 1 }]",
        "layers = [{ vacuum = \"marshak\" }]",
    ] {
        let text = edited(PLANT, "layers = [\"r.msgpack\"]", bad);
        assert!(Plant::parse(&text, base()).is_err(), "{bad}");
    }
    // 誤りのメッセージは、どの値が悪いかを示す
    let text = edited(
        PLANT,
        "layers = [\"r.msgpack\"]",
        "layers = [{ file = \"r.msgpack\", vacuum = \"marshal\" }]",
    );
    let err = Plant::parse(&text, base()).unwrap_err().to_string();
    assert!(err.contains("marshal") && err.contains("marshak"), "{err}");
}
