use super::*;

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
    let m = layout.materials[layout.node_material[2]];
    assert_eq!(
        (layout.files[m.file].clone(), m.turns),
        (PathBuf::from("dir/r.msgpack"), 2)
    );
}

#[test]
fn plant_outline_waits_for_issue_0034() {
    let err =
        Plant::parse(&edited(PLANT, "\"R  *  R\",\n]", "\"-  *  R\",\n]"), base()).unwrap_err();
    assert!(err.to_string().contains("#0034"), "{err}");
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

#[test]
fn split_nodes_wait_for_issue_0033() {
    let err = Case::parse(
        &edited(CASE, "radial-nodes = 1", "radial-nodes = 2"),
        base(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("#0033"), "{err}");
}

#[test]
fn layout_places_layers_and_rotations() {
    let layout = layout_of(&case(), &cycle(), &plant()).unwrap();
    assert_eq!(layout.mesh.dx, vec![21.5; 3]);
    assert_eq!(layout.mesh.dy, vec![21.5; 2]);
    assert_eq!(layout.mesh.dz, vec![4.0, 3.0, 3.0]);
    // ノードの並びは x が最も速く、次に y（下から）、z。図の上の行が +y
    let at = |x: usize, y: usize, z: usize| {
        let m = layout.materials[layout.node_material[x + 3 * (y + 2 * z)]];
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
        k_inf: if fissile { 1.0 } else { 0.0 },
        total: vec![0.5, 1.0],
        absorption: vec![0.01, 0.08],
        nu_fission: vec![nu, nu],
        fission: vec![nu / 2.4, nu / 2.4],
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4, 0.4],
        scatter: vec![vec![0.47, 0.02], vec![0.0, 0.92]],
        // 辺（左、右、上、下）
        adf: (1..=4).map(|s| vec![s as f64, 10.0 * s as f64]).collect(),
        vacuum_ratio: vec![None, Some(vec![0.4, 0.3]), None, None],
        normalized_adf: vec![true, false, false, false],
    }
}

#[test]
fn rotation_turns_sides_counterclockwise() {
    let m = sided(true);
    let r = rotated(&m, 1);
    // 90° 回すと、右の辺が上に、上が左に、左が下に、下が右に来る
    let order = |v: &[Vec<f64>]| v.iter().map(|g| g[0]).collect::<Vec<_>>();
    assert_eq!(order(&r.adf), vec![3.0, 4.0, 2.0, 1.0]);
    assert_eq!(r.vacuum_ratio, vec![None, None, Some(vec![0.4, 0.3]), None]);
    assert_eq!(r.normalized_adf, vec![false, false, false, true]);
    assert_eq!(rotated(&m, 0), m);
    assert_eq!(rotated(&rotated(&m, 1), 3), m);
    assert_eq!(rotated(&rotated(&m, 2), 2), m);
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
    let problem = layout.problem(&constants).unwrap();
    assert_eq!(problem.mesh, layout.mesh);
    assert_eq!(problem.node_material, layout.node_material);
    assert_eq!(problem.boundaries, case().boundaries);
    for (material, m) in problem.materials.iter().zip(&layout.materials) {
        assert_eq!(*material, rotated(&constants[m.file], m.turns));
    }
    // 反射体の層が核分裂すればエラー（燃料は ID を付けて台帳に書く）
    let fissile: Vec<FewGroup> = layout.files.iter().map(|_| sided(true)).collect();
    assert!(layout.problem(&fissile).is_err());
    // 定数の数がファイルの数と違う
    assert!(layout.problem(&constants[1..]).is_err());
}
