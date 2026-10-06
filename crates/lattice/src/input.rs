//! 格子計算の入力ファイル。
//!
//! 1行に `キーワード 値...` を書く。`#` から行末まではコメント。`material`, `cell`, `row` は
//! 何度でも書け、それ以外は1回だけ書く。
//!
//! ```text
//! library m1.lib                 # 多群ライブラリ（cross-sections と排他）
//! # cross-sections c5g7.xs       # 巨視的断面積のファイル（このときは material を書かない）
//! output pincell.fg              # 少数群定数（library のときだけ。few-group-bounds が必要）
//! few-group-bounds 1e7 0.625 1e-5
//! pin-powers pins.txt            # セルごとの核分裂率（平均が 1）
//! pitch 1.26                     # セルのピッチ (cm)
//! azimuthal 32                   # (0, π) の方位角の数
//! spacing 0.05                   # トラック間隔 (cm)
//! boundary right vacuum bottom vacuum    # 外周の辺の境界条件（書かない辺は鏡面反射）
//! material fuel temperature 900 U235 7.0e-4 U238 2.2e-2 O16 4.6e-2   # 核種と個数密度 (1/(barn cm))
//! material water temperature 600 H1 4.9e-2 O16 2.5e-2
//! cell pin fuel 0.4096 3 water - 2       # 物質 外半径 分割数 ...（最後の帯の外半径は - でセルの外周まで）
//! cell box water - 1 grid 4             # 均質なセルは小さな正方形に分けられる（4 × 4）
//! row pin box                            # セルの配置を上の行から
//! row box pin
//! ```
//!
//! セルの帯は外半径の内側を分割数の領域に分ける。外半径のある帯は等体積に、最後の帯は
//! セルの内接円までを半径方向に等間隔に分ける。

use std::path::PathBuf;

use crate::geometry::{Lattice, PinCell};
use crate::tracking::{Boundaries, Boundary};

/// 断面積の出どころ。
#[derive(Debug, Clone, PartialEq)]
pub enum CrossSections {
    /// 多群ライブラリ（`material` で組成を与える）
    Library(PathBuf),
    /// 巨視的断面積のファイル（物質はファイルの名前で参照する）
    Macroscopic(PathBuf),
}

/// ライブラリから作る物質。
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialInput {
    pub name: String,
    /// 温度 (K)。ライブラリにある温度であること
    pub temperature: f64,
    /// 核種と個数密度 (1/(barn cm))
    pub composition: Vec<(String, f64)>,
}

/// セルの同心円の帯。
#[derive(Debug, Clone, PartialEq)]
pub struct Zone {
    pub material: String,
    /// 外半径 (cm)。`None` はセルの外周まで
    pub outer: Option<f64>,
    /// 帯を分ける領域の数
    pub rings: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CellInput {
    pub name: String,
    pub zones: Vec<Zone>,
    /// 均質なセルを分ける小さな正方形の1辺の数（`grid M`、書かなければ 1）
    pub grid: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LatticeInput {
    pub cross_sections: CrossSections,
    /// 少数群定数の出力と、その群境界 (eV)
    pub output: Option<(PathBuf, Vec<f64>)>,
    pub pin_powers: Option<PathBuf>,
    pub pitch: f64,
    pub azimuthal: usize,
    pub spacing: f64,
    pub boundaries: Boundaries,
    pub materials: Vec<MaterialInput>,
    pub cells: Vec<CellInput>,
    /// 上の行から順の、セルの名前の配置
    pub rows: Vec<Vec<String>>,
}

/// 1回だけ書くキーワード
const SINGLE: [&str; 9] = [
    "library",
    "cross-sections",
    "output",
    "few-group-bounds",
    "pin-powers",
    "pitch",
    "azimuthal",
    "spacing",
    "boundary",
];

impl LatticeInput {
    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut single: Vec<(&str, usize, Vec<&str>)> = Vec::new();
        let mut materials = Vec::new();
        let mut cells: Vec<CellInput> = Vec::new();
        let mut rows = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line = index + 1;
            let content = raw.split('#').next().unwrap_or("");
            let tokens: Vec<&str> = content.split_whitespace().collect();
            let Some((&key, values)) = tokens.split_first() else {
                continue;
            };
            let at = |e: String| format!("line {line}: {e}");
            match key {
                "material" => {
                    let material = parse_material(values).map_err(at)?;
                    if materials
                        .iter()
                        .any(|m: &MaterialInput| m.name == material.name)
                    {
                        return Err(at(format!("material `{}` is given twice", material.name)));
                    }
                    materials.push(material);
                }
                "cell" => {
                    let cell = parse_cell(values).map_err(at)?;
                    if cells.iter().any(|c| c.name == cell.name) {
                        return Err(at(format!("cell `{}` is given twice", cell.name)));
                    }
                    cells.push(cell);
                }
                "row" => {
                    if values.is_empty() {
                        return Err(at("`row` needs cell names".into()));
                    }
                    rows.push(values.iter().map(|v| v.to_string()).collect());
                }
                _ if SINGLE.contains(&key) => {
                    if single.iter().any(|(k, ..)| *k == key) {
                        return Err(at(format!("`{key}` is given twice")));
                    }
                    single.push((key, line, values.to_vec()));
                }
                _ => {
                    return Err(at(format!("unknown keyword `{key}`")));
                }
            }
        }

        let get = |key: &str| single.iter().find(|(k, ..)| *k == key);
        let one = |key: &str| -> Result<Option<(usize, &str)>, String> {
            match get(key) {
                None => Ok(None),
                Some((_, line, values)) => match values[..] {
                    [value] => Ok(Some((*line, value))),
                    _ => Err(format!("line {line}: `{key}` takes one value")),
                },
            }
        };
        let required = |key: &str| -> Result<(usize, &str), String> {
            one(key)?.ok_or_else(|| format!("`{key}` is missing"))
        };
        let positive = |key: &str| -> Result<f64, String> {
            let (line, value) = required(key)?;
            parse_positive(value)
                .ok_or_else(|| format!("line {line}: `{key}` must be a positive number"))
        };
        let cross_sections = match (one("library")?, one("cross-sections")?) {
            (Some((_, path)), None) => CrossSections::Library(PathBuf::from(path)),
            (None, Some((_, path))) => CrossSections::Macroscopic(PathBuf::from(path)),
            (Some((line, _)), Some(_)) => {
                return Err(format!(
                    "line {line}: give either `library` or `cross-sections`, not both"
                ));
            }
            (None, None) => return Err("`library` or `cross-sections` is missing".into()),
        };
        let macroscopic = matches!(cross_sections, CrossSections::Macroscopic(_));
        if macroscopic && let Some(m) = materials.first() {
            return Err(format!(
                "material `{}`: materials come from the cross-section file",
                m.name
            ));
        }
        let output = match (one("output")?, get("few-group-bounds")) {
            (None, None) => None,
            (Some((line, _)), _) if macroscopic => {
                return Err(format!(
                    "line {line}: `output` needs a library (group energies)"
                ));
            }
            (Some((_, path)), Some((_, line, bounds))) => {
                let bounds = bounds
                    .iter()
                    .map(|b| b.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| format!("line {line}: invalid few-group bound"))?;
                Some((PathBuf::from(path), bounds))
            }
            (Some((line, _)), None) => {
                return Err(format!("line {line}: `output` needs `few-group-bounds`"));
            }
            (None, Some((_, line, _))) => {
                return Err(format!("line {line}: `few-group-bounds` needs `output`"));
            }
        };
        let azimuthal = {
            let (line, value) = required("azimuthal")?;
            value
                .parse::<usize>()
                .map_err(|_| format!("line {line}: `azimuthal` must be a non-negative integer"))?
        };
        let boundaries = match get("boundary") {
            None => Boundaries::REFLECTIVE,
            Some((_, line, values)) => {
                parse_boundaries(values).map_err(|e| format!("line {line}: {e}"))?
            }
        };
        Ok(Self {
            cross_sections,
            output,
            pin_powers: one("pin-powers")?.map(|(_, path)| PathBuf::from(path)),
            pitch: positive("pitch")?,
            azimuthal,
            spacing: positive("spacing")?,
            boundaries,
            materials,
            cells,
            rows,
        })
    }

    /// 形状と、領域ごとの物質の名前。
    pub fn geometry(&self) -> Result<(Lattice, Vec<String>), String> {
        let half = 0.5 * self.pitch;
        let mut types = Vec::with_capacity(self.cells.len());
        let mut region_materials = Vec::with_capacity(self.cells.len());
        for cell in &self.cells {
            let mut radii = Vec::new();
            let mut names = Vec::new();
            let mut inner = 0.0_f64;
            for zone in &cell.zones {
                if let CrossSections::Library(_) = self.cross_sections
                    && !self.materials.iter().any(|m| m.name == zone.material)
                {
                    return Err(format!(
                        "cell `{}`: material `{}` is not defined",
                        cell.name, zone.material
                    ));
                }
                let n = zone.rings as f64;
                match zone.outer {
                    // 等体積に分ける
                    Some(r) => {
                        radii.extend((1..=zone.rings).map(|k| {
                            (inner * inner + (r * r - inner * inner) * k as f64 / n).sqrt()
                        }))
                    }
                    // セルの内接円までを等間隔に分ける（最後の領域は内接円の外の角を含む）
                    None => {
                        radii.extend((1..zone.rings).map(|k| inner + (half - inner) * k as f64 / n))
                    }
                }
                inner = zone.outer.unwrap_or(half);
                names.extend(std::iter::repeat_n(zone.material.clone(), zone.rings));
            }
            let pin = if cell.grid > 1 {
                names = vec![names[0].clone(); cell.grid * cell.grid];
                PinCell::grid(self.pitch, cell.grid)
            } else {
                PinCell::new(self.pitch, radii)
            }
            .map_err(|e| format!("cell `{}`: {e}", cell.name))?;
            types.push(pin);
            region_materials.push(names);
        }
        let layout = self
            .rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|name| {
                        self.cells
                            .iter()
                            .position(|c| &c.name == name)
                            .ok_or_else(|| format!("cell `{name}` is not defined"))
                    })
                    .collect::<Result<Vec<usize>, String>>()
            })
            .collect::<Result<Vec<_>, String>>()?;
        let materials = layout
            .iter()
            .flatten()
            .flat_map(|&t| region_materials[t].iter().cloned())
            .collect();
        Ok((Lattice::new(types, layout)?, materials))
    }
}

fn parse_positive(value: &str) -> Option<f64> {
    value
        .parse::<f64>()
        .ok()
        .filter(|v| *v > 0.0 && v.is_finite())
}

/// `NAME temperature T (NUCLIDE DENSITY)...`
fn parse_material(values: &[&str]) -> Result<MaterialInput, String> {
    let [name, "temperature", temperature, pairs @ ..] = values else {
        return Err(
            "`material` takes `NAME temperature T` and pairs of nuclide and density".into(),
        );
    };
    let temperature = parse_positive(temperature)
        .ok_or_else(|| format!("material `{name}`: invalid temperature `{temperature}`"))?;
    if pairs.is_empty() || !pairs.len().is_multiple_of(2) {
        return Err(format!(
            "material `{name}`: give pairs of nuclide and density"
        ));
    }
    let composition = pairs
        .chunks(2)
        .map(|pair| {
            let density = pair[1]
                .parse::<f64>()
                .ok()
                .filter(|v| *v >= 0.0 && v.is_finite())
                .ok_or_else(|| format!("material `{name}`: invalid density `{}`", pair[1]))?;
            Ok((pair[0].to_string(), density))
        })
        .collect::<Result<_, String>>()?;
    Ok(MaterialInput {
        name: name.to_string(),
        temperature,
        composition,
    })
}

/// `NAME (MATERIAL OUTER RINGS)... [grid M]`（最後の帯の OUTER は `-`）
fn parse_cell(values: &[&str]) -> Result<CellInput, String> {
    let Some((name, rest)) = values.split_first() else {
        return Err("`cell` needs a name".into());
    };
    let (triples, grid) = match rest {
        [zones @ .., "grid", m] => {
            let m = m
                .parse::<usize>()
                .ok()
                .filter(|&m| m > 0)
                .ok_or_else(|| format!("cell `{name}`: invalid grid `{m}`"))?;
            if !matches!(zones, [_, "-", "1"]) {
                return Err(format!(
                    "cell `{name}`: `grid` needs a homogeneous cell (`MATERIAL - 1`)"
                ));
            }
            (zones, m)
        }
        _ => (rest, 1),
    };
    if triples.is_empty() || !triples.len().is_multiple_of(3) {
        return Err(format!(
            "cell `{name}`: give zones as MATERIAL OUTER-RADIUS RINGS"
        ));
    }
    let count = triples.len() / 3;
    let mut zones = Vec::with_capacity(count);
    let mut inner = 0.0;
    for (k, zone) in triples.chunks(3).enumerate() {
        let last = k + 1 == count;
        let outer = match (zone[1], last) {
            ("-", true) => None,
            ("-", false) => {
                return Err(format!(
                    "cell `{name}`: only the last zone extends to the cell boundary"
                ));
            }
            (_, true) => {
                return Err(format!(
                    "cell `{name}`: the last zone must extend to the cell boundary (`-`)"
                ));
            }
            (r, false) => {
                let r = parse_positive(r)
                    .filter(|&r| r > inner)
                    .ok_or_else(|| format!("cell `{name}`: radii must increase, got `{r}`"))?;
                inner = r;
                Some(r)
            }
        };
        let rings = zone[2]
            .parse::<usize>()
            .ok()
            .filter(|&n| n > 0)
            .ok_or_else(|| format!("cell `{name}`: invalid number of rings `{}`", zone[2]))?;
        zones.push(Zone {
            material: zone[0].to_string(),
            outer,
            rings,
        });
    }
    Ok(CellInput {
        name: name.to_string(),
        zones,
        grid,
    })
}

/// `(SIDE CONDITION)...`
fn parse_boundaries(values: &[&str]) -> Result<Boundaries, String> {
    if values.is_empty() || !values.len().is_multiple_of(2) {
        return Err("`boundary` takes pairs of side and condition".into());
    }
    let mut boundaries = Boundaries::REFLECTIVE;
    for pair in values.chunks(2) {
        let condition = match pair[1] {
            "reflective" => Boundary::Reflective,
            "vacuum" => Boundary::Vacuum,
            other => return Err(format!("unknown boundary condition `{other}`")),
        };
        let side = match pair[0] {
            "left" => &mut boundaries.left,
            "right" => &mut boundaries.right,
            "top" => &mut boundaries.top,
            "bottom" => &mut boundaries.bottom,
            other => return Err(format!("unknown side `{other}`")),
        };
        *side = condition;
    }
    Ok(boundaries)
}

#[cfg(test)]
mod tests {
    use super::*;

    const PIN_CELL: &str = "
# pin cell
library m1.lib
output pincell.fg
few-group-bounds 1e7 0.625 1e-5
pitch 1.26
azimuthal 16
spacing 0.05
material fuel temperature 900 U235 7.0e-4 U238 2.2e-2 O16 4.6e-2   # UO2
material water temperature 600 H1 4.9e-2 O16 2.5e-2
cell pin fuel 0.4 2 water - 2
row pin
";

    const LATTICE: &str = "
cross-sections c5g7.xs
pin-powers pins.txt
pitch 2.0
azimuthal 8
spacing 0.1
boundary right vacuum bottom vacuum
cell pin UO2 0.5 1 clad 0.6 1 water - 1
cell box water - 1
row pin box
row box pin
";

    #[test]
    fn parses_pin_cell() {
        let input = LatticeInput::parse(PIN_CELL).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Library(PathBuf::from("m1.lib"))
        );
        assert_eq!(
            input.output,
            Some((PathBuf::from("pincell.fg"), vec![1e7, 0.625, 1e-5]))
        );
        assert_eq!(input.pin_powers, None);
        assert_eq!(
            (input.pitch, input.azimuthal, input.spacing),
            (1.26, 16, 0.05)
        );
        assert_eq!(input.boundaries, Boundaries::REFLECTIVE);
        assert_eq!(input.materials.len(), 2);
        assert_eq!(input.materials[0].temperature, 900.0);
        assert_eq!(
            input.materials[0].composition[1],
            ("U238".to_string(), 2.2e-2)
        );
        assert_eq!(
            input.cells[0].zones,
            vec![
                Zone {
                    material: "fuel".into(),
                    outer: Some(0.4),
                    rings: 2
                },
                Zone {
                    material: "water".into(),
                    outer: None,
                    rings: 2
                },
            ]
        );
        assert_eq!(input.rows, vec![vec!["pin".to_string()]]);
    }

    #[test]
    fn pin_cell_geometry_splits_fuel_by_volume_and_water_by_radius() {
        let input = LatticeInput::parse(PIN_CELL).unwrap();
        let (lattice, materials) = input.geometry().unwrap();
        assert_eq!(materials, vec!["fuel", "fuel", "water", "water"]);
        // 燃料は等体積の2つのリング、水はペレットとセルの内接円の中点で分ける
        let v = lattice.volumes();
        assert!((v[0] - v[1]).abs() < 1e-14);
        let r = (0.4 + 0.63) / 2.0;
        assert!((v[2] - std::f64::consts::PI * (r * r - 0.16)).abs() < 1e-14);
    }

    #[test]
    fn parses_lattice() {
        let input = LatticeInput::parse(LATTICE).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Macroscopic(PathBuf::from("c5g7.xs"))
        );
        assert_eq!(input.output, None);
        assert_eq!(input.pin_powers, Some(PathBuf::from("pins.txt")));
        assert_eq!(
            input.boundaries,
            Boundaries {
                right: Boundary::Vacuum,
                bottom: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            }
        );
        let (lattice, materials) = input.geometry().unwrap();
        assert_eq!(lattice.width(), 4.0);
        assert_eq!(
            materials,
            vec![
                "UO2", "clad", "water", "water", "water", "UO2", "clad", "water"
            ]
        );
    }

    fn line_of(text: &str, prefix: &str) -> usize {
        text.lines().position(|l| l.starts_with(prefix)).unwrap() + 1
    }

    #[test]
    fn errors_report_lines() {
        let text = PIN_CELL.replace("azimuthal 16", "azimuthal sixteen");
        let err = LatticeInput::parse(&text).unwrap_err();
        assert!(
            err.starts_with(&format!("line {}:", line_of(&text, "azimuthal"))),
            "{err}"
        );
        let text = PIN_CELL.replace("water - 2", "water 0.3 2");
        let err = LatticeInput::parse(&text).unwrap_err();
        assert!(
            err.starts_with(&format!("line {}:", line_of(&text, "cell"))),
            "{err}"
        );
    }

    #[test]
    fn invalid_inputs_are_errors() {
        let parse = LatticeInput::parse;
        // 断面積の出どころは library と cross-sections のどちらか一方
        assert!(parse(&format!("{PIN_CELL}\ncross-sections a.xs")).is_err());
        assert!(parse(&PIN_CELL.replace("library m1.lib", "")).is_err());
        // 巨視的断面積では material を書かない、少数群定数を出さない
        assert!(
            parse(&format!(
                "{LATTICE}\nmaterial UO2 temperature 600 U235 1e-3"
            ))
            .is_err()
        );
        assert!(
            parse(&format!(
                "{LATTICE}\noutput a.fg\nfew-group-bounds 1e7 1e-5"
            ))
            .is_err()
        );
        // output と few-group-bounds は対
        assert!(parse(&PIN_CELL.replace("few-group-bounds 1e7 0.625 1e-5", "")).is_err());
        // 知らないキーワード、2回書いたキーワード、知らない境界の辺と条件
        assert!(parse(&format!("{PIN_CELL}\nalbedo 1")).is_err());
        assert!(parse(&format!("{PIN_CELL}\npitch 2.0")).is_err());
        assert!(parse(&LATTICE.replace("right vacuum", "east vacuum")).is_err());
        assert!(parse(&LATTICE.replace("right vacuum", "right white")).is_err());
        // 帯: 外半径が昇順でない、最後の帯がセルの外周まででない、分割数が 0
        assert!(parse(&LATTICE.replace("clad 0.6 1", "clad 0.4 1")).is_err());
        assert!(parse(&LATTICE.replace("water - 1\ncell box", "water 0.9 1\ncell box")).is_err());
        assert!(parse(&LATTICE.replace("UO2 0.5 1", "UO2 0.5 0")).is_err());
        // 同じ名前のセル・物質
        assert!(parse(&format!("{LATTICE}\ncell box water - 2")).is_err());
        assert!(
            parse(&format!(
                "{PIN_CELL}\nmaterial fuel temperature 600 U238 1e-2"
            ))
            .is_err()
        );
    }

    #[test]
    fn geometry_checks_names_and_layout() {
        let geometry = |text: &str| LatticeInput::parse(text).unwrap().geometry();
        // 配置に知らないセル、行の長さが行の数と違う
        assert!(geometry(&LATTICE.replace("row box pin", "row box rod")).is_err());
        assert!(geometry(&LATTICE.replace("row box pin", "row box pin pin")).is_err());
        // ライブラリのときは、セルの物質が material で定義されていること
        assert!(geometry(&PIN_CELL.replace("water - 2", "steam - 2")).is_err());
        // 帯の外半径がセルの内接円を超える
        assert!(geometry(&LATTICE.replace("clad 0.6 1", "clad 1.2 1")).is_err());
    }

    #[test]
    fn homogeneous_cells_can_be_split_into_a_grid() {
        let text = LATTICE.replace("cell box water - 1", "cell box water - 1 grid 3");
        let input = LatticeInput::parse(&text).unwrap();
        assert_eq!(input.cells[1].grid, 3);
        assert_eq!(input.cells[0].grid, 1);
        let (lattice, materials) = input.geometry().unwrap();
        // ピン 3 領域 + 3 × 3 の箱が2つ + ピン 3 領域
        assert_eq!(lattice.regions(), 24);
        assert_eq!(materials.iter().filter(|m| *m == "water").count(), 2 + 18);
        // grid は同心円のない均質なセル（帯が1つで分割数 1）だけ
        let parse = LatticeInput::parse;
        assert!(
            parse(&LATTICE.replace("water - 1\ncell box", "water - 1 grid 2\ncell box")).is_err()
        );
        assert!(
            parse(&LATTICE.replace("cell box water - 1", "cell box water - 2 grid 2")).is_err()
        );
        assert!(
            parse(&LATTICE.replace("cell box water - 1", "cell box water - 1 grid 0")).is_err()
        );
    }
}
