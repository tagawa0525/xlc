//! 格子計算の入力ファイル（TOML）。
//!
//! ```toml
//! library = "m1.msgpack"          # 多群ライブラリ（cross-sections と排他）
//! # cross-sections = "c5g7.toml"  # 巨視的断面積のファイル（このときは material を書かない）
//! few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }   # 少数群定数（library のとき）
//! pin-powers = "pins.msgpack"     # セルごとの核分裂率（平均が 1）
//! pitch = 1.26                    # セルのピッチ (cm)
//! azimuthal = 32                  # (0, π) の方位角の数
//! spacing = 0.05                  # トラック間隔 (cm)
//! boundary = { right = "vacuum", bottom = "vacuum" }   # 外周の辺の境界条件（書かない辺は鏡面反射）
//! rows = ["pin box", "box pin"]   # セルの配置を上の行から
//!
//! [material.fuel]
//! temperature = 900               # K。ライブラリにある温度であること
//! nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }   # 個数密度 (1/(barn cm))
//!
//! [cell.pin]                      # 帯を内側から。最後の帯は radius を省略してセルの外周まで
//! zones = [{ material = "fuel", radius = 0.4096, rings = 3 }, { material = "water", rings = 2 }]
//! sectors = 8                     # 各リングを 45° ごとの扇形に分ける（1, 4, 8。省略すると 1）
//!
//! [cell.box]
//! zones = [{ material = "water", rings = 1 }]
//! grid = 4                        # 均質なセルを 4 × 4 の正方形に分ける（省略すると 1）
//! ```
//!
//! セルの帯は外半径の内側を分割数の領域に分ける。外半径のある帯は等体積に、最後の帯は
//! セルの内接円までを半径方向に等間隔に分ける。

use std::path::PathBuf;

use std::collections::BTreeMap;

use serde::Deserialize;
use xlc_xs::library::is_nuclide_name;

use crate::geometry::{GeometryError, Lattice, PinCell};
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
    /// 各リングを分ける扇形の数（`sectors S`、書かなければ 1）
    pub sectors: usize,
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

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
    #[error("cell `{cell}`: material `{material}` is not defined")]
    UndefinedMaterial { cell: String, material: String },
    #[error("cell `{0}` is not defined")]
    UndefinedCell(String),
    #[error("cell `{cell}`")]
    Cell { cell: String, source: GeometryError },
    #[error(transparent)]
    Geometry(#[from] GeometryError),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    library: Option<PathBuf>,
    cross_sections: Option<PathBuf>,
    /// 省略すると少数群定数を出さない
    few_group: Option<FewGroupEntry>,
    /// 省略するとピン出力を出さない
    pin_powers: Option<PathBuf>,
    pitch: f64,
    azimuthal: usize,
    spacing: f64,
    /// 省略した辺は鏡面反射
    boundary: Option<BoundaryEntry>,
    rows: Vec<String>,
    /// 巨視的断面積のときは書かない
    material: Option<BTreeMap<String, MaterialEntry>>,
    cell: BTreeMap<String, CellEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FewGroupEntry {
    output: PathBuf,
    bounds: Vec<f64>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
enum Condition {
    Reflective,
    Vacuum,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryEntry {
    left: Option<Condition>,
    right: Option<Condition>,
    top: Option<Condition>,
    bottom: Option<Condition>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MaterialEntry {
    temperature: f64,
    nuclides: BTreeMap<String, f64>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CellEntry {
    zones: Vec<ZoneEntry>,
    /// 省略すると 1（均質なセルだけ）
    grid: Option<usize>,
    /// 省略すると 1
    sectors: Option<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ZoneEntry {
    material: String,
    /// 最後の帯だけ省略する（セルの外周まで）
    radius: Option<f64>,
    rings: usize,
}

fn positive(v: f64) -> bool {
    v > 0.0 && v.is_finite()
}

impl LatticeInput {
    pub fn parse(text: &str) -> Result<Self, InputError> {
        let file: InputFile = toml::from_str(text)?;
        let invalid = |message: String| InputError::Invalid(message);
        let cross_sections = match (file.library, file.cross_sections) {
            (Some(path), None) => CrossSections::Library(path),
            (None, Some(path)) => CrossSections::Macroscopic(path),
            _ => {
                return Err(invalid(
                    "give exactly one of `library` and `cross-sections`".into(),
                ));
            }
        };
        let macroscopic = matches!(cross_sections, CrossSections::Macroscopic(_));
        if macroscopic && file.material.is_some() {
            return Err(invalid("materials come from the cross-section file".into()));
        }
        if macroscopic && file.few_group.is_some() {
            return Err(invalid(
                "`few-group` needs a library (group energies)".into(),
            ));
        }
        for (key, value) in [("pitch", file.pitch), ("spacing", file.spacing)] {
            if !positive(value) {
                return Err(invalid(format!("`{key}` must be positive, got {value}")));
            }
        }
        if file.azimuthal == 0 {
            return Err(invalid("`azimuthal` must be positive".into()));
        }
        let boundaries = match file.boundary {
            None => Boundaries::REFLECTIVE,
            Some(b) => {
                let side = |c: Option<Condition>| match c {
                    None | Some(Condition::Reflective) => Boundary::Reflective,
                    Some(Condition::Vacuum) => Boundary::Vacuum,
                };
                Boundaries {
                    left: side(b.left),
                    right: side(b.right),
                    top: side(b.top),
                    bottom: side(b.bottom),
                }
            }
        };
        let materials = file
            .material
            .unwrap_or_default()
            .into_iter()
            .map(|(name, entry)| entry.into_input(name))
            .collect::<Result<Vec<_>, String>>()
            .map_err(invalid)?;
        let cells = file
            .cell
            .into_iter()
            .map(|(name, entry)| entry.into_input(name))
            .collect::<Result<Vec<_>, String>>()
            .map_err(invalid)?;
        let rows = file
            .rows
            .iter()
            .map(|row| row.split_whitespace().map(String::from).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        if rows.is_empty() || rows.iter().any(Vec::is_empty) {
            return Err(invalid("`rows` needs cell names in every row".into()));
        }
        Ok(Self {
            cross_sections,
            output: file.few_group.map(|f| (f.output, f.bounds)),
            pin_powers: file.pin_powers,
            pitch: file.pitch,
            azimuthal: file.azimuthal,
            spacing: file.spacing,
            boundaries,
            materials,
            cells,
            rows,
        })
    }

    /// 形状と、領域ごとの物質の名前。
    pub fn geometry(&self) -> Result<(Lattice, Vec<String>), InputError> {
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
                    return Err(InputError::UndefinedMaterial {
                        cell: cell.name.clone(),
                        material: zone.material.clone(),
                    });
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
                names = names
                    .into_iter()
                    .flat_map(|n| std::iter::repeat_n(n, cell.sectors))
                    .collect();
                PinCell::with_sectors(self.pitch, radii, cell.sectors)
            }
            .map_err(|source| InputError::Cell {
                cell: cell.name.clone(),
                source,
            })?;
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
                            .ok_or_else(|| InputError::UndefinedCell(name.clone()))
                    })
                    .collect::<Result<Vec<usize>, InputError>>()
            })
            .collect::<Result<Vec<_>, InputError>>()?;
        let materials = layout
            .iter()
            .flatten()
            .flat_map(|&t| region_materials[t].iter().cloned())
            .collect();
        Ok((Lattice::new(types, layout)?, materials))
    }
}

impl MaterialEntry {
    fn into_input(self, name: String) -> Result<MaterialInput, String> {
        if !positive(self.temperature) {
            return Err(format!(
                "material `{name}`: invalid temperature {}",
                self.temperature
            ));
        }
        if self.nuclides.is_empty() {
            return Err(format!("material `{name}` has no nuclides"));
        }
        for (nuclide, &density) in &self.nuclides {
            if !is_nuclide_name(nuclide) {
                return Err(format!(
                    "material `{name}`: nuclide name `{nuclide}` is not of the form U-235, Am-242m"
                ));
            }
            if !(density >= 0.0 && density.is_finite()) {
                return Err(format!(
                    "material `{name}`: invalid density {density} of {nuclide}"
                ));
            }
        }
        Ok(MaterialInput {
            name,
            temperature: self.temperature,
            composition: self.nuclides.into_iter().collect(),
        })
    }
}

impl CellEntry {
    fn into_input(self, name: String) -> Result<CellInput, String> {
        let error = |message: &str| Err(format!("cell `{name}`: {message}"));
        let Some((last, inner)) = self.zones.split_last() else {
            return error("needs at least one zone");
        };
        if last.radius.is_some() {
            return error("the last zone extends to the cell boundary (no radius)");
        }
        let mut previous = 0.0;
        for zone in inner {
            match zone.radius {
                Some(r) if positive(r) && r > previous => previous = r,
                Some(r) => return error(&format!("radii must increase, got {r}")),
                None => return error("only the last zone may omit its radius"),
            }
        }
        if self.zones.iter().any(|z| z.rings == 0) {
            return error("every zone needs at least one ring");
        }
        let (grid, sectors) = match (self.grid, self.sectors) {
            (Some(_), Some(_)) => return error("`grid` and `sectors` exclude each other"),
            (Some(0), _) | (_, Some(0)) => return error("`grid` and `sectors` must be positive"),
            (Some(m), None) => {
                if self.zones.len() != 1 || self.zones[0].rings != 1 {
                    return error("`grid` needs a homogeneous cell (one zone with one ring)");
                }
                (m, 1)
            }
            (None, s) => (1, s.unwrap_or(1)),
        };
        Ok(CellInput {
            name,
            zones: self
                .zones
                .into_iter()
                .map(|z| Zone {
                    material: z.material,
                    outer: z.radius,
                    rings: z.rings,
                })
                .collect(),
            grid,
            sectors,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PIN_CELL: &str = r#"
# pin cell
library = "m1.msgpack"
few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }
pitch = 1.26
azimuthal = 16
spacing = 0.05
rows = ["pin"]

[material.fuel]   # UO2
temperature = 900
nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }

[material.water]
temperature = 600
nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }

[cell.pin]
zones = [{ material = "fuel", radius = 0.4, rings = 2 }, { material = "water", rings = 2 }]
"#;

    const LATTICE: &str = r#"
cross-sections = "c5g7.toml"
pin-powers = "pins.msgpack"
pitch = 2.0
azimuthal = 8
spacing = 0.1
boundary = { right = "vacuum", bottom = "vacuum" }
rows = ["pin box", "box pin"]

[cell.pin]
zones = [
  { material = "UO2", radius = 0.5, rings = 1 },
  { material = "clad", radius = 0.6, rings = 1 },
  { material = "water", rings = 1 },
]

[cell.box]
zones = [{ material = "water", rings = 1 }]
"#;

    fn replaced(text: &str, from: &str, to: &str) -> Result<LatticeInput, InputError> {
        assert!(text.contains(from), "{from}");
        LatticeInput::parse(&text.replace(from, to))
    }

    #[test]
    fn parses_pin_cell() {
        let input = LatticeInput::parse(PIN_CELL).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Library(PathBuf::from("m1.msgpack"))
        );
        assert_eq!(
            input.output,
            Some((PathBuf::from("pincell.fg.msgpack"), vec![1e7, 0.625, 1e-5]))
        );
        assert_eq!(input.pin_powers, None);
        assert_eq!(
            (input.pitch, input.azimuthal, input.spacing),
            (1.26, 16, 0.05)
        );
        assert_eq!(input.boundaries, Boundaries::REFLECTIVE);
        // 物質と核種は名前の順
        assert_eq!(input.materials.len(), 2);
        assert_eq!(input.materials[0].name, "fuel");
        assert_eq!(input.materials[0].temperature, 900.0);
        assert_eq!(
            input.materials[0].composition[1],
            ("U-235".to_string(), 7.0e-4)
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
        assert_eq!((input.cells[0].grid, input.cells[0].sectors), (1, 1));
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
    fn zones_may_reach_beyond_the_inscribed_circle() {
        let geometry = |zones: &str| {
            let text = PIN_CELL.replace(
                "zones = [{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"water\", rings = 2 }]",
                zones,
            );
            LatticeInput::parse(&text).and_then(|i| i.geometry())
        };
        // 内接円（0.63）の外の 0.8 までの水の帯は、セルの外周で切り取った2つのリングになる
        let (lattice, materials) = geometry(
            r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 1 },
]"#,
        )
        .unwrap();
        assert_eq!(materials.len(), 5);
        let v = lattice.volumes();
        assert!((v.iter().sum::<f64>() - 1.26 * 1.26).abs() < 1e-14);
        // 最後の帯は内接円の外から始まるので、等間隔には分けられない
        assert!(
            geometry(
                r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 2 },
]"#
            )
            .is_err()
        );
    }

    #[test]
    fn parses_lattice() {
        let input = LatticeInput::parse(LATTICE).unwrap();
        assert_eq!(
            input.cross_sections,
            CrossSections::Macroscopic(PathBuf::from("c5g7.toml"))
        );
        assert_eq!(input.output, None);
        assert_eq!(input.pin_powers, Some(PathBuf::from("pins.msgpack")));
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
        // セルは配置の順（上の行から）
        assert_eq!(
            materials,
            vec![
                "UO2", "clad", "water", "water", "water", "UO2", "clad", "water"
            ]
        );
    }

    #[test]
    fn grid_and_sectors_split_cells() {
        let text = LATTICE.replace(
            "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]",
            "[cell.box]\nzones = [{ material = \"water\", rings = 1 }]\ngrid = 3",
        );
        let input = LatticeInput::parse(&text).unwrap();
        let (lattice, _) = input.geometry().unwrap();
        assert_eq!(lattice.regions(), 3 + 9 + 9 + 3);
        let text = LATTICE.replace(
            "  { material = \"water\", rings = 1 },\n]",
            "  { material = \"water\", rings = 1 },\n]\nsectors = 8",
        );
        let (lattice, _) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
        assert_eq!(lattice.regions(), 2 * 24 + 2);
    }

    #[test]
    fn invalid_inputs_are_errors() {
        // 断面積の出どころは library と cross-sections のどちらか一方
        assert!(LatticeInput::parse(&format!("cross-sections = \"a.toml\"\n{PIN_CELL}")).is_err());
        assert!(replaced(PIN_CELL, "library = \"m1.msgpack\"\n", "").is_err());
        // 巨視的断面積では material を書かない、少数群定数を出さない
        assert!(
            LatticeInput::parse(&format!(
                "{LATTICE}\n[material.UO2]\ntemperature = 600\nnuclides = {{ U-235 = 1e-3 }}"
            ))
            .is_err()
        );
        assert!(
            replaced(
                LATTICE,
                "pitch = 2.0",
                "pitch = 2.0\nfew-group = { output = \"a\", bounds = [1e7, 1e-5] }"
            )
            .is_err()
        );
        // 知らないキー、知らない辺と条件
        assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 1.26\nalbedo = 1").is_err());
        assert!(replaced(LATTICE, "right = \"vacuum\"", "east = \"vacuum\"").is_err());
        assert!(replaced(LATTICE, "right = \"vacuum\"", "right = \"white\"").is_err());
        // 正でない値、温度、個数密度、核種名
        assert!(replaced(PIN_CELL, "pitch = 1.26", "pitch = 0").is_err());
        assert!(replaced(PIN_CELL, "spacing = 0.05", "spacing = -0.05").is_err());
        assert!(replaced(PIN_CELL, "azimuthal = 16", "azimuthal = 0").is_err());
        assert!(replaced(PIN_CELL, "temperature = 900", "temperature = -900").is_err());
        assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U-238 = -2.2e-2").is_err());
        assert!(replaced(PIN_CELL, "U-238 = 2.2e-2", "U238 = 2.2e-2").is_err());
        assert!(
            replaced(
                PIN_CELL,
                "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
                "nuclides = {}"
            )
            .is_err()
        );
        // 帯: 外半径が昇順でない、最後の帯に外半径がある、途中の帯に外半径がない、分割数が 0
        assert!(replaced(LATTICE, "radius = 0.6", "radius = 0.4").is_err());
        assert!(
            replaced(
                PIN_CELL,
                "{ material = \"water\", rings = 2 }",
                "{ material = \"water\", radius = 0.6, rings = 2 }"
            )
            .is_err()
        );
        assert!(
            replaced(
                PIN_CELL,
                "{ material = \"fuel\", radius = 0.4, rings = 2 }",
                "{ material = \"fuel\", rings = 2 }"
            )
            .is_err()
        );
        assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 0 }]").is_err());
        // grid は均質なセル（帯が1つで分割数 1）だけ、grid と sectors は同時に使わない
        assert!(replaced(PIN_CELL, "rings = 2 }]", "rings = 2 }]\ngrid = 2").is_err());
        assert!(
            replaced(
                LATTICE,
                "zones = [{ material = \"water\", rings = 1 }]",
                "zones = [{ material = \"water\", rings = 1 }]\ngrid = 2\nsectors = 4"
            )
            .is_err()
        );
        assert!(
            replaced(
                LATTICE,
                "zones = [{ material = \"water\", rings = 1 }]",
                "zones = [{ material = \"water\", rings = 1 }]\ngrid = 0"
            )
            .is_err()
        );
    }

    #[test]
    fn geometry_checks_names_and_layout() {
        let geometry = |text: String| LatticeInput::parse(&text).and_then(|i| i.geometry());
        // 配置に知らないセル、行の長さが行の数と違う
        assert!(geometry(LATTICE.replace("\"box pin\"", "\"box rod\"")).is_err());
        assert!(geometry(LATTICE.replace("\"box pin\"", "\"box pin pin\"")).is_err());
        // ライブラリのときは、セルの物質が material で定義されていること
        assert!(
            geometry(PIN_CELL.replace(
                "{ material = \"water\", rings = 2 }",
                "{ material = \"steam\", rings = 2 }"
            ))
            .is_err()
        );
        // 帯の外半径がセルの外接円（ピッチ 2 で √2）を超える、扇形の数が 1, 4, 8 でない
        assert!(geometry(LATTICE.replace("radius = 0.6", "radius = 1.5")).is_err());
        assert!(
            geometry(LATTICE.replace(
                "  { material = \"water\", rings = 1 },\n]",
                "  { material = \"water\", rings = 1 },\n]\nsectors = 3"
            ))
            .is_err()
        );
    }

    #[test]
    fn cell_errors_keep_the_cause_in_the_chain_only() {
        // 原因の文は連鎖（source）にだけ置く。自分の文にも埋め込むと、連鎖を表示するときに2回出る
        let text = LATTICE.replace("radius = 0.6", "radius = 1.2");
        let err = LatticeInput::parse(&text).unwrap().geometry().unwrap_err();
        let cause = std::error::Error::source(&err).expect("the geometry error is the cause");
        assert!(!err.to_string().contains(&cause.to_string()), "{err}");
    }

    #[test]
    fn toml_errors_are_not_repeated_in_the_chain() {
        // 包んだ TOML のエラーは、自分の文と原因の連鎖のどちらか一方にだけ出す
        let err = LatticeInput::parse("pitch = ").unwrap_err();
        let chain = format!("{:#}", anyhow::Error::from(err));
        assert_eq!(chain.matches("TOML parse error").count(), 1, "{chain}");
    }
}
