//! 格子計算の入力ファイル（TOML）。
//!
//! ```toml
//! library = "m1.msgpack"          # 多群ライブラリ（cross-sections と排他）
//! self-shielding = false          # 燃料の自己遮蔽（library のとき。省略すると true）
//! # cross-sections = "c5g7.toml"  # 巨視的断面積のファイル（このときは material を書かない）
//! few-group = { output = "pincell.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }   # 少数群定数（library のとき）
//! pin-powers = "pins.msgpack"     # セルごとの核分裂率（平均が 1）
//! pitch = 1.26                    # セルのピッチ (cm)
//! azimuthal = 32                  # (0, π) の方位角の数
//! spacing = 0.05                  # トラック間隔 (cm)
//! scattering-order = 1            # 散乱の非等方性のルジャンドルの次数（library のとき。省略するとライブラリにあるすべての次数）
//! boundary = { right = "vacuum", bottom = "vacuum" }   # 外周の辺の境界条件（書かない辺は鏡面反射）
//! rows = ["pin box", "box pin"]   # セルの配置を上の行から
//! gap = { width = 0.04, material = "water" }   # 格子の周りの余白（集合体の間の水のギャップ。省略すると無し）
//!
//! [material.fuel]
//! temperature = 900               # K。ライブラリにある温度であること
//! nuclides = { U-235 = 7.0e-4, U-238 = 2.2e-2, O-16 = 4.6e-2 }   # 個数密度 (1/(barn cm))
//!
//! [cell.pin]                      # 帯を内側から。最後の帯は radius を省略してセルの外周まで
//! zones = [{ material = "fuel", radius = 0.4096, rings = 3 }, { material = "void", radius = 0.418, rings = 1 },
//!          { material = "water", rings = 2 }]   # void は真空（定義しない予約の名前）
//! sectors = 8                     # 各リングを 45° ごとの扇形に分ける（1, 4, 8。省略すると 1）
//! frame = { thickness = 0.0141, material = "zircaloy" }   # 外周の正方形の枠（支持格子。省略すると無し）
//!
//! [cell.box]
//! zones = [{ material = "water", rings = 1 }]
//! grid = 4                        # 均質なセルを 4 × 4 の正方形に分ける（省略すると 1）
//! ```
//!
//! セルの帯は外半径の内側を分割数の領域に分ける。外半径のある帯は等体積に（切り取る前の円の
//! 面積で）、最後の帯はセルの内接円までを半径方向に等間隔に分ける。外半径は外接円（ピッチ/√2）
//! まで書け、内接円より大きい円はセルの外周で切り取る。四隅の減速材を半径方向に分けるには、
//! 内接円より大きい外半径の帯を置き、最後の帯は分割数 1 にする。

use std::path::PathBuf;

use std::collections::BTreeMap;

use serde::Deserialize;
use xlc_format::library::is_nuclide_name;

use crate::geometry::{Cell, GeometryError, Lattice, PinCell};
use crate::tracking::{Boundaries, Boundary};

/// 真空の帯の物質の名前（予約語。この名前の物質は定義できない）
pub const VOID: &str = "void";

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
    /// 外周の正方形の枠（支持格子を体積を保って薄い枠にしたもの）
    pub frame: Option<Frame>,
}

/// セルの外周の枠。
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    /// 厚さ (cm)
    pub thickness: f64,
    pub material: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LatticeInput {
    pub cross_sections: CrossSections,
    /// 燃料（核分裂する核種を含む物質）をピンごとに自己遮蔽する（ライブラリのとき）
    pub self_shielding: bool,
    /// 少数群定数の出力と、その群境界 (eV)
    pub output: Option<(PathBuf, Vec<f64>)>,
    pub pin_powers: Option<PathBuf>,
    pub pitch: f64,
    pub azimuthal: usize,
    pub spacing: f64,
    /// 散乱の非等方性のルジャンドルの次数（0 なら等方散乱。`None` ならライブラリにあるすべての次数）
    pub scattering_order: Option<usize>,
    pub boundaries: Boundaries,
    pub materials: Vec<MaterialInput>,
    pub cells: Vec<CellInput>,
    /// 上の行から順の、セルの名前の配置
    pub rows: Vec<Vec<String>>,
    /// 列の幅（左から, cm）。`None` ならどの列もピッチ
    pub widths: Option<Vec<f64>>,
    /// 行の高さ（上から, cm）。`None` ならどの行もピッチ
    pub heights: Option<Vec<f64>>,
    /// 格子の周りの余白
    pub gap: Option<Gap>,
}

/// 格子の周りの余白（集合体の間の水のギャップなど）。
#[derive(Debug, Clone, PartialEq)]
pub struct Gap {
    /// 幅 (cm)
    pub width: f64,
    pub material: String,
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
    #[error(
        "cell `{cell}`: the last zone starts at {start} cm, outside the inscribed circle, \
         so it cannot be split into rings; give the zone an outer radius instead"
    )]
    LastZoneOutside { cell: String, start: f64 },
    #[error("cell `{cell}`")]
    Cell { cell: String, source: GeometryError },
    #[error(transparent)]
    Geometry(#[from] GeometryError),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    library: Option<PathBuf>,
    /// 省略すると、ライブラリのときは true
    self_shielding: Option<bool>,
    cross_sections: Option<PathBuf>,
    /// 省略すると少数群定数を出さない
    few_group: Option<FewGroupEntry>,
    /// 省略するとピン出力を出さない
    pin_powers: Option<PathBuf>,
    pitch: f64,
    azimuthal: usize,
    spacing: f64,
    /// 省略すると、ライブラリにあるすべての次数（巨視的断面積では 0）
    scattering_order: Option<usize>,
    /// 省略した辺は鏡面反射
    boundary: Option<BoundaryEntry>,
    rows: Vec<String>,
    /// 省略すると、どの列もピッチの幅
    widths: Option<Vec<SizeEntry>>,
    /// 省略すると、どの行もピッチの高さ
    heights: Option<Vec<SizeEntry>>,
    /// 省略すると余白なし
    gap: Option<GapEntry>,
    /// 巨視的断面積のときは書かない
    material: Option<BTreeMap<String, MaterialEntry>>,
    cell: BTreeMap<String, CellEntry>,
}

/// 同じ大きさの列（行）の並び
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SizeEntry {
    size: f64,
    /// 省略すると 1
    count: Option<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GapEntry {
    width: f64,
    material: String,
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
    frame: Option<FrameEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FrameEntry {
    thickness: f64,
    material: String,
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
        if macroscopic && file.self_shielding.is_some() {
            return Err(invalid(
                "`self-shielding` needs a library (nuclide data)".into(),
            ));
        }
        if macroscopic && file.few_group.is_some() {
            return Err(invalid(
                "`few-group` needs a library (group energies)".into(),
            ));
        }
        let scattering_order = file.scattering_order;
        if macroscopic && scattering_order.is_some_and(|order| order > 0) {
            return Err(invalid(
                "`scattering-order` needs a library (macroscopic cross sections are P0 only)"
                    .into(),
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
        if file.material.as_ref().is_some_and(|m| m.contains_key(VOID)) {
            return Err(invalid(format!(
                "`{VOID}` is reserved for void zones and cannot be defined"
            )));
        }
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
        if let Some(gap) = &file.gap {
            if !positive(gap.width) {
                return Err(invalid(format!(
                    "the gap width must be positive, got {}",
                    gap.width
                )));
            }
            if gap.material == VOID {
                return Err(invalid(format!(
                    "the gap must be a material, not `{VOID}`; leave out `gap` for no gap"
                )));
            }
        }
        // 列の幅と行の高さ（同じ大きさの並びを展開する）
        let sizes = |entries: Option<Vec<SizeEntry>>, what: &str| -> Result<_, InputError> {
            let Some(entries) = entries else {
                return Ok(None);
            };
            let mut sizes = Vec::new();
            for e in entries {
                let count = e.count.unwrap_or(1);
                if !positive(e.size) || count == 0 {
                    return Err(invalid(format!(
                        "{what} need positive sizes and counts, got {} x {count}",
                        e.size
                    )));
                }
                sizes.extend(std::iter::repeat_n(e.size, count));
            }
            Ok(Some(sizes))
        };
        let widths = sizes(file.widths, "widths")?;
        let heights = sizes(file.heights, "heights")?;
        if (widths.is_some() || heights.is_some()) && file.gap.is_some() {
            return Err(invalid(
                "`gap` and `widths`/`heights` exclude each other; give the gap as columns and rows"
                    .into(),
            ));
        }
        let gap = file.gap.map(|g| Gap {
            width: g.width,
            material: g.material,
        });
        Ok(Self {
            self_shielding: !macroscopic && file.self_shielding.unwrap_or(true),
            cross_sections,
            output: file.few_group.map(|f| (f.output, f.bounds)),
            pin_powers: file.pin_powers,
            pitch: file.pitch,
            azimuthal: file.azimuthal,
            spacing: file.spacing,
            scattering_order,
            boundaries,
            materials,
            cells,
            rows,
            widths,
            heights,
            gap,
        })
    }

    /// 形状と、領域ごとの物質の名前。
    pub fn geometry(&self) -> Result<(Lattice, Vec<String>), InputError> {
        let mut types = Vec::with_capacity(self.cells.len());
        let mut region_materials = Vec::with_capacity(self.cells.len());
        // [セルの種類] の、ピッチに合わないマスに置くときの長方形のセル（均質なセルだけ）
        let mut blocks: Vec<Option<Cell>> = Vec::with_capacity(self.cells.len());
        for cell in &self.cells {
            // 枠の厚さは、それを使って帯を分ける前に確かめる（帯の誤りとして報告しないため）
            if let Some(frame) = &cell.frame
                && !(frame.thickness > 0.0 && frame.thickness < 0.5 * self.pitch)
            {
                return Err(InputError::Cell {
                    cell: cell.name.clone(),
                    source: GeometryError::Frame(frame.thickness),
                });
            }
            // 最後の帯を等間隔に分ける内接円（枠があれば枠の内側の正方形の内接円）の半径
            let half = 0.5 * self.pitch - cell.frame.as_ref().map_or(0.0, |f| f.thickness);
            let materials = cell
                .zones
                .iter()
                .map(|z| &z.material)
                .chain(cell.frame.as_ref().map(|f| &f.material));
            let mut radii = Vec::new();
            let mut names = Vec::new();
            let mut inner = 0.0_f64;
            for material in materials {
                if let CrossSections::Library(_) = self.cross_sections
                    && material != VOID
                    && !self.materials.iter().any(|m| &m.name == material)
                {
                    return Err(InputError::UndefinedMaterial {
                        cell: cell.name.clone(),
                        material: material.clone(),
                    });
                }
            }
            for zone in &cell.zones {
                let n = zone.rings as f64;
                match zone.outer {
                    // 等体積に分ける
                    Some(r) => {
                        radii.extend((1..=zone.rings).map(|k| {
                            (inner * inner + (r * r - inner * inner) * k as f64 / n).sqrt()
                        }))
                    }
                    // セルの内接円までを等間隔に分ける（最後の領域は内接円の外の角を含む）。
                    // 内接円の外から始まる帯は、分ける幅がない
                    None if zone.rings > 1 && inner >= half => {
                        return Err(InputError::LastZoneOutside {
                            cell: cell.name.clone(),
                            start: inner,
                        });
                    }
                    None => {
                        radii.extend((1..zone.rings).map(|k| inner + (half - inner) * k as f64 / n))
                    }
                }
                inner = zone.outer.unwrap_or(half);
                names.extend(std::iter::repeat_n(zone.material.clone(), zone.rings));
            }
            // 均質なセル（帯が1つで分割1、扇形と支持格子なし）は、ピッチの正方形のマスではピンセル（同心円なし）、
            // ほかのマスでは長方形のセルにする（種類はマスで決める。下の `blocks`）
            let homogeneous = cell.zones.len() == 1
                && cell.zones[0].rings == 1
                && cell.sectors == 1
                && cell.frame.is_none();
            blocks.push(homogeneous.then_some(Cell::Block {
                nx: cell.grid,
                ny: cell.grid,
            }));
            let pin = if cell.grid > 1 {
                names = vec![names[0].clone(); cell.grid * cell.grid];
                PinCell::grid(self.pitch, cell.grid).map(Cell::Pin)
            } else {
                names = names
                    .into_iter()
                    .flat_map(|n| std::iter::repeat_n(n, cell.sectors))
                    .collect();
                let pin = PinCell::with_sectors(self.pitch, radii, cell.sectors);
                match &cell.frame {
                    Some(frame) => {
                        names.push(frame.material.clone());
                        pin.and_then(|p| p.with_frame(frame.thickness))
                    }
                    None => pin,
                }
                .map(Cell::Pin)
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
        let columns = layout.first().map_or(0, Vec::len);
        let pitches = |n: usize| vec![self.pitch; n];
        let (mut widths, mut heights) = (
            self.widths.clone().unwrap_or_else(|| pitches(columns)),
            self.heights
                .clone()
                .unwrap_or_else(|| pitches(layout.len())),
        );
        // ピッチの正方形でないマスの均質なセルは、長方形のセルの種類に置き換える
        let mut layout = layout;
        let mut block_type = vec![None; types.len()];
        for (r, row) in layout.iter_mut().enumerate() {
            for (c, t) in row.iter_mut().enumerate() {
                let slot = (widths.get(c), heights.get(r));
                if slot == (Some(&self.pitch), Some(&self.pitch)) {
                    continue;
                }
                if let Some(block) = &blocks[*t] {
                    let index = *block_type[*t].get_or_insert_with(|| {
                        types.push(block.clone());
                        region_materials.push(region_materials[*t].clone());
                        types.len() - 1
                    });
                    *t = index;
                }
            }
        }
        if let Some(gap) = &self.gap {
            if let CrossSections::Library(_) = self.cross_sections
                && !self.materials.iter().any(|m| m.name == gap.material)
            {
                return Err(InputError::UndefinedMaterial {
                    cell: "gap".into(),
                    material: gap.material.clone(),
                });
            }
            // 余白は周りの薄い行と列の均質なセル
            let border = types.len();
            types.push(Cell::BLOCK);
            region_materials.push(vec![gap.material.clone()]);
            let edge = vec![border; columns + 2];
            let mut bordered = vec![edge.clone()];
            bordered.extend(layout.into_iter().map(|row| {
                std::iter::once(border)
                    .chain(row)
                    .chain(std::iter::once(border))
                    .collect()
            }));
            bordered.push(edge);
            layout = bordered;
            for sizes in [&mut widths, &mut heights] {
                sizes.insert(0, gap.width);
                sizes.push(gap.width);
            }
        }
        let materials: Vec<String> = layout
            .iter()
            .flatten()
            .flat_map(|&t| region_materials[t].iter().cloned())
            .collect();
        let lattice = Lattice::rectilinear(widths, heights, types, layout)?;
        Ok((lattice, materials))
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
        if self.nuclides.values().all(|&density| density == 0.0) {
            return Err(format!(
                "material `{name}` has no nuclide with positive density; write `{VOID}` for void zones"
            ));
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
        if self.frame.is_some() && grid > 1 {
            return error("`frame` needs a cell with rings, not `grid`");
        }
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
            frame: self.frame.map(|f| Frame {
                thickness: f.thickness,
                material: f.material,
            }),
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
    fn void_zones_need_no_material() {
        let input = replaced(
            PIN_CELL,
            "{ material = \"fuel\", radius = 0.4, rings = 2 }",
            "{ material = \"fuel\", radius = 0.4, rings = 2 }, { material = \"void\", radius = 0.418, rings = 1 }",
        )
        .unwrap();
        // 燃料 2、真空 1、水 2 の領域
        let (_, regions) = input.geometry().unwrap();
        assert_eq!(regions[2], VOID);
    }

    #[test]
    fn materials_need_some_nuclide() {
        let result = replaced(
            PIN_CELL,
            "nuclides = { H-1 = 4.9e-2, O-16 = 2.5e-2 }",
            "nuclides = { H-1 = 0.0, O-16 = 0.0 }",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn void_cannot_be_defined() {
        let result = replaced(PIN_CELL, "[material.water]", "[material.void]");
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn self_shielding_is_on_unless_turned_off() {
        assert!(LatticeInput::parse(PIN_CELL).unwrap().self_shielding);
        let off = replaced(
            PIN_CELL,
            "library = \"m1.msgpack\"",
            "library = \"m1.msgpack\"\nself-shielding = false",
        )
        .unwrap();
        assert!(!off.self_shielding);
    }

    #[test]
    fn self_shielding_needs_a_library() {
        let result = replaced(
            LATTICE,
            "cross-sections = \"c5g7.toml\"",
            "cross-sections = \"c5g7.toml\"\nself-shielding = true",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
    }

    #[test]
    fn scattering_order_is_optional() {
        // 省略すると、ライブラリにある次数をすべて使う（決めるのは物質を作った後）
        assert_eq!(
            LatticeInput::parse(PIN_CELL).unwrap().scattering_order,
            None
        );
        let p1 = replaced(
            PIN_CELL,
            "pitch = 1.26",
            "pitch = 1.26\nscattering-order = 1",
        )
        .unwrap();
        assert_eq!(p1.scattering_order, Some(1));
        // 巨視的断面積は P0 だけ
        let result = replaced(
            LATTICE,
            "cross-sections = \"c5g7.toml\"",
            "cross-sections = \"c5g7.toml\"\nscattering-order = 1",
        );
        assert!(matches!(result, Err(InputError::Invalid(_))), "{result:?}");
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
        assert!(matches!(
            geometry(
                r#"zones = [
  { material = "fuel", radius = 0.4, rings = 2 },
  { material = "water", radius = 0.8, rings = 2 },
  { material = "water", rings = 2 },
]"#
            ),
            Err(InputError::LastZoneOutside { .. })
        ));
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
    fn columns_and_rows_may_have_their_own_sizes() {
        // ピンの格子の右に、厚さ 0.5 のバッフルと 1.5 の水の列を足す（反射体）
        let reflector = |rows: &str, extra: &str| {
            let text = LATTICE
                .replace(
                    "rows = [\"pin box\", \"box pin\"]",
                    &format!("widths = [{{ size = 2.0, count = 2 }}, {{ size = 0.5 }}, {{ size = 1.5 }}]\n{extra}rows = {rows}"),
                )
                .replace(
                    "[cell.box]",
                    "[cell.baffle]\nzones = [{ material = \"clad\", rings = 1 }]\n\n[cell.box]",
                );
            LatticeInput::parse(&text).and_then(|input| input.geometry())
        };
        let (lattice, materials) =
            reflector("[\"pin box baffle box\", \"box pin baffle box\"]", "").unwrap();
        assert_eq!(lattice.widths(), &[2.0, 2.0, 0.5, 1.5]);
        assert_eq!(lattice.heights(), &[2.0, 2.0]);
        // ピン 3 + 水 1 + バッフル 1 + 水 1 の行が 2 つ
        assert_eq!(lattice.regions(), 12);
        assert_eq!(materials[4], "clad");
        // 均質なセルは、ピッチのマスではピンセル、ほかのマスでは長方形のセル
        assert!(lattice.cell_type(1).pin().is_some());
        assert!(lattice.cell_type(3).pin().is_none());
        assert_eq!(lattice.pin_rows_and_columns(), (vec![0, 1], vec![0, 1]));
        // ピンはピッチの幅の列にだけ置ける
        assert!(reflector("[\"pin box pin box\", \"box pin baffle box\"]", "").is_err());
        // 列の数がそろわない、余白と同時には使えない
        assert!(reflector("[\"pin box baffle\", \"box pin baffle\"]", "").is_err());
        assert!(
            reflector(
                "[\"pin box baffle box\", \"box pin baffle box\"]",
                "gap = { width = 0.1, material = \"water\" }\n"
            )
            .is_err()
        );
    }

    #[test]
    fn sides_may_have_an_albedo() {
        let input = replaced(
            LATTICE,
            "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
            "boundary = { right = \"vacuum\", top = { albedo = 0.3 } }",
        )
        .unwrap();
        assert_eq!(input.boundaries.top, Boundary::Albedo(0.3));
        assert_eq!(input.boundaries.right, Boundary::Vacuum);
        assert_eq!(input.boundaries.bottom, Boundary::Reflective);
        assert!(
            replaced(
                LATTICE,
                "boundary = { right = \"vacuum\", bottom = \"vacuum\" }",
                "boundary = { top = { albedo = 0.3, extra = 1 } }",
            )
            .is_err()
        );
    }

    #[test]
    fn gap_surrounds_the_lattice() {
        let text = LATTICE.replace(
            "rows = [",
            "gap = { width = 0.25, material = \"water\" }\nrows = [",
        );
        let (lattice, materials) = LatticeInput::parse(&text).unwrap().geometry().unwrap();
        assert_eq!(lattice.width(), 4.5);
        assert_eq!(lattice.widths(), &[0.25, 2.0, 2.0, 0.25]);
        // 余白は周りの薄い行と列のセル（各 1 領域）で、領域はセルの順に並ぶ
        assert_eq!(materials.len(), lattice.regions());
        assert_eq!(materials.len(), 8 + 12);
        assert!(materials[..4].iter().all(|m| m == "water"));
        assert_eq!(materials[4..9], ["water", "UO2", "clad", "water", "water"]);
        assert_eq!(materials.iter().filter(|m| *m == "UO2").count(), 2);
        // 幅が正でない、物質がない、真空（巨視的断面積でも）
        for gap in [
            "{ width = 0.0, material = \"water\" }",
            "{ width = 0.1 }",
            "{ width = 0.1, material = \"void\" }",
        ] {
            let text = LATTICE.replace("rows = [", &format!("gap = {gap}\nrows = ["));
            let result = LatticeInput::parse(&text).and_then(|input| input.geometry().map(|_| ()));
            assert!(result.is_err(), "{gap}");
        }
    }

    #[test]
    fn frames_are_cell_edges() {
        let with_frame = |frame: &str| {
            replaced(
                LATTICE,
                "  { material = \"water\", rings = 1 },\n]",
                &format!("  {{ material = \"water\", rings = 2 }},\n]\nframe = {frame}"),
            )
        };
        let input = with_frame("{ thickness = 0.05, material = \"clad\" }").unwrap();
        let (lattice, regions) = input.geometry().unwrap();
        // ピンセル 4 + 枠 1 が 2 つ、均質なセル 1 が 2 つ
        assert_eq!(lattice.regions(), 2 * 5 + 2);
        assert_eq!(regions[4], "clad");
        // 厚すぎる枠は、帯の配置ではなく枠の誤りとして報告する
        for thickness in ["1.0", "-0.1", "nan"] {
            let frame = format!("{{ thickness = {thickness}, material = \"clad\" }}");
            let error = with_frame(&frame).and_then(|input| input.geometry().map(|_| ()));
            assert!(
                matches!(
                    error,
                    Err(InputError::Cell {
                        source: GeometryError::Frame(_),
                        ..
                    })
                ),
                "{thickness}: {error:?}"
            );
        }
        assert!(with_frame("{ thickness = 0.05 }").is_err());
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
        let text = LATTICE.replace("radius = 0.6", "radius = 1.5");
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
