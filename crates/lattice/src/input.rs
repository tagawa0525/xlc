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
//! cell box water - 1
//! row pin box                            # セルの配置を上の行から
//! row box pin
//! ```
//!
//! セルの帯は外半径の内側を分割数の領域に分ける。外半径のある帯は等体積に、最後の帯は
//! セルの内接円までを半径方向に等間隔に分ける。

use std::path::PathBuf;

use crate::geometry::Lattice;
use crate::tracking::Boundaries;

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

impl LatticeInput {
    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        todo!("{text}")
    }

    /// 形状と、領域ごとの物質の名前。
    pub fn geometry(&self) -> Result<(Lattice, Vec<String>), String> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tracking::Boundary;

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
}
