//! ピンセルの入力ファイル。
//!
//! 1行に `キーワード 値...` を書く。`#` から行末まではコメント。
//!
//! ```text
//! library wims69.lib
//! output pincell.fg
//! pitch 1.26                 # cm
//! fuel-radius 0.4096         # cm
//! fuel-rings 3               # 燃料を等体積のリングに分ける数
//! moderator-rings 1          # 減速材に置く円の数（半径方向に等間隔）
//! azimuthal 32               # (0, π) の方位角の数
//! spacing 0.05               # トラック間隔 (cm)
//! fuel U235 7.0e-4 U238 2.2e-2 O16 4.6e-2      # 核種と個数密度 (1/(barn cm))
//! moderator H1 4.9e-2 O16 2.5e-2
//! few-group-bounds 1e7 0.625 1e-5              # 少数群の境界 (eV)
//! ```

use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq)]
pub struct PinCellInput {
    pub library: PathBuf,
    pub output: PathBuf,
    pub pitch: f64,
    pub fuel_radius: f64,
    pub fuel_rings: usize,
    pub moderator_rings: usize,
    pub azimuthal: usize,
    pub spacing: f64,
    pub fuel: Vec<(String, f64)>,
    pub moderator: Vec<(String, f64)>,
    pub few_group_bounds: Vec<f64>,
}

impl PinCellInput {
    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        todo!()
    }

    /// 領域を区切る円の半径。燃料は等体積のリング、減速材は半径方向に等間隔。
    /// 最初の `fuel_rings` 個の領域が燃料、残りが減速材。
    pub fn radii(&self) -> Vec<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "
# pin cell
library wims69.lib
output pincell.fg
pitch 1.26
fuel-radius 0.4
fuel-rings 2
moderator-rings 1
azimuthal 16
spacing 0.05
fuel U235 7.0e-4 U238 2.2e-2 O16 4.6e-2   # UO2
moderator H1 4.9e-2 O16 2.5e-2
few-group-bounds 1e7 0.625 1e-5
";

    #[test]
    fn parses_sample() {
        let input = PinCellInput::parse(SAMPLE).unwrap();
        assert_eq!(input.library, PathBuf::from("wims69.lib"));
        assert_eq!(input.fuel_rings, 2);
        assert_eq!(input.fuel[1], ("U238".to_string(), 2.2e-2));
        assert_eq!(input.moderator.len(), 2);
        assert_eq!(input.few_group_bounds, vec![1e7, 0.625, 1e-5]);
    }

    #[test]
    fn radii_split_fuel_by_volume_and_moderator_by_radius() {
        let input = PinCellInput::parse(SAMPLE).unwrap();
        let r = input.radii();
        assert_eq!(r.len(), 3);
        assert!((r[0] - 0.4 * 0.5_f64.sqrt()).abs() < 1e-15);
        assert_eq!(r[1], 0.4);
        assert!((r[2] - (0.4 + 0.63) / 2.0).abs() < 1e-15);
    }

    fn without(key: &str) -> String {
        SAMPLE
            .lines()
            .filter(|l| !l.starts_with(key))
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn missing_key_is_error() {
        let err = PinCellInput::parse(&without("pitch")).unwrap_err();
        assert!(err.contains("pitch"), "{err}");
    }

    #[test]
    fn bad_values_report_line() {
        let text = SAMPLE.replace("fuel-rings 2", "fuel-rings two");
        let line = text
            .lines()
            .position(|l| l.starts_with("fuel-rings"))
            .unwrap()
            + 1;
        let err = PinCellInput::parse(&text).unwrap_err();
        assert!(err.starts_with(&format!("line {line}:")), "{err}");
        assert!(PinCellInput::parse(&SAMPLE.replace("fuel U235 7.0e-4", "fuel U235")).is_err());
        assert!(PinCellInput::parse(&format!("{SAMPLE}\nunknown 1")).is_err());
        assert!(PinCellInput::parse(&format!("{SAMPLE}\npitch 2.0")).is_err());
    }
}
