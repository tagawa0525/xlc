//! 多群ライブラリ（`xs` から `lattice` に渡すファイル）。
//!
//! フィールド名付きの map 形式の MessagePack（全体計画の「ファイルの約束事」）。f64 は 8 バイトのまま持つので、
//! 書いて読んだ値は最下位のビットまで一致する。核種ごとに温度の表を持ち、温度ごとに σ0（背景断面積）ごとの
//! 群定数を持つ。σ0 は全核種で共通で、降順に並べ、最初の値を無限希釈とみなす。
//!
//! ```text
//! {format: "xlc-multigroup",
//!  group_bounds: [群境界 (eV)、降順],
//!  sigma0: [σ0 (barn)、降順],
//!  nuclides: [{name: "U-235", za, awr,
//!              potential: ポテンシャル散乱 (barn)、分離共鳴を持つ核種だけ（なければ nil）,
//!              temperatures: [{temperature (K), chi: [G 個],
//!                              scatter: [[無限希釈の σ(g→g')、行が入射の群]],
//!                              shielded: [{total, elastic, fission, nu_fission}（σ0 の順）]}]}]}
//! ```

use serde::{Deserialize, Serialize};

use crate::msgpack::{self, MsgpackError};
use crate::multigroup::GroupStructure;

/// 核種名が約束事（全体計画の「ファイルの約束事」）の形か: `U-235`, `H-1`, `Am-242m`, `Hf-178m2`。
pub fn is_nuclide_name(name: &str) -> bool {
    // 数字の並び（先頭は 0 でない）
    let number =
        |s: &str| !s.is_empty() && !s.starts_with('0') && s.bytes().all(|b| b.is_ascii_digit());
    let Some((symbol, rest)) = name.split_once('-') else {
        return false;
    };
    let mut chars = symbol.chars();
    let element = matches!(
        (chars.next(), chars.next(), chars.next()),
        (Some(a), None, None) | (Some(a), Some('a'..='z'), None) if a.is_ascii_uppercase()
    );
    let (mass, isomer) = match rest.split_once('m') {
        // 最初の準安定状態は `m` だけで書く（`m1` を許すと同じ状態に名前が2つできる）
        Some((mass, index)) => (mass, index.is_empty() || (number(index) && index != "1")),
        None => (rest, true),
    };
    element && number(mass) && isomer
}

/// 1つの σ0 での群定数 (barn)。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Shielded {
    pub total: Vec<f64>,
    pub elastic: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
}

/// 1つの温度での群定数。`shielded` は `Library::sigma0` の順。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemperatureData {
    /// 処理温度 (K)
    pub temperature: f64,
    pub chi: Vec<f64>,
    /// 無限希釈の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    pub shielded: Vec<Shielded>,
}

/// 1核種の、温度ごとの群定数。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NuclideData {
    pub name: String,
    pub za: f64,
    pub awr: f64,
    /// ポテンシャル散乱断面積 (barn)。分離共鳴を持つ核種だけ（ESSM で使う）
    pub potential: Option<f64>,
    pub temperatures: Vec<TemperatureData>,
}

/// 1核種の、1つの温度での無限希釈の群定数 (barn)。核分裂のない核種では fission, nu-fission, chi は 0。
#[derive(Debug, Clone, PartialEq)]
pub struct Nuclide {
    pub name: String,
    pub za: f64,
    pub awr: f64,
    /// 処理温度 (K)
    pub temperature: f64,
    pub total: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub chi: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Library {
    pub groups: GroupStructure,
    /// 背景断面積 (barn)、降順。最初の値を無限希釈とみなす
    pub sigma0: Vec<f64>,
    pub nuclides: Vec<NuclideData>,
}

/// 形式の名前
const FORMAT: &str = "xlc-multigroup";

/// ライブラリから核種を引くときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum LookupError {
    #[error("nuclide {0} is not in the library")]
    Nuclide(String),
    #[error("{name} has no data at {temperature} K (available: {available:?})")]
    Temperature {
        name: String,
        temperature: f64,
        available: Vec<f64>,
    },
    #[error("{name}: sigma0 {sigma0} b in group {group} is outside the table {min}..{max} b")]
    Sigma0 {
        name: String,
        group: usize,
        sigma0: f64,
        min: f64,
        max: f64,
    },
    #[error("the library has no sigma0 values")]
    NoSigma0,
    #[error("{name}: {given} sigma0 values for {groups} groups")]
    Sigma0Groups {
        name: String,
        given: usize,
        groups: usize,
    },
}

/// ライブラリの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum LibraryError {
    #[error("cannot encode as MessagePack: {0}")]
    Encode(String),
    #[error("invalid MessagePack: {0}")]
    Decode(String),
    #[error("{0} bytes follow the library")]
    TrailingBytes(usize),
    #[error("structs must be MessagePack maps with string field names")]
    NotMap,
    #[error("expected format {FORMAT}, found {0}")]
    Format(String),
    #[error("{0}")]
    Invalid(String),
}

/// ファイルの中身（フィールド名付きの map 形式で書く）。
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LibraryFile {
    format: String,
    group_bounds: Vec<f64>,
    sigma0: Vec<f64>,
    nuclides: Vec<NuclideData>,
}

impl Library {
    /// MessagePack（フィールド名付きの map 形式）で書く。
    /// MessagePack で書く。読むときと同じ検査（配列の長さ、有限な値）を書く前にもして、壊れたライブラリを書かない。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, LibraryError> {
        let ng = self.groups.len();
        for n in &self.nuclides {
            n.check(ng, self.sigma0.len())
                .map_err(LibraryError::Invalid)?;
        }
        let file = LibraryFile {
            format: FORMAT.into(),
            group_bounds: self.groups.bounds().to_vec(),
            sigma0: self.sigma0.clone(),
            nuclides: self.nuclides.clone(),
        };
        Ok(msgpack::to_named(&file)?)
    }

    /// MessagePack から読み、群境界・σ0 の順序と配列の長さを確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, LibraryError> {
        let value = msgpack::read_value(bytes)?;
        require_maps(&value)?;
        let file: LibraryFile = msgpack::from_value(value)?;
        if file.format != FORMAT {
            return Err(LibraryError::Format(file.format));
        }
        let invalid = LibraryError::Invalid;
        let descending = |v: &[f64]| {
            v.windows(2).all(|w| w[0] > w[1]) && v.iter().all(|&x| x > 0.0 && x.is_finite())
        };
        let groups = GroupStructure::new(file.group_bounds).map_err(|e| invalid(e.to_string()))?;
        if file.sigma0.is_empty() || !descending(&file.sigma0) {
            return Err(invalid(format!(
                "sigma0 must be positive and descending, got {:?}",
                file.sigma0
            )));
        }
        let ng = groups.len();
        for (k, n) in file.nuclides.iter().enumerate() {
            if file.nuclides[..k].iter().any(|m| m.name == n.name) {
                return Err(invalid(format!("nuclide {} appears twice", n.name)));
            }
            n.check(ng, file.sigma0.len()).map_err(invalid)?;
        }
        Ok(Self {
            groups,
            sigma0: file.sigma0,
            nuclides: file.nuclides,
        })
    }

    /// 核種 `name` の、温度 `temperature` (K) での無限希釈の群定数。温度は表の値と完全に一致すること（補間しない）。
    pub fn nuclide(&self, name: &str, temperature: f64) -> Result<Nuclide, LookupError> {
        let (data, t) = self.lookup(name, temperature)?;
        let dilute = &t.shielded[0];
        Ok(Nuclide {
            name: data.name.clone(),
            za: data.za,
            awr: data.awr,
            temperature: t.temperature,
            total: dilute.total.clone(),
            fission: dilute.fission.clone(),
            nu_fission: dilute.nu_fission.clone(),
            chi: t.chi.clone(),
            scatter: t.scatter.clone(),
        })
    }
}

impl Library {
    /// 核種 `name` の温度 `temperature` (K) の表（温度は表の値と完全に一致すること）。
    fn lookup(
        &self,
        name: &str,
        temperature: f64,
    ) -> Result<(&NuclideData, &TemperatureData), LookupError> {
        let data = self
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .ok_or_else(|| LookupError::Nuclide(name.into()))?;
        let t = data
            .temperatures
            .iter()
            .find(|t| t.temperature == temperature)
            .ok_or_else(|| LookupError::Temperature {
                name: name.into(),
                temperature,
                available: data.temperatures.iter().map(|t| t.temperature).collect(),
            })?;
        Ok((data, t))
    }

    /// 核種 `name` の、群ごとの背景断面積 `sigma0` (barn) での遮蔽断面積。表の σ0 の間は ln σ0 に線形に
    /// 補間し、格子点では表の値をそのまま使う。表の範囲の外はエラーにする。
    pub fn shielded(
        &self,
        name: &str,
        temperature: f64,
        sigma0: &[f64],
    ) -> Result<Shielded, LookupError> {
        let (_, t) = self.lookup(name, temperature)?;
        let groups = self.groups.len();
        if sigma0.len() != groups {
            return Err(LookupError::Sigma0Groups {
                name: name.into(),
                given: sigma0.len(),
                groups,
            });
        }
        let (Some(&max), Some(&min)) = (self.sigma0.first(), self.sigma0.last()) else {
            return Err(LookupError::NoSigma0);
        };
        let mut result = Shielded {
            total: Vec::with_capacity(groups),
            elastic: Vec::with_capacity(groups),
            fission: Vec::with_capacity(groups),
            nu_fission: Vec::with_capacity(groups),
        };
        for (g, &s) in sigma0.iter().enumerate() {
            if !(min..=max).contains(&s) {
                return Err(LookupError::Sigma0 {
                    name: name.into(),
                    group: g,
                    sigma0: s,
                    min,
                    max,
                });
            }
            // σ0 は降順。k は s 以下の最初の格子点
            let k = self.sigma0.partition_point(|&x| x > s);
            let value = |of: fn(&Shielded) -> &Vec<f64>| -> f64 {
                let b = of(&t.shielded[k])[g];
                if self.sigma0[k] == s {
                    return b;
                }
                let a = of(&t.shielded[k - 1])[g];
                let u = (s / self.sigma0[k - 1]).ln() / (self.sigma0[k] / self.sigma0[k - 1]).ln();
                a + u * (b - a)
            };
            result.total.push(value(|x| &x.total));
            result.elastic.push(value(|x| &x.elastic));
            result.fission.push(value(|x| &x.fission));
            result.nu_fission.push(value(|x| &x.nu_fission));
        }
        Ok(result)
    }
}

/// 構造体の位置（最上位、核種、温度、遮蔽断面積）がすべて map であることを確かめる。
/// 項目の有無と型は、そのあとの serde の変換で確かめる。
fn require_maps(value: &rmpv::Value) -> Result<(), MsgpackError> {
    use msgpack::{array, map};
    for nuclide in array(map(value)?, "nuclides") {
        for temperature in array(map(nuclide)?, "temperatures") {
            for shielded in array(map(temperature)?, "shielded") {
                map(shielded)?;
            }
        }
    }
    Ok(())
}

impl From<MsgpackError> for LibraryError {
    fn from(error: MsgpackError) -> Self {
        match error {
            MsgpackError::Encode(message) => Self::Encode(message),
            MsgpackError::Decode(message) => Self::Decode(message),
            MsgpackError::TrailingBytes(n) => Self::TrailingBytes(n),
            MsgpackError::NotMap => Self::NotMap,
        }
    }
}

impl NuclideData {
    /// 配列の長さ（群の数 `ng`、σ0 の数 `n_sigma0`）と温度の表があることを確かめる。
    fn check(&self, ng: usize, n_sigma0: usize) -> Result<(), String> {
        let name = &self.name;
        if !is_nuclide_name(name) {
            return Err(format!(
                "nuclide name `{name}` is not of the form U-235, Am-242m"
            ));
        }
        if self.temperatures.is_empty() {
            return Err(format!("{name}: no temperatures"));
        }
        if !(self.za.is_finite() && self.awr.is_finite() && self.awr > 0.0) {
            return Err(format!(
                "{name}: invalid za {} or awr {}",
                self.za, self.awr
            ));
        }
        for (k, t) in self.temperatures.iter().enumerate() {
            let at = |what: String| format!("{name} at {} K: {what}", t.temperature);
            if !(t.temperature > 0.0 && t.temperature.is_finite()) {
                return Err(at("the temperature must be positive and finite".into()));
            }
            if self.temperatures[..k]
                .iter()
                .any(|u| u.temperature == t.temperature)
            {
                return Err(at("the temperature appears twice".into()));
            }
            let sized = |v: &[f64]| v.len() == ng;
            if !sized(&t.chi) {
                return Err(at(format!("chi must have {ng} groups")));
            }
            if t.scatter.len() != ng || !t.scatter.iter().all(|row| sized(row)) {
                return Err(at(format!("scatter must be {ng} x {ng}")));
            }
            if t.shielded.len() != n_sigma0 {
                return Err(at(format!(
                    "{} shielded sets for {n_sigma0} sigma0 values",
                    t.shielded.len()
                )));
            }
            let all_sized = |s: &Shielded| {
                [&s.total, &s.elastic, &s.fission, &s.nu_fission]
                    .iter()
                    .all(|v| sized(v))
            };
            if !t.shielded.iter().all(all_sized) {
                return Err(at(format!("shielded cross sections must have {ng} groups")));
            }
            let values = t.chi.iter().chain(t.scatter.iter().flatten()).chain(
                t.shielded
                    .iter()
                    .flat_map(|s| [&s.total, &s.elastic, &s.fission, &s.nu_fission])
                    .flatten(),
            );
            if let Some(bad) = values.copied().find(|v| !v.is_finite()) {
                return Err(at(format!("non-finite value {bad}")));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn shielded(scale: f64, fissile: bool) -> Shielded {
        Shielded {
            // 1.7929089665572825e-9 は、10 進表記を経ると最下位のビットがずれやすい値
            total: vec![11.5 * scale, 1.7929089665572825e-9],
            elastic: vec![4.0 * scale, 15.0],
            fission: if fissile {
                vec![1.2 * scale, 580.5]
            } else {
                vec![0.0, 0.0]
            },
            nu_fission: if fissile {
                vec![3.0 * scale, 1414.1]
            } else {
                vec![0.0, 0.0]
            },
        }
    }

    fn temperature(t: f64, fissile: bool) -> TemperatureData {
        TemperatureData {
            temperature: t,
            chi: if fissile {
                vec![1.0, 0.0]
            } else {
                vec![0.0, 0.0]
            },
            scatter: vec![vec![4.0, 0.01], vec![1.0e-7, 15.0]],
            shielded: vec![shielded(1.0, fissile), shielded(0.5, fissile)],
        }
    }

    pub(crate) fn sample() -> Library {
        let nuclide = |name: &str, fissile: bool, temperatures: &[f64]| NuclideData {
            name: name.into(),
            za: 92235.0,
            awr: 233.0248,
            potential: None,
            temperatures: temperatures
                .iter()
                .map(|&t| temperature(t, fissile))
                .collect(),
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
            sigma0: vec![1.0e10, 10.0],
            nuclides: vec![
                nuclide("U-235", true, &[293.6, 900.0]),
                nuclide("H-1", false, &[293.6]),
            ],
        }
    }

    #[test]
    fn round_trips_bit_exactly() {
        let library = sample();
        assert_eq!(
            Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap(),
            library
        );
    }

    #[test]
    fn keeps_the_potential_scattering() {
        let mut library = sample();
        library.nuclides[0].potential = Some(11.29);
        let read = Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap();
        assert_eq!(read.nuclides[0].potential, Some(11.29));
        assert_eq!(read.nuclides[1].potential, None);
    }

    #[test]
    fn gives_infinite_dilution_at_a_tabulated_temperature() {
        let library = sample();
        let u235 = library.nuclide("U-235", 900.0).unwrap();
        assert_eq!((u235.name.as_str(), u235.temperature), ("U-235", 900.0));
        assert_eq!(u235.total, vec![11.5, 1.7929089665572825e-9]);
        assert_eq!(u235.chi, vec![1.0, 0.0]);
        assert!(library.nuclide("U-235", 600.0).is_err());
        assert!(library.nuclide("O-16", 293.6).is_err());
    }

    #[test]
    fn interpolates_shielded_cross_sections_in_ln_sigma0() {
        let library = sample();
        // 格子点では表の値
        let at_points = library.shielded("U-235", 900.0, &[1.0e10, 10.0]).unwrap();
        assert_eq!(at_points.total, vec![11.5, 1.7929089665572825e-9]);
        assert_eq!(at_points.fission, vec![1.2, 580.5]);
        // ln σ0 の中点では、両端の値の平均
        let middle = (1.0e10_f64 * 10.0).sqrt();
        let between = library.shielded("U-235", 900.0, &[middle, middle]).unwrap();
        assert!((between.total[0] - 0.75 * 11.5).abs() < 1e-12);
        assert!((between.elastic[0] - 0.75 * 4.0).abs() < 1e-12);
        assert!((between.nu_fission[0] - 0.75 * 3.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_sigma0_outside_the_table() {
        let library = sample();
        for sigma0 in [5.0, 1.0e11, f64::NAN] {
            let error = library.shielded("U-235", 900.0, &[sigma0, 100.0]);
            assert!(
                matches!(error, Err(LookupError::Sigma0 { group: 0, .. })),
                "{sigma0}: {error:?}"
            );
        }
        // σ0 は群の数だけ要る
        assert!(matches!(
            library.shielded("U-235", 900.0, &[100.0]),
            Err(LookupError::Sigma0Groups { .. })
        ));
    }

    fn sample_file() -> LibraryFile {
        let library = sample();
        LibraryFile {
            format: "xlc-multigroup".into(),
            group_bounds: library.groups.bounds().to_vec(),
            sigma0: library.sigma0,
            nuclides: library.nuclides,
        }
    }

    /// サンプルを書き換え、フィールド名付きの map 形式で書いて読む。
    fn read_modified(change: impl FnOnce(&mut LibraryFile)) -> Result<Library, LibraryError> {
        let mut file = sample_file();
        change(&mut file);
        Library::from_msgpack(&rmp_serde::to_vec_named(&file).unwrap())
    }

    #[test]
    fn refuses_to_write_non_finite_values() {
        let mut library = sample();
        library.nuclides[0].temperatures[0].shielded[1].total[0] = f64::NAN;
        assert!(matches!(
            library.to_msgpack(),
            Err(LibraryError::Invalid(_))
        ));
    }

    #[test]
    fn written_as_a_map_with_the_format_name() {
        let bytes = sample().to_msgpack().unwrap();
        let value: std::collections::BTreeMap<String, serde::de::IgnoredAny> =
            rmp_serde::from_slice(&bytes).unwrap();
        let keys: Vec<&str> = value.keys().map(String::as_str).collect();
        assert_eq!(keys, vec!["format", "group_bounds", "nuclides", "sigma0"]);
        #[derive(Deserialize)]
        struct Head {
            format: String,
        }
        let head: Head = rmp_serde::from_slice(&bytes).unwrap();
        assert_eq!(head.format, "xlc-multigroup");
    }

    #[test]
    fn array_form_is_rejected() {
        // 既定の配列形式ではフィールド名がなく、順番の変わったファイルを黙って誤読しうる
        let bytes = rmp_serde::to_vec(&sample_file()).unwrap();
        assert!(matches!(
            Library::from_msgpack(&bytes),
            Err(LibraryError::NotMap)
        ));
    }

    #[test]
    fn unknown_and_missing_fields_are_rejected() {
        #[derive(Serialize)]
        struct Extra {
            format: &'static str,
            group_bounds: Vec<f64>,
            sigma0: Vec<f64>,
            nuclides: Vec<NuclideData>,
            comment: &'static str,
        }
        #[derive(Serialize)]
        struct Missing {
            format: &'static str,
            group_bounds: Vec<f64>,
            nuclides: Vec<NuclideData>,
        }
        let file = sample_file();
        let extra = Extra {
            format: "xlc-multigroup",
            group_bounds: file.group_bounds.clone(),
            sigma0: file.sigma0.clone(),
            nuclides: file.nuclides.clone(),
            comment: "x",
        };
        let missing = Missing {
            format: "xlc-multigroup",
            group_bounds: file.group_bounds.clone(),
            nuclides: file.nuclides.clone(),
        };
        for bytes in [
            rmp_serde::to_vec_named(&extra).unwrap(),
            rmp_serde::to_vec_named(&missing).unwrap(),
        ] {
            assert!(matches!(
                Library::from_msgpack(&bytes),
                Err(LibraryError::Decode(_))
            ));
        }
    }

    #[test]
    fn other_formats_are_rejected() {
        let err = read_modified(|f| f.format = "xlc-fewgroup".into()).unwrap_err();
        assert!(matches!(err, LibraryError::Format(_)), "{err}");
    }

    #[test]
    fn wrong_number_of_values_is_error() {
        let err = read_modified(|f| {
            f.nuclides[0].temperatures[0].shielded[0].total.pop();
        })
        .unwrap_err();
        assert!(err.to_string().contains("U-235"), "{err}");
        assert!(
            read_modified(|f| f.nuclides[1].temperatures[0].scatter[1]
                .pop()
                .map(|_| ())
                .unwrap())
            .is_err()
        );
        assert!(read_modified(|f| f.nuclides[1].temperatures[0].chi.clear()).is_err());
    }

    #[test]
    fn shielded_sets_must_match_sigma0() {
        assert!(
            read_modified(|f| {
                f.nuclides[0].temperatures[0].shielded.pop();
            })
            .is_err()
        );
    }

    #[test]
    fn bounds_and_sigma0_must_be_descending() {
        assert!(read_modified(|f| f.sigma0 = vec![10.0, 1.0e10]).is_err());
        assert!(read_modified(|f| f.group_bounds = vec![1.0e-5, 0.625, 1.0e7]).is_err());
        assert!(read_modified(|f| f.group_bounds = vec![1.0e7]).is_err());
    }

    #[test]
    fn nuclides_need_temperatures() {
        assert!(read_modified(|f| f.nuclides[1].temperatures.clear()).is_err());
    }

    #[test]
    fn broken_bytes_are_rejected() {
        assert!(matches!(
            Library::from_msgpack(&[0x81, 0xc1]),
            Err(LibraryError::Decode(_))
        ));
        assert!(Library::from_msgpack(&[]).is_err());
    }

    #[test]
    fn nested_structs_in_array_form_are_rejected() {
        // 最上位だけ map で、核種を配列形式で書いたファイル
        #[derive(Serialize)]
        struct Mixed {
            format: &'static str,
            group_bounds: Vec<f64>,
            sigma0: Vec<f64>,
            nuclides: embedded::Nested,
        }
        let file = sample_file();
        let nested = rmp_serde::to_vec(&file.nuclides).unwrap();
        let mixed = Mixed {
            format: "xlc-multigroup",
            group_bounds: file.group_bounds,
            sigma0: file.sigma0,
            nuclides: embedded::Nested(nested),
        };
        let bytes = rmp_serde::to_vec_named(&mixed).unwrap();
        assert!(matches!(
            Library::from_msgpack(&bytes),
            Err(LibraryError::NotMap)
        ));
    }

    /// 書き出し済みの MessagePack をそのまま埋め込むための型。
    mod embedded {
        pub struct Nested(pub Vec<u8>);

        impl serde::Serialize for Nested {
            fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                let value = rmpv::decode::read_value(&mut self.0.as_slice()).unwrap();
                value.serialize(s)
            }
        }
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let mut bytes = sample().to_msgpack().unwrap();
        bytes.extend([0xc1, 0xff, 0x00]);
        assert!(Library::from_msgpack(&bytes).is_err());
    }

    #[test]
    fn non_finite_values_are_rejected() {
        assert!(read_modified(|f| f.nuclides[0].temperatures[0].temperature = f64::NAN).is_err());
        assert!(
            read_modified(|f| f.nuclides[0].temperatures[1].shielded[1].total[0] = f64::NAN)
                .is_err()
        );
        assert!(
            read_modified(|f| f.nuclides[1].temperatures[0].scatter[0][1] = f64::INFINITY).is_err()
        );
        assert!(read_modified(|f| f.nuclides[1].awr = f64::NAN).is_err());
        assert!(read_modified(|f| f.group_bounds[0] = f64::INFINITY).is_err());
        assert!(read_modified(|f| f.sigma0[0] = f64::INFINITY).is_err());
    }

    #[test]
    fn duplicate_nuclides_and_temperatures_are_rejected() {
        assert!(read_modified(|f| f.nuclides[1].name = "U-235".into()).is_err());
        assert!(read_modified(|f| f.nuclides[0].temperatures[1].temperature = 293.6).is_err());
    }

    #[test]
    fn map_keys_must_be_strings() {
        // serde はバイナリのキーもフィールド名として受け付けるので、構造の検査をすり抜けさせない
        let file = sample_file();
        let nuclides = rmpv::ext::to_value(&file.nuclides).unwrap();
        let value = rmpv::Value::Map(vec![
            ("format".into(), "xlc-multigroup".into()),
            (
                "group_bounds".into(),
                rmpv::ext::to_value(&file.group_bounds).unwrap(),
            ),
            ("sigma0".into(), rmpv::ext::to_value(&file.sigma0).unwrap()),
            (rmpv::Value::Binary(b"nuclides".to_vec()), nuclides),
        ]);
        let mut bytes = Vec::new();
        rmpv::encode::write_value(&mut bytes, &value).unwrap();
        assert!(matches!(
            Library::from_msgpack(&bytes),
            Err(LibraryError::NotMap)
        ));
    }

    #[test]
    fn temperatures_are_looked_up_exactly() {
        // 重複の判定と同じく完全一致で引く。近い温度を黙って受け入れない
        let library = sample();
        assert!(library.nuclide("U-235", 900.0 * (1.0 + 1.0e-12)).is_err());
    }

    #[test]
    fn nuclide_names_follow_the_convention() {
        for name in ["U-235", "H-1", "Am-242m", "Hf-178m2", "Zr-90", "B-10"] {
            assert!(is_nuclide_name(name), "{name}");
        }
        for name in [
            "", "U235", "u-235", "U-", "U-235M", "Am-242m0", "Am-242m1", "U-0235", "U 235",
            "Uuo-294", "H-1x", "-1",
        ] {
            assert!(!is_nuclide_name(name), "{name}");
        }
    }

    #[test]
    fn library_nuclide_names_are_checked() {
        assert!(read_modified(|f| f.nuclides[1].name = "H1".into()).is_err());
    }
}
