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
//!              temperatures: [{temperature (K), chi: [G 個],
//!                              scatter: [[無限希釈の σ(g→g')、行が入射の群]],
//!                              shielded: [{total, elastic, fission, nu_fission}（σ0 の順）]}]}]}
//! ```

use serde::{Deserialize, Serialize};

use crate::multigroup::GroupStructure;

/// 核種名が約束事（全体計画の「ファイルの約束事」）の形か: `U-235`, `H-1`, `Am-242m`, `Hf-178m2`。
pub fn is_nuclide_name(name: &str) -> bool {
    todo!("{name}")
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

/// ライブラリの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum LibraryError {
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
    pub fn to_msgpack(&self) -> Vec<u8> {
        let file = LibraryFile {
            format: FORMAT.into(),
            group_bounds: self.groups.bounds().to_vec(),
            sigma0: self.sigma0.clone(),
            nuclides: self.nuclides.clone(),
        };
        rmp_serde::to_vec_named(&file).expect("the library has only numbers and strings")
    }

    /// MessagePack から読み、群境界・σ0 の順序と配列の長さを確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, LibraryError> {
        let mut rest = bytes;
        let value =
            rmpv::decode::read_value(&mut rest).map_err(|e| LibraryError::Decode(e.to_string()))?;
        if !rest.is_empty() {
            return Err(LibraryError::TrailingBytes(rest.len()));
        }
        // serde は構造体を配列からも読むので、フィールド名なしで書かれた構造体をここで拒否する
        require_maps(&value)?;
        let file: LibraryFile =
            rmpv::ext::from_value(value).map_err(|e| LibraryError::Decode(e.to_string()))?;
        if file.format != FORMAT {
            return Err(LibraryError::Format(file.format));
        }
        let invalid = LibraryError::Invalid;
        let descending = |v: &[f64]| {
            v.windows(2).all(|w| w[0] > w[1]) && v.iter().all(|&x| x > 0.0 && x.is_finite())
        };
        let bounds = file.group_bounds;
        if bounds.len() < 2 || !descending(&bounds) {
            return Err(invalid(format!(
                "group bounds must be at least two positive descending values, got {bounds:?}"
            )));
        }
        if file.sigma0.is_empty() || !descending(&file.sigma0) {
            return Err(invalid(format!(
                "sigma0 must be positive and descending, got {:?}",
                file.sigma0
            )));
        }
        let ng = bounds.len() - 1;
        for (k, n) in file.nuclides.iter().enumerate() {
            if file.nuclides[..k].iter().any(|m| m.name == n.name) {
                return Err(invalid(format!("nuclide {} appears twice", n.name)));
            }
            n.check(ng, file.sigma0.len()).map_err(invalid)?;
        }
        Ok(Self {
            groups: GroupStructure::new(bounds),
            sigma0: file.sigma0,
            nuclides: file.nuclides,
        })
    }

    /// 核種 `name` の、温度 `temperature` (K) での無限希釈の群定数。温度は表の値と完全に一致すること（補間しない）。
    pub fn nuclide(&self, name: &str, temperature: f64) -> Result<Nuclide, String> {
        let data = self
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .ok_or_else(|| format!("nuclide {name} is not in the library"))?;
        let t = data
            .temperatures
            .iter()
            .find(|t| t.temperature == temperature)
            .ok_or_else(|| {
                let available: Vec<f64> = data.temperatures.iter().map(|t| t.temperature).collect();
                format!("{name} has no data at {temperature} K (available: {available:?})")
            })?;
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

/// 構造体の位置（最上位、核種、温度、遮蔽断面積）がすべて map であることを確かめる。
/// 項目の有無と型は、そのあとの serde の変換で確かめる。
fn require_maps(value: &rmpv::Value) -> Result<(), LibraryError> {
    type Entries = [(rmpv::Value, rmpv::Value)];
    // キーは文字列だけを認める（serde はバイナリのキーもフィールド名として受け付けるため）
    fn map(v: &rmpv::Value) -> Result<&Entries, LibraryError> {
        v.as_map()
            .filter(|entries| entries.iter().all(|(k, _)| k.is_str()))
            .map(Vec::as_slice)
            .ok_or(LibraryError::NotMap)
    }
    fn array<'a>(entries: &'a Entries, key: &str) -> &'a [rmpv::Value] {
        entries
            .iter()
            .find(|(k, _)| k.as_str() == Some(key))
            .and_then(|(_, v)| v.as_array())
            .map_or(&[], Vec::as_slice)
    }
    for nuclide in array(map(value)?, "nuclides") {
        for temperature in array(map(nuclide)?, "temperatures") {
            for shielded in array(map(temperature)?, "shielded") {
                map(shielded)?;
            }
        }
    }
    Ok(())
}

impl NuclideData {
    /// 配列の長さ（群の数 `ng`、σ0 の数 `n_sigma0`）と温度の表があることを確かめる。
    fn check(&self, ng: usize, n_sigma0: usize) -> Result<(), String> {
        let name = &self.name;
        if name.is_empty() {
            return Err("a nuclide has no name".into());
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
            temperatures: temperatures
                .iter()
                .map(|&t| temperature(t, fissile))
                .collect(),
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]),
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
            Library::from_msgpack(&library.to_msgpack()).unwrap(),
            library
        );
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
    fn written_as_a_map_with_the_format_name() {
        let bytes = sample().to_msgpack();
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
        let mut bytes = sample().to_msgpack();
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
            "", "U235", "u-235", "U-", "U-235M", "Am-242m0", "U-0235", "U 235", "Uuo-294", "H-1x",
            "-1",
        ] {
            assert!(!is_nuclide_name(name), "{name}");
        }
    }

    #[test]
    fn library_nuclide_names_are_checked() {
        assert!(read_modified(|f| f.nuclides[1].name = "H1".into()).is_err());
    }
}
