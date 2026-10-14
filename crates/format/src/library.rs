//! 多群ライブラリ（`xs` から `lattice` に渡すファイル）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）。f64 は 8 バイトのまま持つので、
//! 書いて読んだ値は最下位のビットまで一致する。核種ごとに温度の表を持ち、温度ごとに σ0（背景断面積）ごとの
//! 群定数を持つ。σ0 は全核種で共通で、降順に並べ、最初の値を無限希釈とみなす。
//!
//! ```text
//! {format: "xlc-multigroup",
//!  group_bounds: [群境界 (eV)、降順],
//!  sigma0: [σ0 (barn)、降順],
//!  nuclides: [{name: "U-235", za, awr,
//!              resolved: {top: 分離共鳴の上端 (eV), potential: ポテンシャル散乱 (barn)}、分離共鳴を持つ核種だけ（なければ nil）,
//!              temperatures: [{temperature (K), chi: [G 個],
//!                              scatter: [[[無限希釈の σ_l(g→g')、行が入射の群]]（ルジャンドルの次数 l = 0, 1, … の順）],
//!                              shielded: [{total, elastic, fission, nu_fission, capture}（σ0 の順）],
//!                              pointwise: {energy, total, elastic, fission, nu_fission, capture},
//!                              reactions: [{mt: 鎖の反応の代表の MT, channels: [{final_state: nil か {mf, lfs（nil は残り）},
//!                                           cross_section: [無限希釈の G 個]}]}]}]}]}
//! ```

use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::group::GroupStructure;
use crate::msgpack::{self, MsgpackError};

/// 核種名が約束事（ADR-0015）の形か: `U-235`, `H-1`, `Am-242m`, `Hf-178m2`。
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
    /// 捕獲（MT102）
    pub capture: Vec<f64>,
}

/// 燃焼鎖の反応（ADR-0164）の、代表の MT と、それにまとめる MF3 の MT の範囲（両端を含む）。代表の MT が MF3 になければ、
/// 範囲の部分の MT の和を反応の断面積にする（ADR-0172）。鎖（xlc-chain）と xs が同じ規則を使う。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainReaction {
    pub mt: u32,
    pub mts: &'static [(u32, u32)],
}

impl ChainReaction {
    /// MF3 の `mt` がこの反応のものか。
    pub fn includes(&self, mt: u32) -> bool {
        self.mts
            .iter()
            .any(|&(low, high)| (low..=high).contains(&mt))
    }
}

/// 燃焼鎖の反応（核分裂を除く）。OpenMC の `Chain.from_endf` の既定と同じ反応（ADR-0164）
pub const CHAIN_REACTIONS: [ChainReaction; 6] = [
    ChainReaction {
        mt: 16,
        mts: &[(16, 16), (875, 891)],
    },
    ChainReaction {
        mt: 17,
        mts: &[(17, 17)],
    },
    ChainReaction {
        mt: 37,
        mts: &[(37, 37)],
    },
    ChainReaction {
        mt: 102,
        mts: &[(102, 102)],
    },
    ChainReaction {
        mt: 103,
        mts: &[(103, 103), (600, 649)],
    },
    ChainReaction {
        mt: 107,
        mts: &[(107, 107), (800, 849)],
    },
];

/// 燃焼鎖の 1 つの反応（ADR-0164 の反応。代表の MT が 102、16、17、37、103、107）の、生成物の状態ごとの経路の無限希釈の群の
/// 断面積（ADR-0172）。経路の並びと `final_state` は、鎖の反応の経路と同じ。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GroupReaction {
    pub mt: u32,
    pub channels: Vec<GroupChannel>,
}

/// 反応の 1 つの経路。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GroupChannel {
    /// 状態ごとの割合を与えるファイルと状態。状態に分けない反応は nil（項目は省けない）
    #[serde(deserialize_with = "nil_or")]
    pub final_state: Option<ChannelState>,
    /// 無限希釈の群の断面積 (barn)
    pub cross_section: Vec<f64>,
}

/// 経路の状態（MF8 の LMF と LFS）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChannelState {
    /// 割合を与える ENDF のファイル（MF9 か MF10）
    pub mf: u32,
    /// 状態（LFS）。nil は、MF8 が基底状態を挙げない反応の、挙げた状態を除いた残り（基底状態へ）
    #[serde(deserialize_with = "nil_or")]
    pub lfs: Option<u32>,
}

fn nil_or<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    d: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(d)
}

/// 1つの温度での、ドップラー拡がりのポイントワイズの断面積 (barn)。`Library::pointwise_range` の範囲で、エネルギーの昇順
/// （範囲がなければ空）
/// （不連続点は同じエネルギーを2回持ち、左と右の極限の値を並べる）。点の間は lin-lin。
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pointwise {
    /// エネルギー (eV)
    pub energy: Vec<f64>,
    pub total: Vec<f64>,
    pub elastic: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub capture: Vec<f64>,
}

/// 1つの温度での群定数。`shielded` は `Library::sigma0` の順。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemperatureData {
    /// 処理温度 (K)
    pub temperature: f64,
    pub chi: Vec<f64>,
    /// 無限希釈の散乱行列のルジャンドルのモーメント [次数 l][入射の群][出射の群]（少なくとも P0）
    pub scatter: Vec<Vec<Vec<f64>>>,
    pub shielded: Vec<Shielded>,
    /// 超微細群の減速計算に使うポイントワイズの断面積
    pub pointwise: Pointwise,
    /// 燃焼鎖の反応の群の断面積（代表の MT の昇順）
    pub reactions: Vec<GroupReaction>,
}

/// 分離共鳴の範囲（分離共鳴を持つ核種だけが持つ）。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Resolved {
    /// 分離共鳴の上端 (eV)。セルの超微細群の減速計算の範囲を決めるのに使う
    pub top: f64,
    /// ポテンシャル散乱断面積 4π a² (barn)（ESSM で使う）
    pub potential: f64,
}

fn nil_or_resolved<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<Resolved>, D::Error> {
    Option::deserialize(d)
}

/// 1核種の、温度ごとの群定数。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NuclideData {
    pub name: String,
    pub za: f64,
    pub awr: f64,
    /// 分離共鳴を持たない核種は nil を書く。nil でも項目は省けない（`deserialize_with` を付けると、serde は
    /// 省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "nil_or_resolved")]
    pub resolved: Option<Resolved>,
    /// 熱中性子散乱を S(α, β) で作った。S(α, β) は評価済みデータの温度にしかないので、その間の温度は
    /// 群定数を温度について線形に内挿する（それ以外の核種は、使う温度そのもので処理しておく）
    pub thermal_scattering: bool,
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
    /// P0 の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// P1 以上の散乱行列 [l - 1][入射の群][出射の群]（l = 1..=L）
    pub anisotropic: Vec<Vec<Vec<f64>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Library {
    pub groups: GroupStructure,
    /// 背景断面積 (barn)、降順。最初の値を無限希釈とみなす
    pub sigma0: Vec<f64>,
    /// ポイントワイズの断面積を持つエネルギーの範囲 [下限, 上限] (eV)。分離共鳴を持つ核種がなければ、超微細群の計算に
    /// 使うポイントワイズの断面積は要らないので `None`（どの核種のポイントワイズの断面積も空）
    pub pointwise_range: Option<[f64; 2]>,
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
    #[error("{name}: cannot interpolate to {temperature} K: {reason}")]
    Interpolation {
        name: String,
        temperature: f64,
        reason: String,
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
    /// 範囲がなければ空の配列（項目は欠かさない）
    pointwise_range: Vec<f64>,
    nuclides: Vec<NuclideData>,
}

impl Library {
    /// MessagePack（フィールド名付きの map 形式）で書く。読むときと同じ検査（配列の長さ、有限な値）を
    /// 書く前にもして、壊れたライブラリを書かない。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, LibraryError> {
        if let Some(range) = self.pointwise_range {
            check_range(range).map_err(LibraryError::Invalid)?;
        }
        let ng = self.groups.len();
        for n in &self.nuclides {
            n.check(ng, self.sigma0.len(), self.pointwise_range)
                .map_err(LibraryError::Invalid)?;
        }
        let file = LibraryFile {
            format: FORMAT.into(),
            group_bounds: self.groups.bounds().to_vec(),
            sigma0: self.sigma0.clone(),
            pointwise_range: self.pointwise_range.map(Vec::from).unwrap_or_default(),
            nuclides: self.nuclides.clone(),
        };
        Ok(msgpack::to_named(&file)?)
    }

    /// 同じ群構造、σ0、ポイントワイズの範囲を持つライブラリを、核種を並べて 1 つにする（核種の順は `parts` の順）。
    /// 核種ごとに処理した結果を組み立てるのに使う。
    pub fn merge(parts: Vec<Library>) -> Result<Library, LibraryError> {
        let mut parts = parts.into_iter();
        let mut merged = parts
            .next()
            .ok_or_else(|| LibraryError::Invalid("no library to merge".into()))?;
        for part in parts {
            let differs =
                |what: &str| LibraryError::Invalid(format!("the libraries differ in {what}"));
            if part.groups != merged.groups {
                return Err(differs("the group structure"));
            }
            if part.sigma0 != merged.sigma0 {
                return Err(differs("sigma0"));
            }
            if part.pointwise_range != merged.pointwise_range {
                return Err(differs("the pointwise range"));
            }
            for nuclide in part.nuclides {
                if merged.nuclides.iter().any(|n| n.name == nuclide.name) {
                    return Err(LibraryError::Invalid(format!(
                        "nuclide {} appears twice",
                        nuclide.name
                    )));
                }
                merged.nuclides.push(nuclide);
            }
        }
        Ok(merged)
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
        let pointwise_range = match file.pointwise_range[..] {
            [] => None,
            [low, high] => {
                check_range([low, high]).map_err(invalid)?;
                Some([low, high])
            }
            ref other => {
                return Err(invalid(format!(
                    "pointwise range must be empty or [low, high], got {other:?}"
                )));
            }
        };
        let ng = groups.len();
        for (k, n) in file.nuclides.iter().enumerate() {
            if file.nuclides[..k].iter().any(|m| m.name == n.name) {
                return Err(invalid(format!("nuclide {} appears twice", n.name)));
            }
            n.check(ng, file.sigma0.len(), pointwise_range)
                .map_err(invalid)?;
        }
        Ok(Self {
            groups,
            sigma0: file.sigma0,
            pointwise_range,
            nuclides: file.nuclides,
        })
    }

    /// 核種 `name` の、温度 `temperature` (K) での無限希釈の群定数。温度の扱いは [`Library::lookup`] と同じ。
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
            scatter: t.scatter[0].clone(),
            anisotropic: t.scatter[1..].to_vec(),
        })
    }
}

impl Library {
    /// 核種 `name` の温度 `temperature` (K) の表。温度は表の値と一致すること。ただし熱中性子散乱を S(α, β) で
    /// 作った核種（`thermal_scattering`）は、表の2つの温度の間なら群定数を温度について線形に内挿する。
    pub fn lookup(
        &self,
        name: &str,
        temperature: f64,
    ) -> Result<(&NuclideData, Cow<'_, TemperatureData>), LookupError> {
        let data = self
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .ok_or_else(|| LookupError::Nuclide(name.into()))?;
        let missing = || LookupError::Temperature {
            name: name.into(),
            temperature,
            available: data.temperatures.iter().map(|t| t.temperature).collect(),
        };
        if let Some(t) = data
            .temperatures
            .iter()
            .find(|t| t.temperature == temperature)
        {
            return Ok((data, Cow::Borrowed(t)));
        }
        if !data.thermal_scattering {
            return Err(missing());
        }
        let refuse = |reason: String| LookupError::Interpolation {
            name: name.into(),
            temperature,
            reason,
        };
        // 共鳴の遮蔽とドップラー拡がりは温度で線形にならないので、内挿するのは共鳴のない核種だけ
        if data.resolved.is_some() {
            return Err(refuse(
                "nuclides with resolved resonances need their exact temperature".into(),
            ));
        }
        // 表の温度のうち、温度を挟む最も近い2つ
        let below = data
            .temperatures
            .iter()
            .filter(|t| t.temperature < temperature)
            .max_by(|a, b| a.temperature.total_cmp(&b.temperature));
        let above = data
            .temperatures
            .iter()
            .filter(|t| t.temperature > temperature)
            .min_by(|a, b| a.temperature.total_cmp(&b.temperature));
        let (Some(low), Some(high)) = (below, above) else {
            return Err(missing());
        };
        let table = interpolate(low, high, temperature).map_err(refuse)?;
        Ok((data, Cow::Owned(table)))
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
            capture: Vec::with_capacity(groups),
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
            result.capture.push(value(|x| &x.capture));
        }
        Ok(result)
    }
}

/// 2つの温度の表 `low`, `high` から、温度 `temperature` の表を温度について線形に内挿して作る。
/// ポイントワイズの断面積は両方のエネルギーの点を合わせた格子で、それぞれを lin-lin で評価してから内挿する
/// （不連続点のあるポイントワイズの断面積は、左右の極限の対応があいまいになるので内挿しない）。
fn interpolate(
    low: &TemperatureData,
    high: &TemperatureData,
    temperature: f64,
) -> Result<TemperatureData, String> {
    let f = (temperature - low.temperature) / (high.temperature - low.temperature);
    // 形のそろわない表を zip で黙って切り詰めない
    let same = low.chi.len() == high.chi.len()
        && low.scatter.len() == high.scatter.len()
        && low.shielded.len() == high.shielded.len();
    if !same {
        return Err("the two temperatures have tables of different shapes".into());
    }
    let mix = |a: &[f64], b: &[f64]| -> Vec<f64> {
        a.iter().zip(b).map(|(x, y)| x + f * (y - x)).collect()
    };
    let scatter = low
        .scatter
        .iter()
        .zip(&high.scatter)
        .map(|(a, b)| a.iter().zip(b).map(|(x, y)| mix(x, y)).collect())
        .collect();
    let shielded = low
        .shielded
        .iter()
        .zip(&high.shielded)
        .map(|(a, b)| Shielded {
            total: mix(&a.total, &b.total),
            elastic: mix(&a.elastic, &b.elastic),
            fission: mix(&a.fission, &b.fission),
            nu_fission: mix(&a.nu_fission, &b.nu_fission),
            capture: mix(&a.capture, &b.capture),
        })
        .collect();
    let same_channels = |a: &GroupReaction, b: &GroupReaction| {
        a.mt == b.mt
            && a.channels.len() == b.channels.len()
            && a.channels
                .iter()
                .zip(&b.channels)
                .all(|(x, y)| x.final_state == y.final_state)
    };
    if low.reactions.len() != high.reactions.len()
        || !low
            .reactions
            .iter()
            .zip(&high.reactions)
            .all(|(a, b)| same_channels(a, b))
    {
        return Err("the two temperatures have different chain reactions".into());
    }
    let reactions = low
        .reactions
        .iter()
        .zip(&high.reactions)
        .map(|(a, b)| GroupReaction {
            mt: a.mt,
            channels: a
                .channels
                .iter()
                .zip(&b.channels)
                .map(|(x, y)| GroupChannel {
                    final_state: x.final_state,
                    cross_section: mix(&x.cross_section, &y.cross_section),
                })
                .collect(),
        })
        .collect();
    let (p, q) = (&low.pointwise, &high.pointwise);
    let repeated = |e: &[f64]| e.windows(2).any(|w| w[0] == w[1]);
    if repeated(&p.energy) || repeated(&q.energy) {
        return Err("pointwise data have discontinuities".into());
    }
    // 両方のポイントワイズの断面積が同じ範囲を覆うこと（範囲の外を端の値で埋めない）
    if p.energy.first() != q.energy.first() || p.energy.last() != q.energy.last() {
        return Err("pointwise data cover different ranges".into());
    }
    let mut energy: Vec<f64> = p.energy.iter().chain(&q.energy).copied().collect();
    energy.sort_by(f64::total_cmp);
    energy.dedup();
    let column = |a: &[f64], b: &[f64]| -> Vec<f64> {
        energy
            .iter()
            .map(|&e| {
                let (x, y) = (lin_lin(&p.energy, a, e), lin_lin(&q.energy, b, e));
                x + f * (y - x)
            })
            .collect()
    };
    let pointwise = Pointwise {
        total: column(&p.total, &q.total),
        elastic: column(&p.elastic, &q.elastic),
        fission: column(&p.fission, &q.fission),
        nu_fission: column(&p.nu_fission, &q.nu_fission),
        capture: column(&p.capture, &q.capture),
        energy,
    };
    Ok(TemperatureData {
        temperature,
        chi: mix(&low.chi, &high.chi),
        scatter,
        shielded,
        pointwise,
        reactions,
    })
}

/// 昇順の格子 `x` の値 `y` を、`e` で lin-lin に評価する（範囲の外は端の値）
fn lin_lin(x: &[f64], y: &[f64], e: f64) -> f64 {
    let n = x.len().min(y.len());
    if n == 0 {
        return 0.0;
    }
    let k = x[..n].partition_point(|&v| v <= e);
    if k == 0 {
        return y[0];
    }
    if k == n {
        return y[n - 1];
    }
    let (x0, x1, y0, y1) = (x[k - 1], x[k], y[k - 1], y[k]);
    y0 + (y1 - y0) * (e - x0) / (x1 - x0)
}

/// ポイントワイズの断面積の範囲が正の有限な値で、下限 < 上限であることを確かめる。
fn check_range([low, high]: [f64; 2]) -> Result<(), String> {
    if low > 0.0 && high > low && high.is_finite() {
        Ok(())
    } else {
        Err(format!(
            "pointwise range [{low}, {high}] must be positive, finite and increasing"
        ))
    }
}

/// 構造体の位置（最上位、核種、温度、遮蔽断面積）がすべて map であることを確かめる。
/// 項目の有無と型は、そのあとの serde の変換で確かめる。
fn require_maps(value: &rmpv::Value) -> Result<(), MsgpackError> {
    use msgpack::{array, map};
    for nuclide in array(map(value)?, "nuclides") {
        for temperature in array(map(nuclide)?, "temperatures") {
            let temperature = map(temperature)?;
            for shielded in array(temperature, "shielded") {
                map(shielded)?;
            }
            for reaction in array(temperature, "reactions") {
                for channel in array(map(reaction)?, "channels") {
                    let channel = map(channel)?;
                    if let Some((_, state)) = channel
                        .iter()
                        .find(|(k, _)| k.as_str() == Some("final_state"))
                        && !state.is_nil()
                    {
                        map(state)?;
                    }
                }
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
    fn check(&self, ng: usize, n_sigma0: usize, range: Option<[f64; 2]>) -> Result<(), String> {
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
        if let Some(r) = &self.resolved
            && !(r.top > 0.0 && r.top.is_finite() && r.potential >= 0.0 && r.potential.is_finite())
        {
            return Err(format!(
                "{name}: invalid resolved range top {} eV or potential scattering {} b",
                r.top, r.potential
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
            let square = |m: &Vec<Vec<f64>>| m.len() == ng && m.iter().all(|row| sized(row));
            if t.scatter.is_empty() || !t.scatter.iter().all(square) {
                return Err(at(format!(
                    "scatter must be one or more {ng} x {ng} Legendre moments"
                )));
            }
            if t.shielded.len() != n_sigma0 {
                return Err(at(format!(
                    "{} shielded sets for {n_sigma0} sigma0 values",
                    t.shielded.len()
                )));
            }
            let all_sized = |s: &Shielded| {
                [&s.total, &s.elastic, &s.fission, &s.nu_fission, &s.capture]
                    .iter()
                    .all(|v| sized(v))
            };
            if !t.shielded.iter().all(all_sized) {
                return Err(at(format!("shielded cross sections must have {ng} groups")));
            }
            if t.reactions.windows(2).any(|w| w[1].mt <= w[0].mt) {
                return Err(at("chain reactions must be in increasing MT".into()));
            }
            for r in &t.reactions {
                if r.channels.is_empty() {
                    return Err(at(format!("chain reaction MT{} has no channels", r.mt)));
                }
                if !CHAIN_REACTIONS.iter().any(|c| c.mt == r.mt) {
                    return Err(at(format!("MT{} is not a chain reaction", r.mt)));
                }
                // 鎖の経路と同じ並び: 状態に分けないならただ 1 つの経路、分けるなら MF9 か MF10 の状態で、残り（LFS が nil）は最後に 1 つ
                let states: Vec<Option<ChannelState>> =
                    r.channels.iter().map(|c| c.final_state).collect();
                let layout = match states[..] {
                    [None] => true,
                    _ => {
                        let last = states.len() - 1;
                        states.iter().enumerate().all(|(i, s)| match s {
                            Some(ChannelState { mf: 9 | 10, lfs }) => lfs.is_some() || i == last,
                            _ => false,
                        })
                    }
                };
                if !layout {
                    return Err(at(format!(
                        "chain reaction MT{} channels must be one unsplit channel or MF9/MF10 states \
                         with the remainder last",
                        r.mt
                    )));
                }
                if !r.channels.iter().all(|c| sized(&c.cross_section)) {
                    return Err(at(format!(
                        "chain reaction MT{} channels must have {ng} groups",
                        r.mt
                    )));
                }
            }
            let p = &t.pointwise;
            let columns = [&p.total, &p.elastic, &p.fission, &p.nu_fission, &p.capture];
            if columns.iter().any(|c| c.len() != p.energy.len()) {
                return Err(at("pointwise columns must match the energies".into()));
            }
            if p.energy.windows(2).any(|w| w[1] < w[0]) {
                return Err(at("pointwise energies must not decrease".into()));
            }
            // 範囲をちょうど覆うこと（超微細群の減速計算が範囲の外を引かないように）。範囲がなければ空で、
            // 分離共鳴も持たない（範囲は分離共鳴を持つ核種がないときだけない）
            if range.is_none() && self.resolved.is_some() {
                return Err(format!(
                    "{name}: resolved resonances need a pointwise range"
                ));
            }
            match range {
                Some(range) => {
                    if p.energy.first() != Some(&range[0]) || p.energy.last() != Some(&range[1]) {
                        return Err(at(format!("pointwise data must cover the range {range:?}")));
                    }
                }
                None => {
                    if !p.energy.is_empty() {
                        return Err(at("pointwise data without a pointwise range".into()));
                    }
                }
            }
            let values = t
                .chi
                .iter()
                .chain(t.scatter.iter().flatten().flatten())
                .chain(
                    t.shielded
                        .iter()
                        .flat_map(|s| [&s.total, &s.elastic, &s.fission, &s.nu_fission, &s.capture])
                        .flatten(),
                )
                .chain(
                    t.reactions
                        .iter()
                        .flat_map(|r| &r.channels)
                        .flat_map(|c| &c.cross_section),
                )
                .chain(p.energy.iter())
                .chain(columns.into_iter().flatten());
            if let Some(bad) = values.copied().find(|v| !v.is_finite()) {
                return Err(at(format!("non-finite value {bad}")));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
pub(crate) mod tests;
