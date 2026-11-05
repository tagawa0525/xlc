//! 熱水力のフィードバック（ADR-0229）: 集合体の位置ごとの閉じた 1 次元の流路と、燃料棒の径方向の定常の熱伝導。流路の
//! 減速材の状態はノードの中点のエンタルピーの IF97 の温度、燃料の温度はノードの線出力の燃料棒の実効温度。

use std::collections::BTreeMap;

use xlc_format::table::Table;

use crate::input::Layout;

pub mod rod;
pub mod water;

pub use rod::{Effective, Rod, UO2_THEORETICAL_DENSITY};

/// 出力との繰り返しの許容の既定（ADR-0229 の案 I1）: 核分裂のあるノードの出力の変化の最大 / ノードの出力の平均
pub const TOLERANCE: f64 = 1.0e-5;
/// 繰り返しの回数の上限の既定（ADR-0229 の案 I1）
pub const MAX_ITERATIONS: usize = 50;

/// ケースの `[thermal-hydraulics]`（ADR-0229）
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalHydraulics {
    /// 炉心の熱出力 (MW)
    pub power: f64,
    /// 燃料の集合体を流れる流量の和 (kg/s)。燃料の集合体に等しく分ける
    pub flow: f64,
    /// MPa
    pub pressure: f64,
    /// ノードの出力のうち燃料で発生する割合
    pub fuel_heat_fraction: f64,
    pub effective: Effective,
    pub tolerance: f64,
    pub max_iterations: usize,
    pub rod: Rod,
    /// 集合体あたりの燃料棒の本数
    pub rods: usize,
    /// 集合体あたりの流路の面積 (cm²)
    pub flow_area: f64,
    /// 集合体あたりの濡れ縁 (cm)
    pub wetted_perimeter: f64,
}

/// ノードの減速材の温度と、テーブルの燃料の温度の軸に渡す燃料の温度 (K)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodeTemperatures {
    pub moderator: f64,
    pub fuel: f64,
}

/// テーブルが核分裂するか
pub fn fissile(table: &Table) -> bool {
    table
        .states
        .iter()
        .any(|s| s.nu_fission.iter().any(|&v| v > 0.0))
}

/// 並びの流路: 集合体の位置（反射体の位置を除く）ごとに 1 本
#[derive(Debug, Clone, PartialEq)]
pub struct Channels {
    /// 位置ごとの流路の番号（流路のない位置は `None`）
    channel: Vec<Option<usize>>,
    /// 流路ごとの位置（位置の番号の順、つまり下の層から）
    sites: Vec<Vec<usize>>,
    /// 流路ごとの集合体の位置（x、y の番号。下の行から）
    positions: Vec<[usize; 2]>,
    /// 位置ごとに、核分裂のあるテーブルを使うか
    fissile: Vec<bool>,
    /// 位置ごとに、ノードがあるか
    present: Vec<bool>,
    /// 位置ごとの (層の番号, 高さ (cm), 体積 (cm³))
    layers: Vec<(usize, f64, f64)>,
    /// 集合体の位置を各方向に分けたノードの数
    split: usize,
}

/// 流路の解
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelReport {
    /// 集合体の位置（x、y の番号。下の行から）
    pub position: [usize; 2],
    /// 流路の位置（下の層から）
    pub sites: Vec<usize>,
    /// 入口から出口のエンタルピーの上昇 (kJ/kg)
    pub rise: f64,
    /// 出口の温度 (K)
    pub outlet: f64,
}

/// 出力と熱水力の繰り返しの解（状態の解 `case::Solved` の `thermal`）
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalReport {
    /// 固有値計算の回数
    pub iterations: usize,
    /// 最後の回の、核分裂のあるノードの出力の変化の最大 / ノードの出力の平均
    pub change: f64,
    /// 位置ごとのノードの温度（最後の固有値計算に使ったもの）
    pub nodes: Vec<Option<NodeTemperatures>>,
    pub channels: Vec<ChannelReport>,
    /// 位置ごとの出力 (W)（`nodes` を求めた出力）
    pub powers: Vec<f64>,
    /// 流路の核分裂のあるノードの燃料の温度の体積平均と最大 (K)
    pub fuel_average: f64,
    pub fuel_max: f64,
    /// 流路のノードの、テーブルの減速材の軸の格子点の密度を温度で内挿した密度と IF97 の密度の差の最大 (g/cc)
    pub density_error: f64,
}

/// 流路と燃料棒の解
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelSolution {
    /// 位置ごとのノードの温度（ノードのない位置は `None`）
    pub nodes: Vec<Option<NodeTemperatures>>,
    pub channels: Vec<ChannelReport>,
    /// 流路の核分裂のあるノードの燃料の温度の体積平均と最大 (K)
    pub fuel_average: f64,
    pub fuel_max: f64,
}

impl Channels {
    /// `layout` の流路。`tables` は `layout.files` の順
    pub fn new(layout: &Layout, tables: &[Table]) -> Result<Self, String> {
        if tables.len() != layout.files.len() {
            return Err(format!(
                "{} tables for {} few-group constant files",
                tables.len(),
                layout.files.len()
            ));
        }
        let split = layout.split;
        let n = layout.node_material.len();
        let mut channel = vec![None; n];
        let mut fissile_site = vec![false; n];
        let mut present = vec![false; n];
        let mut layers = vec![(0, 0.0, 0.0); n];
        let mut index: BTreeMap<[usize; 2], usize> = BTreeMap::new();
        let mut sites: Vec<Vec<usize>> = Vec::new();
        let mut positions = Vec::new();
        for (site, material) in layout.node_material.iter().enumerate() {
            let [x, y, z] = layout.mesh.position(site);
            layers[site] = (
                z,
                layout.mesh.dz[z],
                layout.mesh.dx[x] * layout.mesh.dy[y] * layout.mesh.dz[z],
            );
            let Some(m) = *material else { continue };
            present[site] = true;
            let (reflector, files): (bool, Vec<usize>) = match layout.materials.get(m) {
                Some(material) => (material.reflector, vec![material.file]),
                None => {
                    let partial = layout
                        .partial
                        .get(m - layout.materials.len())
                        .ok_or_else(|| format!("site {site} has no material {m}"))?;
                    (
                        false,
                        partial
                            .segments
                            .iter()
                            .map(|&(t, _)| layout.tables[t].file)
                            .collect(),
                    )
                }
            };
            fissile_site[site] = files.iter().any(|&f| fissile(&tables[f]));
            if reflector {
                continue;
            }
            let key = [x / split, y / split];
            let c = *index.entry(key).or_insert_with(|| {
                sites.push(Vec::new());
                positions.push(key);
                sites.len() - 1
            });
            sites[c].push(site);
            channel[site] = Some(c);
        }
        if sites.is_empty() {
            return Err("the core has no assembly position for a thermal-hydraulic channel".into());
        }
        Ok(Self {
            channel,
            sites,
            positions,
            fissile: fissile_site,
            present,
            layers,
            split,
        })
    }

    /// 流路の数（燃料の集合体の数）
    pub fn count(&self) -> usize {
        self.sites.len()
    }

    /// 位置 `site` の流路の番号
    pub fn channel_of(&self, site: usize) -> Option<usize> {
        self.channel.get(site).copied().flatten()
    }

    /// 位置 `site` が核分裂のあるテーブルを使うか
    pub fn fissile(&self, site: usize) -> bool {
        self.fissile.get(site).copied().unwrap_or(false)
    }

    /// 位置 `site` の体積 (cm³)
    pub fn volume(&self, site: usize) -> f64 {
        self.layers.get(site).map_or(0.0, |l| l.2)
    }

    /// 位置 `site` にノードがあるか
    pub fn present(&self, site: usize) -> bool {
        self.present.get(site).copied().unwrap_or(false)
    }

    /// 位置ごとの出力 `powers` (W) と入口の温度 `inlet` (K) の、流路の減速材と燃料棒の解。燃料棒の燃焼度は `burnup`（MWd/kg。
    /// なければ 0）。流路のない位置（径方向の反射体）は入口の温度。出力は流路の核分裂のある位置だけが持てる
    pub fn solve(
        &self,
        th: &ThermalHydraulics,
        inlet: f64,
        powers: &[f64],
        burnup: Option<f64>,
    ) -> Result<ChannelSolution, String> {
        if powers.len() != self.channel.len() {
            return Err(format!(
                "{} node powers for {} positions",
                powers.len(),
                self.channel.len()
            ));
        }
        for (site, &p) in powers.iter().enumerate() {
            if !(p >= 0.0 && p.is_finite()) {
                return Err(format!(
                    "site {site}: the node power {p} W must be non-negative"
                ));
            }
            if p > 0.0 && !(self.fissile(site) && self.channel_of(site).is_some()) {
                return Err(format!(
                    "site {site}: the node power {p} W is on a position without a fissile channel node"
                ));
            }
        }
        let p = th.pressure;
        let flow = th.flow / self.count() as f64;
        let h_in = water::region_one(p, inlet)
            .map_err(|e| format!("the inlet: {e}"))?
            .enthalpy;
        let area = th.flow_area * 1e-4;
        let mass_flux = flow / area;
        let diameter = 4.0 * th.flow_area / th.wetted_perimeter * 1e-2;
        let rods = th.rods as f64 / (self.split * self.split) as f64;
        let burnup = burnup.unwrap_or(0.0);
        let mut nodes: Vec<Option<NodeTemperatures>> = self
            .present
            .iter()
            .map(|&here| {
                here.then_some(NodeTemperatures {
                    moderator: inlet,
                    fuel: inlet,
                })
            })
            .collect();
        let mut channels = Vec::with_capacity(self.count());
        let (mut fuel_sum, mut volume_sum, mut fuel_max) = (0.0, 0.0, f64::NEG_INFINITY);
        for (sites, &position) in self.sites.iter().zip(&self.positions) {
            // 下の層から、層ごとの出力の和でエンタルピーを進める
            let mut h = h_in;
            let mut start = 0;
            while start < sites.len() {
                let layer = self.layers[sites[start]].0;
                let end = start
                    + sites[start..]
                        .iter()
                        .take_while(|&&s| self.layers[s].0 == layer)
                        .count();
                let layer_sites = &sites[start..end];
                let layer_power: f64 = layer_sites.iter().map(|&s| powers[s]).sum();
                // kJ/kg
                let gain = layer_power / flow * 1e-3;
                let moderator = water::temperature(p, h + gain / 2.0).map_err(|e| {
                    format!("channel at assembly position {position:?}, layer {layer}: {e}")
                })?;
                for &site in layer_sites {
                    let height = self.layers[site].1 * 1e-2;
                    let fuel = if self.fissile(site) {
                        let q = powers[site] * th.fuel_heat_fraction / (rods * height);
                        let t = th
                            .rod
                            .temperatures(
                                q,
                                moderator,
                                p,
                                mass_flux,
                                diameter,
                                burnup,
                                th.effective,
                            )
                            .map_err(|e| format!("site {site}: the fuel rod: {e}"))?;
                        let volume = self.layers[site].2;
                        fuel_sum += t.effective * volume;
                        volume_sum += volume;
                        fuel_max = fuel_max.max(t.effective);
                        t.effective
                    } else {
                        moderator
                    };
                    nodes[site] = Some(NodeTemperatures { moderator, fuel });
                }
                h += gain;
                start = end;
            }
            let outlet = water::temperature(p, h)
                .map_err(|e| format!("the outlet at assembly position {position:?}: {e}"))?;
            channels.push(ChannelReport {
                position,
                sites: sites.clone(),
                rise: h - h_in,
                outlet,
            });
        }
        Ok(ChannelSolution {
            nodes,
            channels,
            fuel_average: fuel_sum / volume_sum,
            fuel_max,
        })
    }
}

#[cfg(test)]
mod tests;
