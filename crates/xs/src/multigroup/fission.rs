//! 核分裂スペクトル χ。

use xlc_format::library::DelayedGroups;

use super::{GroupStructure, GroupedLaw, GroupedOutgoing, Weight, group_integrals_multi};
use crate::endf::{Frame, Mf5, OutgoingDistributions, Tab1};

/// 核分裂で生まれる中性子の 1 つの成分（即発、遅発の全族、または区別しない全数）の、1 核分裂あたりの数 ν(E) と
/// 出射エネルギーの分布。
#[derive(Debug, Clone)]
pub struct FissionNeutrons {
    pub nubar: Tab1,
    pub spectrum: FissionSpectrum,
}

/// 核分裂の中性子の出射エネルギーの分布。
#[derive(Debug, Clone)]
pub enum FissionSpectrum {
    /// MF5 の部分分布（割合 p_k(E) と分布）
    Mf5(Mf5),
    /// MF6 MT18 の中性子の表の分布（JPN=0 の平均の分布。ENDF-102 の 6.3.9.3）。`multiplicity` は y(E)（1 を読み込みで
    /// 確かめる）、`awr` は重心系から実験室系への変換に使う標的の質量比
    Mf6 {
        multiplicity: Tab1,
        distributions: OutgoingDistributions,
        frame: Frame,
        awr: f64,
    },
}

impl FissionSpectrum {
    /// 部分分布ごとの割合（MF6 は y(E)）と、群に振り分けた分布。MF6 の重心系の分布は、MF6 のほかの反応と同じく重心系で
    /// 等方として実験室系に変換する（`GroupedOutgoing::center_of_mass`）
    fn grouped(&self, groups: &GroupStructure) -> Vec<(&Tab1, GroupedLaw)> {
        match self {
            Self::Mf5(mf5) => mf5
                .subsections
                .iter()
                .map(|(p, law)| (p, GroupedLaw::new(groups, law)))
                .collect(),
            Self::Mf6 {
                multiplicity,
                distributions,
                frame,
                awr,
            } => {
                let grouped = match frame {
                    Frame::Laboratory => GroupedOutgoing::new(groups, distributions),
                    Frame::CenterOfMass => {
                        GroupedOutgoing::center_of_mass(groups, distributions, *awr)
                    }
                };
                vec![(multiplicity, GroupedLaw::Table(grouped))]
            }
        }
    }
}

/// 核分裂源で重み付けした核分裂スペクトル χ_g（和は 1）。成分 c の中性子は ν_c(E)σf(E)w(E) の重みで、その成分の
/// 分布に生まれる。`fission` は σf（MT18）。
pub fn fission_spectrum(
    groups: &GroupStructure,
    weight: &Weight,
    fission: &Tab1,
    parts: &[FissionNeutrons],
) -> Vec<f64> {
    let parts: Vec<(&Tab1, Vec<(&Tab1, GroupedLaw)>)> = parts
        .iter()
        .map(|part| (&part.nubar, part.spectrum.grouped(groups)))
        .collect();
    let mut breakpoints: Vec<f64> = fission.x.clone();
    for (nubar, laws) in &parts {
        breakpoints.extend(&nubar.x);
        for (p, grouped) in laws {
            breakpoints.extend(&p.x);
            breakpoints.extend(grouped.breakpoints());
        }
    }
    let ng = groups.len();
    let rows = group_integrals_multi(groups, weight, &breakpoints, ng, |e, out| {
        let sigma = fission.eval(e).unwrap_or(0.0);
        if sigma == 0.0 {
            return;
        }
        for (nubar, laws) in &parts {
            let source = nubar.eval(e).unwrap_or(0.0) * sigma;
            if source == 0.0 {
                continue;
            }
            for (p, grouped) in laws {
                let pk = p.eval(e).unwrap_or(0.0);
                if pk != 0.0 {
                    grouped.accumulate(e, source * pk, out);
                }
            }
        }
    });
    let mut chi = vec![0.0; ng];
    for row in rows {
        for (c, v) in chi.iter_mut().zip(row) {
            *c += v;
        }
    }
    let total: f64 = chi.iter().sum();
    chi.iter().map(|c| c / total).collect()
}

#[cfg(test)]
mod tests;

/// 遅発中性子の族ごとの群の量（ADR-0238）。族 k は `spectrum`（MF5 MT455）の小区分 k で、その割合 p_k(E) と分布を持つ。
/// 群の遅発中性子数は ν_d,k,g = ∫ν_d p_k σf w dE / ∫σf w dE（`fission` は σf、`nubar` は ν_d）、族のスペクトルは族の中性子
/// （ν_d p_k σf w）で重みを付けた出射の分布で、和は 1（[`fission_spectrum`] に族の小区分だけを渡したもの）。
pub fn delayed_groups(
    groups: &GroupStructure,
    weight: &Weight,
    fission: &Tab1,
    nubar: &Tab1,
    spectrum: &Mf5,
) -> DelayedGroups {
    let families = spectrum.subsections.len();
    let mut breakpoints: Vec<f64> = fission.x.iter().chain(&nubar.x).copied().collect();
    for (p, _) in &spectrum.subsections {
        breakpoints.extend(&p.x);
    }
    // 成分 0 は ∫σf w、成分 1 + k は ∫ν_d p_k σf w
    let rows = group_integrals_multi(groups, weight, &breakpoints, families + 1, |e, out| {
        let sigma = fission.eval(e).unwrap_or(0.0);
        if sigma == 0.0 {
            return;
        }
        out[0] = sigma;
        let nu = nubar.eval(e).unwrap_or(0.0) * sigma;
        for (o, (p, _)) in out[1..].iter_mut().zip(&spectrum.subsections) {
            *o = nu * p.eval(e).unwrap_or(0.0);
        }
    });
    let nu = (0..families)
        .map(|k| {
            rows.iter()
                .map(|row| {
                    if row[0] > 0.0 {
                        row[k + 1] / row[0]
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let chi = spectrum
        .subsections
        .iter()
        .map(|subsection| {
            fission_spectrum(
                groups,
                weight,
                fission,
                &[FissionNeutrons {
                    nubar: nubar.clone(),
                    spectrum: FissionSpectrum::Mf5(Mf5 {
                        subsections: vec![subsection.clone()],
                    }),
                }],
            )
        })
        .collect();
    DelayedGroups { nu, chi }
}
