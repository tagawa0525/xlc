//! 核分裂スペクトル χ。

use super::{GroupStructure, GroupedLaw, Weight, group_integrals_multi};
use crate::endf::{Mf5, Tab1};

/// 核分裂で生まれる中性子の 1 つの成分（即発、遅発の全族、または区別しない全数）の、1 核分裂あたりの数 ν(E) と
/// 出射エネルギーの分布。
#[derive(Debug, Clone)]
pub struct FissionNeutrons {
    pub nubar: Tab1,
    pub spectrum: Mf5,
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
        .map(|part| {
            let laws = part
                .spectrum
                .subsections
                .iter()
                .map(|(p, law)| (p, GroupedLaw::new(groups, law)))
                .collect();
            (&part.nubar, laws)
        })
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
