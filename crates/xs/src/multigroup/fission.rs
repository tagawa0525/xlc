//! 核分裂スペクトル χ。

use super::{GroupStructure, GroupedLaw, Weight, group_integrals_multi};
use crate::endf::{Mf5, Tab1};

/// 核分裂源 ν(E)σf(E)w(E) で重み付けした核分裂スペクトル χ_g（和は 1）。
/// `fission` は σf（MT18）。
pub fn fission_spectrum(
    groups: &GroupStructure,
    weight: &Weight,
    fission: &Tab1,
    nubar: &Tab1,
    mf5: &Mf5,
) -> Vec<f64> {
    let parts: Vec<(&Tab1, GroupedLaw)> = mf5
        .subsections
        .iter()
        .map(|(p, law)| (p, GroupedLaw::new(groups, law)))
        .collect();
    let mut breakpoints: Vec<f64> = fission.x.iter().chain(&nubar.x).copied().collect();
    for (p, grouped) in &parts {
        breakpoints.extend(&p.x);
        breakpoints.extend(grouped.breakpoints());
    }
    let ng = groups.len();
    let rows = group_integrals_multi(groups, weight, &breakpoints, ng, |e, out| {
        let source = nubar.eval(e).unwrap_or(0.0) * fission.eval(e).unwrap_or(0.0);
        if source == 0.0 {
            return;
        }
        for (p, grouped) in &parts {
            let pk = p.eval(e).unwrap_or(0.0);
            if pk != 0.0 {
                grouped.accumulate(e, source * pk, out);
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
