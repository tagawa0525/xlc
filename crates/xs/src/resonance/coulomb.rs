//! 荷電粒子のチャネルの Coulomb の透過因子 P_l(η, ρ) = ρ / (F_l² + G_l²)（ENDF-102 の式 D.84）。

/// 軌道角運動量 l、Coulomb のパラメータ η ≥ 0、ρ > 0 の透過因子 P_l(η, ρ)。
pub(super) fn penetrability(l: u32, eta: f64, rho: f64) -> f64 {
    let _ = (l, eta, rho);
    todo!("the Coulomb penetrability")
}

#[cfg(test)]
mod tests;
