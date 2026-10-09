//! ノードの中の2群の均質な中性子束（ADR-0096）。
//!
//! −D∇²φ + Aφ = q（A は除去、群の間の散乱、核分裂 χνΣf/k。q = −軸方向の漏れはノードの中で一様）を、M = D⁻¹A の
//! 固有ベクトル v_m で φ = Σ_m v_m ψ_m に分けると、モードごとに ∇²ψ_m − μ_m ψ_m = r_m（r = V⁻¹ D⁻¹ × 軸方向の漏れ）
//! になる。ψ_m は 8 つの関数（`Basis`）と特解 r_m P の和で、係数はモードごとに、群の条件を V⁻¹ で移した
//!
//! - 4 つの角の値
//! - ノード平均
//! - 面の平均の外向きの法線微分の 3 つの組み合わせ（x の両面の差、y の両面の差、b(−x と +x の和) − a(−y と +y の和)）
//!
//! で決める。面の微分の重み付きの和は ∮∂ψ/∂n = μ∫ψ で平均と結びつき、μ → 0 では平均を決めなくなるので、和の代わりに
//! 平均を使う。NEM の解はノードの釣り合いを満たすので、外した和も再現する。

use xlc_format::fewgroup::FewGroup;
use xlc_numeric::lu::Lu;

use super::basis::Basis;

/// ノードの条件。群は 2、角は（-x -y、+x -y、+x +y、-x +y）、面は（-x、+x、-y、+y）の順。
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Conditions {
    /// ノードの半幅 [a, b] (cm)。ノードの中心を原点とする
    pub half: [f64; 2],
    /// [群] のノード平均の均質な中性子束
    pub average: [f64; 2],
    /// [角][群] の均質な中性子束
    pub corners: [[f64; 2]; 4],
    /// [面][群] の面の平均の正味の中性子流（座標の正の向きを正）
    pub currents: [[f64; 2]; 4],
    /// [群] の軸方向の正味の漏れ (J+z − J−z)/hz（単位体積あたり、ノードの中で一様とする）
    pub axial_leakage: [f64; 2],
}

/// ノードの中の均質な中性子束。
#[derive(Debug, Clone)]
pub(crate) struct NodeFlux {
    modes: [Basis; 2],
    coefficients: [[f64; 8]; 2],
    /// [モード] の特解の係数 r_m
    sources: [f64; 2],
    /// [モード][群] の固有ベクトル
    vectors: [[f64; 2]; 2],
}

impl NodeFlux {
    pub(crate) fn new(
        m: &FewGroup,
        k: f64,
        c: &Conditions,
        rule: &[(f64, f64)],
    ) -> Result<Self, String> {
        if m.groups() != 2 {
            return Err(format!(
                "pin power reconstruction needs 2 groups, got {}",
                m.groups()
            ));
        }
        let a = |g: usize, h: usize| {
            let diagonal = if g == h { m.total[g] } else { 0.0 };
            diagonal - m.scatter[h][g] - m.chi[g] * m.nu_fission[h] / k
        };
        let d = [m.diffusion[0], m.diffusion[1]];
        let mm = [
            [a(0, 0) / d[0], a(0, 1) / d[0]],
            [a(1, 0) / d[1], a(1, 1) / d[1]],
        ];
        let (half_trace, det) = (
            0.5 * (mm[0][0] + mm[1][1]),
            mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0],
        );
        let disc = half_trace * half_trace - det;
        if disc.is_nan() || disc <= 0.0 {
            return Err(format!(
                "D⁻¹A has complex or equal eigenvalues (discriminant {disc}); \
                 the two-group modes are not separable"
            ));
        }
        let root = disc.sqrt();
        let mus = [half_trace - root, half_trace + root];
        // (M − μ) v = 0 の解を、2 つの行のうち大きい方から取る
        let vector = |mu: f64| {
            let (p, q) = ([mm[0][1], mu - mm[0][0]], [mu - mm[1][1], mm[1][0]]);
            let norm = |v: [f64; 2]| v[0].hypot(v[1]);
            if norm(p) >= norm(q) { p } else { q }
        };
        let vectors = mus.map(vector);
        // V = [v_0 v_1]（列が固有ベクトル）の逆行列
        let det_v = vectors[0][0] * vectors[1][1] - vectors[1][0] * vectors[0][1];
        if det_v.is_nan() || det_v == 0.0 {
            return Err("the two-group modes are not independent".into());
        }
        let to_mode = |f: [f64; 2]| -> [f64; 2] {
            [
                (vectors[1][1] * f[0] - vectors[1][0] * f[1]) / det_v,
                (-vectors[0][1] * f[0] + vectors[0][0] * f[1]) / det_v,
            ]
        };
        let sources = to_mode([c.axial_leakage[0] / d[0], c.axial_leakage[1] / d[1]]);
        let [ha, hb] = c.half;
        // 群の外向きの法線微分 ∂φ/∂n = −J_外向き / D。-x と -y の面では外向きの流れは −J
        let outward = |face: usize| -> [f64; 2] {
            let sign = if face % 2 == 0 { 1.0 } else { -1.0 };
            [
                sign * c.currents[face][0] / d[0],
                sign * c.currents[face][1] / d[1],
            ]
        };
        let corner_points = [(-ha, -hb), (ha, -hb), (ha, hb), (-ha, hb)];
        let modes = mus.map(|mu| Basis::new(mu, ha, hb));
        let mut coefficients = [[0.0; 8]; 2];
        for (mode, basis) in modes.iter().enumerate() {
            let r = sources[mode];
            let mut rows: Vec<[f64; 8]> = Vec::with_capacity(8);
            let mut rhs = Vec::with_capacity(8);
            for (corner, &(x, y)) in corner_points.iter().enumerate() {
                rows.push(basis.values(x, y));
                rhs.push(to_mode(c.corners[corner])[mode] - r * basis.particular(x, y));
            }
            rows.push(basis.averages(-ha, ha, -hb, hb, rule));
            rhs.push(
                to_mode(c.average)[mode] - r * basis.particular_average(-ha, ha, -hb, hb, rule),
            );
            let f: Vec<[f64; 8]> = (0..4)
                .map(|face| basis.face_derivatives(face, rule))
                .collect();
            let data: Vec<f64> = (0..4)
                .map(|face| {
                    to_mode(outward(face))[mode] - r * basis.particular_face_derivative(face, rule)
                })
                .collect();
            rows.push(std::array::from_fn(|k| f[1][k] - f[0][k]));
            rhs.push(data[1] - data[0]);
            rows.push(std::array::from_fn(|k| f[3][k] - f[2][k]));
            rhs.push(data[3] - data[2]);
            rows.push(std::array::from_fn(|k| {
                hb * (f[0][k] + f[1][k]) - ha * (f[2][k] + f[3][k])
            }));
            rhs.push(hb * (data[0] + data[1]) - ha * (data[2] + data[3]));
            // 列ごとに大きさをそろえて解く
            let scales: [f64; 8] =
                std::array::from_fn(|k| rows.iter().map(|r| r[k].abs()).fold(0.0, f64::max));
            if scales.iter().any(|s| s.is_nan() || *s == 0.0) {
                return Err(format!(
                    "mode {mode} (μ = {}) has a vanishing function",
                    mus[mode]
                ));
            }
            let matrix: Vec<f64> = rows
                .iter()
                .flat_map(|r| (0..8).map(move |k| r[k] / scales[k]))
                .collect();
            let lu = Lu::new(matrix, 8)
                .ok_or_else(|| format!("mode {mode} (μ = {}): singular conditions", mus[mode]))?;
            lu.solve(&mut rhs);
            coefficients[mode] = std::array::from_fn(|k| rhs[k] / scales[k]);
        }
        Ok(Self {
            modes,
            coefficients,
            sources,
            vectors,
        })
    }

    /// [モード] の ψ_m(x, y)
    fn modes_at(&self, x: f64, y: f64) -> [f64; 2] {
        std::array::from_fn(|m| {
            let v = self.modes[m].values(x, y);
            let sum: f64 = v
                .iter()
                .zip(&self.coefficients[m])
                .map(|(b, c)| b * c)
                .sum();
            sum + self.sources[m] * self.modes[m].particular(x, y)
        })
    }

    /// (x, y) での [群] の均質な中性子束
    pub(crate) fn flux(&self, x: f64, y: f64) -> [f64; 2] {
        let psi = self.modes_at(x, y);
        std::array::from_fn(|g| self.vectors[0][g] * psi[0] + self.vectors[1][g] * psi[1])
    }

    /// (x, y) での [群][∂/∂x, ∂/∂y]
    pub(crate) fn gradient(&self, x: f64, y: f64) -> [[f64; 2]; 2] {
        let per_mode: [[f64; 2]; 2] = std::array::from_fn(|m| {
            let g = self.modes[m].gradients(x, y);
            let p = self.modes[m].particular_gradient(x, y);
            std::array::from_fn(|axis| {
                let sum: f64 = g
                    .iter()
                    .zip(&self.coefficients[m])
                    .map(|(b, c)| b[axis] * c)
                    .sum();
                sum + self.sources[m] * p[axis]
            })
        });
        std::array::from_fn(|g| {
            std::array::from_fn(|axis| {
                self.vectors[0][g] * per_mode[0][axis] + self.vectors[1][g] * per_mode[1][axis]
            })
        })
    }

    /// 長方形 [x0, x1] × [y0, y1] の上の [群] の平均
    pub(crate) fn average(
        &self,
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        rule: &[(f64, f64)],
    ) -> [f64; 2] {
        let psi: [f64; 2] = std::array::from_fn(|m| {
            let v = self.modes[m].averages(x0, x1, y0, y1, rule);
            let sum: f64 = v
                .iter()
                .zip(&self.coefficients[m])
                .map(|(b, c)| b * c)
                .sum();
            sum + self.sources[m] * self.modes[m].particular_average(x0, x1, y0, y1, rule)
        });
        std::array::from_fn(|g| self.vectors[0][g] * psi[0] + self.vectors[1][g] * psi[1])
    }
}

#[cfg(test)]
mod tests;
