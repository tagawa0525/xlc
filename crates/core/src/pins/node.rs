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
//! 平均を使う。NEM の解はノードのバランスを満たすので、外した和も再現する。

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

/// ノードの定数と大きさで決まる部分（モード、条件の行列の LU、特解の条件の値）。条件の値を変えて何度も解く。
#[derive(Debug, Clone)]
pub(crate) struct NodeModel {
    modes: [Basis; 2],
    lus: [Lu; 2],
    /// [モード][関数] の列の大きさ
    scales: [[f64; 8]; 2],
    /// [モード][群] の固有ベクトル
    vectors: [[f64; 2]; 2],
    /// V⁻¹（行がモード）
    inverse: [[f64; 2]; 2],
    diffusion: [f64; 2],
    half: [f64; 2],
    /// [モード] の特解 P の条件の値（角 4 つ、平均、面の微分 4 つ）
    particular: [[f64; 9]; 2],
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

impl NodeModel {
    pub(crate) fn new(
        m: &FewGroup,
        k: f64,
        half: [f64; 2],
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
        let diffusion = [m.diffusion[0], m.diffusion[1]];
        let mm = [
            [a(0, 0) / diffusion[0], a(0, 1) / diffusion[0]],
            [a(1, 0) / diffusion[1], a(1, 1) / diffusion[1]],
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
        let inverse = [
            [vectors[1][1] / det_v, -vectors[1][0] / det_v],
            [-vectors[0][1] / det_v, vectors[0][0] / det_v],
        ];
        let [ha, hb] = half;
        let corner_points = [(-ha, -hb), (ha, -hb), (ha, hb), (-ha, hb)];
        let modes = mus.map(|mu| Basis::new(mu, ha, hb));
        let mut lus = Vec::with_capacity(2);
        let mut scales = [[0.0; 8]; 2];
        let mut particular = [[0.0; 9]; 2];
        for (mode, basis) in modes.iter().enumerate() {
            let mut rows: Vec<[f64; 8]> = Vec::with_capacity(8);
            for (c, &(x, y)) in corner_points.iter().enumerate() {
                rows.push(basis.values(x, y));
                particular[mode][c] = basis.particular(x, y);
            }
            rows.push(basis.averages(-ha, ha, -hb, hb, rule));
            particular[mode][4] = basis.particular_average(-ha, ha, -hb, hb, rule);
            let f: Vec<[f64; 8]> = (0..4)
                .map(|face| basis.face_derivatives(face, rule))
                .collect();
            for face in 0..4 {
                particular[mode][5 + face] = basis.particular_face_derivative(face, rule);
            }
            rows.push(std::array::from_fn(|k| f[1][k] - f[0][k]));
            rows.push(std::array::from_fn(|k| f[3][k] - f[2][k]));
            rows.push(std::array::from_fn(|k| {
                hb * (f[0][k] + f[1][k]) - ha * (f[2][k] + f[3][k])
            }));
            // 列ごとに大きさをそろえる
            scales[mode] =
                std::array::from_fn(|k| rows.iter().map(|r| r[k].abs()).fold(0.0, f64::max));
            if scales[mode].iter().any(|s| s.is_nan() || *s == 0.0) {
                return Err(format!(
                    "mode {mode} (μ = {}) has a vanishing function",
                    mus[mode]
                ));
            }
            let matrix: Vec<f64> = rows
                .iter()
                .flat_map(|r| (0..8).map(move |k| r[k] / scales[mode][k]))
                .collect();
            lus.push(
                Lu::new(matrix, 8).ok_or_else(|| {
                    format!("mode {mode} (μ = {}): singular conditions", mus[mode])
                })?,
            );
        }
        let lus: [Lu; 2] = lus
            .try_into()
            .map_err(|_| "two modes are expected".to_string())?;
        Ok(Self {
            modes,
            lus,
            scales,
            vectors,
            inverse,
            diffusion,
            half,
            particular,
        })
    }

    /// 条件の値からノードの中の中性子束を求める。`c.half` はこのモデルの大きさに等しいこと
    pub(crate) fn solve(&self, c: &Conditions) -> NodeFlux {
        let to_mode = |f: [f64; 2]| -> [f64; 2] {
            std::array::from_fn(|m| self.inverse[m][0] * f[0] + self.inverse[m][1] * f[1])
        };
        let d = self.diffusion;
        let sources = to_mode([c.axial_leakage[0] / d[0], c.axial_leakage[1] / d[1]]);
        // 群の外向きの法線微分 ∂φ/∂n = −J_外向き / D。-x と -y の面では外向きの流れは −J
        let outward = |face: usize| -> [f64; 2] {
            let sign = if face.is_multiple_of(2) { 1.0 } else { -1.0 };
            [
                sign * c.currents[face][0] / d[0],
                sign * c.currents[face][1] / d[1],
            ]
        };
        let [ha, hb] = self.half;
        let mut coefficients = [[0.0; 8]; 2];
        for mode in 0..2 {
            let (r, p) = (sources[mode], &self.particular[mode]);
            let mut rhs = Vec::with_capacity(8);
            for (corner, pc) in c.corners.iter().zip(p) {
                rhs.push(to_mode(*corner)[mode] - r * pc);
            }
            rhs.push(to_mode(c.average)[mode] - r * p[4]);
            let data: [f64; 4] =
                std::array::from_fn(|face| to_mode(outward(face))[mode] - r * p[5 + face]);
            rhs.push(data[1] - data[0]);
            rhs.push(data[3] - data[2]);
            rhs.push(hb * (data[0] + data[1]) - ha * (data[2] + data[3]));
            self.lus[mode].solve(&mut rhs);
            coefficients[mode] = std::array::from_fn(|k| rhs[k] / self.scales[mode][k]);
        }
        NodeFlux {
            modes: self.modes.clone(),
            coefficients,
            sources,
            vectors: self.vectors,
        }
    }
}

impl NodeFlux {
    #[cfg(test)]
    pub(crate) fn new(
        m: &FewGroup,
        k: f64,
        c: &Conditions,
        rule: &[(f64, f64)],
    ) -> Result<Self, String> {
        Ok(NodeModel::new(m, k, c.half, rule)?.solve(c))
    }

    /// [モード] の ψ_m(x, y)
    #[cfg(test)]
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
    #[cfg(test)]
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
