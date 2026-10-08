//! ノード展開法（NEM）の1ノードの計算と、ADF の1ノードの問題。core と、lattice の ADF（一般化等価理論の
//! 1ノードの問題）が共有する（ADR-0103、ADR-0003）。核分裂のないノードは NEM の代わりに `anm` の ANM で解く
//! （ADR-0102。`is_fissile` で選ぶ）。
//!
//! 横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む2次式で近似する。
//! 面では部分中性子流を使い、入射の部分中性子流 i を、出ていく部分中性子流 o から i = β o + 外からの入射
//! の形で与える。

use crate::anm::{AnmError, AnmModes};

/// νΣf が負か有限でない（核の種類を選べない）。
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[error("nu-fission must be finite and non-negative, got {0}")]
pub struct NuFissionError(pub f64);

/// 核分裂のあるノード（νΣf が正の群がある）か。核分裂のないノード（どの群も 0）は NEM の代わりに ANM で
/// 解く（ADR-0102）。負や有限でない νΣf は、推測せずにエラーにする。
pub fn is_fissile(nu_fission: &[f64]) -> Result<bool, NuFissionError> {
    if let Some(&bad) = nu_fission.iter().find(|v| !(v.is_finite() && **v >= 0.0)) {
        return Err(NuFissionError(bad));
    }
    Ok(nu_fission.iter().any(|&v| v > 0.0))
}

/// 面の入射の部分中性子流 i = β o + `external`（o は同じ面から出ていく部分中性子流）。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Inflow {
    pub beta: f64,
    pub external: f64,
}

/// 1つのノードの1つの方向・1つの群の入力。
#[derive(Debug, Clone, Copy)]
pub struct Line<'a> {
    pub diffusion: f64,
    /// その方向のノードの幅 (cm)
    pub width: f64,
    /// 除去断面積（全断面積 - 自群の散乱）
    pub removal: f64,
    /// [群] の群 `group` への源の結合係数 c(g'→g)（自群の核分裂を含み、自群の散乱を除く）
    pub coupling: &'a [f64],
    /// [群] のこのノード・方向の展開係数 a1〜a4（前の反復の値）
    pub coefficients: &'a [[f64; 4]],
    pub group: usize,
    /// 横方向の漏れの f1, f2 の係数（`transverse_moments`）
    pub transverse: (f64, f64),
    /// [+ 面, - 面] の入射の条件
    pub inflow: [Inflow; 2],
}

/// 1方向の途中結果。`out` と `higher` はノード平均 φ̄ = 0 のときの値で、実際の値は `slope` と
/// `higher_slope` に φ̄ をかけて足したもの。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Outflow {
    /// [+ 面, - 面] の出ていく部分中性子流
    pub out: [f64; 2],
    pub slope: [f64; 2],
    /// [+ 面, - 面] の i = β o + external の β
    pub beta: [f64; 2],
    pub external: [f64; 2],
    /// 3次・4次の係数 [a3, a4]
    pub higher: [f64; 2],
    pub higher_slope: [f64; 2],
}

/// ノード平均を決めたあとの面の値。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Faces {
    /// [+ 面, - 面] の出ていく部分中性子流
    pub out: [f64; 2],
    /// [+ 面, - 面] の入射の部分中性子流
    pub incoming: [f64; 2],
    /// 展開係数 a1〜a4
    pub coefficients: [f64; 4],
}

impl Outflow {
    /// ノード平均 `average` での面の部分中性子流と展開係数。面の中性子束は 2(o + i)。
    pub fn at(&self, average: f64) -> Faces {
        let out: [f64; 2] = std::array::from_fn(|f| self.out[f] + self.slope[f] * average);
        let incoming: [f64; 2] = std::array::from_fn(|f| self.beta[f] * out[f] + self.external[f]);
        let (phi_r, phi_l) = (2.0 * (out[0] + incoming[0]), 2.0 * (out[1] + incoming[1]));
        let [a3, a4]: [f64; 2] =
            std::array::from_fn(|i| self.higher[i] + self.higher_slope[i] * average);
        Faces {
            out,
            incoming,
            coefficients: [phi_r - phi_l, phi_r + phi_l - 2.0 * average, a3, a4],
        }
    }
}

/// 出ていく部分中性子流と 3次・4次の係数の、ノード平均についての線形の式
/// o = out + slope φ̄、(a3, a4) = higher + higher_slope φ̄。
///
/// 面の条件（入射 i = β o + 外からの入射）、a1 と a2（面の中性子束で決まる）、f1, f2 で重み付けした残差の
/// 式の a3, a4 を、このノードのこの群の中で同時に解く。前の反復の値を使うのは、ほかの群の展開係数
/// （源のモーメント）と、`inflow` と `transverse` に入っている隣のノードの値だけである。a3, a4 を前の反復の
/// a1, a2 から求めると、両側の面の正味の流れが決まるノード（1つのノードで両側が反射の方向）では
/// 正味の流れが a1 + a3/2 と 3a2 + a4/5 を決め、D/h² ≪ Σr なら a3 ≈ 10 a1、a4 ≈ 35 a2 なので、誤差が
/// 反復ごとに約 -5 倍と -2.3 倍になって発散する（#0025）。
pub fn outflow(line: &Line) -> Outflow {
    let Line {
        diffusion,
        width: h,
        removal,
        coupling,
        coefficients,
        group: g,
        transverse: (l1, l2),
        inflow,
    } = *line;
    let alpha = diffusion / h;
    // 自群の項を除いた源のモーメント（自群の項は未知の係数とともに左辺に置く）
    let own = coupling[g];
    let mut s = source_moments(coupling, coefficients, (l1, l2));
    for (si, ai) in s.iter_mut().zip(coefficients[g]) {
        *si -= own * ai;
    }
    let [w3, w4] = higher_order_weights(diffusion, h, removal);
    let net = removal - own;

    // 面の中性子束と正味の中性子流の連続から、出ていく部分中性子流 o = (o_R, o_L) は
    //   M o = A i + c,  M = [[1+8α, 4α], [4α, 1+8α]],  A = [[1-8α, -4α], [-4α, 1-8α]]
    // を満たす（c はノード平均と a3, a4 の項）。入射は i = β o + 外からの入射なので、それも同時に解く
    let [
        Inflow {
            beta: beta_r,
            external: ext_r,
        },
        Inflow {
            beta: beta_l,
            external: ext_l,
        },
    ] = inflow;
    let (p, q) = (1.0 + 8.0 * alpha, 1.0 - 8.0 * alpha);
    let w = 4.0 * alpha;
    // K = M - A B（B = diag(β_R, β_L)）
    let k = [
        [p - q * beta_r, w + w * beta_l],
        [w + w * beta_r, p - q * beta_l],
    ];
    let det = k[0][0] * k[1][1] - k[0][1] * k[1][0];
    let solve2 = |r: [f64; 2]| {
        [
            (k[1][1] * r[0] - k[0][1] * r[1]) / det,
            (k[0][0] * r[1] - k[1][0] * r[0]) / det,
        ]
    };
    // o = o0 + o3 a3 + o4 a4 + oφ φ̄
    let o0 = solve2([q * ext_r - w * ext_l, q * ext_l - w * ext_r]);
    let o3 = solve2([-alpha / 2.0, alpha / 2.0]);
    let o4 = solve2([-alpha / 5.0, -alpha / 5.0]);
    let o_average = solve2([6.0 * alpha, 6.0 * alpha]);
    // 面の中性子束 2(o + i) = 2((1 + β) o + 外からの入射) から a1 = φ_R - φ_L、a2 = φ_R + φ_L - 2φ̄
    // （φ̄ の項は下で足す）
    let moments = |o: [f64; 2], external: [f64; 2]| {
        let phi_r = 2.0 * ((1.0 + beta_r) * o[0] + external[0]);
        let phi_l = 2.0 * ((1.0 + beta_l) * o[1] + external[1]);
        [phi_r - phi_l, phi_r + phi_l]
    };
    let constant = moments(o0, [ext_r, ext_l]);
    let by_a3 = moments(o3, [0.0; 2]);
    let by_a4 = moments(o4, [0.0; 2]);
    let mut by_average = moments(o_average, [0.0; 2]);
    by_average[1] -= 2.0;
    // 重み付けした残差 (w3 - c_gg/120) a3 = (Σr - c_gg) a1/12 - s0/12 + s2/120、
    // (w4 - c_gg/700) a4 = (Σr - c_gg) a2/20 - s1/20 + s3/700 に a1, a2 を入れた 2×2 の式
    let (f3, f4) = (net / 12.0, net / 20.0);
    let m = [
        [w3 - own / 120.0 - f3 * by_a3[0], -f3 * by_a4[0]],
        [-f4 * by_a3[1], w4 - own / 700.0 - f4 * by_a4[1]],
    ];
    let det_m = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let solve_m = |r: [f64; 2]| {
        [
            (m[1][1] * r[0] - m[0][1] * r[1]) / det_m,
            (m[0][0] * r[1] - m[1][0] * r[0]) / det_m,
        ]
    };
    let higher = solve_m([
        -s[0] / 12.0 + s[2] / 120.0 + f3 * constant[0],
        -s[1] / 20.0 + s[3] / 700.0 + f4 * constant[1],
    ]);
    let higher_slope = solve_m([f3 * by_average[0], f4 * by_average[1]]);
    let combine = |base: [f64; 2], a: [f64; 2]| -> [f64; 2] {
        std::array::from_fn(|f| base[f] + o3[f] * a[0] + o4[f] * a[1])
    };
    Outflow {
        out: combine(o0, higher),
        slope: combine(o_average, higher_slope),
        beta: [beta_r, beta_l],
        external: [ext_r, ext_l],
        higher,
        higher_slope,
    }
}

/// 源のモーメント s0〜s3: 結合係数 c(g'→g) × 展開係数の和から、横方向の漏れの f1, f2 の係数を引いたもの。
fn source_moments(coupling: &[f64], coefficients: &[[f64; 4]], transverse: (f64, f64)) -> [f64; 4] {
    let mut s = [0.0; 4];
    for (c, a) in coupling.iter().zip(coefficients) {
        for (si, ai) in s.iter_mut().zip(a) {
            *si += c * ai;
        }
    }
    s[0] -= transverse.0;
    s[1] -= transverse.1;
    s
}

/// f1, f2 で重み付けした残差の式
///   w3 a3 = Σr a1/12 - s0/12 + s2/120,  w4 a4 = Σr a2/20 - s1/20 + s3/700
/// の左辺の係数 [w3, w4]。
fn higher_order_weights(diffusion: f64, width: f64, removal: f64) -> [f64; 2] {
    let dh2 = diffusion / (width * width);
    [dh2 / 2.0 + removal / 120.0, dh2 / 5.0 + removal / 700.0]
}

/// 横方向の漏れの2次式で、ノードの片側に置く条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Beside {
    /// 隣（または鏡像）のノード。幅と、その平均の漏れ
    Node { width: f64, leakage: f64 },
    /// 面の上の漏れの値
    Face(f64),
}

/// 方向の幅 `width` のノードの横方向の漏れ（平均 `leakage`）を、両側の条件を含む2次式で近似し、
/// f1, f2 の係数を返す。`minus` は - 側、`plus` は + 側。
pub fn transverse_moments(width: f64, minus: Beside, leakage: f64, plus: Beside) -> (f64, f64) {
    let h = width;
    // ノードの中心を原点とした q(x) = c0 + c1 x + c2 x² への条件 (行, 右辺)
    let average =
        |a: f64, b: f64, value: f64| ([1.0, 0.5 * (a + b), (a * a + a * b + b * b) / 3.0], value);
    let side = |sign: f64, beside: Beside| {
        let near = sign * 0.5 * h;
        match beside {
            Beside::Node { width, leakage } => {
                let far = near + sign * width;
                average(near.min(far), near.max(far), leakage)
            }
            Beside::Face(value) => ([1.0, near, near * near], value),
        }
    };
    let rows = [
        side(-1.0, minus),
        average(-0.5 * h, 0.5 * h, leakage),
        side(1.0, plus),
    ];
    let c = solve3(
        [rows[0].0, rows[1].0, rows[2].0],
        [rows[0].1, rows[1].1, rows[2].1],
    );
    (c[1] * h, c[2] * h * h / 3.0)
}

/// 3×3 の連立一次方程式（Cramer の公式）。
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    std::array::from_fn(|j| {
        let mut m = a;
        for (row, value) in m.iter_mut().zip(b) {
            row[j] = value;
        }
        det(m) / d
    })
}

/// 1ノードの GET の問題で、参照の格子計算の外周の境界。面の条件は参照の正味の中性子流で決まり、
/// これは横方向の漏れの2次式の外側の扱いだけを決める（core の `Reflective` と `Vacuum` と同じ）。
///
/// 複数のノードの格子（ADR-0095）の内側の辺は、core の隣のノードと同じく、隣のノードの幅と参照の
/// 横方向の漏れで2次式を作る。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Edge<'a> {
    /// 鏡面反射: 外に鏡像のノードがあるとみなす
    Reflective,
    /// 真空: 面の漏れをノードの漏れ × 面の中性子束 / ノード平均とする
    Vacuum,
    /// 隣のノード: その方向の幅と、[群] の参照の横方向の漏れ（単位体積あたりの正味の漏れ。その方向に
    /// 垂直な方向の分）
    Neighbor { width: f64, leakage: &'a [f64] },
}

/// 一般化等価理論（GET）の1ノードの問題: 均質化した定数、参照の k、参照のノード平均の中性子束と
/// 面の正味の中性子流から、均質な面の中性子束を求める。2次元（x, y）のノードで、z の漏れは 0。
#[derive(Debug, Clone, Copy)]
pub struct OneNode<'a> {
    /// [群] の拡散係数
    pub diffusion: &'a [f64],
    /// [群] の全断面積
    pub total: &'a [f64],
    /// [出る群][入る群] の散乱
    pub scatter: &'a [Vec<f64>],
    pub nu_fission: &'a [f64],
    pub chi: &'a [f64],
    pub k: f64,
    /// [群] の参照のノード平均の中性子束
    pub flux: &'a [f64],
    /// x, y の幅 (cm)
    pub widths: [f64; 2],
    /// [面（-x, +x, -y, +y）][群] の、面の平均の外向きの正味の中性子流（中性子束と同じ規格化）
    pub currents: [&'a [f64]; 4],
    /// [面（-x, +x, -y, +y）] の参照の境界か隣のノード
    pub edges: [Edge<'a>; 4],
}

/// 1ノードの問題のエラー。
#[derive(Debug, thiserror::Error)]
pub enum OneNodeError {
    #[error("{0}")]
    Invalid(String),
    #[error("the homogeneous flux on face {face} in group {group} is {flux}; it must be positive")]
    NonPositiveSurface {
        face: usize,
        group: usize,
        flux: f64,
    },
    #[error("the one-node NEM did not converge in {0} iterations")]
    NotConverged(usize),
    #[error("analytic nodal method")]
    Anm(#[from] AnmError),
    #[error("the node has no net current, so the direction of the equivalence is undefined")]
    NoNetCurrent,
    #[error(
        "no diffusion coefficient in group {group} gives equal discontinuity factors on both faces"
    )]
    NoEquivalentDiffusion { group: usize },
    #[error("the equivalent diffusion coefficients did not converge in {0} iterations")]
    EquivalenceNotConverged(usize),
}

/// 1ノードの問題の反復の上限
const MAX_ITERATIONS: usize = 10_000;

/// 1ノードの問題の均質な面の中性子束 [面（-x, +x, -y, +y）][群]（ADR-0103）。
///
/// ノード平均は参照の値に固定し（参照の k ではノードのバランス式がそれを満たす）、方向・群ごとに、core と
/// 同じノードの核を、面の正味の中性子流が参照の値になるように解く。核分裂のあるノードは NEM の式（4次の展開、
/// f1, f2 で重み付けした残差）、核分裂のないノードは2群の ANM（`anm::AnmModes`。ADR-0102）で、どちらも2次の
/// 横方向の漏れを使う。横方向の漏れは、もう一方の方向の参照の漏れで、2次式の外側は `edges` で決める。
///
/// 面の条件は、core の入射の形で β = 1、外からの入射 −J（J は参照の外向きの正味の流れ）とし、core と同じ
/// `outflow` か `AnmModes::outflow` で解く（NEM は a1〜a4 をノードの中で同時に解き、群の間の結合と真空の面の
/// 横方向の漏れだけを前の反復の値とする）。面の中性子束の相対変化が丸め誤差の水準になるまで繰り返す。
pub fn surface_fluxes(node: &OneNode) -> Result<[Vec<f64>; 4], OneNodeError> {
    let surface = homogeneous_faces(node)?;
    for (face, values) in surface.iter().enumerate() {
        if let Some((group, &value)) = values
            .iter()
            .enumerate()
            .find(|(_, v)| !(**v > 0.0 && v.is_finite()))
        {
            return Err(OneNodeError::NonPositiveSurface {
                face,
                group,
                flux: value,
            });
        }
    }
    Ok(surface)
}

/// `surface_fluxes` の均質な面の中性子束（正であることを確かめない。等価な拡散係数の根の探索で、符号の変化を見る）
pub(crate) fn homogeneous_faces(node: &OneNode) -> Result<[Vec<f64>; 4], OneNodeError> {
    validate(node)?;
    let OneNode {
        diffusion,
        total,
        scatter,
        nu_fission,
        chi,
        k,
        flux,
        widths,
        currents,
        edges,
    } = *node;
    let ng = flux.len();
    let fissile = is_fissile(nu_fission).map_err(|e| OneNodeError::Invalid(e.to_string()))?;
    // 核分裂のないノードの、方向ごとの ANM のモード（ADR-0102）
    let modes = if fissile {
        None
    } else {
        Some([
            AnmModes::new(diffusion, total, scatter, widths[0])?,
            AnmModes::new(diffusion, total, scatter, widths[1])?,
        ])
    };
    // [方向][群] の単位体積あたりの正味の漏れ
    let leakage: [Vec<f64>; 2] = std::array::from_fn(|d| {
        (0..ng)
            .map(|g| (currents[2 * d][g] + currents[2 * d + 1][g]) / widths[d])
            .collect()
    });
    // [方向][群] の展開係数
    let mut coefficients = [vec![[0.0; 4]; ng], vec![[0.0; 4]; ng]];
    for _ in 0..MAX_ITERATIONS {
        let mut change: f64 = 0.0;
        for g in 0..ng {
            let removal = total[g] - scatter[g][g];
            let coupling: Vec<f64> = (0..ng)
                .map(|from| {
                    let s = if from == g { 0.0 } else { scatter[from][g] };
                    s + chi[g] * nu_fission[from] / k
                })
                .collect();
            for d in 0..2 {
                let h = widths[d];
                // 群 `group` の横方向（z の漏れは 0）の f1, f2 の係数。真空の面では、漏れ × 面の中性子束 /
                // ノード平均（core と同じ）
                let moments_of = |group: usize| {
                    let transverse = leakage[1 - d][group];
                    let [minus, plus] = face_fluxes(&coefficients[d][group], flux[group]);
                    let beside = |edge: Edge, surface: f64| match edge {
                        Edge::Reflective => Beside::Node {
                            width: h,
                            leakage: transverse,
                        },
                        Edge::Vacuum => Beside::Face(transverse * surface / flux[group]),
                        Edge::Neighbor { width, leakage } => Beside::Node {
                            width,
                            leakage: leakage[group],
                        },
                    };
                    transverse_moments(
                        h,
                        beside(edges[2 * d], minus),
                        transverse,
                        beside(edges[2 * d + 1], plus),
                    )
                };
                // 両方の面の外向きの正味の流れ J を参照の値にする: i = o - J（β = 1、外からの入射 -J）
                let inflow_of = |group: usize| {
                    [
                        Inflow {
                            beta: 1.0,
                            external: -currents[2 * d + 1][group],
                        },
                        Inflow {
                            beta: 1.0,
                            external: -currents[2 * d][group],
                        },
                    ]
                };
                let [minus, plus] = face_fluxes(&coefficients[d][g], flux[g]);
                let line = match &modes {
                    None => outflow(&Line {
                        diffusion: diffusion[g],
                        width: h,
                        removal,
                        coupling: &coupling,
                        coefficients: &coefficients[d],
                        group: g,
                        transverse: moments_of(g),
                        inflow: inflow_of(g),
                    }),
                    // 核分裂のないノードは ANM。ほかの群も参照のノード平均と正味の流れで決まる
                    Some(modes) => {
                        let transverse: Vec<(f64, f64)> = (0..ng).map(moments_of).collect();
                        let inflow: Vec<[Inflow; 2]> = (0..ng).map(inflow_of).collect();
                        modes[d].outflow(&transverse, &inflow, flux, g)?
                    }
                };
                let updated = line.at(flux[g]).coefficients;
                for (new, old) in face_fluxes(&updated, flux[g]).iter().zip([minus, plus]) {
                    let relative = (new - old).abs() / flux[g];
                    // NaN も収束していないとみなす
                    change = if relative <= change { change } else { relative };
                }
                coefficients[d][g] = updated;
            }
        }
        if change <= CONVERGED {
            let mut surface = [vec![0.0; ng], vec![0.0; ng], vec![0.0; ng], vec![0.0; ng]];
            for d in 0..2 {
                for g in 0..ng {
                    let [minus, plus] = face_fluxes(&coefficients[d][g], flux[g]);
                    surface[2 * d][g] = minus;
                    surface[2 * d + 1][g] = plus;
                }
            }
            return Ok(surface);
        }
    }
    Err(OneNodeError::NotConverged(MAX_ITERATIONS))
}

/// 1ノードの問題の収束の判定: 面の中性子束の、ノード平均に対する変化
const CONVERGED: f64 = 1.0e-14;

/// 展開係数とノード平均から、[- 面, + 面] の中性子束
fn face_fluxes(a: &[f64; 4], average: f64) -> [f64; 2] {
    [average + 0.5 * (a[1] - a[0]), average + 0.5 * (a[1] + a[0])]
}

fn validate(node: &OneNode) -> Result<(), OneNodeError> {
    let ng = node.flux.len();
    let invalid = |what: &str| Err(OneNodeError::Invalid(what.into()));
    let lengths = [
        node.diffusion.len(),
        node.total.len(),
        node.scatter.len(),
        node.nu_fission.len(),
        node.chi.len(),
    ];
    if ng == 0
        || lengths.iter().any(|&l| l != ng)
        || node.scatter.iter().any(|row| row.len() != ng)
        || node.currents.iter().any(|c| c.len() != ng)
    {
        return invalid("one-node constants, fluxes and currents must have the same groups");
    }
    let finite = |values: &[f64]| values.iter().all(|v| v.is_finite());
    if !(finite(node.total)
        && finite(node.nu_fission)
        && finite(node.chi)
        && node.scatter.iter().all(|r| finite(r))
        && node.currents.iter().all(|c| finite(c)))
    {
        return invalid("one-node cross sections and currents must be finite");
    }
    if !node.flux.iter().all(|&f| f > 0.0 && f.is_finite()) {
        return invalid("one-node average fluxes must be positive");
    }
    if !node.diffusion.iter().all(|&d| d > 0.0 && d.is_finite()) {
        return invalid("one-node diffusion coefficients must be positive");
    }
    if !(0..ng).all(|g| node.total[g] - node.scatter[g][g] > 0.0) {
        return invalid("one-node removal cross sections must be positive");
    }
    if !(node.k > 0.0 && node.k.is_finite()) {
        return invalid("one-node k must be positive");
    }
    if !node.widths.iter().all(|&w| w > 0.0 && w.is_finite()) {
        return invalid("one-node widths must be positive");
    }
    for edge in &node.edges {
        if let Edge::Neighbor { width, leakage } = edge
            && !(*width > 0.0 && width.is_finite() && leakage.len() == ng && finite(leakage))
        {
            return invalid(
                "one-node neighbor edges need a positive width and finite leakages for every group",
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 拡散係数、全断面積、散乱、νΣf、χ
    type Constants = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>);

    /// 2群の定数（χ = (1, 0)、上方散乱あり）
    fn constants() -> Constants {
        (
            vec![1.4, 0.4],
            vec![0.24, 0.83],
            vec![vec![0.21, 0.017], vec![0.001, 0.75]],
            vec![0.006, 0.11],
            vec![1.0, 0.0],
        )
    }

    #[test]
    fn zero_currents_give_the_node_average() {
        // 正味の流れが 0 で外周が反射なら、均質な解は平らで、面の中性子束はノード平均
        let (diffusion, total, scatter, nu_fission, chi) = constants();
        let flux = [1.3, 0.4];
        let zero = [0.0, 0.0];
        let node = OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &nu_fission,
            chi: &chi,
            k: 1.05,
            flux: &flux,
            widths: [21.5, 20.0],
            currents: [&zero; 4],
            edges: [Edge::Reflective; 4],
        };
        let surface = surface_fluxes(&node).unwrap();
        for face in &surface {
            for (s, f) in face.iter().zip(&flux) {
                assert!((s - f).abs() < 1e-14 * f, "{surface:?}");
            }
        }
    }

    #[test]
    fn outward_current_lowers_the_face_flux() {
        // +x から外へ流れ出ると、均質な解は +x に向かって下がる。-x と ±y は反射
        let (diffusion, total, scatter, nu_fission, chi) = constants();
        let flux = [1.3, 0.4];
        let zero = [0.0, 0.0];
        let out = [0.01, 0.002];
        let node = OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &nu_fission,
            chi: &chi,
            k: 0.98,
            flux: &flux,
            widths: [21.5, 21.5],
            currents: [&zero, &out, &zero, &zero],
            edges: [
                Edge::Reflective,
                Edge::Vacuum,
                Edge::Reflective,
                Edge::Reflective,
            ],
        };
        let surface = surface_fluxes(&node).unwrap();
        for g in 0..2 {
            assert!(
                surface[1][g] < flux[g] && flux[g] < surface[0][g],
                "{surface:?}"
            );
            // 面の中性子束の差 φ+ - φ- は a1。正味の流れ J+ = -(D/h)(a1 + 3a2 + a3/2 + a4/5) の
            // 主な項から、a1 はおよそ -J h / D の大きさ
            let a1 = surface[1][g] - surface[0][g];
            let scale = out[g] * 21.5 / diffusion[g];
            assert!(-a1 > 0.2 * scale && -a1 < 2.0 * scale, "{a1} vs {scale}");
        }
    }

    #[test]
    fn neighbor_edges_shape_the_transverse_leakage() {
        // 隣のノードの辺（ADR-0095 の複数のノードの格子）は、core と同じく、隣の幅とその参照の横方向の
        // 漏れで2次式を作る。自分と同じ幅と漏れの隣は鏡像と同じ（鏡面反射の辺と一致する）
        let (diffusion, total, scatter, nu_fission, chi) = constants();
        let flux = [1.3, 0.4];
        let zero = [0.0, 0.0];
        // y の漏れ（x の方向の横方向の漏れ）が正: ±y の面から外へ流れ出る
        let out_y = [0.004, 0.001];
        let out_x = [0.01, 0.002];
        let width = [21.5, 21.5];
        let solve_with = |plus: Edge| {
            surface_fluxes(&OneNode {
                diffusion: &diffusion,
                total: &total,
                scatter: &scatter,
                nu_fission: &nu_fission,
                chi: &chi,
                k: 0.98,
                flux: &flux,
                widths: width,
                currents: [&zero, &out_x, &out_y, &out_y],
                edges: [Edge::Reflective, plus, Edge::Reflective, Edge::Reflective],
            })
            .unwrap()
        };
        let own: Vec<f64> = (0..2).map(|g| 2.0 * out_y[g] / width[1]).collect();
        let mirror = solve_with(Edge::Reflective);
        let same = solve_with(Edge::Neighbor {
            width: width[0],
            leakage: &own,
        });
        for (a, b) in mirror.iter().flatten().zip(same.iter().flatten()) {
            assert!((a - b).abs() < 1e-13 * b, "{mirror:?} vs {same:?}");
        }
        // 隣の横方向の漏れが大きいと、2次式は +x に向かって大きくなり、+x の面の中性子束は下がる
        let larger: Vec<f64> = own.iter().map(|l| 3.0 * l).collect();
        let raised = solve_with(Edge::Neighbor {
            width: width[0],
            leakage: &larger,
        });
        for g in 0..2 {
            assert!(raised[1][g] < mirror[1][g], "{raised:?} vs {mirror:?}");
        }
    }

    #[test]
    fn invalid_one_node_problems_are_errors() {
        let (diffusion, total, scatter, nu_fission, chi) = constants();
        let zero = [0.0, 0.0];
        let base = |flux: &'static [f64]| (flux, [21.5, 21.5]);
        for (flux, widths) in [
            base(&[0.0, 0.4]),
            base(&[f64::NAN, 0.4]),
            base(&[1.0]),
            (&[1.3, 0.4][..], [0.0, 21.5]),
        ] {
            let node = OneNode {
                diffusion: &diffusion,
                total: &total,
                scatter: &scatter,
                nu_fission: &nu_fission,
                chi: &chi,
                k: 1.0,
                flux,
                widths,
                currents: [&zero; 4],
                edges: [Edge::Reflective; 4],
            };
            assert!(
                matches!(surface_fluxes(&node), Err(OneNodeError::Invalid(_))),
                "{flux:?} {widths:?}"
            );
        }
    }

    #[test]
    fn non_positive_surface_flux_is_error() {
        // 外へ出る流れがノード平均に比べて大きすぎると、均質な面の中性子束が負になり ADF を作れない
        let (diffusion, total, scatter, nu_fission, chi) = constants();
        let flux = [1.3, 0.4];
        let zero = [0.0, 0.0];
        let out = [1.0, 0.5];
        let node = OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &nu_fission,
            chi: &chi,
            k: 0.9,
            flux: &flux,
            widths: [21.5, 21.5],
            currents: [&zero, &out, &zero, &zero],
            edges: [Edge::Reflective; 4],
        };
        assert!(matches!(
            surface_fluxes(&node),
            Err(OneNodeError::NonPositiveSurface { face: 1, .. })
        ));
    }

    #[test]
    fn fissile_nodes_have_positive_nu_fission() {
        // νΣf が正の群があれば核分裂のあるノード（NEM）、どの群も 0 なら核分裂のないノード（ANM。ADR-0102）。
        // 負や有限でない値はエラー
        assert_eq!(is_fissile(&[0.0, 0.0]), Ok(false));
        assert_eq!(is_fissile(&[0.0, 1.0e-3]), Ok(true));
        for bad in [-1.0e-3, f64::NAN, f64::INFINITY] {
            assert!(is_fissile(&[0.0, bad]).is_err(), "{bad}");
        }
        let (diffusion, total, scatter, _, chi) = constants();
        let zero = [0.0, 0.0];
        let node = OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &[0.006, -0.1],
            chi: &chi,
            k: 1.0,
            flux: &[1.3, 0.4],
            widths: [21.5, 21.5],
            currents: [&zero; 4],
            edges: [Edge::Reflective; 4],
        };
        assert!(matches!(
            surface_fluxes(&node),
            Err(OneNodeError::Invalid(_))
        ));
    }
}
