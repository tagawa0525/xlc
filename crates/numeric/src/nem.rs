//! ノード展開法（NEM）の1ノードの計算。core の NEM と、lattice の ADF（一般化等価理論の1ノードの問題）が
//! 共有する（ADR-0094、ADR-0003）。
//!
//! 横方向に積分した1次元の拡散方程式を、ノード平均と4次までの多項式
//! f1 = ξ, f2 = 3ξ² - 1/4, f3 = ξ³ - ξ/4, f4 = (ξ² - 1/20)(ξ² - 1/4)（ξ ∈ [-1/2, 1/2]）で解く。
//! 3次・4次の係数は f1, f2 で重み付けした残差から、横方向の漏れは隣のノードを含む2次式で近似する。
//! 面では部分中性子流を使い、入射の部分中性子流 i を、出ていく部分中性子流 o から i = β o + 外からの入射
//! の形で与える。

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

/// 1方向の途中結果。`out` はノード平均 φ̄ = 0 のときの出ていく部分中性子流で、
/// 実際の値は `slope` × φ̄ を足したもの。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Outflow {
    /// [+ 面, - 面]
    pub out: [f64; 2],
    pub slope: [f64; 2],
    /// [+ 面, - 面] の i = β o + external の β
    pub beta: [f64; 2],
    pub external: [f64; 2],
    pub a3: f64,
    pub a4: f64,
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
        Faces {
            out,
            incoming,
            coefficients: [
                phi_r - phi_l,
                phi_r + phi_l - 2.0 * average,
                self.a3,
                self.a4,
            ],
        }
    }
}

/// 3次・4次の係数と、出ていく部分中性子流のノード平均についての線形の式 o = out + slope φ̄。
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
    // 源のモーメント（他の群・核分裂の展開係数から）
    let mut s = [0.0; 4];
    for (c, a) in coupling.iter().zip(coefficients) {
        for (si, ai) in s.iter_mut().zip(a) {
            *si += c * ai;
        }
    }
    s[0] -= l1;
    s[1] -= l2;

    let [a1, a2, _, _] = coefficients[g];
    let dh2 = diffusion / (h * h);
    // f1, f2 で重み付けした残差から 3次・4次の係数
    let a3 = (removal * a1 / 12.0 - s[0] / 12.0 + s[2] / 120.0) / (dh2 / 2.0 + removal / 120.0);
    let a4 = (removal * a2 / 20.0 - s[1] / 20.0 + s[3] / 700.0) / (dh2 / 5.0 + removal / 700.0);

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
    // φ̄ = 0 のときの値と、φ̄ に対する傾き
    let base = solve2([
        q * ext_r - w * ext_l - alpha * a3 / 2.0 - alpha * a4 / 5.0,
        q * ext_l - w * ext_r + alpha * a3 / 2.0 - alpha * a4 / 5.0,
    ]);
    let slope = solve2([6.0 * alpha, 6.0 * alpha]);
    Outflow {
        out: base,
        slope,
        beta: [beta_r, beta_l],
        external: [ext_r, ext_l],
        a3,
        a4,
    }
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
