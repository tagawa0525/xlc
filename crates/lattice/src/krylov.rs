//! 再始動つきの GMRES（Saad と Schultz の、Givens の回転で最小二乗を解く形）。一般化随伴の固定源の問題（`moc::generalized`）と、
//! その粗メッシュの補正の固定源の問題（`cmfd::Coarse::solve_projected`）が使う

use rayon::prelude::*;

/// GMRES が進めなくなった（`Krylov::receive`）
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Breakdown;

/// 線形の問題 A x = b の、再始動つきの GMRES の過程。作用素 A を掛ける値を `request` で出し、掛けた値を `receive` で受ける（作用素を
/// 掛ける側が、スイープを一緒にするなど掛け方を決める）。内積は値の和
pub(crate) struct Krylov {
    /// 解の近似
    x: Vec<f64>,
    /// 右辺
    b: Vec<f64>,
    /// 止める残差のノルム（判定 × 右辺のノルム）
    target: f64,
    restart: usize,
    /// 正規直交の基底
    basis: Vec<Vec<f64>>,
    /// Givens の回転で上三角にした Hessenberg 行列の列
    columns: Vec<Vec<f64>>,
    rotations: Vec<(f64, f64)>,
    /// 最小二乗の右辺
    g: Vec<f64>,
    /// 次に掛ける値が解の近似 x（再始動の残差）か
    residual: bool,
    done: bool,
    /// 判定を満たしたときの、作用素を掛けた回数（呼ぶ側の数え方）
    sweeps: usize,
}

impl Krylov {
    /// 右辺 `b`、判定 `tolerance`（残差の 2 ノルムが右辺のノルムのこの倍で止める）、`restart` 回ごとに再始動
    pub(crate) fn new(b: Vec<f64>, tolerance: f64, restart: usize) -> Self {
        let target = tolerance * norm(&b);
        let mut k = Self {
            x: vec![0.0; b.len()],
            b,
            target,
            restart,
            basis: Vec::new(),
            columns: Vec::new(),
            rotations: Vec::new(),
            g: Vec::new(),
            residual: false,
            done: false,
            sweeps: 1,
        };
        let r = k.b.clone();
        k.begin(r);
        k
    }

    /// 残差 `r` から Arnoldi の過程を始める（判定を満たしていれば終える）
    fn begin(&mut self, r: Vec<f64>) {
        let beta = norm(&r);
        self.residual = false;
        if beta <= self.target {
            self.done = true;
            return;
        }
        self.basis = vec![divided(&r, beta)];
        self.columns.clear();
        self.rotations.clear();
        self.g = vec![beta];
    }

    /// 判定を満たしたか
    pub(crate) fn done(&self) -> bool {
        self.done
    }

    /// 解の近似
    pub(crate) fn solution(&self) -> &[f64] {
        &self.x
    }

    /// 判定を満たしたときの、作用素を掛けた回数
    pub(crate) fn steps(&self) -> usize {
        self.sweeps
    }

    /// 作用素を掛ける値（終えていれば `None`）
    pub(crate) fn request(&self) -> Option<&[f64]> {
        if self.done {
            None
        } else if self.residual {
            Some(&self.x)
        } else {
            self.basis.last().map(Vec::as_slice)
        }
    }

    /// `request` の値に作用素を掛けた値 `applied` を受ける（`sweeps` はそこまでに作用素を掛けた回数）。Givens の回転の半径が 0（作用素が
    /// 基底の値を、すでに張った空間の 0 に移す。作用素が特異）か、解の近似が有限でなければ `Breakdown`
    pub(crate) fn receive(&mut self, applied: Vec<f64>, sweeps: usize) -> Result<(), Breakdown> {
        if self.residual {
            let r = difference(&self.b, &applied);
            self.begin(r);
            self.sweeps = sweeps;
            return Ok(());
        }
        let j = self.columns.len();
        let mut w = applied;
        let mut h = Vec::with_capacity(j + 2);
        for v in &self.basis {
            let c = inner(&w, v);
            subtract_scaled(&mut w, c, v);
            h.push(c);
        }
        let next = norm(&w);
        h.push(next);
        for (i, &(c, s)) in self.rotations.iter().enumerate() {
            let (a, b) = (h[i], h[i + 1]);
            h[i] = c * a + s * b;
            h[i + 1] = -s * a + c * b;
        }
        let (a, b) = (h[j], h[j + 1]);
        let radius = a.hypot(b);
        if !(radius > 0.0 && radius.is_finite()) {
            return Err(Breakdown);
        }
        let (c, s) = (a / radius, b / radius);
        h[j] = radius;
        h.truncate(j + 1);
        self.rotations.push((c, s));
        self.g.push(-s * self.g[j]);
        self.g[j] *= c;
        self.columns.push(h);
        let residual = self.g[j + 1].abs();
        tracing::debug!(
            sweeps,
            residual,
            target = self.target,
            "generalized adjoint GMRES"
        );
        if residual <= self.target || next == 0.0 || j + 1 == self.restart {
            // 上三角の系を後退代入で解き、x に基底の和を足す
            let n = self.columns.len();
            let mut y = vec![0.0; n];
            for i in (0..n).rev() {
                let sum: f64 = ((i + 1)..n).map(|l| self.columns[l][i] * y[l]).sum();
                y[i] = (self.g[i] - sum) / self.columns[i][i];
            }
            for (v, coefficient) in self.basis.iter().zip(&y) {
                add_scaled(&mut self.x, *coefficient, v);
            }
            if self.x.iter().any(|v| !v.is_finite()) {
                return Err(Breakdown);
            }
            self.sweeps = sweeps;
            if residual <= self.target || next == 0.0 {
                self.done = true;
            } else {
                // 再始動の残差は、Arnoldi の関係から基底で組み立てられるが、作用素を 1 回掛け直した b − Ax で丸めの積み重ねを
                // 除く
                self.residual = true;
            }
            return Ok(());
        }
        self.basis.push(divided(&w, next));
        Ok(())
    }
}

fn norm(v: &[f64]) -> f64 {
    inner(v, v).sqrt()
}

/// 並列に計算する値の長さの下限。粗メッシュの補正の問題（2A で約 4,000）の内積と和は、スレッドに分けて待つ時間が計算より
/// 長い（2A の補正 1 回 69 回の GMRES で 0.23 秒、1 回あたり 3.3 ms。ADR-0247）ので、この長さより短ければ 1 つのスレッドで
/// 計算する。一般化随伴の状態（2A で約 400 万）は並列のまま
const PARALLEL: usize = 1 << 16;

fn inner(a: &[f64], b: &[f64]) -> f64 {
    if a.len() < PARALLEL {
        a.iter().zip(b).map(|(x, y)| x * y).sum()
    } else {
        a.par_iter().zip(b).map(|(x, y)| x * y).sum()
    }
}

/// `v` を `d` で割った値
fn divided(v: &[f64], d: f64) -> Vec<f64> {
    if v.len() < PARALLEL {
        v.iter().map(|x| x / d).collect()
    } else {
        v.par_iter().map(|x| x / d).collect()
    }
}

/// `a` − `b`
fn difference(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.len() < PARALLEL {
        a.iter().zip(b).map(|(x, y)| x - y).collect()
    } else {
        a.par_iter().zip(b).map(|(x, y)| x - y).collect()
    }
}

/// `x` から `v` の `factor` 倍を引く
fn subtract_scaled(x: &mut [f64], factor: f64, v: &[f64]) {
    if x.len() < PARALLEL {
        x.iter_mut().zip(v).for_each(|(a, b)| *a -= factor * b);
    } else {
        x.par_iter_mut().zip(v).for_each(|(a, b)| *a -= factor * b);
    }
}

/// `x` に `v` の `factor` 倍を足す
fn add_scaled(x: &mut [f64], factor: f64, v: &[f64]) {
    if x.len() < PARALLEL {
        x.iter_mut().zip(v).for_each(|(a, b)| *a += factor * b);
    } else {
        x.par_iter_mut().zip(v).for_each(|(a, b)| *a += factor * b);
    }
}

#[cfg(test)]
mod tests;
