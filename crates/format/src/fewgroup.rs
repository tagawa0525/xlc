//! 均質化した少数群定数（1つの状態）。ファイルでは分岐のテーブル（`table`。ADR-0097）の状態として受け渡す。
//!
//! 項目: `bounds`（群境界 (eV)、降順）, `k_inf`, `total`, `absorption`, `nu_fission`, `fission`, `chi`, `diffusion`
//! （G 個ずつ）, `scatter`（[入射の群][出射の群]）, `adf`（[辺（左、右、上、下）][群]）, `vacuum_ratio`（[辺][群]。
//! 真空の辺だけ値を持ち、ほかの辺は nil。ADR-0101）、`normalized_adf`（[辺]。`adf` が燃料の ADF で割った値 f* の辺。
//! ADR-0095）。

/// セル全体で均質化した少数群定数。断面積は 1/cm、拡散係数は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroup {
    /// 群境界 (eV)、降順
    pub bounds: Vec<f64>,
    pub k_inf: f64,
    pub total: Vec<f64>,
    pub absorption: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    pub chi: Vec<f64>,
    pub diffusion: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// [辺（左、右、上、下）][群] の不連続因子（辺の上の中性子束 / 均質化した中性子束）
    pub adf: Vec<Vec<f64>>,
    /// [辺（左、右、上、下）][群] の、格子計算の真空の辺の外向きの正味の中性子流 / 辺の上の不均質な中性子束
    /// （辺にわたる平均の比）。真空でない辺は `None`（ADR-0101）
    pub vacuum_ratio: Vec<Option<Vec<f64>>>,
    /// [辺（左、右、上、下）] の、`adf` が格子計算で隣の燃料のノードの ADF で割った値
    /// f* = f_反射体 / f_燃料 である辺（反射体の定数の、燃料と接する辺。ADR-0095）。core は、その辺の ADF を
    /// f* × 隣のノードの ADF とする
    pub normalized_adf: Vec<bool>,
}

/// 少数群定数の値の誤り。
#[derive(Debug, thiserror::Error)]
pub enum FewGroupError {
    #[error("{0}")]
    Invalid(String),
}

impl FewGroup {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// 群境界の順序、配列の長さ、値が有限であること、不連続因子が正であること、真空の辺の比が
    /// 0 < r ≤ 1 であること、f* の辺が真空の辺でないことを確かめる。
    /// ファイルから読んだ値にも、ファイルを経ずに作った値にも使う。
    pub fn validate(&self) -> Result<(), FewGroupError> {
        let invalid = |message: String| Err(FewGroupError::Invalid(message));
        let b = &self.bounds;
        if b.len() < 2
            || !b.windows(2).all(|w| w[0] > w[1])
            || !b.iter().all(|&x| x > 0.0 && x.is_finite())
        {
            return invalid(format!("bounds must be positive and descending, got {b:?}"));
        }
        let ng = b.len() - 1;
        let vectors = [
            ("total", &self.total),
            ("absorption", &self.absorption),
            ("nu_fission", &self.nu_fission),
            ("fission", &self.fission),
            ("chi", &self.chi),
            ("diffusion", &self.diffusion),
        ];
        for (name, v) in vectors
            .into_iter()
            .chain(self.scatter.iter().map(|row| ("scatter", row)))
        {
            if v.len() != ng {
                return invalid(format!("`{name}` has {} values for {ng} groups", v.len()));
            }
            if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
                return invalid(format!("`{name}` has a non-finite value {bad}"));
            }
        }
        if self.scatter.len() != ng {
            return invalid(format!(
                "`scatter` has {} rows for {ng} groups",
                self.scatter.len()
            ));
        }
        if self.adf.len() != 4 {
            return invalid(format!(
                "`adf` (discontinuity factors) has {} sides, expected 4",
                self.adf.len()
            ));
        }
        for side in &self.adf {
            if side.len() != ng || !side.iter().all(|v| *v > 0.0 && v.is_finite()) {
                return invalid(format!(
                    "`adf` (discontinuity factors) needs {ng} positive values per side"
                ));
            }
        }
        if self.vacuum_ratio.len() != 4 {
            return invalid(format!(
                "`vacuum_ratio` has {} sides, expected 4",
                self.vacuum_ratio.len()
            ));
        }
        // 真空の辺では入射がないので、J/φ は出ていく角度中性子束の μ の平均
        for side in self.vacuum_ratio.iter().flatten() {
            if side.len() != ng || !side.iter().all(|r| *r > 0.0 && *r <= 1.0) {
                return invalid(format!(
                    "`vacuum_ratio` needs {ng} values in (0, 1] on a vacuum side, got {side:?}"
                ));
            }
        }
        if self.normalized_adf.len() != 4 {
            return invalid(format!(
                "`normalized_adf` has {} sides, expected 4",
                self.normalized_adf.len()
            ));
        }
        // f* は燃料と接する辺の値で、外周の真空の辺は燃料と接しない
        if let Some(side) =
            (0..4).find(|&s| self.normalized_adf[s] && self.vacuum_ratio[s].is_some())
        {
            return invalid(format!(
                "`normalized_adf` marks side {side}, which is a vacuum side"
            ));
        }
        if !(self.k_inf > 0.0 && self.k_inf.is_finite()) {
            return invalid(format!("invalid k-inf {}", self.k_inf));
        }
        Ok(())
    }
}

/// テストの少数群定数（2群、真空の辺と f* の辺を1つずつ持つ）
#[cfg(test)]
pub(crate) fn sample() -> FewGroup {
    // 10 進で短く書けない値にして、ビット単位の往復を確かめる
    let x = |v: f64| v / 3.0;
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: x(3.5),
        total: vec![x(1.6), x(2.2)],
        absorption: vec![x(0.03), x(0.3)],
        nu_fission: vec![x(0.02), x(0.4)],
        fission: vec![x(0.008), x(0.16)],
        chi: vec![1.0, 0.0],
        diffusion: vec![x(4.1), x(1.2)],
        scatter: vec![vec![x(1.5), x(0.05)], vec![x(0.002), x(1.9)]],
        adf: vec![vec![x(3.0), x(3.3)]; 4],
        vacuum_ratio: vec![None, Some(vec![x(1.6), x(1.7)]), None, None],
        normalized_adf: vec![true, false, false, false],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_checks_values_in_memory() {
        // ファイルを経ずに作った値（core の Problem など）も同じ検査を通す
        assert!(sample().validate().is_ok());
        let mut fg = sample();
        fg.adf.pop();
        assert!(fg.validate().is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut fg = sample();
            fg.adf[1][0] = bad;
            let err = fg.validate().unwrap_err().to_string();
            assert!(err.contains("discontinuity"), "{err}");
        }
        let mut fg = sample();
        fg.diffusion.pop();
        assert!(fg.validate().is_err());
        let mut fg = sample();
        fg.k_inf = 0.0;
        assert!(fg.validate().is_err());
        // 群境界が降順でない、散乱の行が短い
        let mut fg = sample();
        fg.bounds.reverse();
        assert!(fg.validate().is_err());
        let mut fg = sample();
        fg.scatter[0].pop();
        assert!(fg.validate().is_err());
    }

    #[test]
    fn vacuum_ratios_are_checked() {
        // 真空の辺の J/φ は、出ていく角度中性子束の μ の平均なので 0 < r ≤ 1（ADR-0101）
        let mut fg = sample();
        fg.vacuum_ratio[1] = Some(vec![1.0, 0.5]);
        assert!(fg.validate().is_ok());
        for bad in [0.0, -0.1, 1.0 + 1e-12, f64::NAN, f64::INFINITY] {
            let mut fg = sample();
            fg.vacuum_ratio[1] = Some(vec![bad, 0.5]);
            let err = fg.validate().unwrap_err().to_string();
            assert!(err.contains("vacuum_ratio"), "{err}");
        }
        let mut fg = sample();
        fg.vacuum_ratio.pop();
        assert!(fg.validate().is_err());
        let mut fg = sample();
        fg.vacuum_ratio[0] = Some(vec![0.5]);
        assert!(fg.validate().is_err());
        // 真空の辺がなければ、どの辺も持たない
        let mut fg = sample();
        fg.vacuum_ratio = vec![None; 4];
        assert!(fg.validate().is_ok());
    }

    #[test]
    fn normalized_sides_are_checked() {
        // 反射体の定数は、燃料と接する辺の ADF を燃料の ADF で割った f* で持ち、その辺に印を付ける
        // （ADR-0095）。印は4辺ぶんで、外周の真空の辺（燃料と接しない）には付けない
        assert!(sample().validate().is_ok());
        let mut fg = sample();
        fg.normalized_adf.pop();
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains("normalized_adf"), "{err}");
        let mut fg = sample();
        fg.normalized_adf[1] = true;
        let err = fg.validate().unwrap_err().to_string();
        assert!(err.contains("normalized_adf"), "{err}");
        let mut fg = sample();
        fg.normalized_adf.push(false);
        assert!(fg.validate().is_err());
    }
}
