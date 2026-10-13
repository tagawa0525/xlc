//! 均質化した少数群定数（1つの状態）。ファイルでは分岐のテーブル（`table`。ADR-0097）の状態として受け渡す。
//!
//! 項目: `bounds`（群境界 (eV)、降順）, `k_inf`, `total`, `absorption`, `nu_fission`, `fission`, `chi`, `diffusion`
//! （G 個ずつ）, `scatter`（[入射の群][出射の群]）, `adf`（[辺（左、右、上、下、−z、+z）][群]）, `vacuum_ratio`（[辺][群]。
//! 真空の辺だけ値を持ち、ほかの辺は nil。ADR-0101）、`normalized_adf`（[辺]。`adf` が燃料の ADF で割った値 f* の辺。
//! ADR-0095）、`height`（z の厚さ (cm)。軸方向の定数だけが持ち、ほかは nil。ADR-0128）、`pins`（ピン出力の再構成の
//! 形状関数と角の不連続因子。ないときは nil。ADR-0096）。

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
    /// [辺（左、右、上、下、−z、+z）][群] の不連続因子（辺の上の中性子束 / 均質化した中性子束）。2次元の格子計算の
    /// 定数の z の辺は 1（ADR-0128）
    pub adf: Vec<Vec<f64>>,
    /// [辺（左、右、上、下、−z、+z）][群] の、格子計算の真空の辺の外向きの正味の中性子流 / 辺の上の不均質な中性子束
    /// （辺にわたる平均の比）。真空でない辺は `None`（ADR-0101）
    pub vacuum_ratio: Vec<Option<Vec<f64>>>,
    /// [辺（左、右、上、下、−z、+z）] の、`adf` が格子計算で隣の燃料のノードの ADF で割った値
    /// f* = f_反射体 / f_燃料 である辺（反射体の定数の、燃料と接する辺。ADR-0095）。core は、その辺の ADF を
    /// f* × 隣のノードの ADF とする
    pub normalized_adf: Vec<bool>,
    /// z の厚さ (cm)。軸方向の定数（格子の x を炉心の z とした定数）は、その格子のノードの幅を持ち、core は
    /// それを置くノードの z の幅と比べる。2次元の格子計算の定数は軸方向に一様で、`None`（ADR-0128）
    pub height: Option<f64>,
    /// ピン出力の再構成の形状関数と角の不連続因子（ADR-0096）。外周がすべて鏡面反射の単一の格子でだけ作る
    pub pins: Option<Pins>,
}

/// 格子のセルごとの形状関数と、ノードの角の不連続因子（ADR-0096）。
#[derive(Debug, Clone, PartialEq)]
pub struct Pins {
    /// 列の幅 (cm)。-x の辺から順に。和はノードの幅
    pub dx: Vec<f64>,
    /// 行の高さ (cm)。-y の辺から順に。和はノードの高さ
    pub dy: Vec<f64>,
    /// [群][行（-y から）][列（-x から）] の形状関数（セルの群の核分裂率の平均 / ノードの群の核分裂率の平均）。
    /// 核分裂のないセルは 0
    pub form: Vec<Vec<Vec<f64>>>,
    /// [角（-x -y、+x -y、+x +y、-x +y の順）][群] の角の不連続因子（角の不均質な中性子束 / 角の均質な中性子束）
    pub corners: Vec<Vec<f64>>,
    /// [角] の、`corners` が角の点に接する燃料（核分裂のある）のノードの CDF の調和平均で割った値 C* である角
    /// （反射体の定数の、燃料に接する角。ADR-0131）。core は、その角の CDF を C* × 炉心で角の点に接する燃料のノードの
    /// CDF の調和平均とする
    pub normalized_corners: Vec<bool>,
    /// [辺（左、右、上、下）] の、その辺に接するセルの列（行）が、隣の部分ノードと分け合った1つのセルの半分か
    /// （集合体を分けた部分ノード。ADR-0108）
    pub halves: Vec<bool>,
}

/// 辺（左、右、上、下、−z、+z）の数
pub const SIDES: usize = 6;

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
    /// 0 < r ≤ 1 であること、f* の辺が真空の辺でないこと、z の厚さが正であることを確かめる。
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
        if self.adf.len() != SIDES {
            return invalid(format!(
                "`adf` (discontinuity factors) has {} sides, expected {SIDES}",
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
        if self.vacuum_ratio.len() != SIDES {
            return invalid(format!(
                "`vacuum_ratio` has {} sides, expected {SIDES}",
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
        if self.normalized_adf.len() != SIDES {
            return invalid(format!(
                "`normalized_adf` has {} sides, expected {SIDES}",
                self.normalized_adf.len()
            ));
        }
        // f* は燃料と接する辺の値で、外周の真空の辺は燃料と接しない
        if let Some(side) =
            (0..SIDES).find(|&s| self.normalized_adf[s] && self.vacuum_ratio[s].is_some())
        {
            return invalid(format!(
                "`normalized_adf` marks side {side}, which is a vacuum side"
            ));
        }
        if let Some(h) = self.height
            && !(h > 0.0 && h.is_finite())
        {
            return invalid(format!("`height` must be positive, got {h}"));
        }
        if !(self.k_inf > 0.0 && self.k_inf.is_finite()) {
            return invalid(format!("invalid k-inf {}", self.k_inf));
        }
        if let Some(pins) = &self.pins {
            pins.validate(ng)?;
        }
        Ok(())
    }
}

impl Pins {
    /// 列の幅と行の高さが正で、形状関数が [群][行][列] でセルの数だけあって負でなく、角の不連続因子が
    /// 4 つの角に群の数だけあって正であることを確かめる
    fn validate(&self, ng: usize) -> Result<(), FewGroupError> {
        let invalid = |message: String| Err(FewGroupError::Invalid(message));
        for (name, sizes) in [("dx", &self.dx), ("dy", &self.dy)] {
            if sizes.is_empty() || !sizes.iter().all(|s| *s > 0.0 && s.is_finite()) {
                return invalid(format!("pins: `{name}` must be positive, got {sizes:?}"));
            }
        }
        let (rows, columns) = (self.dy.len(), self.dx.len());
        if self.form.len() != ng
            || self.form.iter().any(|g| {
                g.len() != rows
                    || g.iter().any(|row| {
                        row.len() != columns || !row.iter().all(|v| *v >= 0.0 && v.is_finite())
                    })
            })
        {
            return invalid(format!(
                "pins: `form` needs {ng} groups of {rows} rows of {columns} non-negative values"
            ));
        }
        if self.corners.len() != 4
            || self
                .corners
                .iter()
                .any(|c| c.len() != ng || !c.iter().all(|v| *v > 0.0 && v.is_finite()))
        {
            return invalid(format!(
                "pins: `corners` needs 4 corners of {ng} positive values"
            ));
        }
        if self.normalized_corners.len() != 4 {
            return invalid(format!(
                "pins: `normalized_corners` has {} corners, expected 4",
                self.normalized_corners.len()
            ));
        }
        if self.halves.len() != 4 {
            return invalid(format!(
                "pins: `halves` has {} sides, expected 4",
                self.halves.len()
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
pub(crate) mod tests;
