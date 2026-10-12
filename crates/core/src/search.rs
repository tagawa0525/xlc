//! 目標の値を挟む範囲を囲い込んだまま、Illinois 法で根を探す（ADR-0145）。

/// 探索で見つけた点。
#[derive(Debug, Clone, PartialEq)]
pub struct Found<T> {
    /// 根の位置
    pub x: f64,
    /// その位置の値（目標から許容の中）
    pub value: f64,
    /// 評価の回数（範囲の両端を含む）
    pub solves: usize,
    /// その位置の評価が返したもの
    pub at: T,
}

/// 探索のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SearchError<E> {
    /// 評価そのもののエラー
    #[error("{0}")]
    Solve(E),
    #[error(
        "the values {at_low} at {low} and {at_high} at {high} do not bracket the target {target}"
    )]
    NotBracketed {
        low: f64,
        high: f64,
        at_low: f64,
        at_high: f64,
        target: f64,
    },
    #[error(
        "no value within {tolerance} of the target {target} in {solves} solves; the last bracket \
         is {low} to {high} with values {at_low} and {at_high}"
    )]
    NotConverged {
        tolerance: f64,
        target: f64,
        solves: usize,
        low: f64,
        high: f64,
        at_low: f64,
        at_high: f64,
    },
}

/// `bracket` の両端で `evaluate` の値が `target` を挟むとき、|値 − `target`| ≤ `tolerance` になる点を探す。
/// 評価は `max_solves` 回まで（両端を含む）。
pub fn find<T, E>(
    bracket: [f64; 2],
    target: f64,
    tolerance: f64,
    max_solves: usize,
    evaluate: impl FnMut(f64) -> Result<(f64, T), E>,
) -> Result<Found<T>, SearchError<E>> {
    let _ = (bracket, target, tolerance, max_solves, evaluate);
    todo!()
}

#[cfg(test)]
mod tests;
