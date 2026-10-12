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
    #[error("the value at {x} is {value}, not a finite number")]
    NotFinite { x: f64, value: f64 },
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
    mut evaluate: impl FnMut(f64) -> Result<(f64, T), E>,
) -> Result<Found<T>, SearchError<E>> {
    let [mut low, mut high] = bracket;
    let mut solves = 0;
    // 目標からの差 f と、その点の評価
    let mut at = |x: f64, solves: &mut usize| {
        *solves += 1;
        evaluate(x)
            .map(|(value, t)| (value - target, value, t))
            .map_err(SearchError::Solve)
    };
    let found = |x, value, solves, at| Found {
        x,
        value,
        solves,
        at,
    };
    let (mut f_low, value, t) = at(low, &mut solves)?;
    if f_low.abs() <= tolerance {
        return Ok(found(low, value, solves, t));
    }
    // エラーで示す端の値は評価した値（Illinois 法が半分にした f ではない）
    let not_converged = |solves, low, high, at_low, at_high| {
        Err(SearchError::NotConverged {
            tolerance,
            target,
            solves,
            low,
            high,
            at_low,
            at_high,
        })
    };
    let mut value_low = value;
    if solves >= max_solves {
        return not_converged(solves, low, high, value_low, f64::NAN);
    }
    let (mut f_high, mut value_high, t) = at(high, &mut solves)?;
    if f_high.abs() <= tolerance {
        return Ok(found(high, value_high, solves, t));
    }
    if f_low.signum() == f_high.signum() {
        return Err(SearchError::NotBracketed {
            low,
            high,
            at_low: f_low + target,
            at_high: f_high + target,
            target,
        });
    }
    // `high` は最後に評価した点、`low` は符号の違う残った端。残った端が2回続けて残れば、その f を半分にして
    // 割線を残った端に寄せる（Dowell and Jarratt, BIT 11, 168 (1971)）。f を半分にしても符号は変わらないので、
    // 新しい点はいつも2つの端の間にある
    loop {
        if solves >= max_solves {
            return not_converged(solves, low, high, value_low, value_high);
        }
        let x = high - f_high * (high - low) / (f_high - f_low);
        let (f, value, t) = at(x, &mut solves)?;
        if f.abs() <= tolerance {
            return Ok(found(x, value, solves, t));
        }
        if f.signum() == f_high.signum() {
            f_low /= 2.0;
        } else {
            low = high;
            f_low = f_high;
            value_low = value_high;
        }
        high = x;
        f_high = f;
        value_high = value;
    }
}

#[cfg(test)]
mod tests;
