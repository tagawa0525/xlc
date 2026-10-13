//! 適応 Gauss-Kronrod 積分（ベクトル値の関数）。

use crate::numeric::NumericError;

/// 適応積分の区間の数の上限
const MAX_INTERVALS: usize = 10_000;

/// [a, b] でベクトル値の関数を適応 Gauss-Kronrod（7-15 点）積分する。誤差の見積もりが最も大きい
/// 区間を順に二分し、全成分で誤差の合計が積分値の `tolerance` 倍以下になったら止める。
/// 被積分関数のエラーはそのまま返し、区間の数が上限に達したら収束しないエラーにする。
pub(super) fn gauss_kronrod<E: From<NumericError>>(
    f: impl Fn(f64, &mut [f64]) -> Result<(), E>,
    n: usize,
    a: f64,
    b: f64,
    tolerance: f64,
) -> Result<Vec<f64>, E> {
    let mut intervals = vec![kronrod_interval(&f, n, a, b)?];
    loop {
        let mut total = vec![0.0; n];
        let mut error = vec![0.0; n];
        for iv in &intervals {
            for i in 0..n {
                total[i] += iv.value[i];
                error[i] += iv.error[i];
            }
        }
        let converged = (0..n).all(|i| error[i] <= tolerance * total[i].abs());
        if converged {
            return Ok(total);
        }
        if intervals.len() >= MAX_INTERVALS {
            return Err(NumericError::NotConverged(format!(
                "Gauss-Kronrod integration over [{a}, {b}] in {MAX_INTERVALS} intervals"
            ))
            .into());
        }
        // 相対誤差が最大の区間を分ける
        let score = |iv: &Interval| {
            (0..n)
                .map(|i| iv.error[i] / total[i].abs().max(f64::MIN_POSITIVE))
                .fold(0.0, f64::max)
        };
        let (worst, _) = intervals.iter().enumerate().skip(1).fold(
            (0, score(&intervals[0])),
            |best, (p, iv)| {
                // max_by と同じく、同じ値なら後の区間を選ぶ
                let s = score(iv);
                if s.total_cmp(&best.1).is_ge() {
                    (p, s)
                } else {
                    best
                }
            },
        );
        let iv = intervals.swap_remove(worst);
        let mid = 0.5 * (iv.a + iv.b);
        intervals.push(kronrod_interval(&f, n, iv.a, mid)?);
        intervals.push(kronrod_interval(&f, n, mid, iv.b)?);
    }
}

/// 1区間の Kronrod 15 点の値と、Gauss 7 点との差による誤差の見積もり。
struct Interval {
    a: f64,
    b: f64,
    value: Vec<f64>,
    error: Vec<f64>,
}

fn kronrod_interval<E>(
    f: &impl Fn(f64, &mut [f64]) -> Result<(), E>,
    n: usize,
    a: f64,
    b: f64,
) -> Result<Interval, E> {
    // Kronrod の節点（正の側）と重み、Gauss 7 点の重み（節点は奇数番目と中心）
    const XK: [f64; 8] = [
        0.991_455_371_120_812_6,
        0.949_107_912_342_758_5,
        0.864_864_423_359_769_1,
        0.741_531_185_599_394_4,
        0.586_087_235_467_691_1,
        0.405_845_151_377_397_2,
        0.207_784_955_007_898_5,
        0.0,
    ];
    const WK: [f64; 8] = [
        0.022_935_322_010_529_22,
        0.063_092_092_629_978_55,
        0.104_790_010_322_250_2,
        0.140_653_259_715_525_9,
        0.169_004_726_639_267_9,
        0.190_350_578_064_785_4,
        0.204_432_940_075_298_9,
        0.209_482_141_084_727_8,
    ];
    const WG: [f64; 4] = [
        0.129_484_966_168_869_7,
        0.279_705_391_489_276_7,
        0.381_830_050_505_118_9,
        0.417_959_183_673_469_4,
    ];
    let (center, half) = (0.5 * (a + b), 0.5 * (b - a));
    let mut kronrod = vec![0.0; n];
    let mut gauss = vec![0.0; n];
    let mut buffer = vec![0.0; n];
    let mut add = |x: f64, wk: f64, wg: f64, kronrod: &mut [f64], gauss: &mut [f64]| {
        f(x, &mut buffer)?;
        for i in 0..n {
            kronrod[i] += wk * buffer[i];
            gauss[i] += wg * buffer[i];
        }
        Ok(())
    };
    for (j, (&x, &wk)) in XK.iter().zip(&WK).enumerate() {
        // 奇数番目（j = 1, 3, 5）と中心（j = 7）は Gauss の節点でもある
        let wg = if j % 2 == 1 { WG[j / 2] } else { 0.0 };
        if x == 0.0 {
            add(center, wk, wg, &mut kronrod, &mut gauss)?;
        } else {
            add(center - half * x, wk, wg, &mut kronrod, &mut gauss)?;
            add(center + half * x, wk, wg, &mut kronrod, &mut gauss)?;
        }
    }
    Ok(Interval {
        a,
        b,
        value: kronrod.iter().map(|v| v * half).collect(),
        error: kronrod
            .iter()
            .zip(&gauss)
            .map(|(k, g)| ((k - g) * half).abs())
            .collect(),
    })
}

#[cfg(test)]
mod tests;
