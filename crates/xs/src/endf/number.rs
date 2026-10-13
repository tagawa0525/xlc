//! ENDF のデータ欄（11桁）の数値の解釈。

/// 実数欄を解釈する。空欄は 0 とする。
///
/// ENDF では `1.234567+5` のように `E` を省いた指数表記が使われるが、
/// `1.234567E+05` のような通常の表記も受け付ける。
pub(super) fn parse_real(field: &str) -> Option<f64> {
    let s = field.trim();
    if s.is_empty() {
        return Some(0.0);
    }
    if s.contains(['e', 'E']) {
        return s.parse().ok();
    }
    // 先頭以外にある符号は指数の符号なので、その前に `e` を補う
    let exponent_sign = s
        .char_indices()
        .skip(1)
        .find(|&(_, c)| c == '+' || c == '-')
        .map(|(i, _)| i);
    match exponent_sign {
        Some(i) => format!("{}e{}", &s[..i], &s[i..]).parse().ok(),
        None => s.parse().ok(),
    }
}

/// 整数欄を解釈する。空欄は 0 とする。
pub(super) fn parse_int(field: &str) -> Option<i64> {
    let s = field.trim();
    if s.is_empty() {
        return Some(0);
    }
    s.parse().ok()
}

#[cfg(test)]
mod tests;
