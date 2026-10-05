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
mod tests {
    use super::*;

    #[test]
    fn real_without_exponent_letter() {
        assert_eq!(parse_real(" 1.234567+5"), Some(1.234567e5));
        assert_eq!(parse_real("-1.234567-5"), Some(-1.234567e-5));
        assert_eq!(parse_real(" 2.000000+7"), Some(2.0e7));
        assert_eq!(parse_real("-1.00000+10"), Some(-1.0e10));
    }

    #[test]
    fn real_with_exponent_letter() {
        assert_eq!(parse_real(" 1.0000E+05"), Some(1.0e5));
        assert_eq!(parse_real("1.0e-3"), Some(1.0e-3));
    }

    #[test]
    fn real_plain() {
        assert_eq!(parse_real(" 9.991673-1"), Some(9.991673e-1));
        assert_eq!(parse_real("   0.5     "), Some(0.5));
        assert_eq!(parse_real("        125"), Some(125.0));
    }

    #[test]
    fn blank_field_is_zero() {
        assert_eq!(parse_real("           "), Some(0.0));
        assert_eq!(parse_int("           "), Some(0));
    }

    #[test]
    fn int_field() {
        assert_eq!(parse_int("        153"), Some(153));
        assert_eq!(parse_int("         -1"), Some(-1));
    }

    #[test]
    fn invalid_fields() {
        assert_eq!(parse_real(" abc       "), None);
        assert_eq!(parse_real(" 1.0+      "), None);
        assert_eq!(parse_int(" 1.5       "), None);
    }
}
