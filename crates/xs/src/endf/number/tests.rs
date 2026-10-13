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
