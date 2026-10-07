//! Dancoff 係数を合わせた円筒の等価なセル。

#[cfg(test)]
mod tests {
    use super::*;

    /// 燃料・ギャップ（真空）・被覆管の外半径と、その外の減速材を含むリングの全断面積（燃料の値は使わない）
    const INNER: [f64; 3] = [0.4096, 0.418, 0.475];
    const TOTALS: [f64; 4] = [0.0, 0.0, 0.3, 1.2];

    fn cell(outer: f64) -> Vec<f64> {
        let mut radii = INNER.to_vec();
        radii.push(outer);
        radii
    }

    #[test]
    fn isolated_pin_is_not_shadowed() {
        // 減速材が十分に厚ければ、隣のピンの影はない
        let c = cell_dancoff(&cell(30.0), &TOTALS).unwrap();
        assert!(c.abs() < 1e-8, "{c}");
    }

    #[test]
    fn dancoff_falls_as_the_cell_grows() {
        let c = [0.55, 0.711, 0.9].map(|r| cell_dancoff(&cell(r), &TOTALS).unwrap());
        assert!(c[0] > c[1] && c[1] > c[2] && c[2] > 0.0, "{c:?}");
    }

    #[test]
    fn equivalent_radius_recovers_the_cell() {
        let c = cell_dancoff(&cell(0.711), &TOTALS).unwrap();
        let r = equivalent_radius(&INNER, &TOTALS, c).unwrap();
        assert!((r - 0.711).abs() < 1e-9, "{r}");
    }

    #[test]
    fn rejects_unreachable_dancoff() {
        // 減速材をなくしても届かない影、0 以下、1 以上
        for c in [0.999, 0.0, -0.1, 1.0, f64::NAN] {
            assert!(equivalent_radius(&INNER, &TOTALS, c).is_err(), "{c}");
        }
        // 減速材の全断面積が 0、リングの数が合わない
        assert!(equivalent_radius(&INNER, &[0.0, 0.0, 0.3, 0.0], 0.3).is_err());
        assert!(cell_dancoff(&cell(0.711), &TOTALS[..3]).is_err());
    }
}
