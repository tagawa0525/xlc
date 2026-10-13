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
    let r = equivalent_radius(&INNER, &TOTALS, None, c).unwrap();
    assert!((r - 0.711).abs() < 1e-9, "{r}");
}

#[test]
fn rejects_unreachable_dancoff() {
    // 減速材をなくしても届かない影、0 以下、1 以上、セルの Dancoff 係数の数値誤差より小さい値
    for c in [0.999, 0.0, -0.1, 1.0, f64::NAN, 1e-12] {
        assert!(equivalent_radius(&INNER, &TOTALS, None, c).is_err(), "{c}");
    }
    // 減速材の全断面積が 0、リングの数が合わない
    assert!(equivalent_radius(&INNER, &[0.0, 0.0, 0.3, 0.0], None, 0.3).is_err());
    assert!(cell_dancoff(&cell(0.711), &TOTALS[..3]).is_err());
}

#[test]
fn outer_ring_keeps_its_area() {
    // 減速材の外に面積の決まったリング（枠）があっても、そのセルの外半径を取り戻す
    let outer = OuterRing {
        area: 0.12,
        total: 0.3,
    };
    let moderator = 0.66;
    let mut radii = cell(moderator);
    radii.push(outer.radius(moderator));
    let mut totals = TOTALS.to_vec();
    totals.push(outer.total);
    let c = cell_dancoff(&radii, &totals).unwrap();
    let r = equivalent_radius(&INNER, &TOTALS, Some(outer), c).unwrap();
    assert!((r - moderator).abs() < 1e-9, "{r}");
    // 枠の影の分だけ、枠のないセルより減速材は薄い
    let without = equivalent_radius(&INNER, &TOTALS, None, c).unwrap();
    assert!(without > r, "{without} {r}");
}

#[test]
fn rejects_bad_outer_rings() {
    let c = cell_dancoff(&cell(0.7), &TOTALS).unwrap();
    for (area, total) in [
        (0.0, 0.3),
        (-0.1, 0.3),
        (0.1, -0.3),
        (f64::NAN, 0.3),
        (0.1, f64::INFINITY),
    ] {
        let outer = Some(OuterRing { area, total });
        assert!(
            equivalent_radius(&INNER, &TOTALS, outer, c).is_err(),
            "{area} {total}"
        );
    }
}
