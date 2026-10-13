use super::*;
use crate::geometry::{Lattice, PinCell};

fn pin(radii: Vec<f64>) -> Tracking {
    let cell = PinCell::new(1.26, radii).unwrap();
    Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * b.abs().max(1.0)
}

#[test]
fn uniform_source_ratio_gives_that_flux() {
    let tracking = pin(vec![0.4]);
    let flux = fixed_source(
        &tracking,
        &[vec![2.0, 0.3], vec![0.5, 1.2]],
        &[vec![6.0, 0.9], vec![1.5, 3.6]],
        &FixedSourceOptions::default(),
    )
    .unwrap();
    for f in flux.iter().flatten() {
        assert!(close(*f, 3.0), "{f}");
    }
}

#[test]
fn void_regions_carry_the_surrounding_flux() {
    let tracking = pin(vec![0.4, 0.45]);
    let flux = fixed_source(
        &tracking,
        &[vec![1.0], vec![0.0], vec![1.0]],
        &[vec![2.0], vec![0.0], vec![2.0]],
        &FixedSourceOptions::default(),
    )
    .unwrap();
    for f in flux.iter().flatten() {
        assert!(close(*f, 2.0), "{f}");
    }
}

#[test]
fn rejects_untracked_regions() {
    // 外接円（0.891 cm）のすぐ内側の円の外は四隅のごく小さな領域で、間隔 0.3 cm のトラックは通らない
    let cell = PinCell::new(1.26, vec![0.89]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 4, 0.3).unwrap();
    let result = fixed_source(
        &tracking,
        &[vec![1.0], vec![1.0]],
        &[vec![1.0], vec![1.0]],
        &FixedSourceOptions::default(),
    );
    assert!(
        matches!(result, Err(FixedSourceError::Untracked(1))),
        "{result:?}"
    );
}

#[test]
fn small_groups_converge_on_their_own_scale() {
    // 黒体の燃料と、収束の速い群（Σt 大）・遅い群（Σt 小、線源はごく小さい）の減速材
    let tracking = pin(vec![0.4]);
    let both = fixed_source(
        &tracking,
        &[vec![f64::INFINITY; 2], vec![10.0, 0.05]],
        &[vec![0.0; 2], vec![10.0, 0.05e-10]],
        &FixedSourceOptions::default(),
    )
    .unwrap();
    let alone = fixed_source(
        &tracking,
        &[vec![f64::INFINITY], vec![0.05]],
        &[vec![0.0], vec![0.05]],
        &FixedSourceOptions::default(),
    )
    .unwrap();
    let (a, b) = (both[1][1], 1.0e-10 * alone[1][0]);
    assert!((a - b).abs() < 1e-8 * b, "{a} vs {b}");
}

#[test]
fn black_regions_absorb_everything() {
    let tracking = pin(vec![0.4]);
    let flux = fixed_source(
        &tracking,
        &[vec![f64::INFINITY], vec![1.0]],
        &[vec![0.0], vec![1.0]],
        &FixedSourceOptions::default(),
    )
    .unwrap();
    assert_eq!(flux[0][0], 0.0);
    assert!(flux[1][0] > 0.0 && flux[1][0] < 1.0, "{}", flux[1][0]);
}
