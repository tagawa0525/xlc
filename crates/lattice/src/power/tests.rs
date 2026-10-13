use super::*;
use crate::geometry::PinCell;
use crate::material::FissionTerm;
use crate::moc::{Options, solve};
use crate::tracking::Tracking;

fn material(total: f64, scatter: f64, nu_fission: f64) -> Material {
    Material {
        total: vec![total],
        scatter: vec![vec![scatter]],
        anisotropic: Vec::new(),
        fission: vec![nu_fission / 2.5],
        nu_fission: vec![nu_fission],
        fission_terms: if nu_fission > 0.0 {
            vec![FissionTerm {
                chi: vec![1.0],
                nu_fission: vec![nu_fission],
            }]
        } else {
            vec![]
        },
    }
}

#[test]
fn symmetric_pins_have_equal_rates() {
    // 対角にピンを置いた 2×2 の格子は、鏡面反射で2本のピンが等価
    let pin = PinCell::new(1.26, vec![0.4]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let (fuel, moderator) = (material(0.5, 0.3, 0.3), material(1.0, 0.99, 0.0));
    let materials = [&fuel, &moderator, &moderator, &moderator, &fuel, &moderator];
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    let solution = solve(&tracking, &materials, &Options::default()).unwrap();
    let rates = cell_fission_rates(&solution, &lattice, &materials).unwrap();
    assert_eq!(rates.len(), 4);
    assert_eq!((rates[1], rates[2]), (0.0, 0.0));
    assert!(
        (rates[0] - 1.0).abs() < 1e-6 && (rates[3] - 1.0).abs() < 1e-6,
        "{rates:?}"
    );
}

#[test]
fn pin_powers_drop_margins() {
    #[derive(serde::Deserialize)]
    struct File {
        rows: Vec<Vec<f64>>,
    }
    let read = |bytes: Vec<u8>| -> Vec<Vec<f64>> {
        let file: File = rmp_serde::from_slice(&bytes).unwrap();
        file.rows
    };
    let pin = PinCell::new(1.26, vec![0.4]).unwrap();
    let lattice = Lattice::new(vec![pin], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let rows = read(to_msgpack(&[1.0, 0.0, 0.5, 1.5], &lattice).unwrap());
    assert_eq!(rows, vec![vec![1.0, 0.0], vec![0.5, 1.5]]);
    // 余白の行と列は書かない
    let pin = PinCell::new(1.26, vec![0.4]).unwrap();
    let margin = Lattice::with_margin(vec![pin], vec![vec![0]], 0.1).unwrap();
    let mut rates = vec![0.0; 9];
    rates[4] = 1.0;
    assert_eq!(read(to_msgpack(&rates, &margin).unwrap()), vec![vec![1.0]]);
}
