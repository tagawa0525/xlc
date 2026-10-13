use super::*;

#[test]
fn pin_powers_are_written_by_rows() {
    #[derive(serde::Deserialize, Debug, PartialEq)]
    #[serde(deny_unknown_fields)]
    struct File {
        format: String,
        rows: Vec<Vec<f64>>,
    }
    let rows = vec![vec![1.0, 0.0], vec![0.5, 1.5]];
    let bytes = PinPowers { rows: rows.clone() }.to_msgpack().unwrap();
    let value: File = msgpack::from_value(msgpack::read_value(&bytes).unwrap()).unwrap();
    assert_eq!(
        value,
        File {
            format: "xlc-pin-powers".into(),
            rows,
        }
    );
}

#[test]
fn axial_pin_powers_are_written_by_layers() {
    #[derive(serde::Deserialize, Debug, PartialEq)]
    #[serde(deny_unknown_fields)]
    struct File {
        format: String,
        bounds: Vec<f64>,
        axial: Vec<f64>,
        layers: Vec<Vec<Vec<f64>>>,
    }
    let powers = AxialPinPowers {
        bounds: vec![10.0, 12.0, 15.0],
        axial: vec![0.8, 1.2],
        layers: vec![vec![vec![0.8]], vec![vec![1.2]]],
    };
    let bytes = powers.to_msgpack().unwrap();
    let value: File = msgpack::from_value(msgpack::read_value(&bytes).unwrap()).unwrap();
    assert_eq!(
        value,
        File {
            format: "xlc-axial-pin-powers".into(),
            bounds: powers.bounds,
            axial: powers.axial,
            layers: powers.layers,
        }
    );
}
