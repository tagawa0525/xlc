use super::*;
use xlc_format::fewgroup::Pins;

/// 2次元の格子の定数（2群）。辺ごとに違う ADF、右の辺が真空、左の辺が f*
fn radial() -> FewGroup {
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: 1.1,
        total: vec![0.5, 1.2],
        absorption: vec![0.01, 0.1],
        nu_fission: vec![0.0, 0.0],
        fission: vec![0.0, 0.0],
        chi: vec![0.0, 0.0],
        diffusion: vec![1.3, 0.4],
        scatter: vec![vec![0.47, 0.02], vec![0.0, 1.1]],
        adf: vec![
            vec![1.1, 0.9],
            vec![1.2, 0.8],
            vec![1.3, 0.7],
            vec![1.4, 0.6],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ],
        vacuum_ratio: vec![None, Some(vec![0.6, 0.5]), None, None, None, None],
        normalized_adf: vec![true, false, false, false, false, false],
        height: None,
        pins: None,
    }
}

#[test]
fn axial_constants_put_the_lattice_x_on_z() {
    // 格子の（左、右）の辺が炉心の（−z、+z）になり、径方向の4辺の ADF は 1（横に一様な平板では厳密に 1）で、
    // 真空の比と f* の印を持たない。z の厚さは格子のノードの幅（ADR-0128）
    let few = radial();
    let a = axial(&few, 6.053).unwrap();
    assert_eq!(
        a.adf,
        vec![
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.1, 0.9],
            vec![1.2, 0.8],
        ]
    );
    assert_eq!(
        a.vacuum_ratio,
        vec![None, None, None, None, None, Some(vec![0.6, 0.5])]
    );
    assert_eq!(
        a.normalized_adf,
        vec![false, false, false, false, true, false]
    );
    assert_eq!(a.height, Some(6.053));
    // 断面積はそのまま
    assert_eq!(
        (a.total.clone(), a.diffusion.clone()),
        (few.total.clone(), few.diffusion.clone())
    );
    a.validate().unwrap();
    // 形状関数は径方向のピンのもので、軸方向の定数は持たない
    let with_pins = FewGroup {
        pins: Some(Pins {
            dx: vec![1.0],
            dy: vec![1.0],
            form: vec![vec![vec![1.0]]; 2],
            corners: vec![vec![1.0, 1.0]; 4],
            normalized_corners: vec![false; 4],
            halves: vec![false; 4],
        }),
        ..radial()
    };
    let err = axial(&with_pins, 6.053).unwrap_err();
    assert!(err.contains("pins"), "{err}");
}

#[test]
fn axial_nodes_span_every_row() {
    // 平板の各ノードは格子のすべての行にわたる（行の向きは炉心の径方向で、一様）
    let node = |columns: std::ops::Range<usize>, rows: std::ops::Range<usize>| Node {
        columns,
        rows,
        reflector: false,
    };
    assert!(check_nodes(2, &[node(0..3, 0..2), node(3..5, 0..2)]).is_ok());
    let err = check_nodes(2, &[node(0..3, 0..2), node(3..5, 0..1), node(3..5, 1..2)]).unwrap_err();
    assert!(err.contains("row"), "{err}");
}

#[test]
fn axial_constants_need_a_laterally_uniform_slab() {
    // 上下の辺が真空の比や f* の印を持つ定数、z の辺の ADF が 1 でない定数は、横に一様な平板のものでないのでエラー
    // （黙って捨てない）
    let mut top_vacuum = radial();
    top_vacuum.vacuum_ratio[2] = Some(vec![0.5, 0.5]);
    let mut bottom_mark = radial();
    bottom_mark.normalized_adf[3] = true;
    let mut z_factor = radial();
    z_factor.adf[5] = vec![1.1, 1.0];
    for few in [top_vacuum, bottom_mark, z_factor] {
        assert!(axial(&few, 6.053).is_err());
    }
}
