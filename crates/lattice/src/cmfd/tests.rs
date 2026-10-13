use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::material::{FissionTerm, Material};

fn two_by_two() -> Tracking {
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    Tracking::new(&lattice, 16, 0.05).unwrap()
}

#[test]
fn face_numbering() {
    let mesh = CoarseMesh::new(&two_by_two());
    assert_eq!((mesh.cells(), mesh.faces()), (4, 12));
    // 左上のセルの右の面は、右上のセルの左の面
    let (right, out) = mesh.cell_faces(0)[1];
    let (left, into) = mesh.cell_faces(1)[0];
    assert_eq!((right, out, into), (left, 1.0, -1.0));
    // 左上のセルの下の面は、左下のセルの上の面
    assert_eq!(mesh.cell_faces(0)[3].0, mesh.cell_faces(2)[2].0);
    // 全セルの面を数えると、内側の面は2回、外周の面は1回
    let mut count = vec![0; mesh.faces()];
    for cell in 0..mesh.cells() {
        for (face, _) in mesh.cell_faces(cell) {
            count[face] += 1;
        }
    }
    assert_eq!(count.iter().filter(|&&c| c == 2).count(), 4);
    assert_eq!(count.iter().filter(|&&c| c == 1).count(), 8);
}

#[test]
fn homogeneous_cell_gives_the_infinite_medium_k() {
    // 4辺が鏡面反射の均質なセルでは面の流れが 0 で、粗メッシュの k は無限媒質の k
    // k∞ = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
    let m = Material {
        total: vec![0.5, 1.2],
        scatter: vec![vec![0.45, 0.02], vec![0.0, 0.9]],
        anisotropic: Vec::new(),
        fission: vec![0.008 / 2.5, 0.25 / 2.5],
        nu_fission: vec![0.008, 0.25],
        fission_terms: vec![FissionTerm {
            chi: vec![1.0, 0.0],
            nu_fission: vec![0.008, 0.25],
        }],
    };
    let a = [[0.5 - 0.45, 0.0], [-0.02, 1.2 - 0.9]];
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let expected = (0.008 * a[1][1] - 0.25 * a[1][0]) / det;
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let mesh = CoarseMesh::new(&tracking);
    let volumes = tracking.volumes(2);
    // 均質なので、領域の中性子束の形によらない
    let flux = vec![vec![1.0, 0.3], vec![0.8, 0.5]];
    let current = vec![0.0; mesh.faces() * 2];
    let orbits = Orbits::identity(mesh.cells());
    let coarse = Coarse::homogenize(
        &mesh,
        &orbits,
        &[&m, &m],
        &volumes,
        &flux,
        &[0, 1],
        &current,
    )
    .unwrap();
    let (k, ratio) = coarse.solve(1.0, &CmfdOptions::default()).unwrap();
    assert!((k - expected).abs() < 1e-9, "{k} vs {expected}");
    assert_eq!(ratio.len(), 2);
}

#[test]
fn artificial_diffusion_follows_the_optically_thick_fit() {
    // Zhu ら (2016): 光学的に薄いと 0（従来の CMFD）、厚いと 0.127 に飽和し、Σt h = 4 で約 0.1
    assert_eq!(artificial_diffusion(0.5), 0.0);
    assert!((artificial_diffusion(4.0) - 0.1).abs() < 0.01);
    assert_eq!(artificial_diffusion(20.0), 0.127);
    // 多項式は飽和値に滑らかにつながる
    assert!((artificial_diffusion(14.0 - 1e-9) - 0.127).abs() < 0.002);
    let samples: Vec<f64> = (2..28)
        .map(|i| artificial_diffusion(0.5 * i as f64))
        .collect();
    assert!(samples.windows(2).all(|w| w[1] >= w[0]), "{samples:?}");
}

/// 面の外向きの重みの和（セル `cell` から出る分が正）
fn outflow(mesh: &CoarseMesh, cell: usize, faces: &[(usize, f64)]) -> f64 {
    mesh.cell_faces(cell)
        .iter()
        .map(|&(face, out)| {
            faces
                .iter()
                .filter(|&&(f, _)| f == face)
                .map(|&(_, w)| out * w)
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn paths_leave_one_cell_and_enter_the_other() {
    let mesh = CoarseMesh::new(&two_by_two());
    for from in 0..4 {
        for to in (0..4).filter(|&to| to != from) {
            let path = mesh.path(from, to).unwrap();
            let faces = path.faces();
            assert_eq!(outflow(&mesh, from, faces), 1.0);
            assert_eq!(outflow(&mesh, to, faces), -1.0);
            // 角を通るときの経由するセルは、入って出るので正味 0
            for other in (0..4).filter(|&c| c != from && c != to) {
                assert_eq!(outflow(&mesh, other, faces), 0.0);
            }
        }
    }
}

#[test]
fn cells_that_do_not_touch_have_no_path() {
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![water], vec![vec![0; 3]; 3]).unwrap();
    let mesh = CoarseMesh::new(&Tracking::new(&lattice, 8, 0.2).unwrap());
    assert!(matches!(
        mesh.path(0, 2),
        Err(CmfdError::NonAdjacent { from: 0, to: 2 })
    ));
}

#[test]
fn faces_know_their_cells() {
    let mesh = CoarseMesh::new(&two_by_two());
    // 左上のセルの左の面は外周（外向きは面の負の向き）、右の面は右上のセルとの間
    let [(left, _), (right, _), (top, _), (bottom, _)] = mesh.cell_faces(0);
    assert_eq!(
        mesh.face_cells(left),
        FaceCells::Outer {
            cell: 0,
            outward: -1.0
        }
    );
    assert_eq!(mesh.face_cells(right), FaceCells::Inner { low: 0, high: 1 });
    assert_eq!(
        mesh.face_cells(top),
        FaceCells::Outer {
            cell: 0,
            outward: -1.0
        }
    );
    assert_eq!(
        mesh.face_cells(bottom),
        FaceCells::Inner { low: 0, high: 2 }
    );
}

#[test]
fn split_track_ends_share_their_flux_between_both_cells() {
    // 端が外周の 2 つの面の境にあるトラックの角度中性子束は、経路と同じく境の両側のセルに 1/2 ずつ属する
    let tracking = two_by_two();
    let mesh = CoarseMesh::new(&tracking);
    let mut split = 0;
    for index in 0..tracking.tracks.len() {
        for (path, exit) in mesh.ends(index).into_iter().zip([false, true]) {
            let cells = mesh.end_cells(index, exit);
            if path.faces().len() == 1 {
                assert_eq!(cells.len(), 1, "track {index}");
                assert_eq!(cells[0].1, 1.0, "track {index}");
                continue;
            }
            split += 1;
            // 経路が外周の面を通るセルが、その 2 つ
            let mut outer: Vec<usize> = path
                .faces()
                .iter()
                .filter_map(|&(face, _)| match mesh.face_cells(face) {
                    FaceCells::Outer { cell, .. } => Some(cell),
                    FaceCells::Inner { .. } => None,
                })
                .collect();
            let mut found: Vec<usize> = cells.iter().map(|&(c, _)| c).collect();
            outer.sort_unstable();
            found.sort_unstable();
            assert_eq!(found, outer, "track {index}");
            assert!(
                cells.iter().all(|&(_, w)| w == 0.5),
                "track {index}: {cells:?}"
            );
        }
    }
    assert!(split > 0, "no track ends where a cell line meets the edge");
}

#[test]
fn track_ends_are_on_the_outer_faces_of_their_cells() {
    let tracking = two_by_two();
    let mesh = CoarseMesh::new(&tracking);
    let mut split = 0;
    for (index, track) in tracking.tracks.iter().enumerate() {
        let [start, end] = mesh.ends(index);
        let first = mesh.cell(track.segments[0].region);
        let last = mesh.cell(track.segments[track.segments.len() - 1].region);
        assert_eq!(outflow(&mesh, first, start.faces()), -1.0, "track {index}");
        assert_eq!(outflow(&mesh, last, end.faces()), 1.0, "track {index}");
        // 端が外周の 2 つの面の境にあって隣のセルを経由する経路でも、端のセル以外のセルの出入りは釣り合う
        for cell in 0..mesh.cells() {
            if cell != first {
                assert_eq!(outflow(&mesh, cell, start.faces()), 0.0, "track {index}");
            }
            if cell != last {
                assert_eq!(outflow(&mesh, cell, end.faces()), 0.0, "track {index}");
            }
        }
        split += usize::from(start.faces().len() > 1) + usize::from(end.faces().len() > 1);
    }
    assert!(split > 0, "no track ends where a cell line meets the edge");
}
