use super::*;

/// 拡散係数、全断面積、散乱、νΣf、χ
type Constants = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>);

/// 2群の定数（χ = (1, 0)、上方散乱あり）
fn constants() -> Constants {
    (
        vec![1.4, 0.4],
        vec![0.24, 0.83],
        vec![vec![0.21, 0.017], vec![0.001, 0.75]],
        vec![0.006, 0.11],
        vec![1.0, 0.0],
    )
}

#[test]
fn zero_currents_give_the_node_average() {
    // 正味の中性子流が 0 で外周が反射なら、均質な解は平らで、面の中性子束はノード平均
    let (diffusion, total, scatter, nu_fission, chi) = constants();
    let flux = [1.3, 0.4];
    let zero = [0.0, 0.0];
    let node = OneNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        nu_fission: &nu_fission,
        chi: &chi,
        k: 1.05,
        flux: &flux,
        widths: [21.5, 20.0],
        currents: [&zero; 4],
        edges: [Edge::Reflective; 4],
    };
    let surface = surface_fluxes(&node).unwrap();
    for face in &surface {
        for (s, f) in face.iter().zip(&flux) {
            assert!((s - f).abs() < 1e-14 * f, "{surface:?}");
        }
    }
}

#[test]
fn outward_current_lowers_the_face_flux() {
    // +x から外へ流れ出ると、均質な解は +x に向かって下がる。-x と ±y は反射
    let (diffusion, total, scatter, nu_fission, chi) = constants();
    let flux = [1.3, 0.4];
    let zero = [0.0, 0.0];
    let out = [0.01, 0.002];
    let node = OneNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        nu_fission: &nu_fission,
        chi: &chi,
        k: 0.98,
        flux: &flux,
        widths: [21.5, 21.5],
        currents: [&zero, &out, &zero, &zero],
        edges: [
            Edge::Reflective,
            Edge::Vacuum,
            Edge::Reflective,
            Edge::Reflective,
        ],
    };
    let surface = surface_fluxes(&node).unwrap();
    for g in 0..2 {
        assert!(
            surface[1][g] < flux[g] && flux[g] < surface[0][g],
            "{surface:?}"
        );
        // 面の中性子束の差 φ+ - φ- は a1。正味の中性子流 J+ = -(D/h)(a1 + 3a2 + a3/2 + a4/5) の
        // 主な項から、a1 はおよそ -J h / D の大きさ
        let a1 = surface[1][g] - surface[0][g];
        let scale = out[g] * 21.5 / diffusion[g];
        assert!(-a1 > 0.2 * scale && -a1 < 2.0 * scale, "{a1} vs {scale}");
    }
}

#[test]
fn neighbor_edges_shape_the_transverse_leakage() {
    // 隣のノードの辺（ADR-0095 の複数のノードの格子）は、core と同じく、隣の幅とその参照の横方向の
    // 漏れで2次式を作る。自分と同じ幅と漏れの隣は鏡像と同じ（鏡面反射の辺と一致する）
    let (diffusion, total, scatter, nu_fission, chi) = constants();
    let flux = [1.3, 0.4];
    let zero = [0.0, 0.0];
    // y の漏れ（x の方向の横方向の漏れ）が正: ±y の面から外へ流れ出る
    let out_y = [0.004, 0.001];
    let out_x = [0.01, 0.002];
    let width = [21.5, 21.5];
    let solve_with = |plus: Edge| {
        surface_fluxes(&OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &nu_fission,
            chi: &chi,
            k: 0.98,
            flux: &flux,
            widths: width,
            currents: [&zero, &out_x, &out_y, &out_y],
            edges: [Edge::Reflective, plus, Edge::Reflective, Edge::Reflective],
        })
        .unwrap()
    };
    let own: Vec<f64> = (0..2).map(|g| 2.0 * out_y[g] / width[1]).collect();
    let mirror = solve_with(Edge::Reflective);
    let same = solve_with(Edge::Neighbor {
        width: width[0],
        leakage: &own,
    });
    for (a, b) in mirror.iter().flatten().zip(same.iter().flatten()) {
        assert!((a - b).abs() < 1e-13 * b, "{mirror:?} vs {same:?}");
    }
    // 隣の横方向の漏れが大きいと、2次式は +x に向かって大きくなり、+x の面の中性子束は下がる
    let larger: Vec<f64> = own.iter().map(|l| 3.0 * l).collect();
    let raised = solve_with(Edge::Neighbor {
        width: width[0],
        leakage: &larger,
    });
    for g in 0..2 {
        assert!(raised[1][g] < mirror[1][g], "{raised:?} vs {mirror:?}");
    }
}

#[test]
fn invalid_one_node_problems_are_errors() {
    let (diffusion, total, scatter, nu_fission, chi) = constants();
    let zero = [0.0, 0.0];
    let base = |flux: &'static [f64]| (flux, [21.5, 21.5]);
    for (flux, widths) in [
        base(&[0.0, 0.4]),
        base(&[f64::NAN, 0.4]),
        base(&[1.0]),
        (&[1.3, 0.4][..], [0.0, 21.5]),
    ] {
        let node = OneNode {
            diffusion: &diffusion,
            total: &total,
            scatter: &scatter,
            nu_fission: &nu_fission,
            chi: &chi,
            k: 1.0,
            flux,
            widths,
            currents: [&zero; 4],
            edges: [Edge::Reflective; 4],
        };
        assert!(
            matches!(surface_fluxes(&node), Err(OneNodeError::Invalid(_))),
            "{flux:?} {widths:?}"
        );
    }
}

#[test]
fn non_positive_surface_flux_is_error() {
    // 外へ出る流れがノード平均に比べて大きすぎると、均質な面の中性子束が負になり ADF を作れない
    let (diffusion, total, scatter, nu_fission, chi) = constants();
    let flux = [1.3, 0.4];
    let zero = [0.0, 0.0];
    let out = [1.0, 0.5];
    let node = OneNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        nu_fission: &nu_fission,
        chi: &chi,
        k: 0.9,
        flux: &flux,
        widths: [21.5, 21.5],
        currents: [&zero, &out, &zero, &zero],
        edges: [Edge::Reflective; 4],
    };
    assert!(matches!(
        surface_fluxes(&node),
        Err(OneNodeError::NonPositiveSurface { face: 1, .. })
    ));
}

#[test]
fn fissile_nodes_have_positive_nu_fission() {
    // νΣf が正の群があれば核分裂のあるノード（NEM）、どの群も 0 なら核分裂のないノード（ANM。ADR-0102）。
    // 負や有限でない値はエラー
    assert_eq!(is_fissile(&[0.0, 0.0]), Ok(false));
    assert_eq!(is_fissile(&[0.0, 1.0e-3]), Ok(true));
    for bad in [-1.0e-3, f64::NAN, f64::INFINITY] {
        assert!(is_fissile(&[0.0, bad]).is_err(), "{bad}");
    }
    let (diffusion, total, scatter, _, chi) = constants();
    let zero = [0.0, 0.0];
    let node = OneNode {
        diffusion: &diffusion,
        total: &total,
        scatter: &scatter,
        nu_fission: &[0.006, -0.1],
        chi: &chi,
        k: 1.0,
        flux: &[1.3, 0.4],
        widths: [21.5, 21.5],
        currents: [&zero; 4],
        edges: [Edge::Reflective; 4],
    };
    assert!(matches!(
        surface_fluxes(&node),
        Err(OneNodeError::Invalid(_))
    ));
}
