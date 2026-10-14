use super::*;
use xlc_format::group::GroupStructure;
use xlc_format::library::{Library, NuclideData, Shielded, TemperatureData};

const T: f64 = 900.0;

/// 群 0 だけが σ0 に依存する核種
fn library() -> Library {
    let shielded = |total: f64, elastic: f64| Shielded {
        total: vec![total, 20.0],
        elastic: vec![elastic, 15.0],
        fission: vec![0.0, 0.0],
        nu_fission: vec![0.0, 0.0],
        capture: vec![0.0, 0.0],
    };
    Library {
        pointwise_range: Some([4.0, 1.0e4]),
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10, 100.0, 1.0],
        nuclides: vec![NuclideData {
            name: "U-238".into(),
            za: 92238.0,
            awr: 236.0,
            resolved: None,
            thermal_scattering: false,
            temperatures: vec![TemperatureData {
                reactions: Vec::new(),
                pointwise: Default::default(),
                temperature: T,
                chi: vec![0.0, 0.0],
                scatter: vec![vec![vec![6.0, 2.0], vec![0.5, 15.0]]],
                shielded: vec![shielded(10.0, 8.0), shielded(8.0, 7.0), shielded(5.0, 5.0)],
            }],
        }],
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

#[test]
fn isolated_rod_gives_carlvik_coefficients() {
    let [a, b] = carlvik(0.0).unwrap();
    assert!(close(a.alpha, 2.0) && close(a.beta, 2.0), "{a:?}");
    assert!(close(b.alpha, 3.0) && close(b.beta, -1.0), "{b:?}");
}

#[test]
fn black_lattice_has_no_escape() {
    // C = 1 では燃料から出た中性子はすべて隣の燃料に入るので、脱出の項がなくなる
    let [a, b] = carlvik(1.0).unwrap();
    assert!(close(a.alpha, 0.0) && close(a.beta, 1.0), "{a:?}");
    assert!(close(b.alpha, 6.0) && close(b.beta, 0.0), "{b:?}");
}

#[test]
fn coefficients_are_ordered_and_betas_sum_to_one() {
    for c in [0.1, 0.5, 0.9] {
        let [a, b] = carlvik(c).unwrap();
        assert!(a.alpha < b.alpha, "{c}");
        assert!(close(a.beta + b.beta, 1.0), "{c}");
    }
}

#[test]
fn rejects_dancoff_factors_outside_zero_to_one() {
    for c in [-0.1, 1.1, f64::NAN] {
        assert!(matches!(carlvik(c), Err(ShieldingError::Dancoff(_))), "{c}");
    }
}

#[test]
fn rejects_non_positive_chord_and_density() {
    let library = library();
    let equivalence = |chord| Equivalence {
        chord,
        terms: vec![carlvik(0.0).unwrap(); 2],
    };
    let shield = |density, chord| {
        effective(
            &library,
            "U-238",
            T,
            density,
            &[0.2; 2],
            &equivalence(chord),
        )
    };
    assert!(matches!(shield(0.02, 0.0), Err(ShieldingError::Chord(_))));
    assert!(matches!(shield(0.0, 0.8), Err(ShieldingError::Density(_))));
    assert!(matches!(
        shield(f64::NAN, 0.8),
        Err(ShieldingError::Density(_))
    ));
}

#[test]
fn without_escape_the_background_alone_sets_sigma0() {
    let library = library();
    let equivalence = Equivalence {
        chord: 0.8192,
        terms: vec![carlvik(1.0).unwrap(); 2],
    };
    // σ0 = 0.2 / 0.02 = 10 b
    let result = effective(&library, "U-238", T, 0.02, &[0.2, 0.2], &equivalence).unwrap();
    let expected = library.shielded("U-238", T, &[10.0, 10.0]).unwrap();
    for (a, b) in result.total.iter().zip(&expected.total) {
        assert!(close(*a, *b), "{a} vs {b}");
    }
}

#[test]
fn two_terms_are_weighted_by_the_narrow_resonance_flux() {
    // NR 近似の燃料の中性子束 (σp + σ0) / (σt + σ0) は、吸収体自身のポテンシャル散乱 σp を線源に含む
    let potential = 11.0;
    let mut library = library();
    library.nuclides[0].resolved = Some(xlc_format::library::Resolved {
        top: 2.0e4,
        potential,
    });
    let (chord, density, background) = (0.8192, 0.02, 0.05);
    let terms = carlvik(0.0).unwrap();
    let equivalence = Equivalence {
        chord,
        terms: vec![terms; 2],
    };
    let result = effective(
        &library,
        "U-238",
        T,
        density,
        &[background; 2],
        &equivalence,
    )
    .unwrap();
    let (mut numerator, mut denominator) = (0.0, 0.0);
    for term in terms {
        let sigma0 = (background + term.alpha / chord) / density;
        let table = library.shielded("U-238", T, &[sigma0; 2]).unwrap();
        let flux = (potential + sigma0) / (table.total[0] + sigma0);
        numerator += term.beta * flux * table.total[0];
        denominator += term.beta * flux;
    }
    assert!(close(result.total[0], numerator / denominator));
    // σ0 に依存しない群はそのまま
    assert!(close(result.total[1], 20.0));
}

mod lattice {
    use super::super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::input::MaterialInput;
    use crate::material::Material;
    use crate::symmetry::{Symmetry, detect};
    use crate::tracking::Tracking;
    use xlc_format::group::GroupStructure;
    use xlc_format::library::{Library, NuclideData, Resolved, Shielded, TemperatureData};

    const T: f64 = 600.0;

    /// 共鳴のある U-238（群 0 が σ0 に依存）、核分裂する U-235、水素
    fn library() -> Library {
        let data = |name: &str, tables: Vec<Shielded>, fissile: bool| NuclideData {
            name: name.into(),
            za: 0.0,
            awr: 1.0,
            resolved: None,
            thermal_scattering: false,
            temperatures: vec![TemperatureData {
                reactions: Vec::new(),
                pointwise: Default::default(),
                temperature: T,
                chi: if fissile {
                    vec![1.0, 0.0]
                } else {
                    vec![0.0, 0.0]
                },
                scatter: vec![vec![vec![6.0, 1.0], vec![0.0, 4.0]]],
                shielded: tables,
            }],
        };
        let table = |total: [f64; 2], elastic: [f64; 2], fission: [f64; 2]| Shielded {
            total: total.to_vec(),
            elastic: elastic.to_vec(),
            fission: fission.to_vec(),
            nu_fission: fission.iter().map(|f| 2.4 * f).collect(),
            capture: fission.to_vec(),
        };
        let u238 = |t: f64, e: f64| table([t, 10.0], [e, 4.0], [0.0, 0.0]);
        let constant = |total, fission| vec![table(total, [6.0, 4.0], fission); 3];
        Library {
            pointwise_range: Some([4.0, 1.0e4]),
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
            sigma0: vec![1.0e10, 100.0, 1.0],
            nuclides: vec![
                data(
                    "U-238",
                    vec![u238(30.0, 8.0), u238(20.0, 7.5), u238(9.0, 7.0)],
                    false,
                ),
                data("U-235", constant([10.0, 50.0], [2.0, 40.0]), true),
                data("H-1", constant([20.0, 30.0], [0.0, 0.0]), false),
            ],
        }
    }

    fn materials() -> Vec<MaterialInput> {
        vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![("U-238".into(), 0.022), ("U-235".into(), 7.0e-4)],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05)],
            },
        ]
    }

    /// ピン（燃料と水）と水のセルを `layout` に並べた格子と、領域ごとの物質の名前
    fn geometry(layout: Vec<Vec<usize>>) -> (Tracking, Vec<String>) {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let names = |t: usize| -> Vec<String> {
            if t == 0 {
                vec!["fuel".into(), "water".into()]
            } else {
                vec!["water".into()]
            }
        };
        let regions = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        (Tracking::new(&lattice, 16, 0.02).unwrap(), regions)
    }

    fn fuel_total(shielded: &FuelShielding, region: usize) -> f64 {
        shielded.materials[shielded.region_material[region]].total[0]
    }

    #[test]
    fn shielding_scales_every_scattering_moment() {
        // 群 0 の弾性散乱が無限希釈の 8 から 7 に減ると、群 0 の散乱の行はどの次数も (7 - 1) / 7 倍
        let mut library = library();
        library.nuclides[0].temperatures[0]
            .scatter
            .push(vec![vec![2.0, 0.5], vec![0.0, 1.0]]);
        let table = library.shielded("U-238", T, &[1.0; 2]).unwrap();
        let nuclide = shielded_nuclide(&library, "U-238", T, &table).unwrap();
        let factor = 6.0 / 7.0;
        assert_eq!(nuclide.scatter[0], vec![6.0 * factor, factor]);
        assert_eq!(nuclide.anisotropic[0][0], vec![2.0 * factor, 0.5 * factor]);
        assert_eq!(nuclide.anisotropic[0][1], vec![0.0, 1.0]);
    }

    #[test]
    fn shields_fuel_and_keeps_moderator_dilute() {
        let library = library();
        let (tracking, regions) = geometry(vec![vec![0]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let dilute = |name: &str| {
            let m = materials().into_iter().find(|m| m.name == name).unwrap();
            Material::mix(&library, &m.composition, T).unwrap()
        };
        assert!(fuel_total(&shielded, 0) < dilute("fuel").total[0]);
        assert_eq!(
            shielded.materials[shielded.region_material[1]],
            dilute("water")
        );
    }

    #[test]
    fn reports_each_pin_with_its_nuclides() {
        let library = library();
        let (tracking, regions) = geometry(vec![vec![0, 0], vec![0, 1]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let cells: Vec<usize> = shielded.pins.iter().map(|p| p.cell).collect();
        assert_eq!(cells, vec![0, 1, 2]);
        let pin = &shielded.pins[0];
        assert_eq!(pin.dancoff.len(), 2);
        let names: Vec<&str> = pin.nuclides.iter().map(|(n, _)| n.name.as_str()).collect();
        assert_eq!(names, vec!["U-238", "U-235"]);
        // 実効断面積（弾性散乱を含む）も核種ごとに持ち、全断面積は核種の群定数と同じ
        assert_eq!(pin.effective.len(), 2);
        assert_eq!(pin.effective[0].total, pin.nuclides[0].0.total);
        assert!(pin.effective[0].elastic[0] > 0.0);
        // 物質の全断面積は、核種の全断面積と個数密度の積の和
        let material = &shielded.materials[shielded.region_material[0]];
        let sum: f64 = pin
            .nuclides
            .iter()
            .map(|(n, density)| density * n.total[0])
            .sum();
        assert!((material.total[0] - sum).abs() < 1e-12 * sum);
    }

    #[test]
    fn pins_next_to_water_are_shielded_less() {
        let library = library();
        // 領域: セル0 (燃料, 水), セル1 (燃料, 水), セル2 (燃料, 水), セル3 (水)
        let (tracking, regions) = geometry(vec![vec![0, 0], vec![0, 1]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let (far, beside, below) = (
            fuel_total(&shielded, 0),
            fuel_total(&shielded, 2),
            fuel_total(&shielded, 4),
        );
        // 水のセルの隣のピンは影が少なく（Dancoff 係数が小さく）、遮蔽が弱い
        assert!(beside > far, "{beside} vs {far}");
        assert!(
            (beside - below).abs() < 1e-3 * beside,
            "{beside} vs {below}"
        );
    }

    /// 分離共鳴を持つ（群 0 が σ0 に依存する）Zr を加えた `library()`
    fn library_with_zirconium() -> Library {
        let mut library = library();
        let mut zr = library.nuclides[2].clone();
        zr.name = "Zr-90".into();
        zr.resolved = Some(Resolved {
            top: 1.0e4,
            potential: 6.0,
        });
        zr.temperatures[0].shielded = [15.0, 10.0, 5.0]
            .iter()
            .map(|&t| Shielded {
                total: vec![t, 8.0],
                elastic: vec![6.0, 4.0],
                fission: vec![0.0, 0.0],
                nu_fission: vec![0.0, 0.0],
                capture: vec![0.0, 0.0],
            })
            .collect();
        library.nuclides.push(zr);
        library
    }

    #[test]
    fn structural_frames_are_shielded_as_lumps() {
        let library = library_with_zirconium();
        let mut materials = materials();
        materials.push(MaterialInput {
            name: "zirc".into(),
            temperature: T,
            composition: vec![("Zr-90".into(), 0.04)],
        });
        // 領域: セルごとに (燃料, 水, 枠)
        let pin = PinCell::new(1.26, vec![0.4096])
            .unwrap()
            .with_frame(0.05)
            .unwrap();
        let lattice = Lattice::new(vec![pin], vec![vec![0, 0], vec![0, 0]]).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let regions: Vec<String> = (0..4)
            .flat_map(|_| ["fuel", "water", "zirc"].map(String::from))
            .collect();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let dilute = |name: &str| {
            let m = materials.iter().find(|m| m.name == name).unwrap();
            Material::mix(&library, &m.composition, T).unwrap()
        };
        let frame = shielded.region_material[2];
        assert!(shielded.materials[frame].total[0] < dilute("zirc").total[0]);
        // 対称に移り合う枠は同じ物質を共有する
        assert!((0..4).all(|c| shielded.region_material[3 * c + 2] == frame));
        // 分離共鳴を持たない水は無限希釈のまま
        assert_eq!(
            shielded.materials[shielded.region_material[1]],
            dilute("water")
        );
    }

    /// `library_with_zirconium()` に、Zr-90 と同じ群定数で分離共鳴を持つ核種 `resonant` と持たない核種
    /// `plain` を加えたもの
    fn library_with_copies(resonant: &str, plain: &str) -> Library {
        let mut library = library_with_zirconium();
        let zr = library.nuclides[3].clone();
        let mut copy = |name: &str, resolved: bool| {
            let mut n = zr.clone();
            n.name = name.into();
            if !resolved {
                n.resolved = None;
            }
            library.nuclides.push(n);
        };
        copy(resonant, true);
        copy(plain, false);
        library
    }

    /// 枠（物質 zirc、厚さ 0.05）付きのセルの種類（半径と、内側からの物質の名前）を `layout` に並べた格子と、
    /// 領域ごとの物質の名前
    fn framed(types: &[(Vec<f64>, Vec<&str>)], layout: Vec<Vec<usize>>) -> (Tracking, Vec<String>) {
        let cells = types
            .iter()
            .map(|(radii, _)| {
                PinCell::new(1.26, radii.clone())
                    .unwrap()
                    .with_frame(0.05)
                    .unwrap()
            })
            .collect();
        let regions = layout
            .iter()
            .flatten()
            .flat_map(|&t| types[t].1.iter().chain(&["zirc"]).map(|n| n.to_string()))
            .collect();
        let lattice = Lattice::new(cells, layout).unwrap();
        (Tracking::new(&lattice, 16, 0.02).unwrap(), regions)
    }

    fn with_material(name: &str, nuclide: &str) -> Vec<MaterialInput> {
        let mut materials = materials();
        for (name, nuclide) in [("zirc", "Zr-90"), (name, nuclide)] {
            materials.push(MaterialInput {
                name: name.into(),
                temperature: T,
                composition: vec![(nuclide.into(), 0.04)],
            });
        }
        materials
    }

    #[test]
    fn reports_each_lump() {
        let library = library_with_zirconium();
        let types = [(vec![0.4096], vec!["fuel", "water"])];
        let (tracking, regions) = framed(&types, vec![vec![0, 0], vec![0, 0]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &with_material("unused", "H-1"),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // 対称に移り合う4つの枠は1つの塊。露出した面は枠の内側だけ
        assert_eq!(shielded.lumps.len(), 1, "{:?}", shielded.lumps);
        let lump = &shielded.lumps[0];
        assert_eq!((lump.cell, lump.material.as_str()), (0, "zirc"));
        assert_eq!(lump.regions, vec![2, 5, 8, 11]);
        let (p, t): (f64, f64) = (1.26, 0.05);
        let chord = (p * p - (p - 2.0 * t).powi(2)) / (p - 2.0 * t);
        assert!((lump.chord - chord).abs() < 1e-2 * chord, "{}", lump.chord);
        assert_eq!(lump.dancoff.len(), 2);
        // 実効断面積は、塊の物質の核種の全断面積
        let frame = &shielded.materials[shielded.region_material[2]];
        assert!((frame.total[0] - 0.04 * lump.effective[0].total[0]).abs() < 1e-12);
    }

    #[test]
    fn touching_pieces_are_one_lump() {
        // 管（外径 0.6）は枠の内側（半辺 0.58）に届き、同じ物質の枠と1つの塊になる
        let library = library_with_zirconium();
        let types = [
            (vec![0.4096], vec!["fuel", "water"]),
            (vec![0.5, 0.6], vec!["water", "zirc", "water"]),
        ];
        let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &with_material("unused", "H-1"),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // セル1の領域: 水, 管, 水, 枠
        assert_eq!(shielded.region_material[4], shielded.region_material[6]);
    }

    #[test]
    fn touching_materials_with_the_same_resonant_nuclides_are_one_body() {
        // 同じ核種で密度の違う物質 zirc2 を zirc が囲む。2つは1つの黒体で、Dancoff 係数と平均弦長を共有する
        let library = library_with_zirconium();
        let types = [
            (vec![0.4096], vec!["fuel", "water"]),
            (
                vec![0.3, 0.35, 0.4, 0.45],
                vec!["water", "zirc", "zirc2", "zirc", "water"],
            ),
        ];
        let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &with_material("zirc2", "Zr-90"),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let lump = |material: &str| {
            shielded
                .lumps
                .iter()
                .find(|l| l.material == material && l.cell == 1)
                .unwrap()
        };
        let (outer, inner) = (lump("zirc"), lump("zirc2"));
        assert_eq!(outer.dancoff, inner.dancoff);
        assert_eq!(outer.chord, inner.chord);
        assert_eq!(outer.regions, vec![4, 6, 10, 12]);
        assert_eq!(inner.regions, vec![5, 11]);
    }

    /// 自己遮蔽が像どうしでまとめた塊の組（`shielded.lumps`）を、組の中のセルごとの塊に分けて `lump_dancoff` で
    /// 別々に求め、どれも組の値と 1e-10 以内で一致すること（ADR-0112。対称な位置の塊の値の食い違い（#0055）が
    /// 再発すれば、組にまとめた値に隠れずにここで見つかる）。影の全断面積は `shield_lumps` と同じく、燃料はピンの
    /// 実効断面積、それ以外は無限希釈
    fn assert_lump_members_agree(
        tracking: &Tracking,
        inputs: &[MaterialInput],
        regions: &[String],
        shielded: &FuelShielding,
        groups: usize,
    ) {
        let lattice = &tracking.lattice;
        let dilute = |r: usize| -> &Material {
            let i = inputs.iter().position(|m| m.name == regions[r]).unwrap();
            &shielded.materials[i]
        };
        let fuel = |r: usize| dilute(r).nu_fission.iter().any(|&v| v > 0.0);
        let mut totals: Vec<Vec<f64>> = (0..regions.len())
            .map(|r| dilute(r).total.clone())
            .collect();
        for pin in &shielded.pins {
            let total = Material::from_nuclides(&pin.nuclides, groups)
                .unwrap()
                .total;
            for r in lattice.cell_regions(pin.cell).filter(|&r| fuel(r)) {
                totals[r].clone_from(&total);
            }
        }
        let lumped: Vec<&str> = shielded.lumps.iter().map(|l| l.material.as_str()).collect();
        let black: Vec<bool> = regions
            .iter()
            .map(|n| lumped.contains(&n.as_str()))
            .collect();
        // [領域] の、組の中の塊（セルの中で接する領域のまとまり。`shield_lumps` と同じく `Cell::adjacent` でつなぐ）の
        // 番号と、その塊の組。同じセルの中で像どうしの塊（対角線の鏡映で移り合う隅など）も別々に求める
        let mut parent: Vec<usize> = (0..regions.len()).collect();
        fn root(parent: &mut [usize], mut r: usize) -> usize {
            while parent[r] != r {
                parent[r] = parent[parent[r]];
                r = parent[r];
            }
            r
        }
        let lump_of: Vec<Option<usize>> = (0..regions.len())
            .map(|r| shielded.lumps.iter().position(|l| l.regions.contains(&r)))
            .collect();
        for cell in 0..lattice.cells() {
            let first = lattice.cell_regions(cell).start;
            for (a, b) in lattice.cell_type(cell).adjacent() {
                let (a, b) = (first + a, first + b);
                if lump_of[a].is_some() && lump_of[a] == lump_of[b] {
                    let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut member = vec![None; regions.len()];
        let mut owner = Vec::new();
        let mut numbered: Vec<(usize, usize)> = Vec::new();
        for r in 0..regions.len() {
            let Some(k) = lump_of[r] else { continue };
            let top = root(&mut parent, r);
            let m = match numbered.iter().find(|&&(t, _)| t == top) {
                Some(&(_, m)) => m,
                None => {
                    owner.push(k);
                    numbered.push((top, owner.len() - 1));
                    owner.len() - 1
                }
            };
            member[r] = Some(m);
        }
        let found = lump_dancoff(
            tracking,
            &totals,
            &black,
            &member,
            owner.len(),
            &DancoffOptions::default(),
        )
        .unwrap();
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-10 * x.abs().max(y.abs());
        for (m, (one, &k)) in found.iter().zip(&owner).enumerate() {
            let group = &shielded.lumps[k];
            assert!(
                close(one.chord, group.chord),
                "member {m} of lump {k}: chord {} vs {}",
                one.chord,
                group.chord
            );
            for (g, (x, y)) in one.dancoff.iter().zip(&group.dancoff).enumerate() {
                assert!(
                    close(*x, *y),
                    "member {m} of lump {k}, group {g}: {x} vs {y}"
                );
            }
        }
        assert!(owner.len() > shielded.lumps.len(), "no lump was grouped");
    }

    #[test]
    fn pieces_split_by_sectors_are_grouped_by_images() {
        // 外のリング（内径 0.65 > 半ピッチ 0.63）は四隅に分かれる。成り立つ操作は恒等、180° 回転、対角線と
        // 反対角線の鏡映で、隅（扇形2つずつ）はその像どうしで組になる（ADR-0112）: 格子の中心の 2 つ、外の角の
        // 2 つ、残りの 4 つ
        let library = library_with_zirconium();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let corner = PinCell::with_sectors(1.26, vec![0.65], 8).unwrap();
        let layout = vec![vec![0, 1], vec![1, 0]];
        let names = |t: usize| -> Vec<String> {
            let names: &[&str] = if t == 0 {
                &["fuel", "water"]
            } else {
                &["water"; 8]
            };
            let mut names: Vec<String> = names.iter().map(|n| n.to_string()).collect();
            if t == 1 {
                names.extend(std::iter::repeat_n("zirc".to_string(), 8));
            }
            names
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, corner], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let mut materials = materials();
        materials.push(MaterialInput {
            name: "zirc".into(),
            temperature: T,
            composition: vec![("Zr-90".into(), 0.04)],
        });
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let mut sizes: Vec<usize> = shielded.lumps.iter().map(|l| l.regions.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![4, 4, 8], "{:?}", shielded.lumps);
        assert_lump_members_agree(
            &tracking,
            &materials,
            &regions,
            &shielded,
            library.groups.len(),
        );
        // 組の領域は1つの物質を共有する
        for lump in &shielded.lumps {
            let m = shielded.region_material[lump.regions[0]];
            assert!(
                lump.regions
                    .iter()
                    .all(|&r| shielded.region_material[r] == m),
                "{lump:?}"
            );
        }
    }

    #[test]
    fn half_turned_lumps_have_equal_dancoff() {
        // 2 × 2 の格子（左上と右下が燃料、右上と左下が、四隅に分かれた Zr のリング）は 180° の回転で
        // 自分自身に移る。右上のセルの隅の塊と、その像の左下のセルの隅の塊（扇形を 4 つ回した領域）は、トラックも
        // 正確に移り合うので、Dancoff 係数と平均弦長が丸めの程度で一致する（#0055）。自己遮蔽は像の隅を1つの黒体に
        // まとめる（ADR-0112）ので、隅ごとの値は、隅を別々の塊にして `lump_dancoff` で求めて比べる
        let library = library_with_zirconium();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let corner = PinCell::with_sectors(1.26, vec![0.65], 8).unwrap();
        let layout = vec![vec![0, 1], vec![1, 0]];
        let names = |t: usize| -> Vec<String> {
            let names: &[&str] = if t == 0 {
                &["fuel", "water"]
            } else {
                &["water"; 8]
            };
            let mut names: Vec<String> = names.iter().map(|n| n.to_string()).collect();
            if t == 1 {
                names.extend(std::iter::repeat_n("zirc".to_string(), 8));
            }
            names
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, corner], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let mut materials = materials();
        materials.push(MaterialInput {
            name: "zirc".into(),
            temperature: T,
            composition: vec![("Zr-90".into(), 0.04)],
        });
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // 隅ごとの塊: Zr のセル（1、2）の外のリングの扇形の対 (0, 1)、(2, 3)、(4, 5)、(6, 7)
        let lattice = &tracking.lattice;
        let corner = |region: usize| -> Option<usize> {
            let cell = lattice.region_cell(region);
            let local = region - lattice.cell_regions(cell).start;
            (!cell.is_multiple_of(3) && local >= 8).then(|| (cell - 1) * 4 + (local % 8) / 2)
        };
        let lump: Vec<Option<usize>> = (0..lattice.regions()).map(corner).collect();
        let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
        let totals: Vec<Vec<f64>> = shielded
            .region_material
            .iter()
            .map(|&m| shielded.materials[m].total.clone())
            .collect();
        let found = lump_dancoff(
            &tracking,
            &totals,
            &black,
            &lump,
            8,
            &DancoffOptions::default(),
        )
        .unwrap();
        // 180° の回転の像: 右上のセル（1）と左下のセル（2）が入れ替わり、扇形は 4 つ回る（隅は 2 つ回る）
        let image = |l: usize| (1 - l / 4) * 4 + (l % 4 + 2) % 4;
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-10 * x.abs().max(y.abs());
        for l in 0..8 {
            let (a, b) = (&found[l], &found[image(l)]);
            assert!(
                close(a.chord, b.chord),
                "chord {l}: {} vs {}",
                a.chord,
                b.chord
            );
            for (g, (x, y)) in a.dancoff.iter().zip(&b.dancoff).enumerate() {
                assert!(close(*x, *y), "group {g}, corner {l}: {x} vs {y}");
            }
        }
    }

    #[test]
    fn other_resonant_materials_do_not_shadow_a_lump() {
        // 管の物質が分離共鳴を持つかどうかで、同じセルの枠（別の物質）の遮蔽は変わらない（群定数は同じ）
        let library = library_with_copies("Fe-56", "Fe-56x");
        let types = [
            (vec![0.4096], vec!["fuel", "water"]),
            (vec![0.45, 0.5], vec!["water", "tube", "water"]),
        ];
        let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
        let frame = |nuclide: &str| {
            let shielded = shield_fuel(
                &library,
                &tracking,
                &with_material("tube", nuclide),
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            // セル1の領域: 水, 管, 水, 枠
            shielded.materials[shielded.region_material[6]].total[0]
        };
        let (resonant, plain) = (frame("Fe-56"), frame("Fe-56x"));
        assert!(
            (resonant - plain).abs() < 1e-12 * plain,
            "{resonant} vs {plain}"
        );
    }

    #[test]
    fn separate_pieces_of_a_material_are_separate_lumps() {
        // 管と枠が同じ物質でも、離れた領域は別の塊として形に合った平均弦長で遮蔽する
        let library = library_with_zirconium();
        let types = [
            (vec![0.4096], vec!["fuel", "water"]),
            (vec![0.45, 0.5], vec!["water", "zirc", "water"]),
        ];
        let (tracking, regions) = framed(&types, vec![vec![0, 1], vec![1, 0]]);
        let shielded = shield_fuel(
            &library,
            &tracking,
            &with_material("unused", "H-1"),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // セル1の領域: 水, 管, 水, 枠
        let (tube, frame) = (shielded.region_material[4], shielded.region_material[6]);
        assert_ne!(tube, frame);
        let total = |m: usize| shielded.materials[m].total[0];
        assert!(total(tube) != total(frame), "{}", total(tube));
    }

    #[test]
    fn void_gap_is_kept_void() {
        let library = library();
        let pin = PinCell::new(1.26, vec![0.4096, 0.418]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let regions: Vec<String> = ["fuel", crate::input::VOID, "water"]
            .map(String::from)
            .to_vec();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        assert!(shielded.materials[shielded.region_material[1]].is_void());
        assert!(fuel_total(&shielded, 0) > 0.0);
    }

    #[test]
    fn nuclides_without_density_stay_dilute() {
        let library = library();
        let (tracking, regions) = geometry(vec![vec![0]]);
        let shield = |materials: &[MaterialInput]| {
            let s = shield_fuel(
                &library,
                &tracking,
                materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap();
            s.materials[s.region_material[0]].clone()
        };
        let mut with_zero = materials();
        with_zero[0].composition.push(("H-1".into(), 0.0));
        assert_eq!(shield(&with_zero), shield(&materials()));
    }

    #[test]
    fn rejects_non_finite_shielded_cross_sections() {
        let mut library = library();
        library.nuclides[0].temperatures[0].shielded[1].total[0] = f64::NAN;
        let (tracking, regions) = geometry(vec![vec![0]]);
        let result = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        );
        assert!(
            matches!(result, Err(SelfShieldingError::NonFinite(0))),
            "{result:?}"
        );
    }

    /// 分離共鳴（200 eV の共鳴）を持つ U-238 と、ポイントワイズの断面積を持つ核種のライブラリ。群 1, 2 が UFG の範囲
    fn resonant_library() -> Library {
        let bounds = vec![1.0e7, 1.0e4, 1.0e2, 4.0, 0.625, 1.0e-5];
        let ng = bounds.len() - 1;
        let resonance = |e: f64| 1.0e4 / (1.0 + ((e - 200.0) / 0.5).powi(2));
        let mut energies: Vec<f64> = (0..=400).map(|k| 180.0 + 0.1 * k as f64).collect();
        energies.splice(0..0, [4.0, 100.0]);
        energies.extend([1.0e3, 1.0e4]);
        // (名前, awr, ポイントワイズの断面積 [全, 弾性散乱, 核分裂], 分離共鳴の上端)
        let nuclide = |name: &str, awr: f64, points: Vec<[f64; 4]>, top: Option<f64>| {
            let column = |k: usize| points.iter().map(|p| p[k]).collect::<Vec<_>>();
            let (total, elastic, fission) = (points[0][1], points[0][2], points[0][3]);
            // σ0 の表: 無限希釈から下がる（Carlvik で引く範囲の外の群に使う）
            let table = |factor: f64| Shielded {
                total: vec![total * factor; ng],
                elastic: vec![elastic; ng],
                fission: vec![fission * factor; ng],
                nu_fission: vec![2.4 * fission * factor; ng],
                // 捕獲は吸収から核分裂を除いたもの
                capture: vec![(total - elastic - fission) * factor; ng],
            };
            NuclideData {
                name: name.into(),
                za: 0.0,
                awr,
                resolved: top.map(|top| xlc_format::library::Resolved {
                    top,
                    potential: 10.0,
                }),
                thermal_scattering: false,
                temperatures: vec![TemperatureData {
                    reactions: Vec::new(),
                    temperature: T,
                    chi: (0..ng)
                        .map(|g| if fission > 0.0 && g == 0 { 1.0 } else { 0.0 })
                        .collect(),
                    // 弾性散乱は自群に 9 割、次の群に 1 割（行の和が弾性散乱。最後の群は自群だけ）
                    scatter: vec![
                        (0..ng)
                            .map(|g| {
                                (0..ng)
                                    .map(|h| match (h == g, h == g + 1, g + 1 == ng) {
                                        (true, _, true) => elastic,
                                        (true, _, false) => 0.9 * elastic,
                                        (_, true, _) => 0.1 * elastic,
                                        _ => 0.0,
                                    })
                                    .collect()
                            })
                            .collect(),
                    ],
                    shielded: vec![table(1.0), table(0.9), table(0.8)],
                    pointwise: xlc_format::library::Pointwise {
                        energy: column(0),
                        total: column(1),
                        elastic: column(2),
                        fission: column(3),
                        nu_fission: column(3).iter().map(|f| 2.4 * f).collect(),
                        capture: (0..points.len())
                            .map(|i| points[i][1] - points[i][2] - points[i][3])
                            .collect(),
                    },
                }],
            }
        };
        let flat = |t: f64, e: f64, f: f64| vec![[4.0, t, e, f], [1.0e4, t, e, f]];
        let u238 = energies
            .iter()
            .map(|&e| [e, 11.0 + resonance(e), 11.0, 0.0])
            .collect();
        Library {
            pointwise_range: Some([4.0, 1.0e4]),
            groups: GroupStructure::new(bounds).unwrap(),
            sigma0: vec![1.0e10, 100.0, 1.0],
            nuclides: vec![
                nuclide("U-238", 236.0, u238, Some(1.0e4)),
                nuclide("U-235", 233.0, flat(10.0, 8.0, 1.0), Some(2.0e3)),
                nuclide("O-16", 15.86, flat(3.8, 3.8, 0.0), None),
                nuclide("Zr-90", 89.0, flat(6.5, 6.5, 0.0), None),
                nuclide("H-1", 0.99917, flat(20.0, 20.0, 0.0), None),
            ],
        }
    }

    /// 被覆管が物質 `clad`、枠（厚さ 0.05）が分離共鳴を持つ Zr-91 の物質 grid の、UFG で遮蔽するピンセル。
    /// 枠の塊の Dancoff 係数と、`black` の領域を黒体にして求め直した Dancoff 係数（燃料はピンの SPH 因子を
    /// 掛けない実効断面積、ほかは無限希釈の全断面積）
    fn frame_dancoff(clad: &str, black: &[usize]) -> (Vec<f64>, Vec<f64>) {
        let mut library = resonant_library();
        let mut zr91 = library.nuclides[3].clone();
        zr91.name = "Zr-91".into();
        zr91.resolved = Some(xlc_format::library::Resolved {
            top: 1.0e3,
            potential: 6.5,
        });
        library.nuclides.push(zr91);
        let material = |name: &str, composition: &[(&str, f64)]| MaterialInput {
            name: name.into(),
            temperature: T,
            composition: composition.iter().map(|&(n, d)| (n.into(), d)).collect(),
        };
        let materials = vec![
            material(
                "fuel",
                &[("U-238", 0.022), ("U-235", 5.0e-4), ("O-16", 0.045)],
            ),
            material("clad", &[("Zr-90", 0.043)]),
            material("grid", &[("Zr-91", 0.043)]),
            material("grid2", &[("Zr-91", 0.03)]),
            material("water", &[("H-1", 0.05), ("O-16", 0.025)]),
        ];
        let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475])
            .unwrap()
            .with_frame(0.05)
            .unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let regions: Vec<String> = [
            "fuel",
            "fuel",
            "fuel",
            crate::input::VOID,
            clad,
            "water",
            "grid",
        ]
        .map(String::from)
        .to_vec();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let pin = &shielded.pins[0];
        assert!(pin.sph.is_some());
        let fuel = Material::from_nuclides(&pin.nuclides, library.groups.len()).unwrap();
        let totals: Vec<Vec<f64>> = regions
            .iter()
            .map(|name| match name.as_str() {
                "fuel" => fuel.total.clone(),
                crate::input::VOID => vec![0.0; library.groups.len()],
                name => {
                    let m = materials.iter().find(|m| m.name == name).unwrap();
                    Material::mix(&library, &m.composition, T).unwrap().total
                }
            })
            .collect();
        let lump: Vec<Option<usize>> = (0..regions.len())
            .map(|r| black.contains(&r).then_some(usize::from(r != 6)))
            .collect();
        let lumps = 1 + usize::from(black.len() > 1);
        let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
        let expected = lump_dancoff(
            &tracking,
            &totals,
            &black,
            &lump,
            lumps,
            &DancoffOptions::default(),
        )
        .unwrap();
        assert_eq!(shielded.lumps.len(), 1, "{:?}", shielded.lumps);
        (
            shielded.lumps[0].dancoff.clone(),
            expected[0].dancoff.clone(),
        )
    }

    #[test]
    fn lumps_see_the_physical_fuel() {
        // 燃料は SPH 因子を掛けない実効断面積で、塊の影の計算に入る
        let (dancoff, expected) = frame_dancoff("clad", &[6]);
        for (c, e) in dancoff.iter().zip(&expected) {
            assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
        }
    }

    #[test]
    fn same_material_in_the_pin_shadows_a_lump() {
        // 被覆管が枠と同じ物質なら、等価なセルで遮蔽した被覆管も、枠と同じ共鳴の影を作る（黒体）
        let (dancoff, expected) = frame_dancoff("grid", &[4, 6]);
        for (c, e) in dancoff.iter().zip(&expected) {
            assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
        }
    }

    #[test]
    fn materials_with_the_same_resonant_nuclides_shadow_a_lump() {
        // 被覆管が枠と別の物質でも、分離共鳴を持つ核種が同じなら同じ共鳴の影を作る（黒体）
        let (dancoff, expected) = frame_dancoff("grid2", &[4, 6]);
        for (c, e) in dancoff.iter().zip(&expected) {
            assert!((c - e).abs() < 1e-9, "{dancoff:?} vs {expected:?}");
        }
    }

    #[test]
    fn pieces_of_a_lump_follow_contact() {
        // 小さな正方形のセルの縦の列と、扇形の角度の 0 と 2π をまたぐ扇形は、それぞれ1つの塊
        let library = library_with_zirconium();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let grid = PinCell::grid(1.26, 3).unwrap();
        let sectored = PinCell::with_sectors(1.26, vec![0.5, 0.6], 4).unwrap();
        let layout = vec![vec![0, 1, 0], vec![2, 0, 0], vec![0, 0, 0]];
        let names = |t: usize| -> Vec<String> {
            let names = match t {
                0 => "fuel water",
                1 => "water zirc water water zirc water water zirc water",
                _ => "water water water water zirc water water zirc water water water water",
            };
            names.split(' ').map(String::from).collect()
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, grid, sectored], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &with_material("unused", "H-1"),
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let m = &shielded.region_material;
        // 正方形のセル（セル1）は領域 2〜10、扇形のセル（セル3）は 13〜24（リング 1 の扇形は 17〜20）
        assert!(m[3] == m[6] && m[6] == m[9], "{m:?}");
        assert_eq!(m[17], m[20], "{m:?}");
        assert_eq!(shielded.lumps.len(), 2, "{:?}", shielded.lumps);
    }

    #[test]
    fn fuel_rings_are_shielded_by_the_cell() {
        let library = resonant_library();
        let materials = vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![
                    ("U-238".into(), 0.022),
                    ("U-235".into(), 5.0e-4),
                    ("O-16".into(), 0.045),
                ],
            },
            MaterialInput {
                name: "clad".into(),
                temperature: T,
                composition: vec![("Zr-90".into(), 0.043)],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
            },
        ];
        let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let regions: Vec<String> = ["fuel", "fuel", "fuel", crate::input::VOID, "clad", "water"]
            .map(String::from)
            .to_vec();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let total =
            |region: usize, g: usize| shielded.materials[shielded.region_material[region]].total[g];
        // UFG の範囲の群（100 eV〜10 keV）では、燃料の内側ほど共鳴が遮蔽される
        assert!(
            total(0, 1) < total(1, 1) && total(1, 1) < total(2, 1),
            "{:?}",
            [total(0, 1), total(1, 1), total(2, 1)]
        );
        // 範囲の外の群は、燃料のリングによらず等価原理の値
        assert_eq!(total(0, 0), total(2, 0));
        // 等価なセルの外径は、正方形のセルと面積の等しい円に近い
        let radius = shielded.pins[0].equivalent.unwrap();
        assert!((0.6..0.8).contains(&radius), "{radius}");
        assert!(shielded.materials[shielded.region_material[3]].is_void());
        // SPH 因子は等価なセルのリング（燃料3、ギャップ、被覆管、減速材）と UFG の範囲の群（1, 2）ごと。
        // 真空のギャップは 1、燃料では 1 と違う
        let sph = &shielded.pins[0].sph.as_ref().unwrap().factors;
        assert_eq!(sph.len(), 6);
        assert!(sph.iter().all(|row| row.len() == 2));
        assert_eq!(sph[3], vec![1.0, 1.0]);
        assert!(sph[0].iter().any(|mu| (mu - 1.0).abs() > 1e-6), "{sph:?}");
    }

    /// 3 つの燃料のリングの物質を内側から順に置いた、UFG で遮蔽するピンセル。物質は `rings` の名前と U-235 の数密度
    /// の組（核種の並びは U-238、U-235、O-16）と、被覆管と水。飛程は格子点ごとに求める（飛程の表は物質ごとの軸の
    /// 直積なので、燃料のリングの物質が 3 つあると大きすぎる）
    fn ring_fuel(rings: [(&str, f64); 3]) -> (Library, FuelShielding, Vec<MaterialInput>) {
        let library = resonant_library();
        let fuel = |name: &str, u235: f64| MaterialInput {
            name: name.into(),
            temperature: T,
            composition: vec![
                ("U-238".into(), 0.022),
                ("U-235".into(), u235),
                ("O-16".into(), 0.045),
            ],
        };
        let mut materials: Vec<MaterialInput> = Vec::new();
        for (name, u235) in rings {
            if !materials.iter().any(|m| m.name == name) {
                materials.push(fuel(name, u235));
            }
        }
        materials.push(MaterialInput {
            name: "clad".into(),
            temperature: T,
            composition: vec![("Zr-90".into(), 0.043)],
        });
        materials.push(MaterialInput {
            name: "water".into(),
            temperature: T,
            composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
        });
        let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let regions: Vec<String> = [
            rings[0].0,
            rings[1].0,
            rings[2].0,
            crate::input::VOID,
            "clad",
            "water",
        ]
        .map(String::from)
        .to_vec();
        let options = ShieldingOptions {
            collision: crate::ufg::Collision::Exact,
            ..ShieldingOptions::default()
        };
        let shielded = shield_fuel(&library, &tracking, &materials, &regions, &options).unwrap();
        (library, shielded, materials)
    }

    #[test]
    fn fuel_rings_of_equal_composition_shield_like_one_material() {
        let (_, one, _) = ring_fuel([("fuel", 5.0e-4), ("fuel", 5.0e-4), ("fuel", 5.0e-4)]);
        let (_, split, _) = ring_fuel([("a", 5.0e-4), ("b", 5.0e-4), ("c", 5.0e-4)]);
        for region in 0..6 {
            let (x, y) = (
                &one.materials[one.region_material[region]],
                &split.materials[split.region_material[region]],
            );
            for (a, b) in x.total.iter().zip(&y.total) {
                assert!(
                    (a - b).abs() <= 1e-12 * a.abs(),
                    "region {region}: {a} vs {b}"
                );
            }
        }
        assert_eq!(one.pins[0].rings, split.pins[0].rings);
    }

    #[test]
    fn fuel_rings_take_their_own_compositions() {
        let (_, shielded, materials) = ring_fuel([("a", 2.0e-4), ("b", 5.0e-4), ("c", 9.0e-4)]);
        let pin = &shielded.pins[0];
        assert_eq!(pin.rings.len(), 3);
        // リングの物質の巨視的断面積は、そのリングの数密度とそのリングの微視的断面積の積の和（SPH 因子のない範囲の外の群）
        for (ring, material) in materials.iter().take(3).enumerate() {
            let total = shielded.materials[shielded.region_material[ring]].total[0];
            let expected: f64 = material
                .composition
                .iter()
                .zip(&pin.rings[ring])
                .map(|((_, n), s)| n * s.total[0])
                .sum();
            assert!(
                (total - expected).abs() <= 1e-12 * expected,
                "ring {ring}: {total} vs {expected}"
            );
        }
        // 等価原理（範囲の外）は体積平均の組成で解くので、どのリングも同じ微視的断面積
        assert_eq!(pin.rings[0][0].total[0], pin.rings[2][0].total[0]);
        // ピンの核種の数密度は、リングの体積の重みの平均（体積は 0.2²、0.35² − 0.2²、0.4096² − 0.35² に比例）
        let volumes = [
            0.04,
            0.35f64.powi(2) - 0.04,
            0.4096f64.powi(2) - 0.35f64.powi(2),
        ];
        let average = (volumes[0] * 2.0e-4 + volumes[1] * 5.0e-4 + volumes[2] * 9.0e-4)
            / volumes.iter().sum::<f64>();
        let u235 = pin.nuclides[1].1;
        assert!(
            (u235 - average).abs() <= 1e-9 * average,
            "{u235} vs {average}"
        );
    }

    #[test]
    fn ring_cross_sections_resolve_the_rim_in_the_ufg_range() {
        let (_, shielded, _) = ring_fuel([("fuel", 5.0e-4), ("fuel", 5.0e-4), ("fuel", 5.0e-4)]);
        let rings = &shielded.pins[0].rings;
        // U-238 の 200 eV の共鳴（群 2）は内側ほど遮蔽され、範囲の外の群はリングによらない
        assert!(rings[0][0].total[1] < rings[1][0].total[1]);
        assert!(rings[1][0].total[1] < rings[2][0].total[1]);
        assert_eq!(rings[0][0].total[0], rings[2][0].total[0]);
    }

    #[test]
    fn rejects_fuel_rings_of_different_nuclides() {
        let library = resonant_library();
        let material = |name: &str, composition: &[(&str, f64)]| MaterialInput {
            name: name.into(),
            temperature: T,
            composition: composition.iter().map(|&(n, d)| (n.into(), d)).collect(),
        };
        let materials = vec![
            material("a", &[("U-238", 0.022), ("U-235", 5.0e-4), ("O-16", 0.045)]),
            material("b", &[("U-235", 5.0e-4), ("U-238", 0.022), ("O-16", 0.045)]),
            material("water", &[("H-1", 0.05), ("O-16", 0.025)]),
        ];
        let pin = PinCell::new(1.26, vec![0.2, 0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let regions: Vec<String> = ["a", "b", "water"].map(String::from).to_vec();
        let result = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        );
        assert!(
            matches!(result, Err(SelfShieldingError::MixedFuel(0))),
            "{result:?}"
        );
    }

    /// 燃料のピンに囲まれた、分離共鳴を持つ U-238 だけの（核分裂しない）棒（中央のセル 4 の内側の 2 リング）の
    /// 3 × 3 の格子を遮蔽する
    fn rod_lattice() -> (Library, Vec<MaterialInput>, FuelShielding) {
        let library = resonant_library();
        let material = |name: &str, composition: &[(&str, f64)]| MaterialInput {
            name: name.into(),
            temperature: T,
            composition: composition.iter().map(|&(n, d)| (n.into(), d)).collect(),
        };
        let materials = vec![
            material(
                "fuel",
                &[("U-238", 0.022), ("U-235", 5.0e-4), ("O-16", 0.045)],
            ),
            material("absorber", &[("U-238", 0.02)]),
            material("clad", &[("Zr-90", 0.043)]),
            material("water", &[("H-1", 0.05), ("O-16", 0.025)]),
        ];
        let pin = PinCell::new(1.26, vec![0.4096, 0.475]).unwrap();
        let rod = PinCell::new(1.26, vec![0.2, 0.38, 0.48]).unwrap();
        let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
        let names = |t: usize| -> Vec<String> {
            let names = match t {
                0 => "fuel clad water",
                _ => "absorber absorber clad water",
            };
            names.split(' ').map(String::from).collect()
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, rod], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        (library, materials, shielded)
    }

    #[test]
    fn absorber_rods_are_shielded_by_slowing_down() {
        // UFG の範囲の群（1, 2）の棒の実効断面積は、塊の Dancoff 係数を合わせた等価なセルの減速計算を、棒の
        // リングの中性子束と体積で平均したもの
        let (library, materials, shielded) = rod_lattice();
        let lump = shielded
            .lumps
            .iter()
            .find(|l| l.material == "absorber")
            .unwrap();
        // 棒のセル（セル 4）の領域は 12〜15 で、棒は 12 と 13
        assert_eq!(lump.regions, vec![12, 13]);
        let dilute = |name: &str| {
            let m = materials.iter().find(|m| m.name == name).unwrap();
            Material::mix(&library, &m.composition, T).unwrap().total
        };
        let (clad, water) = (dilute("clad"), dilute("water"));
        let ufg = [1, 2];
        let radius = ufg
            .iter()
            .map(|&g| {
                crate::equivalent::equivalent_radius(
                    &[0.38, 0.48],
                    &[0.0, clad[g], water[g]],
                    None,
                    lump.dancoff[g].max(MIN_DANCOFF),
                )
                .unwrap()
            })
            .sum::<f64>()
            / ufg.len() as f64;
        let cylinder = Cylinder::new(vec![0.2, 0.38, 0.48, radius]).unwrap();
        let rings: Vec<Option<&MaterialInput>> = [1, 1, 2, 3].map(|i| Some(&materials[i])).to_vec();
        let spectrum = crate::ufg::slowing_down(&library, &cylinder, &rings, [4.0, 1.0e4]).unwrap();
        assert_eq!(spectrum.groups, 1..3);
        let volumes = cylinder.volumes();
        for g in ufg {
            let weights = [0, 1].map(|j| spectrum.flux[j][g - 1] * volumes[j]);
            let expected = (0..2)
                .map(|j| weights[j] * spectrum.nuclides[j][0].total[g - 1])
                .sum::<f64>()
                / weights.iter().sum::<f64>();
            let total = lump.effective[0].total[g];
            assert!(
                (total - expected).abs() < 1e-9 * expected,
                "group {g}: {total} vs {expected}"
            );
            // 捕獲も同じ中性子束と体積で平均する
            let expected = (0..2)
                .map(|j| weights[j] * spectrum.nuclides[j][0].capture[g - 1])
                .sum::<f64>()
                / weights.iter().sum::<f64>();
            let capture = lump.effective[0].capture[g];
            assert!(
                (capture - expected).abs() <= 1e-9 * expected.abs(),
                "group {g}: capture {capture} vs {expected}"
            );
        }
    }

    #[test]
    fn absorber_rod_rings_take_their_own_cross_sections() {
        // 棒のリングは、等価なセルのリングの群定数に SPH 因子を掛けたもの（燃料のピンと同じ）。UFG の範囲の群では
        // 内側のリングほど共鳴が遮蔽される
        let (_, _, shielded) = rod_lattice();
        let total =
            |region: usize, g: usize| shielded.materials[shielded.region_material[region]].total[g];
        assert!(
            total(12, 1) < total(13, 1),
            "{} vs {}",
            total(12, 1),
            total(13, 1)
        );
    }

    #[test]
    fn frames_are_the_outer_ring_of_the_cell() {
        let library = resonant_library();
        let materials = vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![
                    ("U-238".into(), 0.022),
                    ("U-235".into(), 5.0e-4),
                    ("O-16".into(), 0.045),
                ],
            },
            MaterialInput {
                name: "clad".into(),
                temperature: T,
                composition: vec![("Zr-90".into(), 0.043)],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
            },
        ];
        let shield = |frame: bool| {
            let mut pin = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
            let mut regions = vec!["fuel", crate::input::VOID, "clad", "water"];
            if frame {
                pin = pin.with_frame(0.05).unwrap();
                regions.push("clad");
            }
            let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
            let regions: Vec<String> = regions.into_iter().map(String::from).collect();
            shield_fuel(
                &library,
                &tracking,
                &materials,
                &regions,
                &ShieldingOptions::default(),
            )
            .unwrap()
        };
        let (plain, framed) = (shield(false), shield(true));
        // 等価なセルのリング: 燃料、ギャップ、被覆管、減速材と、枠
        let sph = |s: &FuelShielding| s.pins[0].sph.as_ref().unwrap().factors.len();
        assert_eq!(sph(&plain), 4);
        assert_eq!(sph(&framed), 5);
        // 枠の影は枠のリングが作るので、減速材のリングは枠のないセルより小さい
        let radius = |s: &FuelShielding| s.pins[0].equivalent.unwrap();
        assert!(
            radius(&framed) < radius(&plain),
            "{} {}",
            radius(&framed),
            radius(&plain)
        );
        // 枠のリングの断面積は格子に使わないので、SPH 因子は 1 のまま
        let rings = &framed.pins[0].sph.as_ref().unwrap().factors;
        assert!(rings[4].iter().all(|&mu| mu == 1.0), "{rings:?}");
        // 格子の枠は等価なセルのリングの物質で上書きされない（この Zr-90 は分離共鳴を持たないので、
        // 被覆管と同じ物質の無限希釈のまま。被覆管のリングは UFG の物質）
        assert_eq!(framed.region_material[4], 1);
        assert_ne!(framed.region_material[4], framed.region_material[2]);
    }

    /// 燃料、ギャップ、被覆管、水のピンセル（`frame` なら外周に被覆管と同じ物質の枠）を自己遮蔽する（領域ごとの物質の
    /// 名前とともに）。領域は 0 燃料、1 ギャップ、2 被覆管、3 水、4 枠
    fn shielded_pin_cell(frame: bool) -> (Tracking, FuelShielding, Vec<String>) {
        let library = resonant_library();
        let materials = vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![
                    ("U-238".into(), 0.022),
                    ("U-235".into(), 5.0e-4),
                    ("O-16".into(), 0.045),
                ],
            },
            MaterialInput {
                name: "clad".into(),
                temperature: T,
                composition: vec![("Zr-90".into(), 0.043)],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
            },
        ];
        let mut pin = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
        let mut names = vec!["fuel", crate::input::VOID, "clad", "water"];
        if frame {
            pin = pin.with_frame(0.05).unwrap();
            names.push("clad");
        }
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.02).unwrap();
        let names: Vec<String> = names.into_iter().map(String::from).collect();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &names,
            &ShieldingOptions::default(),
        )
        .unwrap();
        (tracking, shielded, names)
    }

    #[test]
    fn sph_factors_name_the_regions_of_their_rings() {
        // 等価なセルのリング（燃料、ギャップ、被覆管、減速材）は、格子の領域 0〜3 の 1 つずつ。枠のリングの物質は格子に
        // 使わないので、枠の領域 4 はどのリングにも入らず、因子の行は枠の分だけ 1 つ多い
        for frame in [false, true] {
            let (tracking, shielded, names) = shielded_pin_cell(frame);
            let pin = &shielded.pins[0];
            let rings = pin.ring_regions(&tracking.lattice).unwrap();
            assert_eq!(
                rings,
                vec![vec![0], vec![1], vec![2], vec![3]],
                "frame {frame}"
            );
            let sph = pin.sph.as_ref().unwrap();
            // UFG の範囲はライブラリの 1、2 群
            assert_eq!(sph.groups, 1..3);
            assert_eq!(sph.factors.len(), 4 + usize::from(frame));
            assert!(sph.factors.iter().all(|row| row.len() == 2));
            // 減速材のリングの領域の物質は水を遮蔽したもので、燃料や被覆管と別
            assert_eq!(names[3], "water");
            let m = |r: usize| shielded.region_material[r];
            assert!(m(3) != m(0) && m(3) != m(2), "frame {frame}");
        }
    }

    #[test]
    fn physical_flux_restores_the_rod_sph_factors() {
        // 棒のリングにも SPH 因子を掛けるので、MOC の中性子束がどこも 1 なら、棒の領域の物理的な中性子束は UFG の
        // 範囲の群（1, 2）で 1 と違い、範囲の外の群で 1
        let (library, _, shielded) = rod_lattice();
        let lattice = &shielded.lattice;
        let flux = vec![vec![1.0; library.groups.len()]; lattice.regions()];
        let physical = shielded.physical_flux(lattice, &flux).unwrap();
        for region in [12, 13] {
            assert_eq!(physical[region][0], 1.0, "region {region}");
            assert!(
                (physical[region][1] - 1.0).abs() > 1e-6,
                "region {region}: {:?}",
                physical[region]
            );
        }
    }

    #[test]
    fn physical_flux_restores_the_sph_factors() {
        // MOC の中性子束がどこも 1 なら、物理的な中性子束は、リングの領域の UFG の範囲の群でそのリングの因子、
        // ほかは 1（枠の領域を含む）
        for frame in [false, true] {
            let (tracking, shielded, _) = shielded_pin_cell(frame);
            let lattice = &tracking.lattice;
            let ng = shielded.materials[0].groups();
            let flux = vec![vec![1.0; ng]; lattice.regions()];
            let physical = shielded.physical_flux(lattice, &flux).unwrap();
            let factors = &shielded.pins[0].sph.as_ref().unwrap().factors;
            for (region, row) in physical.iter().enumerate() {
                for (g, &value) in row.iter().enumerate() {
                    let expected = if region < 4 && (1..3).contains(&g) {
                        factors[region][g - 1]
                    } else {
                        1.0
                    };
                    assert_eq!(value, expected, "frame {frame}, region {region}, group {g}");
                }
            }
            // 燃料では因子が 1 と違う（因子が実際に効いている）
            assert!(
                physical[0].iter().any(|&f| (f - 1.0).abs() > 1e-6),
                "{:?}",
                physical[0]
            );
        }
    }

    #[test]
    fn physical_flux_refuses_mismatched_sph() {
        let refused = |shielded: &FuelShielding, tracking: &Tracking| {
            let ng = shielded.materials[0].groups();
            let flux = vec![vec![1.0; ng]; tracking.lattice.regions()];
            matches!(
                shielded.physical_flux(&tracking.lattice, &flux),
                Err(SelfShieldingError::SphLayout(0))
            )
        };
        let (tracking, shielded, _) = shielded_pin_cell(false);
        let lattice = &tracking.lattice;
        let ng = shielded.materials[0].groups();
        let flux = vec![vec![1.0; ng]; lattice.regions()];
        // 領域の数、群の数の違う中性子束
        assert!(matches!(
            shielded.physical_flux(lattice, &flux[1..]),
            Err(SelfShieldingError::Regions { .. })
        ));
        let short: Vec<Vec<f64>> = flux.iter().map(|f| f[1..].to_vec()).collect();
        assert!(matches!(
            shielded.physical_flux(lattice, &short),
            Err(SelfShieldingError::FluxGroups { .. })
        ));
        let edit = |change: &dyn Fn(&mut Sph)| {
            let mut edited = shielded.clone();
            change(edited.pins[0].sph.as_mut().unwrap());
            edited
        };
        // リングより少ない行、群の数の違う行、ライブラリの群の外の列、枠のないセルの余分な（1 だけの）行
        assert!(refused(
            &edit(&|s| s.factors.truncate(s.factors.len() - 1)),
            &tracking
        ));
        assert!(refused(&edit(&|s| s.factors[0].truncate(1)), &tracking));
        assert!(refused(&edit(&|s| s.groups = ng - 1..ng + 1), &tracking));
        assert!(refused(&edit(&|s| s.factors.push(vec![1.0; 2])), &tracking));
        // 枠のあるセルで、枠の行がないか、枠の行の因子が 1 でない
        let (tracking, shielded, _) = shielded_pin_cell(true);
        let edit = |change: &dyn Fn(&mut Sph)| {
            let mut edited = shielded.clone();
            change(edited.pins[0].sph.as_mut().unwrap());
            edited
        };
        assert!(refused(
            &edit(&|s| s.factors.truncate(s.factors.len() - 1)),
            &tracking
        ));
        assert!(refused(&edit(&|s| s.factors[4][0] = 1.1), &tracking));
    }

    #[test]
    fn physical_flux_refuses_another_lattice() {
        // 領域の数が同じでも、自己遮蔽したのと違う格子（燃料の半径が違う）
        let (_, shielded, _) = shielded_pin_cell(false);
        let other = Lattice::pin_cell(PinCell::new(1.26, vec![0.3, 0.418, 0.475]).unwrap());
        let ng = shielded.materials[0].groups();
        let flux = vec![vec![1.0; ng]; other.regions()];
        assert!(matches!(
            shielded.physical_flux(&other, &flux),
            Err(SelfShieldingError::ForeignLattice)
        ));
    }

    /// 自己遮蔽した格子で成り立つ対称性（物質の参照で確かめる）が `expected` で、その対称性で解いた MOC の解が全体を
    /// 解いた解と一致すること。自己遮蔽が対称な部分を同じ組として解いていなければ、物質の参照が移り合わない
    fn assert_symmetric_solves(tracking: &Tracking, shielded: &FuelShielding, expected: Symmetry) {
        use crate::moc::{Options, solve, solve_symmetric};
        let materials: Vec<&Material> = shielded
            .region_material
            .iter()
            .map(|&m| &shielded.materials[m])
            .collect();
        let identity: Vec<usize> = materials
            .iter()
            .map(|&m| std::ptr::from_ref(m) as usize)
            .collect();
        let detected = detect(tracking, &identity).unwrap();
        assert_eq!(
            detected.symmetry(),
            expected,
            "{}",
            detected.symmetry().name()
        );
        let options = Options::default();
        let full = solve(tracking, &materials, &options).unwrap();
        let reduced = solve_symmetric(tracking, &materials, &options, &detected).unwrap();
        assert!(
            (full.k - reduced.k).abs() < options.k_tolerance,
            "{} vs {}",
            full.k,
            reduced.k
        );
        for (a, b) in full
            .flux
            .iter()
            .flatten()
            .zip(reduced.flux.iter().flatten())
        {
            assert!(
                (a - b).abs() < options.flux_tolerance * a.abs(),
                "{a} vs {b}"
            );
        }
    }

    #[test]
    fn corner_lumps_keep_the_symmetry() {
        // 3 × 3 の、辺の中央の 4 つのセルが、四隅に分かれた構造材（外のリングの内径 0.65 > 半ピッチ 0.63）のセル。
        // 移り合う隅は同じ塊の物質を共有し、1/8 対称のまま解ける（#0054 の (a)）
        let library = library_with_zirconium();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let corner = PinCell::with_sectors(1.26, vec![0.65], 8).unwrap();
        let layout = vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]];
        let names = |t: usize| -> Vec<String> {
            let names: Vec<&str> = if t == 0 {
                vec!["fuel", "water"]
            } else {
                [["water"; 8], ["zirc"; 8]].concat()
            };
            names.into_iter().map(String::from).collect()
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, corner], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let mut materials = materials();
        materials.push(MaterialInput {
            name: "zirc".into(),
            temperature: T,
            composition: vec![("Zr-90".into(), 0.04)],
        });
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        assert_lump_members_agree(
            &tracking,
            &materials,
            &regions,
            &shielded,
            library.groups.len(),
        );
        assert_symmetric_solves(&tracking, &shielded, Symmetry::OCTANT);
    }

    #[test]
    fn sectored_cells_keep_the_symmetry() {
        // 3 × 3 の中央が、扇形 8 で物質の違う水のセル（扇形 0、3、4、7 が water、1、2、5、6 が water2）。左右と
        // 上下の鏡映では移り合うが対角線の鏡映では移り合わないので、ピンは鏡映で組になり、1/4 対称で解ける（#0054 の
        // (b)）
        let library = library();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let sectored = PinCell::with_sectors(1.26, vec![], 8).unwrap();
        let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
        let names = |t: usize| -> Vec<String> {
            let names: Vec<&str> = if t == 0 {
                vec!["fuel", "water"]
            } else {
                (0..8)
                    .map(|k| {
                        if [1, 2, 5, 6].contains(&k) {
                            "water2"
                        } else {
                            "water"
                        }
                    })
                    .collect()
            };
            names.into_iter().map(String::from).collect()
        };
        let regions: Vec<String> = layout.iter().flatten().flat_map(|&t| names(t)).collect();
        let lattice = Lattice::new(vec![pin, sectored], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
        let mut materials = materials();
        materials.push(MaterialInput {
            name: "water2".into(),
            temperature: T,
            composition: vec![("H-1".into(), 0.04)],
        });
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // 角のピンどうし、上下の辺のピンどうしは同じ物質を共有する
        let first = |cell: usize| shielded.region_material[lattice.cell_regions(cell).start];
        assert_eq!(first(0), first(8));
        assert_eq!(first(1), first(7));
        assert_symmetric_solves(&tracking, &shielded, Symmetry::QUARTER);
    }

    #[test]
    fn symmetry_of_another_tracking_is_an_error() {
        // 渡された対称操作は、同じ格子とトラッキングで求めたものでなければならない（トラックの数やセルの数が違えば、
        // 像の対応が成り立たない）。黙って使ったりパニックしたりせず、エラーにする
        let library = library();
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let lattice = |n: usize| Lattice::new(vec![pin.clone()], vec![vec![0; n]; n]).unwrap();
        let names = |n: usize| -> Vec<String> {
            (0..n * n)
                .flat_map(|_| ["fuel".to_string(), "water".to_string()])
                .collect()
        };
        let (small, large) = (lattice(2), lattice(3));
        let tracking = Tracking::new(&large, 16, 0.05).unwrap();
        let regions = names(3);
        // 方位角の数の違うトラッキング（領域とセルは同じ）と、セルの数の違う格子
        let other = Tracking::new(&large, 8, 0.05).unwrap();
        let smaller = Tracking::new(&small, 16, 0.05).unwrap();
        let foreign = [
            detect(&other, &regions).unwrap(),
            detect(&smaller, &names(2)).unwrap(),
        ];
        for symmetry in &foreign {
            let result = shield_fuel_with(
                &library,
                &tracking,
                &materials(),
                &regions,
                symmetry,
                &ShieldingOptions::default(),
            );
            assert!(
                matches!(result, Err(SelfShieldingError::ForeignSymmetry)),
                "{:?}",
                result.err()
            );
        }
    }

    #[test]
    fn symmetric_pins_are_solved_once() {
        // 3×3 の中央が水のセル: 角の 4 本と辺の 4 本がそれぞれ同じ組になり、等価なセルの減速計算は組ごとに1回
        let library = resonant_library();
        let materials = vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![
                    ("U-238".into(), 0.022),
                    ("U-235".into(), 5.0e-4),
                    ("O-16".into(), 0.045),
                ],
            },
            MaterialInput {
                name: "clad".into(),
                temperature: T,
                composition: vec![("Zr-90".into(), 0.043)],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
            },
        ];
        let pin = PinCell::new(1.26, vec![0.2, 0.35, 0.4096, 0.418, 0.475]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
        let lattice = Lattice::new(vec![pin, water], layout).unwrap();
        let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
        let pin_regions = ["fuel", "fuel", "fuel", crate::input::VOID, "clad", "water"];
        let regions: Vec<String> = (0..9)
            .flat_map(|cell| {
                if cell == 4 {
                    vec!["water"]
                } else {
                    pin_regions.to_vec()
                }
            })
            .map(String::from)
            .collect();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        // 希釈した物質（3 と真空）に、2 組 × 真空でない 5 リングの物質
        assert_eq!(shielded.materials.len(), 4 + 2 * 5);
        // 角どうし、辺どうしは同じ物質を共有する。燃料のリングの最初の領域で比べる
        let first = |cell: usize| shielded.region_material[lattice.cell_regions(cell).start];
        assert_eq!(first(0), first(8));
        assert_eq!(first(1), first(7));
        assert_ne!(first(0), first(1));
        // ピンごとの結果は 8 本とも残り、組の中では同じ
        assert_eq!(shielded.pins.len(), 8);
        let radius = |cell: usize| {
            shielded
                .pins
                .iter()
                .find(|p| p.cell == cell)
                .unwrap()
                .equivalent
                .unwrap()
        };
        assert_eq!(radius(0), radius(8));
        assert_eq!(radius(1), radius(3));
        // 報告する Dancoff 係数も、外径を決めた代表の値
        let dancoff = |cell: usize| {
            &shielded
                .pins
                .iter()
                .find(|p| p.cell == cell)
                .unwrap()
                .dancoff
        };
        assert_eq!(dancoff(0), dancoff(8));
    }

    #[test]
    fn isolated_pins_are_shielded() {
        // 水が厚く、Dancoff 係数がセルの数値誤差より小さいピンも、孤立したピンとして解ける
        let library = resonant_library();
        let materials = vec![
            MaterialInput {
                name: "fuel".into(),
                temperature: T,
                composition: vec![
                    ("U-238".into(), 0.022),
                    ("U-235".into(), 5.0e-4),
                    ("O-16".into(), 0.045),
                ],
            },
            MaterialInput {
                name: "water".into(),
                temperature: T,
                composition: vec![("H-1".into(), 0.05), ("O-16".into(), 0.025)],
            },
        ];
        let pin = PinCell::new(30.0, vec![0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(pin), 16, 0.05).unwrap();
        let regions: Vec<String> = ["fuel", "water"].map(String::from).to_vec();
        let shielded = shield_fuel(
            &library,
            &tracking,
            &materials,
            &regions,
            &ShieldingOptions::default(),
        )
        .unwrap();
        let radius = shielded.pins[0].equivalent.unwrap();
        assert!(radius > 3.0, "{radius}");
    }

    #[test]
    fn rejects_fuel_outside_the_inner_rings() {
        let library = library();
        let (tracking, mut regions) = geometry(vec![vec![0]]);
        regions.swap(0, 1);
        let result = shield_fuel(
            &library,
            &tracking,
            &materials(),
            &regions,
            &ShieldingOptions::default(),
        );
        assert!(
            matches!(result, Err(SelfShieldingError::FuelShape(0))),
            "{result:?}"
        );
    }
}
