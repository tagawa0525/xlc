use super::*;
use crate::input::MaterialInput;
use xlc_format::group::GroupStructure;
use xlc_format::library::{NuclideData, Shielded, TemperatureData};

fn library() -> Library {
    let nuclide = |name: &str, fissile: bool| {
        let fission = |values: Vec<f64>| if fissile { values } else { vec![0.0, 0.0] };
        NuclideData {
            name: name.into(),
            za: 0.0,
            awr: 1.0,
            resolved: None,
            thermal_scattering: false,
            // 600 K の値は 293.6 K の 1.5 倍
            temperatures: [(293.6, 1.0), (600.0, 1.5)]
                .map(|(temperature, c)| {
                    let v = |values: Vec<f64>| values.iter().map(|x| c * x).collect();
                    TemperatureData {
                        sab_transfer: None,
                        reactions: Vec::new(),
                        pointwise: Default::default(),
                        temperature,
                        chi: fission(vec![1.0, 0.0]),
                        scatter: vec![
                            vec![v(vec![6.0, 2.0]), v(vec![0.5, 15.0])],
                            vec![v(vec![3.0, 0.4]), v(vec![0.1, 1.0])],
                        ],
                        shielded: vec![Shielded {
                            total: v(vec![10.0, 20.0]),
                            elastic: v(vec![8.0, 15.5]),
                            fission: fission(v(vec![1.0, 2.0])),
                            nu_fission: fission(v(vec![2.5, 5.0])),
                            capture: v(vec![0.5, 1.0]),
                        }],
                    }
                })
                .to_vec(),
        }
    };
    Library {
        pointwise_range: Some([4.0, 1.0e4]),
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10],
        nuclides: vec![nuclide("F", true), nuclide("M", false)],
    }
}

fn mix(library: &Library, composition: &[(String, f64)]) -> Result<Material, MaterialError> {
    Material::mix(library, composition, 293.6)
}

#[test]
fn mixes_number_densities() {
    let m = mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
    assert_eq!(m.groups(), 2);
    assert!((m.total[1] - 0.8).abs() < 1e-15);
    assert!((m.scatter[0][1] - 0.08).abs() < 1e-15);
    // P1 も個数密度で足す
    assert_eq!(m.anisotropic.len(), 1);
    assert!((m.anisotropic[0][0][1] - 0.016).abs() < 1e-15);
    assert!((m.nu_fission[1] - 0.05).abs() < 1e-15);
    assert_eq!(m.fission_terms.len(), 1);
    assert!((m.fission_terms[0].nu_fission[0] - 0.025).abs() < 1e-15);
    let a = m.absorption();
    assert!((a[0] - (0.4 - 0.32)).abs() < 1e-15);
}

#[test]
fn scaled_multiplies_the_reactions_of_each_group() {
    let m = mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
    assert!(m.scaled(&[2.0]).is_err());
    let s = m.scaled(&[2.0, 0.5]).unwrap();
    assert_eq!(s.total, vec![m.total[0] * 2.0, m.total[1] * 0.5]);
    // 散乱は入射の群（行）の係数
    assert_eq!(
        s.scatter[0],
        m.scatter[0].iter().map(|v| v * 2.0).collect::<Vec<_>>()
    );
    assert_eq!(
        s.scatter[1],
        m.scatter[1].iter().map(|v| v * 0.5).collect::<Vec<_>>()
    );
    assert_eq!(
        s.anisotropic[0][1],
        m.anisotropic[0][1]
            .iter()
            .map(|v| v * 0.5)
            .collect::<Vec<_>>()
    );
    assert_eq!(s.nu_fission[1], m.nu_fission[1] * 0.5);
    assert_eq!(s.fission[0], m.fission[0] * 2.0);
    assert_eq!(
        s.fission_terms[0].nu_fission[1],
        m.fission_terms[0].nu_fission[1] * 0.5
    );
    assert_eq!(s.fission_terms[0].chi, m.fission_terms[0].chi);
    s.validate().unwrap();
}

#[test]
fn fission_matrix_sums_each_term() {
    // 核種ごとの νΣf(g) χ(g') の和 [入射の群][出射の群]
    let mut m = mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
    m.fission_terms.push(FissionTerm {
        chi: vec![0.25, 0.75],
        nu_fission: vec![0.1, 0.2],
    });
    let expected = [[0.025 + 0.025, 0.075], [0.05 + 0.05, 0.15]];
    let matrix = m.fission_matrix();
    for (row, e) in matrix.iter().zip(expected) {
        for (v, e) in row.iter().zip(e) {
            assert!((v - e).abs() < 1e-15, "{matrix:?}");
        }
    }
}

#[test]
fn fissile_nuclide_needs_normalized_chi() {
    let mut lib = library();
    lib.nuclides[0].temperatures[0].chi = vec![0.0, 0.0];
    assert!(mix(&lib, &[("F".into(), 0.01)]).is_err());
}

#[test]
fn validate_detects_inconsistent_data() {
    let good = mix(&library(), &[("F".into(), 0.01), ("M".into(), 0.03)]).unwrap();
    assert!(good.validate().is_ok());
    let mut short_row = good.clone();
    short_row.scatter[0].pop();
    assert!(short_row.validate().is_err());
    let mut short_moment = good.clone();
    short_moment.anisotropic[0][1].pop();
    assert!(short_moment.validate().is_err());
    let mut missing_term = good.clone();
    missing_term.fission_terms.clear();
    assert!(missing_term.validate().is_err());
    let mut bad_chi = good.clone();
    bad_chi.fission_terms[0].chi = vec![0.5, 0.0];
    assert!(bad_chi.validate().is_err());
}

#[test]
fn nuclides_must_have_the_same_scattering_order() {
    let mut lib = library();
    lib.nuclides[1].temperatures[0].scatter.pop();
    assert!(mix(&lib, &[("F".into(), 0.01), ("M".into(), 0.03)]).is_err());
    // 核種が1つなら、その次数のまま
    assert!(
        mix(&lib, &[("M".into(), 0.03)])
            .unwrap()
            .anisotropic
            .is_empty()
    );
}

#[test]
fn temperature_must_be_in_the_library() {
    assert!(Material::mix(&library(), &[("F".into(), 0.01)], 900.0).is_err());
}

fn component(name: &str, temperature: f64, composition: &[(&str, f64)]) -> MaterialInput {
    MaterialInput {
        name: name.into(),
        temperature,
        composition: composition
            .iter()
            .map(|&(n, d)| (n.to_string(), d))
            .collect(),
    }
}

#[test]
fn mixtures_add_number_densities_by_volume() {
    // 成分の（核種、温度）ごとの個数密度 × 体積の割合の和（ADR-0128）。同じ温度の同じ核種は個数密度を足す
    let a = component("a", 293.6, &[("F", 0.01), ("M", 0.03)]);
    let b = component("b", 293.6, &[("M", 0.02)]);
    assert_eq!(
        mixed_nuclides(&[(&a, 0.25), (&b, 0.75)]),
        vec![
            ("F".to_string(), 293.6, 0.01 * 0.25),
            ("M".to_string(), 293.6, 0.03 * 0.25 + 0.02 * 0.75)
        ]
    );
    let mixed = Material::mixture(&library(), &[(&a, 0.25), (&b, 0.75)]).unwrap();
    let direct = mix(
        &library(),
        &[
            ("F".into(), 0.01 * 0.25),
            ("M".into(), 0.03 * 0.25 + 0.02 * 0.75),
        ],
    );
    assert_eq!(mixed, direct.unwrap());
    // 違う温度の成分は、それぞれの温度のライブラリの値で別の項のまま足す
    let hot = component("hot", 600.0, &[("M", 0.02)]);
    assert_eq!(
        mixed_nuclides(&[(&a, 0.25), (&hot, 0.75)]),
        vec![
            ("F".to_string(), 293.6, 0.01 * 0.25),
            ("M".to_string(), 293.6, 0.03 * 0.25),
            ("M".to_string(), 600.0, 0.02 * 0.75)
        ]
    );
    let m = Material::mixture(&library(), &[(&a, 0.25), (&hot, 0.75)]).unwrap();
    let expected = 0.0025 * 20.0 + 0.0075 * 20.0 + 0.015 * 30.0;
    assert!((m.total[1] - expected).abs() < 1e-14, "{}", m.total[1]);
    // 真空の成分（呼ぶ側で除く）だけの混合物は核種がない
    assert!(Material::mixture(&library(), &[]).is_err());
}

#[test]
fn mixtures_reproduce_the_vera_nozzle_and_plate_densities() {
    // VERA の仕様書の表 P3-3 の SS304 とホウ酸水（3A の 1300 ppm と 3B の 1066 ppm）を、水の体積の割合（トップ
    // ノズル 0.80853、ボトムノズル 0.72077、炉心板 0.5。表の H-1 の比）で混ぜると、表のノズルと炉心板の H-1、O-16
    // と鋼の核種を再現する（ADR-0128）。鋼の核種の差は、表の鋼と水の割合（Fe-56 と H-1 の比）の和が 1 からずれる
    // （トップノズル 1 + 1.0 × 10⁻⁵、ボトムノズル 1 − 1.1 × 10⁻⁵）ためで、最大 7.3 × 10⁻⁵。表のノズルと炉心板の
    // B-10 と B-11 は、減速材の値 × 水の割合と 0.05〜0.45% 違い（3A の炉心板で −0.43%、3B の炉心板で +0.45%）、
    // 割合の取り方で再現できないので比べない
    let ss304 = component(
        "ss304",
        600.0,
        &[
            ("C", 3.20895e-4),
            ("Si-28", 1.58197e-3),
            ("Si-29", 8.03653e-5),
            ("Si-30", 5.30394e-5),
            ("P-31", 6.99938e-5),
            ("Cr-50", 7.64915e-4),
            ("Cr-52", 1.47506e-2),
            ("Cr-53", 1.67260e-3),
            ("Cr-54", 4.16346e-4),
            ("Mn-55", 1.75387e-3),
            ("Fe-54", 3.44776e-3),
            ("Fe-56", 5.41225e-2),
            ("Fe-57", 1.24992e-3),
            ("Fe-58", 1.66342e-4),
            ("Ni-58", 5.30854e-3),
            ("Ni-60", 2.04484e-3),
            ("Ni-61", 8.88879e-5),
            ("Ni-62", 2.83413e-4),
            ("Ni-64", 7.21770e-5),
        ],
    );
    let water = |h: f64, o: f64| component("water", 600.0, &[("H-1", h), ("O-16", o)]);
    let moderator = [water(4.96224e-2, 2.48112e-2), water(4.96340e-2, 2.48170e-2)];
    // （水の割合、[3A と 3B の (H-1, O-16)]、鋼の核種（ss304 と同じ順））
    let table: [(f64, [[f64; 2]; 2], [f64; 19]); 3] = [
        (
            0.80853,
            [[4.01211e-2, 2.00606e-2], [4.01305e-2, 2.00653e-2]],
            [
                6.14459e-5, 3.02920e-4, 1.53886e-5, 1.01561e-5, 1.34026e-5, 1.46468e-4, 2.82449e-3,
                3.20275e-4, 7.97232e-5, 3.35836e-4, 6.60188e-4, 1.03635e-2, 2.39339e-4, 3.18517e-5,
                1.01650e-3, 3.91552e-4, 1.70205e-5, 5.42688e-5, 1.38207e-5,
            ],
        ),
        (
            0.72077,
            [[3.57661e-2, 1.78830e-2], [3.57744e-2, 1.78872e-2]],
            [
                8.96008e-5, 4.41720e-4, 2.24397e-5, 1.48097e-5, 1.95438e-5, 2.13581e-4, 4.11869e-3,
                4.67027e-4, 1.16253e-4, 4.89719e-4, 9.62690e-4, 1.51122e-2, 3.49006e-4, 4.64463e-5,
                1.48226e-3, 5.70964e-4, 2.48194e-5, 7.91351e-5, 2.01534e-5,
            ],
        ),
        (
            0.5,
            [[2.48112e-2, 1.24056e-2], [2.48171e-2, 1.24085e-2]],
            [
                1.60447e-4, 7.90985e-4, 4.01826e-5, 2.65197e-5, 3.49969e-5, 3.82458e-4, 7.37532e-3,
                8.36302e-4, 2.08173e-4, 8.76936e-4, 1.72388e-3, 2.70613e-2, 6.24963e-4, 8.31710e-5,
                2.65427e-3, 1.02242e-3, 4.44439e-5, 1.41707e-4, 3.60885e-5,
            ],
        ),
    ];
    for (fraction, waters, steel) in table {
        for (water, [h, o]) in moderator.iter().zip(waters) {
            let mixed = mixed_nuclides(&[(&ss304, 1.0 - fraction), (water, fraction)]);
            let density = |name: &str| {
                mixed
                    .iter()
                    .find(|(n, _, _)| n == name)
                    .map(|&(_, _, d)| d)
                    .unwrap()
            };
            // 表は有効数字 6 桁
            for (name, want) in [("H-1", h), ("O-16", o)] {
                let got = density(name);
                assert!((got / want - 1.0).abs() < 1e-5, "{name}: {got} vs {want}");
            }
            for ((name, _), want) in ss304.composition.iter().zip(steel) {
                let got = density(name);
                assert!((got / want - 1.0).abs() < 1e-4, "{name}: {got} vs {want}");
            }
        }
    }
}

#[test]
fn unknown_nuclide_is_error() {
    assert!(mix(&library(), &[("X".into(), 0.01)]).is_err());
    assert!(mix(&library(), &[("F".into(), -0.01)]).is_err());
    assert!(mix(&library(), &[]).is_err());
}
