use super::*;

fn shielded(scale: f64, fissile: bool) -> Shielded {
    Shielded {
        // 1.7929089665572825e-9 は、10 進表記を経ると最下位のビットがずれやすい値
        total: vec![11.5 * scale, 1.7929089665572825e-9],
        elastic: vec![4.0 * scale, 15.0],
        fission: if fissile {
            vec![1.2 * scale, 580.5]
        } else {
            vec![0.0, 0.0]
        },
        nu_fission: if fissile {
            vec![3.0 * scale, 1414.1]
        } else {
            vec![0.0, 0.0]
        },
        capture: vec![0.7 * scale, 98.3],
    }
}

/// チェーンの反応（代表の MT の昇順）: (n,2n) は分けず、捕獲を準安定状態（MF9 の LFS 2）と残り（LFS nil）に分ける
fn reactions() -> Vec<GroupReaction> {
    vec![
        GroupReaction {
            mt: 16,
            channels: vec![GroupChannel {
                final_state: None,
                cross_section: vec![0.002, 0.0],
            }],
        },
        GroupReaction {
            mt: 102,
            channels: vec![
                GroupChannel {
                    final_state: Some(ChannelState {
                        mf: 9,
                        lfs: Some(2),
                    }),
                    cross_section: vec![0.1, 10.0],
                },
                GroupChannel {
                    final_state: Some(ChannelState { mf: 9, lfs: None }),
                    cross_section: vec![0.6, 88.3],
                },
            ],
        },
    ]
}

fn temperature(t: f64, fissile: bool) -> TemperatureData {
    TemperatureData {
        // sample() の範囲 [4, 1e4] eV を覆う
        pointwise: Pointwise {
            energy: vec![4.0, 1.0e4],
            total: vec![11.0, 10.0],
            elastic: vec![9.0, 9.0],
            fission: vec![0.0, 0.0],
            nu_fission: vec![0.0, 0.0],
            capture: vec![0.5, 0.4],
        },
        temperature: t,
        chi: if fissile {
            vec![1.0, 0.0]
        } else {
            vec![0.0, 0.0]
        },
        scatter: vec![vec![vec![4.0, 0.01], vec![1.0e-7, 15.0]]],
        shielded: vec![shielded(1.0, fissile), shielded(0.5, fissile)],
        reactions: reactions(),
    }
}

pub(crate) fn sample() -> Library {
    let nuclide = |name: &str, fissile: bool, temperatures: &[f64]| NuclideData {
        name: name.into(),
        za: 92235.0,
        awr: 233.0248,
        resolved: None,
        thermal_scattering: false,
        temperatures: temperatures
            .iter()
            .map(|&t| temperature(t, fissile))
            .collect(),
    };
    Library {
        pointwise_range: Some([4.0, 1.0e4]),
        groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        sigma0: vec![1.0e10, 10.0],
        nuclides: vec![
            nuclide("U-235", true, &[293.6, 900.0]),
            nuclide("H-1", false, &[293.6]),
        ],
    }
}

/// `sample()` の核種 `index` だけを持つライブラリ
fn part(index: usize) -> Library {
    let mut library = sample();
    library.nuclides = vec![library.nuclides.swap_remove(index)];
    library
}

#[test]
fn merge_joins_the_nuclides_in_order() {
    let merged = Library::merge(vec![part(0), part(1)]).unwrap();
    assert_eq!(merged, sample());
    let reversed = Library::merge(vec![part(1), part(0)]).unwrap();
    let names: Vec<_> = reversed.nuclides.iter().map(|n| n.name.as_str()).collect();
    assert_eq!(names, ["H-1", "U-235"]);
}

#[test]
fn merge_rejects_libraries_that_do_not_share_the_grids() {
    let mut sigma0 = part(1);
    sigma0.sigma0 = vec![1.0e10, 5.0];
    let mut groups = part(1);
    groups.groups = GroupStructure::new(vec![1.0e7, 1.0, 1.0e-5]).unwrap();
    let mut range = part(1);
    range.pointwise_range = Some([4.0, 1.0e3]);
    let mut no_range = part(1);
    no_range.pointwise_range = None;
    for other in [sigma0, groups, range, no_range] {
        assert!(Library::merge(vec![part(0), other]).is_err());
    }
}

#[test]
fn merge_rejects_a_repeated_nuclide_and_an_empty_list() {
    assert!(Library::merge(vec![part(0), part(0)]).is_err());
    assert!(Library::merge(Vec::new()).is_err());
}

#[test]
fn round_trips_bit_exactly() {
    let library = sample();
    assert_eq!(
        Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap(),
        library
    );
}

#[test]
fn resolved_resonances_need_a_pointwise_range() {
    // 分離共鳴を持つ核種があるのに範囲がないライブラリは、書くときにも読むときにも拒む
    let resolved = Some(Resolved {
        top: 1.0e4,
        potential: 10.0,
    });
    let mut library = sample();
    library.pointwise_range = None;
    library.nuclides[0].resolved = resolved.clone();
    for n in &mut library.nuclides {
        for t in &mut n.temperatures {
            t.pointwise = Pointwise::default();
        }
    }
    assert!(library.to_msgpack().is_err());
    let mut file = sample_file();
    file.pointwise_range = Vec::new();
    file.nuclides[0].resolved = resolved;
    for n in &mut file.nuclides {
        for t in &mut n.temperatures {
            t.pointwise = Pointwise::default();
        }
    }
    let bytes = rmp_serde::to_vec_named(&file).unwrap();
    assert!(Library::from_msgpack(&bytes).is_err());
}

#[test]
fn pointwise_range_must_be_empty_or_a_pair() {
    for range in [vec![4.0], vec![4.0, 1.0e3, 1.0e4]] {
        let mut file = sample_file();
        file.pointwise_range = range;
        let bytes = rmp_serde::to_vec_named(&file).unwrap();
        assert!(Library::from_msgpack(&bytes).is_err());
    }
}

/// 温度 550 K と 573.6 K の H-1 様の核種（573.6 K は値を 2 倍にしたもの）
fn two_temperatures(thermal: bool) -> Library {
    let mut library = sample();
    library.nuclides.retain(|n| n.name != "H-1");
    let mut n = library.nuclides[0].clone();
    n.name = "H-1".into();
    n.thermal_scattering = thermal;
    let mut low = n.temperatures[0].clone();
    low.temperature = 550.0;
    let mut high = low.clone();
    high.temperature = 573.6;
    high.scatter[0][0][0] *= 2.0;
    high.shielded[0].total[0] *= 2.0;
    n.temperatures = vec![low, high];
    library.nuclides.push(n);
    library
}

#[test]
fn thermal_scattering_is_interpolated_in_temperature() {
    let library = two_temperatures(true);
    let (_, low) = library.lookup("H-1", 550.0).unwrap();
    let (_, t) = library.lookup("H-1", 565.0).unwrap();
    let f = 1.0 + (565.0 - 550.0) / (573.6 - 550.0);
    assert_eq!(t.temperature, 565.0);
    assert!((t.scatter[0][0][0] - f * low.scatter[0][0][0]).abs() < 1e-12);
    assert!((t.shielded[0].total[0] - f * low.shielded[0].total[0]).abs() < 1e-12);
    // 表の温度ちょうどは表のまま、範囲の外はエラー
    assert!(matches!(
        library.lookup("H-1", 550.0).unwrap().1,
        Cow::Borrowed(_)
    ));
    assert!(library.lookup("H-1", 600.0).is_err());
}

#[test]
fn interpolation_refuses_what_it_cannot_do() {
    let lookup = |change: &dyn Fn(&mut NuclideData)| {
        let mut library = two_temperatures(true);
        let n = library
            .nuclides
            .iter_mut()
            .find(|n| n.name == "H-1")
            .unwrap();
        change(n);
        library.lookup("H-1", 565.0).map(|_| ())
    };
    assert!(lookup(&|_| {}).is_ok());
    // 分離共鳴の遮蔽とドップラー拡がりは温度で内挿できない
    assert!(
        lookup(&|n| {
            n.resolved = Some(Resolved {
                top: 1.0e3,
                potential: 10.0,
            })
        })
        .is_err()
    );
    // 散乱行列の次数が違う
    assert!(
        lookup(&|n| {
            let p0 = n.temperatures[1].scatter[0].clone();
            n.temperatures[1].scatter.push(p0);
        })
        .is_err()
    );
    // ポイントワイズの断面積の範囲が違う
    assert!(lookup(&|n| n.temperatures[1].pointwise = Pointwise::default()).is_err());
}

#[test]
fn other_nuclides_are_not_interpolated() {
    let library = two_temperatures(false);
    assert!(library.lookup("H-1", 565.0).is_err());
}

#[test]
fn library_without_pointwise_range_round_trips() {
    // 分離共鳴を持つ核種がないライブラリ: 範囲はなく、ポイントワイズの断面積は空
    let mut library = sample();
    library.pointwise_range = None;
    for n in &mut library.nuclides {
        n.resolved = None;
        for t in &mut n.temperatures {
            t.pointwise = Pointwise::default();
        }
    }
    let bytes = library.to_msgpack().unwrap();
    assert_eq!(Library::from_msgpack(&bytes).unwrap(), library);
    // 範囲がないのにポイントワイズの断面積がある
    let mut bad = library.clone();
    bad.nuclides[0].temperatures[0].pointwise =
        sample().nuclides[0].temperatures[0].pointwise.clone();
    assert!(bad.to_msgpack().is_err());
    // 読むときも拒む
    let mut file = sample_file();
    file.pointwise_range = Vec::new();
    for n in &mut file.nuclides {
        n.resolved = None;
    }
    let bytes = rmp_serde::to_vec_named(&file).unwrap();
    assert!(Library::from_msgpack(&bytes).is_err());
}

#[test]
fn refuses_to_write_a_bad_pointwise_range() {
    for range in [[0.0, 1.0e4], [1.0e4, 4.0], [4.0, f64::INFINITY]] {
        let mut library = sample();
        library.pointwise_range = Some(range);
        assert!(
            matches!(library.to_msgpack(), Err(LibraryError::Invalid(_))),
            "{range:?}"
        );
    }
}

#[test]
fn pointwise_data_must_cover_the_range() {
    let mut library = sample();
    let pointwise = &mut library.nuclides[0].temperatures[0].pointwise;
    // 上端 1e4 eV に届かない
    *pointwise = Pointwise {
        energy: vec![4.0, 100.0],
        total: vec![11.0, 12.0],
        elastic: vec![9.0, 10.0],
        fission: vec![0.5, 0.4],
        nu_fission: vec![1.2, 1.0],
        capture: vec![0.3, 0.2],
    };
    assert!(matches!(
        library.to_msgpack(),
        Err(LibraryError::Invalid(_))
    ));
    library.nuclides[0].temperatures[0].pointwise = Pointwise::default();
    assert!(matches!(
        library.to_msgpack(),
        Err(LibraryError::Invalid(_))
    ));
}

#[test]
fn keeps_the_pointwise_data() {
    let mut library = sample();
    library.pointwise_range = Some([4.0, 1.0e4]);
    library.nuclides[0].temperatures[0].pointwise = Pointwise {
        energy: vec![4.0, 1.0e4],
        total: vec![11.0, 12.0],
        elastic: vec![9.0, 10.0],
        fission: vec![0.5, 0.4],
        nu_fission: vec![1.2, 1.0],
        capture: vec![0.3, 0.2],
    };
    let read = Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap();
    assert_eq!(read, library);
}

#[test]
fn keeps_the_scattering_moments() {
    // 散乱行列はルジャンドルの次数ごと。P0 は無限希釈の核種の散乱行列になる
    let mut library = sample();
    let p0 = library.nuclides[0].temperatures[0].scatter[0].clone();
    library.nuclides[0].temperatures[0]
        .scatter
        .push(vec![vec![0.1, 0.02], vec![0.0, 0.3]]);
    let read = Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap();
    assert_eq!(read, library);
    let u235 = read.nuclide("U-235", 293.6).unwrap();
    assert_eq!(u235.scatter, p0);
    // P1 以上は次数の順に
    assert_eq!(
        u235.anisotropic,
        vec![vec![vec![0.1, 0.02], vec![0.0, 0.3]]]
    );
    assert!(read.nuclide("H-1", 293.6).unwrap().anisotropic.is_empty());
    // 次数がない、群の数が合わない
    let mut empty = sample();
    empty.nuclides[0].temperatures[0].scatter.clear();
    assert!(empty.to_msgpack().is_err());
    let mut short = sample();
    short.nuclides[0].temperatures[0]
        .scatter
        .push(vec![vec![0.1, 0.02]]);
    assert!(short.to_msgpack().is_err());
}

#[test]
fn keeps_the_resolved_range() {
    let mut library = sample();
    let resolved = Resolved {
        top: 2.0e4,
        potential: 11.29,
    };
    library.nuclides[0].resolved = Some(resolved.clone());
    let read = Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap();
    assert_eq!(read.nuclides[0].resolved, Some(resolved));
    assert_eq!(read.nuclides[1].resolved, None);
}

#[test]
fn refuses_a_bad_resolved_range() {
    for (top, potential) in [(0.0, 11.29), (f64::INFINITY, 11.29), (2.0e4, -1.0)] {
        let mut library = sample();
        library.nuclides[0].resolved = Some(Resolved { top, potential });
        assert!(library.to_msgpack().is_err(), "{top} {potential}");
    }
}

#[test]
fn gives_infinite_dilution_at_a_tabulated_temperature() {
    let library = sample();
    let u235 = library.nuclide("U-235", 900.0).unwrap();
    assert_eq!((u235.name.as_str(), u235.temperature), ("U-235", 900.0));
    assert_eq!(u235.total, vec![11.5, 1.7929089665572825e-9]);
    assert_eq!(u235.chi, vec![1.0, 0.0]);
    assert!(library.nuclide("U-235", 600.0).is_err());
    assert!(library.nuclide("O-16", 293.6).is_err());
}

#[test]
fn interpolates_shielded_cross_sections_in_ln_sigma0() {
    let library = sample();
    // 格子点では表の値
    let at_points = library.shielded("U-235", 900.0, &[1.0e10, 10.0]).unwrap();
    assert_eq!(at_points.total, vec![11.5, 1.7929089665572825e-9]);
    assert_eq!(at_points.fission, vec![1.2, 580.5]);
    // ln σ0 の中点では、両端の値の平均
    let middle = (1.0e10_f64 * 10.0).sqrt();
    let between = library.shielded("U-235", 900.0, &[middle, middle]).unwrap();
    assert!((between.total[0] - 0.75 * 11.5).abs() < 1e-12);
    assert!((between.elastic[0] - 0.75 * 4.0).abs() < 1e-12);
    assert!((between.nu_fission[0] - 0.75 * 3.0).abs() < 1e-12);
    assert!((between.capture[0] - 0.75 * 0.7).abs() < 1e-12);
}

#[test]
fn keeps_the_chain_reactions_and_the_capture() {
    let library = sample();
    let read = Library::from_msgpack(&library.to_msgpack().unwrap()).unwrap();
    let t = &read.nuclides[0].temperatures[0];
    assert_eq!(t.reactions, reactions());
    assert_eq!(t.shielded[1].capture, vec![0.35, 98.3]);
    assert_eq!(t.pointwise.capture, vec![0.5, 0.4]);
}

#[test]
fn reaction_channels_need_one_value_per_group() {
    let mut library = sample();
    library.nuclides[0].temperatures[0].reactions[0].channels[0]
        .cross_section
        .push(1.0);
    assert!(library.to_msgpack().is_err());
}

#[test]
fn reactions_must_not_repeat_and_need_channels() {
    let mut repeated = sample();
    let first = repeated.nuclides[0].temperatures[0].reactions[0].clone();
    repeated.nuclides[0].temperatures[0].reactions.push(first);
    assert!(repeated.to_msgpack().is_err());
    let mut empty = sample();
    empty.nuclides[0].temperatures[0].reactions[1]
        .channels
        .clear();
    assert!(empty.to_msgpack().is_err());
}

#[test]
fn interpolation_mixes_the_capture_and_the_reactions() {
    let mut library = two_temperatures(true);
    let h1 = library.nuclides.last_mut().unwrap();
    h1.temperatures[1].shielded[0].capture[0] *= 2.0;
    h1.temperatures[1].pointwise.capture = vec![1.0, 0.8];
    h1.temperatures[1].reactions[1].channels[0].cross_section[0] *= 2.0;
    let (_, low) = library.lookup("H-1", 550.0).unwrap();
    let low = low.into_owned();
    let (_, t) = library.lookup("H-1", 565.0).unwrap();
    let f = 1.0 + (565.0 - 550.0) / (573.6 - 550.0);
    assert!((t.shielded[0].capture[0] - f * low.shielded[0].capture[0]).abs() < 1e-12);
    assert!((t.pointwise.capture[0] - f * low.pointwise.capture[0]).abs() < 1e-12);
    let r = |d: &TemperatureData| d.reactions[1].channels[0].cross_section[0];
    assert!((r(&t) - f * r(&low)).abs() < 1e-15);
}

#[test]
fn rejects_sigma0_outside_the_table() {
    let library = sample();
    for sigma0 in [5.0, 1.0e11, f64::NAN] {
        let error = library.shielded("U-235", 900.0, &[sigma0, 100.0]);
        assert!(
            matches!(error, Err(LookupError::Sigma0 { group: 0, .. })),
            "{sigma0}: {error:?}"
        );
    }
    // σ0 は群の数だけ要る
    assert!(matches!(
        library.shielded("U-235", 900.0, &[100.0]),
        Err(LookupError::Sigma0Groups { .. })
    ));
}

fn sample_file() -> LibraryFile {
    let library = sample();
    LibraryFile {
        format: "xlc-multigroup".into(),
        group_bounds: library.groups.bounds().to_vec(),
        sigma0: library.sigma0,
        pointwise_range: library.pointwise_range.map(Vec::from).unwrap_or_default(),
        nuclides: library.nuclides,
    }
}

/// サンプルを書き換え、フィールド名付きの map 形式で書いて読む。
fn read_modified(change: impl FnOnce(&mut LibraryFile)) -> Result<Library, LibraryError> {
    let mut file = sample_file();
    change(&mut file);
    Library::from_msgpack(&rmp_serde::to_vec_named(&file).unwrap())
}

#[test]
fn refuses_to_write_non_finite_values() {
    let mut library = sample();
    library.nuclides[0].temperatures[0].shielded[1].total[0] = f64::NAN;
    assert!(matches!(
        library.to_msgpack(),
        Err(LibraryError::Invalid(_))
    ));
}

#[test]
fn written_as_a_map_with_the_format_name() {
    let bytes = sample().to_msgpack().unwrap();
    let value: std::collections::BTreeMap<String, serde::de::IgnoredAny> =
        rmp_serde::from_slice(&bytes).unwrap();
    let keys: Vec<&str> = value.keys().map(String::as_str).collect();
    assert_eq!(
        keys,
        vec![
            "format",
            "group_bounds",
            "nuclides",
            "pointwise_range",
            "sigma0"
        ]
    );
    #[derive(Deserialize)]
    struct Head {
        format: String,
    }
    let head: Head = rmp_serde::from_slice(&bytes).unwrap();
    assert_eq!(head.format, "xlc-multigroup");
}

#[test]
fn array_form_is_rejected() {
    // 既定の配列形式ではフィールド名がなく、順番の変わったファイルを黙って誤読しうる
    let bytes = rmp_serde::to_vec(&sample_file()).unwrap();
    assert!(matches!(
        Library::from_msgpack(&bytes),
        Err(LibraryError::NotMap)
    ));
}

#[test]
fn unknown_and_missing_fields_are_rejected() {
    #[derive(Serialize)]
    struct Extra {
        format: &'static str,
        group_bounds: Vec<f64>,
        sigma0: Vec<f64>,
        nuclides: Vec<NuclideData>,
        comment: &'static str,
    }
    #[derive(Serialize)]
    struct Missing {
        format: &'static str,
        group_bounds: Vec<f64>,
        nuclides: Vec<NuclideData>,
    }
    let file = sample_file();
    let extra = Extra {
        format: "xlc-multigroup",
        group_bounds: file.group_bounds.clone(),
        sigma0: file.sigma0.clone(),
        nuclides: file.nuclides.clone(),
        comment: "x",
    };
    let missing = Missing {
        format: "xlc-multigroup",
        group_bounds: file.group_bounds.clone(),
        nuclides: file.nuclides.clone(),
    };
    for bytes in [
        rmp_serde::to_vec_named(&extra).unwrap(),
        rmp_serde::to_vec_named(&missing).unwrap(),
    ] {
        assert!(matches!(
            Library::from_msgpack(&bytes),
            Err(LibraryError::Decode(_))
        ));
    }
}

#[test]
fn other_formats_are_rejected() {
    let err = read_modified(|f| f.format = "xlc-fewgroup".into()).unwrap_err();
    assert!(matches!(err, LibraryError::Format(_)), "{err}");
}

#[test]
fn wrong_number_of_values_is_error() {
    let err = read_modified(|f| {
        f.nuclides[0].temperatures[0].shielded[0].total.pop();
    })
    .unwrap_err();
    assert!(err.to_string().contains("U-235"), "{err}");
    assert!(
        read_modified(|f| f.nuclides[1].temperatures[0].scatter[0][1]
            .pop()
            .map(|_| ())
            .unwrap())
        .is_err()
    );
    assert!(read_modified(|f| f.nuclides[1].temperatures[0].chi.clear()).is_err());
}

#[test]
fn shielded_sets_must_match_sigma0() {
    assert!(
        read_modified(|f| {
            f.nuclides[0].temperatures[0].shielded.pop();
        })
        .is_err()
    );
}

#[test]
fn bounds_and_sigma0_must_be_descending() {
    assert!(read_modified(|f| f.sigma0 = vec![10.0, 1.0e10]).is_err());
    assert!(read_modified(|f| f.group_bounds = vec![1.0e-5, 0.625, 1.0e7]).is_err());
    assert!(read_modified(|f| f.group_bounds = vec![1.0e7]).is_err());
}

#[test]
fn nuclides_need_temperatures() {
    assert!(read_modified(|f| f.nuclides[1].temperatures.clear()).is_err());
}

#[test]
fn broken_bytes_are_rejected() {
    assert!(matches!(
        Library::from_msgpack(&[0x81, 0xc1]),
        Err(LibraryError::Decode(_))
    ));
    assert!(Library::from_msgpack(&[]).is_err());
}

#[test]
fn nested_structs_in_array_form_are_rejected() {
    // 最上位だけ map で、核種を配列形式で書いたファイル
    #[derive(Serialize)]
    struct Mixed {
        format: &'static str,
        group_bounds: Vec<f64>,
        sigma0: Vec<f64>,
        nuclides: embedded::Nested,
    }
    let file = sample_file();
    let nested = rmp_serde::to_vec(&file.nuclides).unwrap();
    let mixed = Mixed {
        format: "xlc-multigroup",
        group_bounds: file.group_bounds,
        sigma0: file.sigma0,
        nuclides: embedded::Nested(nested),
    };
    let bytes = rmp_serde::to_vec_named(&mixed).unwrap();
    assert!(matches!(
        Library::from_msgpack(&bytes),
        Err(LibraryError::NotMap)
    ));
}

/// 書き出し済みの MessagePack をそのまま埋め込むための型。
mod embedded {
    pub struct Nested(pub Vec<u8>);

    impl serde::Serialize for Nested {
        fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
            let value = rmpv::decode::read_value(&mut self.0.as_slice()).unwrap();
            value.serialize(s)
        }
    }
}

#[test]
fn trailing_bytes_are_rejected() {
    let mut bytes = sample().to_msgpack().unwrap();
    bytes.extend([0xc1, 0xff, 0x00]);
    assert!(Library::from_msgpack(&bytes).is_err());
}

#[test]
fn non_finite_values_are_rejected() {
    assert!(read_modified(|f| f.nuclides[0].temperatures[0].temperature = f64::NAN).is_err());
    assert!(
        read_modified(|f| f.nuclides[0].temperatures[1].shielded[1].total[0] = f64::NAN).is_err()
    );
    assert!(
        read_modified(|f| f.nuclides[1].temperatures[0].scatter[0][0][1] = f64::INFINITY).is_err()
    );
    assert!(read_modified(|f| f.nuclides[1].awr = f64::NAN).is_err());
    assert!(read_modified(|f| f.group_bounds[0] = f64::INFINITY).is_err());
    assert!(read_modified(|f| f.sigma0[0] = f64::INFINITY).is_err());
}

#[test]
fn nuclide_without_the_resolved_field_is_rejected() {
    // 項目を省いた核種を、「分離共鳴なし」として読まない（ADR-0013）。nil は書ける
    let bytes = rmp_serde::to_vec_named(&sample_file()).unwrap();
    let mut value = rmpv::decode::read_value(&mut bytes.as_slice()).unwrap();
    assert!(Library::from_msgpack(&bytes).is_ok());
    let rmpv::Value::Map(file) = &mut value else {
        panic!("a library is a map");
    };
    let (_, nuclides) = file
        .iter_mut()
        .find(|(k, _)| k.as_str() == Some("nuclides"))
        .unwrap();
    let rmpv::Value::Array(nuclides) = nuclides else {
        panic!("nuclides is an array");
    };
    let rmpv::Value::Map(nuclide) = &mut nuclides[0] else {
        panic!("a nuclide is a map");
    };
    nuclide.retain(|(k, _)| k.as_str() != Some("resolved"));
    let mut bytes = Vec::new();
    rmpv::encode::write_value(&mut bytes, &value).unwrap();
    assert!(matches!(
        Library::from_msgpack(&bytes),
        Err(LibraryError::Decode(_))
    ));
}

#[test]
fn duplicate_nuclides_and_temperatures_are_rejected() {
    assert!(read_modified(|f| f.nuclides[1].name = "U-235".into()).is_err());
    assert!(read_modified(|f| f.nuclides[0].temperatures[1].temperature = 293.6).is_err());
}

#[test]
fn map_keys_must_be_strings() {
    // serde はバイナリのキーもフィールド名として受け付けるので、構造の検査をすり抜けさせない
    let file = sample_file();
    let nuclides = rmpv::ext::to_value(&file.nuclides).unwrap();
    let value = rmpv::Value::Map(vec![
        ("format".into(), "xlc-multigroup".into()),
        (
            "group_bounds".into(),
            rmpv::ext::to_value(&file.group_bounds).unwrap(),
        ),
        ("sigma0".into(), rmpv::ext::to_value(&file.sigma0).unwrap()),
        (rmpv::Value::Binary(b"nuclides".to_vec()), nuclides),
    ]);
    let mut bytes = Vec::new();
    rmpv::encode::write_value(&mut bytes, &value).unwrap();
    assert!(matches!(
        Library::from_msgpack(&bytes),
        Err(LibraryError::NotMap)
    ));
}

#[test]
fn temperatures_are_looked_up_exactly() {
    // 重複の判定と同じく完全一致で引く。近い温度を黙って受け入れない
    let library = sample();
    assert!(library.nuclide("U-235", 900.0 * (1.0 + 1.0e-12)).is_err());
}

#[test]
fn nuclide_names_follow_the_convention() {
    for name in ["U-235", "H-1", "Am-242m", "Hf-178m2", "Zr-90", "B-10"] {
        assert!(is_nuclide_name(name), "{name}");
    }
    for name in [
        "", "U235", "u-235", "U-", "U-235M", "Am-242m0", "Am-242m1", "U-0235", "U 235", "Uuo-294",
        "H-1x", "-1",
    ] {
        assert!(!is_nuclide_name(name), "{name}");
    }
}

#[test]
fn library_nuclide_names_are_checked() {
    assert!(read_modified(|f| f.nuclides[1].name = "H1".into()).is_err());
}

#[test]
fn reaction_channels_follow_the_chain_layout() {
    let bad = |change: fn(&mut Vec<GroupReaction>)| {
        let mut library = sample();
        change(&mut library.nuclides[0].temperatures[0].reactions);
        library.to_msgpack().is_err()
    };
    // 状態に分けない経路はただ 1 つ
    assert!(bad(|r| r[1].channels[0].final_state = None));
    // 残り（LFS が nil）は最後
    assert!(bad(|r| r[1].channels.swap(0, 1)));
    // 割合は MF9 か MF10
    assert!(bad(|r| r[1].channels[0].final_state = Some(ChannelState {
        mf: 8,
        lfs: Some(2)
    })));
    // チェーンの反応の MT だけ
    assert!(bad(|r| r[0].mt = 18));
}
