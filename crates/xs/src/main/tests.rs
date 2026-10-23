use xlc_format::library::TemperatureData;

use super::*;

const INPUT: &str = r#"
# M1 のライブラリ
output = "m1.msgpack"
endf-directory = "endf"
thermal-directory = "tsl"
sigma0 = [1e10, 10, 1]

[nuclide.H-1]          # 水素
endf = "h.endf"
temperatures = [293.6, 600]
thermal = "h2o.endf"

[nuclide.U-238]
endf = "u.endf"
temperatures = [293.6, 900]
"#;

#[test]
fn no_resolved_resonances_means_no_pointwise_range() {
    // 分離共鳴を持つ核種がなければ、超微細群の計算に使うポイントワイズの断面積は要らない（エラーにしない）
    assert!(matches!(
        pointwise_range(&[], &[], &GroupStructure::wims69()),
        Ok(None)
    ));
}

#[test]
fn parses_input_relative_to_its_directory() {
    let input = parse_input(INPUT, Path::new("work")).unwrap();
    assert_eq!(input.output, PathBuf::from("work/m1.msgpack"));
    assert_eq!(input.sigma0, vec![1.0e10, 10.0, 1.0]);
    assert_eq!(input.scattering_order, 3);
    assert_eq!(
        input.nuclides,
        vec![
            NuclideInput {
                name: "H-1".into(),
                endf: PathBuf::from("work/endf/h.endf"),
                temperatures: vec![293.6, 600.0],
                thermal: Some(ThermalInput {
                    path: PathBuf::from("work/tsl/h2o.endf"),
                    upper: SabUpper::Evaluated,
                    outside: SabOutside::ShortCollisionTime,
                    free_gas_upper: DEFAULT_FREE_GAS_UPPER,
                }),
            },
            NuclideInput {
                name: "U-238".into(),
                endf: PathBuf::from("work/endf/u.endf"),
                temperatures: vec![293.6, 900.0],
                thermal: None,
            },
        ]
    );
}

#[test]
fn accepts_a_scattering_order() {
    let text = INPUT.replace("sigma0 =", "scattering-order = 1\nsigma0 =");
    let input = parse_input(&text, Path::new(".")).unwrap();
    assert_eq!(input.scattering_order, 1);
    let text = INPUT.replace("sigma0 =", "scattering-order = 0\nsigma0 =");
    assert_eq!(
        parse_input(&text, Path::new(".")).unwrap().scattering_order,
        0
    );
}

#[test]
fn accepts_the_thermal_scattering_limits() {
    let text = INPUT.replace(
        "thermal = \"h2o.endf\"",
        "thermal = \"h2o.endf\"\nsab-upper = 4\nsab-outside = \"zero\"\nfree-gas-upper = 48.052",
    );
    let input = parse_input(&text, Path::new(".")).unwrap();
    let thermal = input.nuclides[0].thermal.as_ref().unwrap();
    assert_eq!(thermal.upper, SabUpper::Energy(4.0));
    assert_eq!(thermal.outside, SabOutside::Zero);
    assert_eq!(thermal.free_gas_upper, 48.052);
    let text = INPUT.replace(
        "thermal = \"h2o.endf\"",
        "thermal = \"h2o.endf\"\nsab-upper = \"evaluated\"\nsab-outside = \"sct\"",
    );
    let input = parse_input(&text, Path::new(".")).unwrap();
    let thermal = input.nuclides[0].thermal.as_ref().unwrap();
    assert_eq!(thermal.upper, SabUpper::Evaluated);
    assert_eq!(thermal.outside, SabOutside::ShortCollisionTime);
}

#[test]
fn rejects_bad_thermal_scattering_limits() {
    let with = |lines: &str| {
        let text = INPUT.replace(
            "thermal = \"h2o.endf\"",
            &format!("thermal = \"h2o.endf\"\n{lines}"),
        );
        parse_input(&text, Path::new("."))
    };
    for lines in [
        "sab-upper = \"top\"",
        "sab-upper = 0",
        "sab-upper = -4",
        "sab-outside = \"free-gas\"",
        "free-gas-upper = 0",
        "free-gas-upper = \"all\"",
    ] {
        assert!(with(lines).is_err(), "{lines}");
    }
    // S(α, β) を使わない核種には書けない
    for key in [
        "sab-upper = 4",
        "sab-outside = \"sct\"",
        "free-gas-upper = 100",
    ] {
        let text = INPUT.replace(
            "temperatures = [293.6, 900]",
            &format!("temperatures = [293.6, 900]\n{key}"),
        );
        assert!(parse_input(&text, Path::new(".")).is_err(), "{key}");
    }
}

#[test]
fn rejects_bad_input() {
    let base = Path::new(".");
    let replaced = |from: &str, to: &str| {
        assert!(INPUT.contains(from), "{from}");
        parse_input(&INPUT.replace(from, to), base)
    };
    // 必須のキーがない、知らないキー
    for line in [
        "output = \"m1.msgpack\"\n",
        "endf-directory = \"endf\"\n",
        "sigma0 = [1e10, 10, 1]\n",
    ] {
        assert!(replaced(line, "").is_err(), "without {line}");
    }
    assert!(replaced("output", "out").is_err());
    assert!(replaced("endf = \"u.endf\"", "endf = \"u.endf\"\nlibrary = \"x\"").is_err());
    // 核種名が約束事の形でない
    assert!(replaced("[nuclide.U-238]", "[nuclide.U238]").is_err());
    assert!(replaced("[nuclide.U-238]", "[nuclide.\"\"]").is_err());
    // 核種がない、同じ名前の核種（TOML の重複キー）
    let header = INPUT.split("[nuclide.H-1]").next().unwrap();
    assert!(parse_input(header, base).is_err());
    assert!(replaced("[nuclide.U-238]", "[nuclide.H-1]").is_err());
    // thermal を使うのに thermal-directory がない
    assert!(replaced("thermal-directory = \"tsl\"\n", "").is_err());
    // σ0 が降順でない、正でない
    assert!(replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1, 10]").is_err());
    assert!(replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1e10, 0]").is_err());
    // 散乱行列の次数が負、整数でない、文字列
    for value in ["-1", "1.5", "\"3\""] {
        let text = format!("scattering-order = {value}\n{INPUT}");
        assert!(parse_input(&text, base).is_err(), "{value}");
    }
    // 温度がない、正でない、重複している
    assert!(replaced("temperatures = [293.6, 900]", "temperatures = []").is_err());
    assert!(replaced("temperatures = [293.6, 900]", "temperatures = [-1]").is_err());
    assert!(replaced("temperatures = [293.6, 900]", "temperatures = [900, 900]").is_err());
    // TOML の文法の誤りは位置つきで報告する
    let err = replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1e10, 10").unwrap_err();
    assert!(err.to_string().contains("line"), "{err}");
}

#[test]
fn the_pointwise_range_file_round_trips() {
    for range in [Some((4.0, 9118.0)), Some((0.1 + 0.2, 1.0e5)), None] {
        let text = range_to_toml(range);
        assert!(text.contains("pointwise-range"), "{text}");
        assert_eq!(range_from_toml(&text).unwrap(), range);
    }
}

#[test]
fn rejects_a_bad_pointwise_range_file() {
    for text in [
        "",
        "pointwise-range = [4.0]",
        "pointwise-range = [4.0, 5.0, 6.0]",
        "pointwise-range = [5.0, 4.0]",
        "pointwise-range = [4.0, 4.0]",
        "pointwise-range = [0.0, 4.0]",
        "pointwise-range = [4.0, 1e5]\nextra = 1",
        "pointwise-range = \"4 1e5\"",
    ] {
        assert!(range_from_toml(text).is_err(), "{text}");
    }
}

fn part_of(name: &str, temperatures: &[f64]) -> Library {
    Library {
        groups: GroupStructure::wims69(),
        sigma0: vec![1.0e10],
        pointwise_range: None,
        nuclides: vec![NuclideData {
            name: name.into(),
            za: 0.0,
            awr: 1.0,
            resolved: None,
            thermal_scattering: false,
            temperatures: temperatures
                .iter()
                .map(|&temperature| TemperatureData {
                    temperature,
                    chi: Vec::new(),
                    scatter: Vec::new(),
                    shielded: Vec::new(),
                    pointwise: Default::default(),
                    reactions: Vec::new(),
                })
                .collect(),
        }],
    }
}

#[test]
fn a_nuclide_result_must_match_its_input() {
    let input = parse_input(INPUT, Path::new(".")).unwrap();
    let uranium = &input.nuclides[1];
    assert!(check_part(uranium, &part_of("U-238", &[293.6, 900.0])).is_ok());
    // 名前、温度の値と順、核種の数が入力と違う
    for part in [
        part_of("U-235", &[293.6, 900.0]),
        part_of("U-238", &[293.6]),
        part_of("U-238", &[900.0, 293.6]),
        part_of("U-238", &[293.6, 600.0]),
    ] {
        assert!(check_part(uranium, &part).is_err());
    }
    let mut two = part_of("U-238", &[293.6, 900.0]);
    two.nuclides.extend(part_of("H-1", &[293.6]).nuclides);
    assert!(check_part(uranium, &two).is_err());
    let mut none = part_of("U-238", &[293.6, 900.0]);
    none.nuclides.clear();
    assert!(check_part(uranium, &none).is_err());
}

#[test]
fn a_nuclide_result_must_share_the_conditions_of_the_input() {
    let input = parse_input(INPUT, Path::new(".")).unwrap();
    let range = Some((4.0, 9118.0));
    let mut part = part_of("U-238", &[293.6, 900.0]);
    part.sigma0 = input.sigma0.clone();
    part.pointwise_range = range.map(|(low, high)| [low, high]);
    for t in &mut part.nuclides[0].temperatures {
        t.scatter = vec![Vec::new(); input.scattering_order + 1];
    }
    assert!(check_shared(&input, range, &part).is_ok());
    let mut sigma0 = part.clone();
    sigma0.sigma0 = vec![1.0e10, 5.0, 1.0];
    let mut groups = part.clone();
    groups.groups = GroupStructure::new(vec![1.0e7, 1.0, 1.0e-5]).unwrap();
    let mut order = part.clone();
    order.nuclides[0].temperatures[1].scatter.pop();
    let mut other_range = part.clone();
    other_range.pointwise_range = Some([4.0, 1.0e3]);
    let mut no_range = part.clone();
    no_range.pointwise_range = None;
    for stale in [sigma0, groups, order, other_range, no_range] {
        assert!(check_shared(&input, range, &stale).is_err());
    }
    assert!(check_shared(&input, None, &part).is_err());
}

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_string()).collect()
}

#[test]
fn reads_the_number_of_threads_before_the_command() {
    let args = strings(&["--threads", "4", "range", "xs.toml", "xs-range.toml"]);
    let (threads, rest) = thread_option(&args).unwrap();
    assert_eq!(threads, NonZeroUsize::new(4));
    assert_eq!(rest, &args[2..]);
}

#[test]
fn leaves_the_number_of_threads_to_the_default_without_the_option() {
    let args = strings(&["xs.toml"]);
    let (threads, rest) = thread_option(&args).unwrap();
    assert_eq!(threads, None);
    assert_eq!(rest, &args[..]);
}

#[test]
fn rejects_a_bad_number_of_threads() {
    for args in [
        &["--threads"][..],
        &["--threads", "0", "xs.toml"],
        &["--threads", "-1", "xs.toml"],
        &["--threads", "four", "xs.toml"],
    ] {
        assert!(thread_option(&strings(args)).is_err(), "{args:?}");
    }
}
