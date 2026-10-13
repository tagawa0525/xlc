use super::*;
use crate::endf::Tape;
use crate::endf::tape::tests::rec;

fn r(f: [&str; 6]) -> String {
    rec(f, 825, 4, 2)
}

fn parse(text: &str) -> Result<Mf4, Error> {
    let tape = Tape::parse(text).unwrap();
    Mf4::parse(&tape.section(825, 4, 2).unwrap())
}

fn close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-12, "{actual:?} vs {expected:?}");
    }
}

/// Legendre 係数（LTT=1）の部分: 入射エネルギー 1e5 と 2e5 eV で a_1 = 0.1 と 0.3、a_2 = 0.02 と 0.04
fn legendre_part() -> String {
    [
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["0.0", "1.0e5", "0", "0", "2", "0"]),
        r(["0.1", "0.02", "", "", "", ""]),
        r(["0.0", "2.0e5", "0", "0", "2", "0"]),
        r(["0.3", "0.04", "", "", "", ""]),
    ]
    .concat()
}

/// 表（LTT=2）の部分: 2e5 eV で f(μ) = 0.5 + 0.3 μ、4e5 eV で等方
fn tabulated_part() -> String {
    [
        r(["0.0", "0.0", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["0.0", "2.0e5", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["-1.0", "0.2", "1.0", "0.8", "", ""]),
        r(["0.0", "4.0e5", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
    ]
    .concat()
}

fn head(ltt: &str, li: &str, lct: &str) -> String {
    [
        r(["8016.0", "15.85751", "0", ltt, "0", "0"]),
        r(["0.0", "15.85751", li, lct, "0", "0"]),
    ]
    .concat()
}

#[test]
fn isotropic_distribution() {
    let mf4 = parse(&head("0", "1", "1")).unwrap();
    assert_eq!(mf4.frame, Frame::Laboratory);
    assert_eq!(mf4.angular, Angular::Isotropic);
    close(&mf4.moments(1.0e6, 2), &[1.0, 0.0, 0.0]);
}

#[test]
fn legendre_coefficients_are_interpolated_between_energies() {
    let mf4 = parse(&(head("1", "0", "2") + &legendre_part())).unwrap();
    assert_eq!(mf4.frame, Frame::CenterOfMass);
    close(&mf4.moments(1.5e5, 3), &[1.0, 0.2, 0.03, 0.0]);
    // 次数を打ち切る
    close(&mf4.moments(1.0e5, 1), &[1.0, 0.1]);
}

#[test]
fn tabulated_distribution_gives_its_moments() {
    // f(μ) = 0.5 + 0.3 μ の a_1 = ∫ f μ dμ = 0.2、a_2 = 0
    let mf4 = parse(&(head("2", "0", "2") + &tabulated_part())).unwrap();
    close(&mf4.moments(2.0e5, 2), &[1.0, 0.2, 0.0]);
    close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
}

#[test]
fn mixed_distribution_switches_at_the_last_legendre_energy() {
    let text = head("3", "0", "2") + &legendre_part() + &tabulated_part();
    let mf4 = parse(&text).unwrap();
    close(&mf4.moments(1.5e5, 2), &[1.0, 0.2, 0.03]);
    close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
}

#[test]
fn log_lin_tables_are_integrated_exactly() {
    // ln f が μ に線形（INT=4）: f(μ) = 0.2 e^{b(μ+1)}、b = ln 3 / 2（f(1) = 0.6）の
    // a_1 = ∫ μ f dμ = 0.2 e^b [e^{bμ}(μ/b − 1/b²)]_{−1}^{1}
    let text = head("2", "0", "2")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "1"]),
            r(["1", "2", "", "", "", ""]),
            r(["0.0", "1.0e5", "0", "0", "1", "2"]),
            r(["2", "4", "", "", "", ""]),
            r(["-1.0", "0.2", "1.0", "0.6", "", ""]),
        ]
        .concat();
    let mf4 = parse(&text).unwrap();
    let b = 3.0_f64.ln() / 2.0;
    let primitive = |mu: f64| (b * mu).exp() * (mu / b - 1.0 / (b * b));
    // ∫ f dμ = 0.2 (e^{2b} − 1)/b で割って規格化する
    let norm = 0.2 * ((2.0 * b).exp() - 1.0) / b;
    let expected = 0.2 * b.exp() * (primitive(1.0) - primitive(-1.0)) / norm;
    let a = mf4.moments(1.0e5, 1);
    assert!((a[1] - expected).abs() < 1e-12, "{a:?} vs {expected}");
}

#[test]
fn repeated_incident_energy_is_a_discontinuity() {
    // 先頭の入射エネルギーが2回ある（Pb-204 など）: 2つ目の分布から使う
    let text = head("1", "0", "2")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "3"]),
            r(["3", "2", "", "", "", ""]),
            r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
            r(["0.5", "", "", "", "", ""]),
            r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
            r(["0.0", "", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "0"]),
            r(["0.2", "", "", "", "", ""]),
        ]
        .concat();
    let mf4 = parse(&text).unwrap();
    // 1e-5 eV の a_1 = 0 と 2e5 eV の 0.2 の間の線形補間
    let expected = 0.2 * (1.0e5 - 1.0e-5) / (2.0e5 - 1.0e-5);
    close(&mf4.moments(1.0e5, 1), &[1.0, expected]);
}

#[test]
fn density_of_each_form() {
    let iso = parse(&head("0", "1", "1")).unwrap();
    assert_eq!(iso.density(1.0e6, 0.3), 0.5);
    // Legendre: f(μ) = 1/2 + 3/2 a_1 μ + 5/2 a_2 P_2(μ)、1.5e5 eV で a_1 = 0.2、a_2 = 0.03
    let legendre = parse(&(head("1", "0", "2") + &legendre_part())).unwrap();
    let mu = 0.4_f64;
    let expected = 0.5 + 1.5 * 0.2 * mu + 2.5 * 0.03 * 0.5 * (3.0 * mu * mu - 1.0);
    assert!((legendre.density(1.5e5, mu) - expected).abs() < 1e-12);
    // 表: 2e5 eV の 0.5 + 0.3 μ と 4e5 eV の 0.5 を 3e5 eV で線形に補間
    let tabulated = parse(&(head("2", "0", "2") + &tabulated_part())).unwrap();
    assert!((tabulated.density(3.0e5, mu) - (0.5 + 0.15 * mu)).abs() < 1e-12);
}

#[test]
fn repeated_first_energy_uses_the_later_distribution_there() {
    let text = head("1", "0", "2")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "3"]),
            r(["3", "2", "", "", "", ""]),
            r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
            r(["0.5", "", "", "", "", ""]),
            r(["0.0", "1.0e-5", "0", "0", "1", "0"]),
            r(["0.0", "", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "0"]),
            r(["0.2", "", "", "", "", ""]),
        ]
        .concat();
    let mf4 = parse(&text).unwrap();
    close(&mf4.moments(1.0e-5, 1), &[1.0, 0.0]);
}

#[test]
fn tabulated_moments_follow_the_interpolated_density() {
    // 入射エネルギーの間を INT=4（ln f が E に線形）で補間する: モーメントは補間した密度のモーメント
    let text = head("2", "0", "2")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "4", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.2", "1.0", "0.8", "", ""]),
            r(["0.0", "4.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
        ]
        .concat();
    let mf4 = parse(&text).unwrap();
    let e = 3.0e5;
    let nodes = xlc_numeric::quadrature::gauss_legendre(64).unwrap();
    let integral = |g: &dyn Fn(f64) -> f64| -> f64 {
        nodes
            .iter()
            .map(|&(x, w)| w * g(x) * mf4.density(e, x))
            .sum()
    };
    // 補間した密度は規格化されていないので、∫ f dμ で割る
    let expected = integral(&|x| x) / integral(&|_| 1.0);
    let a = mf4.moments(e, 1);
    assert!((a[1] - expected).abs() < 1e-12, "{a:?} vs {expected}");
}

#[test]
fn tables_must_cover_all_angles() {
    let text = head("2", "0", "2")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "1"]),
            r(["1", "2", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-0.5", "0.5", "1.0", "0.5", "", ""]),
        ]
        .concat();
    assert!(parse(&text).is_err());
}

#[test]
fn transformation_matrix_is_not_supported() {
    // LVT = 1（HEAD の L1。Legendre 係数の変換行列が続く）は使われていないので扱わない
    let text = [
        r(["8016.0", "15.85751", "1", "0", "0", "0"]),
        r(["0.0", "15.85751", "1", "1", "0", "0"]),
    ]
    .concat();
    assert!(parse(&text).is_err());
}

#[test]
fn tabulated_moments_are_normalized() {
    // f(μ) = 0.98 (0.5 + 0.3 μ)（∫ f = 0.98）でも a_1 = 0.2
    let text = head("2", "0", "1")
        + &[
            r(["0.0", "0.0", "0", "0", "1", "1"]),
            r(["1", "2", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.196", "1.0", "0.784", "", ""]),
        ]
        .concat();
    let mf4 = parse(&text).unwrap();
    close(&mf4.moments(2.0e5, 1), &[1.0, 0.2]);
}

#[test]
fn mixed_parts_must_meet() {
    // Legendre の部分が 2e5 eV で終わり、表の部分が 4e5 eV から始まる
    let tabulated = [
        r(["0.0", "0.0", "0", "0", "1", "1"]),
        r(["1", "2", "", "", "", ""]),
        r(["0.0", "4.0e5", "0", "0", "1", "2"]),
        r(["2", "2", "", "", "", ""]),
        r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
    ]
    .concat();
    let text = head("3", "0", "2") + &legendre_part() + &tabulated;
    assert!(parse(&text).is_err());
}
