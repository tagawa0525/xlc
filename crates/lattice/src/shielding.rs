//! 等価原理による燃料の自己遮蔽（Dancoff 係数を入れた Carlvik の2項有理近似）。
//!
//! 燃料の脱出確率を P(x) = Σ_i β_i x / (x + α_i)（x = Σt l̄）で近似し、核種 k の背景断面積を
//! σ0,i = (Σ_{j≠k} N_j σt,j + α_i / l̄) / N_k とする。実効断面積は、2つの σ0 でのライブラリの値を、
//! NR 近似の群平均の中性子束 φ_i = σ0,i / (σt(σ0,i) + σ0,i) と β_i で重みづけして合わせる
//! （ライブラリの重み w/(σt + σ0) での群平均では、この式が厳密に成り立つ）。

use xlc_xs::library::{Library, LookupError, Shielded};

/// 自己遮蔽の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum ShieldingError {
    #[error("Dancoff factor {0} must be from 0 to 1")]
    Dancoff(f64),
    #[error("{given} rational terms for {groups} groups")]
    Groups { given: usize, groups: usize },
    #[error(transparent)]
    Lookup(#[from] LookupError),
}

/// 有理近似の1つの項。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RationalTerm {
    pub alpha: f64,
    pub beta: f64,
}

/// 1本のピンの燃料の等価原理のパラメータ: ペレットの平均弦長 l̄ (cm) と、群ごとの有理近似の項。
#[derive(Debug, Clone, PartialEq)]
pub struct Equivalence {
    pub chord: f64,
    pub terms: Vec<[RationalTerm; 2]>,
}

/// Dancoff 係数 `dancoff` での Carlvik の2項有理近似の係数。
///
/// A = (1 − C)/C として α₁,₂ = [(5A + 6) ∓ √(A² + 36A + 36)] / (2(A + 1))、
/// β₁ = ((4A + 6)/(A + 1) − α₁)/(α₂ − α₁)、β₂ = 1 − β₁。分子と分母に C を掛けた形
/// α₁,₂ = [(5 + C) ∓ √(1 + 34C + C²)] / 2、(4A + 6)/(A + 1) = 4 + 2C で計算し、C = 0（孤立した棒）でも使えるようにする。
pub fn carlvik(dancoff: f64) -> Result<[RationalTerm; 2], ShieldingError> {
    if !(0.0..=1.0).contains(&dancoff) {
        return Err(ShieldingError::Dancoff(dancoff));
    }
    let c = dancoff;
    let root = (1.0 + 34.0 * c + c * c).sqrt();
    let (alpha1, alpha2) = (0.5 * (5.0 + c - root), 0.5 * (5.0 + c + root));
    let beta1 = (4.0 + 2.0 * c - alpha1) / (alpha2 - alpha1);
    Ok([
        RationalTerm {
            alpha: alpha1,
            beta: beta1,
        },
        RationalTerm {
            alpha: alpha2,
            beta: 1.0 - beta1,
        },
    ])
}

/// 燃料の核種 `name`（個数密度 `density`、1/(barn cm)）の実効断面積。`background` は群ごとの、
/// 同じ燃料の他の核種の全断面積の和 (1/cm)。
pub fn effective(
    library: &Library,
    name: &str,
    temperature: f64,
    density: f64,
    background: &[f64],
    equivalence: &Equivalence,
) -> Result<Shielded, ShieldingError> {
    if equivalence.terms.len() != background.len() {
        return Err(ShieldingError::Groups {
            given: equivalence.terms.len(),
            groups: background.len(),
        });
    }
    let sigma0 = |i: usize| -> Vec<f64> {
        background
            .iter()
            .zip(&equivalence.terms)
            .map(|(b, terms)| (b + terms[i].alpha / equivalence.chord) / density)
            .collect()
    };
    let (sigma0_1, sigma0_2) = (sigma0(0), sigma0(1));
    let tables = [
        library.shielded(name, temperature, &sigma0_1)?,
        library.shielded(name, temperature, &sigma0_2)?,
    ];
    let sigma0 = [sigma0_1, sigma0_2];
    let groups = background.len();
    // 群ごとの重み β_i φ_i
    let weights: Vec<[f64; 2]> = (0..groups)
        .map(|g| {
            std::array::from_fn(|i| {
                let s = sigma0[i][g];
                equivalence.terms[g][i].beta * s / (tables[i].total[g] + s)
            })
        })
        .collect();
    let combine = |of: fn(&Shielded) -> &Vec<f64>| -> Vec<f64> {
        (0..groups)
            .map(|g| {
                let [w1, w2] = weights[g];
                (w1 * of(&tables[0])[g] + w2 * of(&tables[1])[g]) / (w1 + w2)
            })
            .collect()
    };
    Ok(Shielded {
        total: combine(|x| &x.total),
        elastic: combine(|x| &x.elastic),
        fission: combine(|x| &x.fission),
        nu_fission: combine(|x| &x.nu_fission),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use xlc_xs::library::{Library, NuclideData, Shielded, TemperatureData};
    use xlc_xs::multigroup::GroupStructure;

    const T: f64 = 900.0;

    /// 群 0 だけが σ0 に依存する核種
    fn library() -> Library {
        let shielded = |total: f64, elastic: f64| Shielded {
            total: vec![total, 20.0],
            elastic: vec![elastic, 15.0],
            fission: vec![0.0, 0.0],
            nu_fission: vec![0.0, 0.0],
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
            sigma0: vec![1.0e10, 100.0, 1.0],
            nuclides: vec![NuclideData {
                name: "U-238".into(),
                za: 92238.0,
                awr: 236.0,
                temperatures: vec![TemperatureData {
                    temperature: T,
                    chi: vec![0.0, 0.0],
                    scatter: vec![vec![6.0, 2.0], vec![0.5, 15.0]],
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
        let library = library();
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
            let flux = sigma0 / (table.total[0] + sigma0);
            numerator += term.beta * flux * table.total[0];
            denominator += term.beta * flux;
        }
        assert!(close(result.total[0], numerator / denominator));
        // σ0 に依存しない群はそのまま
        assert!(close(result.total[1], 20.0));
    }
}
