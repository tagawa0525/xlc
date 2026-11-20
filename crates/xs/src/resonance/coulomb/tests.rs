use super::*;
use crate::endf::OrbitalL;

mod cases;
use cases::*;

/// P の相対誤差の許容。ln P を積分するので、P の相対誤差は ln P の丸めの誤差（|ln P| に比例）の大きさになる
fn tolerance(p: f64) -> f64 {
    2e-15 * p.ln().abs().max(10.0)
}

#[test]
fn penetrability_matches_multiprecision_values() {
    for &(l, eta, rho, expected) in COULOMB_REFERENCE {
        let actual = penetrability(l, eta, rho).unwrap();
        let error = (actual / expected - 1.0).abs();
        let allowed = tolerance(expected);
        assert!(
            error <= allowed,
            "l={l} eta={eta} rho={rho}: {actual:e} against {expected:e} ({error:.2e} > {allowed:.2e})"
        );
    }
}

#[test]
fn penetrability_underflows_to_zero_far_below_the_barrier() {
    for &(l, eta, rho) in COULOMB_UNDERFLOW {
        assert_eq!(
            penetrability(l, eta, rho).unwrap(),
            0.0,
            "l={l} eta={eta} rho={rho}"
        );
    }
}

#[test]
fn zero_eta_is_the_hard_sphere_penetrability() {
    let ls = [
        OrbitalL::S,
        OrbitalL::P,
        OrbitalL::D,
        OrbitalL::F,
        OrbitalL::G,
    ];
    for l in ls {
        // RHO_MIN より下（積分する）、RHO_MIN、その上（CF2 だけ）
        for rho in [0.01, 0.5, 1.9, 2.0, 3.0, 30.0] {
            let expected = crate::resonance::factors::penetrability(l, rho);
            let actual = penetrability(l.number(), 0.0, rho).unwrap();
            let error = (actual / expected - 1.0).abs();
            assert!(
                error <= tolerance(expected),
                "l={} rho={rho}: {actual:e} against {expected:e}",
                l.number()
            );
        }
    }
}
