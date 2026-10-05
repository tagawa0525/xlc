//! 群平均の重み関数（中性子スペクトルの近似）。

/// 熱 Maxwell 分布 + 1/E + 核分裂スペクトルをつないだ重み関数（NJOY GROUPR の iwt=4）。
///
/// - E ≤ eb: ab E exp(-E/tb)
/// - eb < E ≤ ec: 1/E
/// - E > ec: ac √E exp(-E/tc)
///
/// ab と ac は境界で連続になるように決める。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Weight {
    eb: f64,
    tb: f64,
    ab: f64,
    ec: f64,
    tc: f64,
    ac: f64,
}

impl Weight {
    pub fn thermal_1e_fission(eb: f64, tb: f64, ec: f64, tc: f64) -> Self {
        Self {
            eb,
            tb,
            ab: 1.0 / ((-eb / tb).exp() * eb * eb),
            ec,
            tc,
            ac: 1.0 / ((-ec / tc).exp() * ec.powf(1.5)),
        }
    }

    /// 軽水炉向けの値（eb = 0.1 eV, tb = 0.0253 eV, ec = 820.3 keV, tc = 1.4 MeV）。
    pub fn lwr() -> Self {
        Self::thermal_1e_fission(0.1, 0.0253, 820.3e3, 1.4e6)
    }

    pub fn value(&self, e: f64) -> f64 {
        if e <= self.eb {
            self.ab * e * (-e / self.tb).exp()
        } else if e <= self.ec {
            1.0 / e
        } else {
            self.ac * e.sqrt() * (-e / self.tc).exp()
        }
    }

    /// 関数の形が切り替わるエネルギー（積分の区間の境界に入れる）。
    pub fn breakpoints(&self) -> [f64; 2] {
        [self.eb, self.ec]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn continuous_at_breakpoints() {
        let w = Weight::lwr();
        for e in w.breakpoints() {
            let (below, above) = (w.value(e * (1.0 - 1e-12)), w.value(e * (1.0 + 1e-12)));
            assert!((below - above).abs() <= 1e-9 * above, "at {e}");
        }
        assert_eq!(w.value(1.0), 1.0);
    }

    #[test]
    fn shapes() {
        let w = Weight::lwr();
        let r = w.value(0.02) / w.value(0.01);
        assert!((r - 2.0 * (-0.01_f64 / 0.0253).exp()).abs() < 1e-12);
        let r = w.value(4.0e6) / w.value(2.0e6);
        assert!((r - 2.0_f64.sqrt() * (-2.0e6_f64 / 1.4e6).exp()).abs() < 1e-12);
    }
}
