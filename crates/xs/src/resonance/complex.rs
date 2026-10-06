//! 共鳴の計算に使う複素数。

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct Complex {
    pub(super) re: f64,
    pub(super) im: f64,
}

impl Complex {
    pub(super) fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    pub(super) fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub(super) fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
    pub(super) fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    pub(super) fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
    pub(super) fn recip(self) -> Self {
        let n = self.norm_sqr();
        Self::new(self.re / n, -self.im / n)
    }
}
