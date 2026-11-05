//! 熱水力のフィードバック（ADR-0229）: 集合体ごとの閉じた 1 次元の流路と、燃料棒の径方向の定常の熱伝導。

pub mod rod;
pub mod water;

#[cfg(test)]
mod tests;
