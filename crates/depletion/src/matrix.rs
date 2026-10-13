//! 燃焼行列（dn/dt = A n）の疎な形。

use crate::DepletionError;

/// 燃焼行列 A。
#[derive(Debug, Clone)]
pub struct DepletionMatrix {
    n: usize,
}

impl DepletionMatrix {
    /// (行, 列, 値) の並びから作る。
    pub fn from_triplets(
        n: usize,
        _entries: &[(usize, usize, f64)],
    ) -> Result<Self, DepletionError> {
        Ok(Self { n })
    }

    /// 核種の数
    pub fn len(&self) -> usize {
        self.n
    }

    /// 核種がない
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

#[cfg(test)]
mod tests;
