//! ステージの実行ファイルのログのイベントのマクロ（ADR-0255）。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

#[doc(hidden)]
pub use tracing as __tracing;

/// `tracing::trace!` に渡す。
#[macro_export]
macro_rules! trace {
    ($($arg:tt)+) => {
        $crate::__tracing::trace!($($arg)+)
    };
}

/// `tracing::debug!` に渡す。
#[macro_export]
macro_rules! debug {
    ($($arg:tt)+) => {
        $crate::__tracing::debug!($($arg)+)
    };
}

/// `tracing::info!` に渡す。
#[macro_export]
macro_rules! info {
    ($($arg:tt)+) => {
        $crate::__tracing::info!($($arg)+)
    };
}

/// `tracing::warn!` に渡す。
#[macro_export]
macro_rules! warn {
    ($($arg:tt)+) => {
        $crate::__tracing::warn!($($arg)+)
    };
}

/// `tracing::error!` に渡す。
#[macro_export]
macro_rules! error {
    ($($arg:tt)+) => {
        $crate::__tracing::error!($($arg)+)
    };
}

#[cfg(test)]
mod tests;
