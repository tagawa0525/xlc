//! ステージの実行ファイルのログのイベントのマクロ（ADR-0255）。
//!
//! `tracing` のイベントのマクロ（`info!` など）は、名前を指定しないと `"event <ファイル>:<行>"` にするので、注釈の行を
//! 足すだけで実行ファイルの文字列が変わる。名前の文字列の長さが定数の大きさ（32 バイトなど）に合うと、リンカが
//! まとめる定数の節に入って中身で並びが決まり、行の番号を 0 にしたコピー（scripts/normalize_binary.py）でもコードの
//! 参照の位置の差が残る。ここのマクロは、名前を `"event <ファイル>"` にして `tracing::event!` に渡す。引数は
//! `tracing` のマクロのうち parent を除いた形（target、フィールド、書式）で、腕の分け方は `tracing::info!` の
//! `name:` の腕と同じ。
//!
//! `tracing::info!` などに `name:` を付けて呼ばないのは、clippy の `disallowed_macros` が展開の中の呼び出しも
//! 止め、`#[allow]` を展開の中に置いても効かないため（clippy.toml が `tracing::info` などを止める）。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

#[doc(hidden)]
pub use tracing as __tracing;

/// target を `$target`、水準を `$lvl` にしたイベントを、名前に行の番号を含めずに出す。水準ごとのマクロの中身。
#[doc(hidden)]
#[macro_export]
macro_rules! __event {
    ($target:expr, $lvl:expr, { $($field:tt)* }, $($arg:tt)*) => {
        $crate::__tracing::event!(name: $crate::__name!(), target: $target, $lvl, { $($field)* }, $($arg)*)
    };
    ($target:expr, $lvl:expr, $($k:ident).+ $($field:tt)*) => {
        $crate::__tracing::event!(name: $crate::__name!(), target: $target, $lvl, { $($k).+ $($field)* })
    };
    ($target:expr, $lvl:expr, ?$($k:ident).+ $($field:tt)*) => {
        $crate::__tracing::event!(name: $crate::__name!(), target: $target, $lvl, { ?$($k).+ $($field)* })
    };
    ($target:expr, $lvl:expr, %$($k:ident).+ $($field:tt)*) => {
        $crate::__tracing::event!(name: $crate::__name!(), target: $target, $lvl, { %$($k).+ $($field)* })
    };
    ($target:expr, $lvl:expr, $($arg:tt)+) => {
        $crate::__tracing::event!(name: $crate::__name!(), target: $target, $lvl, {}, $($arg)+)
    };
}

/// イベントの名前。`file!()` は最初のマクロの呼び出し（呼ぶ側のファイル）を指す
#[doc(hidden)]
#[macro_export]
macro_rules! __name {
    () => {
        ::core::concat!("event ", ::core::file!())
    };
}

/// `tracing::trace!` と同じ引数（parent を除く）で、名前に行の番号を含まないイベントを出す。
#[macro_export]
macro_rules! trace {
    (target: $target:expr, $($arg:tt)+) => {
        $crate::__event!($target, $crate::__tracing::Level::TRACE, $($arg)+)
    };
    ($($arg:tt)+) => {
        $crate::__event!(::core::module_path!(), $crate::__tracing::Level::TRACE, $($arg)+)
    };
}

/// `tracing::debug!` と同じ引数（parent を除く）で、名前に行の番号を含まないイベントを出す。
#[macro_export]
macro_rules! debug {
    (target: $target:expr, $($arg:tt)+) => {
        $crate::__event!($target, $crate::__tracing::Level::DEBUG, $($arg)+)
    };
    ($($arg:tt)+) => {
        $crate::__event!(::core::module_path!(), $crate::__tracing::Level::DEBUG, $($arg)+)
    };
}

/// `tracing::info!` と同じ引数（parent を除く）で、名前に行の番号を含まないイベントを出す。
#[macro_export]
macro_rules! info {
    (target: $target:expr, $($arg:tt)+) => {
        $crate::__event!($target, $crate::__tracing::Level::INFO, $($arg)+)
    };
    ($($arg:tt)+) => {
        $crate::__event!(::core::module_path!(), $crate::__tracing::Level::INFO, $($arg)+)
    };
}

/// `tracing::warn!` と同じ引数（parent を除く）で、名前に行の番号を含まないイベントを出す。
#[macro_export]
macro_rules! warn {
    (target: $target:expr, $($arg:tt)+) => {
        $crate::__event!($target, $crate::__tracing::Level::WARN, $($arg)+)
    };
    ($($arg:tt)+) => {
        $crate::__event!(::core::module_path!(), $crate::__tracing::Level::WARN, $($arg)+)
    };
}

/// `tracing::error!` と同じ引数（parent を除く）で、名前に行の番号を含まないイベントを出す。
#[macro_export]
macro_rules! error {
    (target: $target:expr, $($arg:tt)+) => {
        $crate::__event!($target, $crate::__tracing::Level::ERROR, $($arg)+)
    };
    ($($arg:tt)+) => {
        $crate::__event!(::core::module_path!(), $crate::__tracing::Level::ERROR, $($arg)+)
    };
}

#[cfg(test)]
mod tests;
