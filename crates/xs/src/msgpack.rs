//! データファイルの MessagePack の読み書き（全体計画の「ファイルの約束事」）。
//!
//! フィールド名付きの map 形式で書き、読むときは末尾の余り、文字列でないキー、map でない構造体を拒否する。
//! serde は構造体を配列からも読み、バイナリのキーもフィールド名として受け付けるので、型に変換する前に
//! 値の木で確かめる。入れ子の構造体の位置は型ごとに違うので、`map` と `field` で呼び出し側が確かめる。

use rmpv::Value;
use serde::Serialize;
use serde::de::DeserializeOwned;

/// MessagePack の読み書きのエラー。
#[derive(Debug, thiserror::Error)]
pub enum MsgpackError {
    #[error("cannot encode as MessagePack: {0}")]
    Encode(String),
    #[error("invalid MessagePack: {0}")]
    Decode(String),
    #[error("{0} bytes follow the value")]
    TrailingBytes(usize),
    #[error("structs must be MessagePack maps with string field names")]
    NotMap,
}

/// フィールド名付きの map 形式で書く。
pub fn to_named<T: Serialize>(value: &T) -> Result<Vec<u8>, MsgpackError> {
    rmp_serde::to_vec_named(value).map_err(|e| MsgpackError::Encode(e.to_string()))
}

/// バイト列全体を1つの値として読み、最上位が文字列キーの map であることを確かめる。
pub fn read_value(bytes: &[u8]) -> Result<Value, MsgpackError> {
    let mut rest = bytes;
    let value =
        rmpv::decode::read_value(&mut rest).map_err(|e| MsgpackError::Decode(e.to_string()))?;
    if !rest.is_empty() {
        return Err(MsgpackError::TrailingBytes(rest.len()));
    }
    map(&value)?;
    Ok(value)
}

/// 文字列キーの map の中身。
pub fn map(value: &Value) -> Result<&[(Value, Value)], MsgpackError> {
    value
        .as_map()
        .filter(|entries| entries.iter().all(|(k, _)| k.is_str()))
        .map(Vec::as_slice)
        .ok_or(MsgpackError::NotMap)
}

/// map の中の、キーが `key` の配列（なければ空）。項目の有無と型は、型への変換で確かめる。
pub fn array<'a>(entries: &'a [(Value, Value)], key: &str) -> &'a [Value] {
    entries
        .iter()
        .find(|(k, _)| k.as_str() == Some(key))
        .and_then(|(_, v)| v.as_array())
        .map_or(&[], Vec::as_slice)
}

/// 確かめた値を型に変換する。
pub fn from_value<T: DeserializeOwned>(value: Value) -> Result<T, MsgpackError> {
    rmpv::ext::from_value(value).map_err(|e| MsgpackError::Decode(e.to_string()))
}
