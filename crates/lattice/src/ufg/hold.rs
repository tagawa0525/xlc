//! 自由ガスの散乱源の区間を、減速計算のスイープの間で持つ量の上限（ADR-0213）。
//!
//! 自己遮蔽の全体で持つ量は、入力（`ufg-hold-memory`）が無ければ、自己遮蔽を呼ぶたびに使用可能メモリから求め（[`hold_budget`]）、
//! 同時に解きうる減速計算の数で分ける（[`held_cells`]）。

use std::path::{Path, PathBuf};

use super::free_gas::HELD_CELL_BYTES;

/// 使用可能メモリのうち、OS とほかの処理のために残す量 (B)
const RESERVE: u64 = 2 << 30;

/// 使用可能メモリが読めない
#[derive(Debug, thiserror::Error)]
pub enum MemoryError {
    #[error("cannot read {path}: {source} (give `ufg-hold-memory` to set the held memory)")]
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error(
        "cannot parse the available memory: {0} (give `ufg-hold-memory` to set the held memory)"
    )]
    Parse(String),
    #[error(
        "no cgroup v2 entry in /proc/self/cgroup (give `ufg-hold-memory` to set the held memory)"
    )]
    NoCgroup,
}

/// 使用可能メモリ `available` (B) のうち、自己遮蔽の全体で持つ量 (B): [`RESERVE`] を除いた 8 割
pub fn hold_budget(available: u64) -> u64 {
    (u128::from(available.saturating_sub(RESERVE)) * 4 / 5) as u64
}

/// 自己遮蔽の全体で持つ量 (B): 入力の量 `given`（`ufg-hold-memory`）か、無ければ今の使用可能メモリからの [`hold_budget`]
pub fn hold_total(given: Option<u64>) -> Result<u64, MemoryError> {
    match given {
        Some(bytes) => Ok(bytes),
        None => Ok(hold_budget(available_memory()?)),
    }
}

/// 全体で持つ量 `total` (B) を、同時に解きうる `at_once` 個の減速計算で分けたときの、1 回の減速計算が持つ区間の数
pub fn held_cells(total: u64, at_once: usize) -> usize {
    let cells = total / at_once.max(1) as u64 / HELD_CELL_BYTES as u64;
    usize::try_from(cells).unwrap_or(usize::MAX)
}

/// このプロセスの使用可能メモリ (B): /proc/meminfo の MemAvailable と、このプロセスの cgroup（v2）とその祖先の memory.max の
/// 残り（memory.max − memory.current）の小さい方。memory.max の無い階層（メモリの制御が無い）は制限しない
pub fn available_memory() -> Result<u64, MemoryError> {
    let read = |path: &Path| {
        std::fs::read_to_string(path).map_err(|source| MemoryError::Read {
            path: path.to_path_buf(),
            source,
        })
    };
    let meminfo = read(Path::new("/proc/meminfo"))?;
    let cgroup = read(Path::new("/proc/self/cgroup"))?;
    let relative = cgroup
        .lines()
        .find_map(|line| line.strip_prefix("0::"))
        .ok_or(MemoryError::NoCgroup)?;
    let mut limits = Vec::new();
    let mut dir = Path::new("/sys/fs/cgroup").join(relative.trim_start_matches('/'));
    loop {
        let max = dir.join("memory.max");
        if max.exists() {
            limits.push((read(&max)?, read(&dir.join("memory.current"))?));
        }
        if dir == Path::new("/sys/fs/cgroup") || !dir.pop() {
            break;
        }
    }
    let limits: Vec<(&str, &str)> = limits
        .iter()
        .map(|(m, c)| (m.as_str(), c.as_str()))
        .collect();
    available_from(&meminfo, &limits)
}

/// /proc/meminfo の内容 `meminfo` と、cgroup の（memory.max、memory.current）の内容 `limits` から求めた使用可能メモリ (B)
pub(super) fn available_from(meminfo: &str, limits: &[(&str, &str)]) -> Result<u64, MemoryError> {
    let line = meminfo
        .lines()
        .find_map(|line| line.strip_prefix("MemAvailable:"))
        .ok_or_else(|| MemoryError::Parse("no MemAvailable in /proc/meminfo".into()))?;
    let mut words = line.split_whitespace();
    let (Some(value), Some("kB"), None) = (words.next(), words.next(), words.next()) else {
        return Err(MemoryError::Parse(format!("MemAvailable:{line}")));
    };
    let number = |text: &str| {
        text.trim()
            .parse::<u64>()
            .map_err(|_| MemoryError::Parse(format!("{text:?} is not a byte count")))
    };
    let mut available = number(value)?.saturating_mul(1024);
    for &(max, current) in limits {
        if max.trim() == "max" {
            continue;
        }
        available = available.min(number(max)?.saturating_sub(number(current)?));
    }
    Ok(available)
}
