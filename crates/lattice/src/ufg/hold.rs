//! 自由ガスの散乱源の区間を、減速計算のスイープの間で持つ量の上限（ADR-0213）。
//!
//! 自己遮蔽の全体で持つ量は、入力（`ufg-hold-memory`）が無ければ、自己遮蔽を呼ぶたびに使用可能メモリから求め（[`hold_budget`]）、
//! 自由ガスの核の表（ADR-0248）の分を除いて、同時に解きうる減速計算の数で分ける（[`hold_total`]、[`held_cells`]）。1 回の減速計算の量は、その減速計算の作業の配列と持つ区間を合わせたもので、
//! 作業の配列の分を除いた残りの区間を持つ（[`cells_beside_work`]）。掃き直す格子点の飛程（`HeldTracks`）は、収まるときだけ作業の
//! 配列に足して区間より先に除く。

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
    #[error("{0} is not in a cgroup v2 tree (give `ufg-hold-memory` to set the held memory)")]
    NotCgroup2(PathBuf),
}

/// 使用可能メモリ `available` (B) のうち、自己遮蔽の全体で持つ量 (B): [`RESERVE`] を除いた 8 割
pub fn hold_budget(available: u64) -> u64 {
    (u128::from(available.saturating_sub(RESERVE)) * 4 / 5) as u64
}

/// 自己遮蔽の全体で持つ量 (B) のうち、前の自己遮蔽から残った自由ガスの核の表 `held` (B。ADR-0248) を除いて新たに確保できる量:
/// 入力の量 `given`（`ufg-hold-memory`）− `held` か、無ければ今の使用可能メモリに `held` を足し戻した [`hold_budget`] − `held`
/// （表は使用可能メモリを減らしているので、足し戻して全体を求め、その中に表を入れる）
pub fn hold_total(given: Option<u64>, held: u64) -> Result<u64, MemoryError> {
    let whole = match given {
        Some(bytes) => bytes,
        None => hold_budget(available_memory()?.saturating_add(held)),
    };
    Ok(whole.saturating_sub(held))
}

/// 自由ガスの核の表（ADR-0248）の上限 (B): 残った表 `held` を含む全体（新たに確保できる量 `total` + `held`）を、同時に解く
/// `at_once` 個の減速計算と表で分けた 1 つ分。ピンと塊には、`total` から表がまだ増えうる量（上限 − `held`、または 0）を除いて渡す
pub fn kernel_cap(total: u64, held: u64, at_once: usize) -> u64 {
    total.saturating_add(held) / (at_once as u64 + 1)
}

/// 全体で持つ量 `total` (B) を、同時に解きうる `at_once` 個の減速計算で分けたときの、1 回の減速計算の量（1 区間の大きさで数えた
/// 数。作業の配列を含む。[`cells_beside_work`]）
pub fn held_cells(total: u64, at_once: usize) -> usize {
    let cells = total / at_once.max(1) as u64 / HELD_CELL_BYTES as u64;
    usize::try_from(cells).unwrap_or(usize::MAX)
}

/// 1 回の減速計算の量 `hold`（[`held_cells`]）から、その減速計算の作業の配列 `work` (B) の分（1 区間の大きさで切り上げ）を除いた、
/// 持てる区間の数。作業の配列が量を超えれば 0
pub(super) fn cells_beside_work(hold: usize, work: u64) -> usize {
    let work = work.div_ceil(HELD_CELL_BYTES as u64);
    usize::try_from(work).map_or(0, |work| hold.saturating_sub(work))
}

/// このプロセスの使用可能メモリ (B): /proc/meminfo の MemAvailable と、このプロセスの cgroup（v2）とその祖先の memory.max と
/// memory.high の残り（上限 − memory.current）の小さい方（[`cgroup_limits`]）
pub fn available_memory() -> Result<u64, MemoryError> {
    let meminfo = read(Path::new("/proc/meminfo"))?;
    let cgroup = read(Path::new("/proc/self/cgroup"))?;
    let relative = cgroup
        .lines()
        .find_map(|line| line.strip_prefix("0::"))
        .ok_or(MemoryError::NoCgroup)?;
    let limits = cgroup_limits(Path::new("/sys/fs/cgroup"), relative)?;
    let limits: Vec<(&str, &str)> = limits
        .iter()
        .map(|(m, c)| (m.as_str(), c.as_str()))
        .collect();
    available_from(&meminfo, &limits)
}

fn read(path: &Path) -> Result<String, MemoryError> {
    std::fs::read_to_string(path).map_err(|source| MemoryError::Read {
        path: path.to_path_buf(),
        source,
    })
}

/// cgroup（v2）の根 `root` の下の `relative`（/proc/self/cgroup の `0::` の後）から根までの階層の、上限（memory.max と memory.high）
/// ごとの（上限、memory.current）の内容。上限のファイルの無い階層（メモリの制御が無い）は読まない。`root` が cgroup v2 の根
/// （cgroup.controllers のあるところ）でないとき、`relative` が根の外を指す（..）とき、ディレクトリを調べられないとき、上限のファイルの
/// 有無を確かめられないときは、制限を読み落とさないように止める
pub(super) fn cgroup_limits(
    root: &Path,
    relative: &str,
) -> Result<Vec<(String, String)>, MemoryError> {
    let relative = Path::new(relative.trim_start_matches('/'));
    let inside = relative
        .components()
        .all(|c| matches!(c, std::path::Component::Normal(_)));
    let exists = |path: &Path| {
        path.try_exists().map_err(|source| MemoryError::Read {
            path: path.to_path_buf(),
            source,
        })
    };
    if !inside || !exists(&root.join("cgroup.controllers"))? {
        return Err(MemoryError::NotCgroup2(root.join(relative)));
    }
    let mut dir = root.join(relative);
    let metadata = std::fs::metadata(&dir).map_err(|source| MemoryError::Read {
        path: dir.clone(),
        source,
    })?;
    if !metadata.is_dir() {
        return Err(MemoryError::NotCgroup2(dir));
    }
    let mut limits = Vec::new();
    loop {
        let mut found = Vec::new();
        for name in ["memory.max", "memory.high"] {
            let limit = dir.join(name);
            if exists(&limit)? {
                found.push(read(&limit)?);
            }
        }
        // 同じ階層の上限は、1 度に読んだ同じ使用量と比べる
        if !found.is_empty() {
            let current = read(&dir.join("memory.current"))?;
            limits.extend(found.into_iter().map(|limit| (limit, current.clone())));
        }
        if dir == root || !dir.pop() {
            break;
        }
    }
    Ok(limits)
}

/// /proc/meminfo の内容 `meminfo` と、cgroup の（上限、memory.current）の内容 `limits` から求めた使用可能メモリ (B)
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
