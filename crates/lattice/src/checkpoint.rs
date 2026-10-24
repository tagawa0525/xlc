//! 長い計算（燃焼の点、分岐の状態）を途中で止めても、書き終えた点と状態から再開するためのチェックポイント（ADR-0196）。
//!
//! フィールド名付きの map 形式の MessagePack: `format` ("xlc-lattice-checkpoint"), `fingerprint`（`Fingerprint`）と、燃焼なら
//! `depletion: {history（`xlc-depletion-history` のバイト列）, next（次の点の組成 [燃焼の領域][チェーンの核種]）, start（直前の
//! 輸送の MOC の収束した解。ADR-0195）}`、分岐なら `states: {states: [{k, outputs（[出力][部分ノード] の `xlc-table` の状態 1 つの
//! バイト列）, medium（細かい群の媒質か nil）}], starts: {state（解いた状態の数）, solutions（MOC の解。同じ解は 1 つ）, slots（[軸] の解の番号か nil）}}`。

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, anyhow, bail};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use xlc_format::fewgroup::FewGroup;
use xlc_format::medium::MediumState;
use xlc_format::msgpack;
use xlc_format::table::Table;

use crate::branches::Starts;
use crate::depletion::{History, Partial};
use crate::moc::Start;

const FORMAT: &str = "xlc-lattice-checkpoint";

/// 計算の結果を決めるもの（実行ファイル、入力、読むファイル、スレッドの数）の SHA-256
#[derive(Default)]
pub struct Fingerprint(Sha256);

impl Fingerprint {
    /// `label` の部分として `bytes` を足す。名前と長さを前置して、境界をずらした別の組と同じ hash にしない
    pub fn add(&mut self, label: &str, bytes: &[u8]) {
        for part in [label.as_bytes(), bytes] {
            self.0.update((part.len() as u64).to_le_bytes());
            self.0.update(part);
        }
    }

    /// 16 進の文字列
    pub fn finish(self) -> String {
        self.0
            .finalize()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    }
}

/// `bytes` の SHA-256（大きなファイルを指紋に足すときの要約。読むスレッドで求められる）
pub fn digest(bytes: &[u8]) -> Vec<u8> {
    Sha256::digest(bytes).to_vec()
}

/// 分岐の 1 つの状態の結果
#[derive(Debug, Clone, PartialEq)]
pub struct SavedState {
    pub k: f64,
    /// [出力][部分ノード] の少数群定数（少数群定数の出力がなければ空）
    pub outputs: Vec<Vec<FewGroup>>,
    /// 細かい群の媒質の出力（ADR-0128）があるときの媒質
    pub medium: Option<MediumState>,
}

/// チェックポイントのファイルと、今の実行の指紋
pub struct Checkpoint {
    path: PathBuf,
    fingerprint: String,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CheckpointFile {
    format: String,
    fingerprint: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    depletion: Option<DepletionFile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    states: Option<StatesFile>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DepletionFile {
    history: rmpv::Value,
    next: Vec<Vec<f64>>,
    start: Start,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StatesFile {
    states: Vec<StateFile>,
    starts: StartsFile,
}

/// 分岐の MOC の反復を始める解（`branches::Starts`）
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StartsFile {
    state: usize,
    /// 残した解（同じ解は 1 つ）
    solutions: Vec<rmpv::Value>,
    /// [軸] の解の番号か nil
    slots: Vec<Option<usize>>,
}

/// 書く分岐のチェックポイント（解は大きいので複製せず、参照のまま書く）
#[derive(Serialize)]
struct StatesFileOut<'a, T> {
    format: &'static str,
    fingerprint: &'a str,
    states: StatesOut<'a, T>,
}

#[derive(Serialize)]
struct StatesOut<'a, T> {
    states: Vec<StateFile>,
    starts: StartsOut<'a, T>,
}

#[derive(Serialize)]
struct StartsOut<'a, T> {
    state: usize,
    solutions: Vec<&'a T>,
    slots: Vec<Option<usize>>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StateFile {
    k: f64,
    outputs: Vec<Vec<rmpv::Value>>,
    /// nil を書いても、キーは省けない（`deserialize_with` を付けると、serde は省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "xlc_format::msgpack::nil_or")]
    medium: Option<MediumState>,
}

impl Checkpoint {
    pub fn new(path: PathBuf, fingerprint: String) -> Self {
        Self { path, fingerprint }
    }

    /// 入力ファイルのチェックポイントのパス: 入力の拡張子を `checkpoint` にしたもの
    pub fn beside(input: &Path) -> PathBuf {
        input.with_extension("checkpoint")
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// 前の実行が残したファイルを読む。なければ `None`。指紋が今の実行と違えばエラー
    fn read(&self) -> anyhow::Result<Option<CheckpointFile>> {
        let bytes = match std::fs::read(&self.path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(error).with_context(|| self.shown()),
        };
        let file: CheckpointFile = msgpack::read_value(&bytes)
            .and_then(msgpack::from_value)
            .with_context(|| format!("{}: unreadable; remove it to start over", self.shown()))?;
        if file.format != FORMAT {
            bail!(
                "{}: expected format {FORMAT}, found {}",
                self.shown(),
                file.format
            );
        }
        if file.fingerprint != self.fingerprint {
            bail!(
                "{}: written by another executable, input, input file or thread count; remove it to start over",
                self.shown()
            );
        }
        match (&file.depletion, &file.states) {
            (Some(_), Some(_)) => {
                bail!("{}: holds both a depletion and branch states", self.shown())
            }
            (None, None) => bail!(
                "{}: holds neither a depletion nor branch states",
                self.shown()
            ),
            _ => Ok(Some(file)),
        }
    }

    /// 一時ファイルに書いて `fsync` し、rename で置き換える（殺されても、ノードが止まっても、前の版か新しい版が残る）
    fn write(&self, file: &impl Serialize) -> anyhow::Result<()> {
        let start = std::time::Instant::now();
        let bytes = msgpack::to_named(file)?;
        let temporary = self.temporary();
        let shown = || temporary.display().to_string();
        let mut out = File::create(&temporary).with_context(shown)?;
        out.write_all(&bytes).with_context(shown)?;
        out.sync_all().with_context(shown)?;
        std::fs::rename(&temporary, &self.path).with_context(|| self.shown())?;
        sync_directory(&self.path)?;
        tracing::info!(
            bytes = bytes.len(),
            seconds = start.elapsed().as_secs_f64(),
            "checkpoint"
        );
        Ok(())
    }

    fn temporary(&self) -> PathBuf {
        let mut name = self.path.as_os_str().to_owned();
        name.push(".tmp");
        PathBuf::from(name)
    }

    fn shown(&self) -> String {
        format!("checkpoint {}", self.path.display())
    }

    /// 燃焼の途中の履歴。なければ `None`
    pub fn depletion(&self) -> anyhow::Result<Option<Partial>> {
        let Some(file) = self.read()? else {
            return Ok(None);
        };
        let Some(DepletionFile {
            history,
            next,
            start,
        }) = file.depletion
        else {
            bail!("{}: holds branch states, not a depletion", self.shown());
        };
        let rmpv::Value::Binary(bytes) = history else {
            bail!("{}: the history must be binary", self.shown());
        };
        let history = History::from_msgpack(&bytes).with_context(|| self.shown())?;
        Ok(Some(Partial {
            history,
            next,
            start,
        }))
    }

    /// 燃焼の解き終えた点の履歴 `history`、次の点の組成 `next`、直前の輸送の MOC の解 `start` を書く
    pub fn save_depletion(
        &self,
        history: &History,
        next: &[Vec<f64>],
        start: &Start,
    ) -> anyhow::Result<()> {
        // `Start` は大きいので複製せず、参照のまま書く
        #[derive(Serialize)]
        struct DepletionOut<'a> {
            history: rmpv::Value,
            next: &'a [Vec<f64>],
            start: &'a Start,
        }
        #[derive(Serialize)]
        struct FileOut<'a> {
            format: &'static str,
            fingerprint: &'a str,
            depletion: DepletionOut<'a>,
        }
        self.write(&FileOut {
            format: FORMAT,
            fingerprint: &self.fingerprint,
            depletion: DepletionOut {
                history: rmpv::Value::Binary(history.to_msgpack()?),
                next,
                start,
            },
        })
    }

    /// 分岐の解き終えた状態（解く状態の順）と、次に解く状態の MOC の反復を始める解（`counts` は [軸] の格子点の数、`present` は
    /// 解く状態の、軸の直積の順の番号）。なければ空と最初の状態から
    pub fn states<T: DeserializeOwned>(
        &self,
        counts: &[usize],
        present: &[usize],
    ) -> anyhow::Result<(Vec<SavedState>, Starts<T>)> {
        let Some(file) = self.read()? else {
            return Ok((Vec::new(), Starts::over(counts, present.to_vec())));
        };
        let Some(StatesFile { states, starts }) = file.states else {
            bail!("{}: holds a depletion, not branch states", self.shown());
        };
        if starts.state != states.len() {
            bail!(
                "{}: {} states with the starts of state {}",
                self.shown(),
                states.len(),
                starts.state
            );
        }
        let solutions = starts
            .solutions
            .into_iter()
            .map(rmpv::ext::from_value)
            .collect::<Result<Vec<T>, _>>()
            .map_err(|e| anyhow!("{}: starts: {e}", self.shown()))?;
        let starts = Starts::restored(
            counts,
            present.to_vec(),
            starts.state,
            solutions,
            &starts.slots,
        )
        .map_err(|e| anyhow!("{}: {e}; remove it to start over", self.shown()))?;
        let states = states
            .into_iter()
            .enumerate()
            .map(|(i, s)| -> anyhow::Result<SavedState> {
                let outputs = s
                    .outputs
                    .into_iter()
                    .map(|subs| subs.into_iter().map(few_group).collect())
                    .collect::<anyhow::Result<Vec<Vec<FewGroup>>>>();
                let outputs = outputs.with_context(|| format!("{}: state {i}", self.shown()))?;
                Ok(SavedState {
                    k: s.k,
                    outputs,
                    medium: s.medium,
                })
            })
            .collect::<anyhow::Result<_>>()?;
        Ok((states, starts))
    }

    /// 分岐の解き終えた状態 `states`（状態 0 から）と、次に解く状態の MOC の反復を始める解 `starts` を書く
    pub fn save_states<T: Serialize>(
        &self,
        states: &[SavedState],
        starts: &Starts<T>,
    ) -> anyhow::Result<()> {
        let states = states
            .iter()
            .map(|s| -> anyhow::Result<StateFile> {
                let outputs = s
                    .outputs
                    .iter()
                    .map(|subs| {
                        subs.iter()
                            .map(|few| {
                                Ok(rmpv::Value::Binary(
                                    Table::single(few.clone()).to_msgpack()?,
                                ))
                            })
                            .collect::<anyhow::Result<_>>()
                    })
                    .collect::<anyhow::Result<_>>()?;
                Ok(StateFile {
                    k: s.k,
                    outputs,
                    medium: s.medium.clone(),
                })
            })
            .collect::<anyhow::Result<_>>()?;
        let (state, solutions, slots) = starts.saved();
        self.write(&StatesFileOut {
            format: FORMAT,
            fingerprint: &self.fingerprint,
            states: StatesOut {
                states,
                starts: StartsOut {
                    state,
                    solutions,
                    slots,
                },
            },
        })
    }

    /// 計算を終えて出力を書いたあとに消す（一時ファイルも。どちらもなくてよい）
    pub fn remove(&self) -> anyhow::Result<()> {
        for path in [self.path.clone(), self.temporary()] {
            match std::fs::remove_file(&path) {
                Err(error) if error.kind() != std::io::ErrorKind::NotFound => {
                    return Err(error).with_context(|| path.display().to_string());
                }
                _ => {}
            }
        }
        Ok(())
    }
}

/// `path` のあるディレクトリを `fsync` する（rename や unlink を、ノードが止まっても残す）
fn sync_directory(path: &Path) -> anyhow::Result<()> {
    let directory = match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    };
    File::open(directory)
        .and_then(|d| d.sync_all())
        .with_context(|| directory.display().to_string())
}

/// 出力を書いて `fsync` する。チェックポイントを消す前に書く出力は、ノードが止まっても残す（消した後に中身のない出力が残ると、
/// 計算をすべて失う）
pub fn write_synced(path: &Path, bytes: &[u8]) -> anyhow::Result<()> {
    let shown = || path.display().to_string();
    let mut out = File::create(path).with_context(shown)?;
    out.write_all(bytes).with_context(shown)?;
    out.sync_all().with_context(shown)?;
    sync_directory(path)
}

/// 状態 1 つの `xlc-table` のバイト列の少数群定数
fn few_group(value: rmpv::Value) -> anyhow::Result<FewGroup> {
    let rmpv::Value::Binary(bytes) = value else {
        bail!("the few-group constants must be binary");
    };
    let table = Table::from_msgpack(&bytes)?;
    match (table.axes.is_empty(), table.states.as_slice()) {
        (true, [one]) => Ok(one.clone()),
        _ => Err(anyhow!(
            "the few-group constants of a state must be one state"
        )),
    }
}

/// 分岐の状態のうち `present`（[軸] の格子点の数 `counts` の直積の順の番号。増える順）を順に解く。`solve` は状態を、MOC の反復を始める解
/// （`branches::Starts`。ADR-0195）から解き、結果と収束した解を返す。`checkpoint` があれば、書いてある状態を `restored` に渡して
/// （今の実行の出力に合わなければ `restored` がエラーを返し、残りの状態を解く前に止まる）使い、書いてある始める解から続きの
/// 状態を解いて、状態を解き終えるたびに書く
pub fn solve_states<T: Serialize + DeserializeOwned>(
    present: &[usize],
    counts: &[usize],
    checkpoint: Option<&Checkpoint>,
    mut restored: impl FnMut(usize, &SavedState) -> anyhow::Result<()>,
    mut solve: impl FnMut(usize, Option<&T>) -> anyhow::Result<(SavedState, T)>,
) -> anyhow::Result<Vec<SavedState>> {
    let count = present.len();
    let (mut states, mut starts) = match checkpoint {
        Some(c) => c.states(counts, present)?,
        None => (Vec::new(), Starts::over(counts, present.to_vec())),
    };
    if states.len() > count {
        bail!(
            "{}: {} states for a run of {count}",
            checkpoint.map_or_else(String::new, Checkpoint::shown),
            states.len()
        );
    }
    if let (Some(c), false) = (checkpoint, states.is_empty()) {
        tracing::info!(
            checkpoint = %c.path.display(),
            states = states.len(),
            "resuming after the saved states"
        );
    }
    for (&state, saved) in present.iter().zip(&states) {
        restored(state, saved).with_context(|| {
            format!(
                "{}: state {state}",
                checkpoint.map_or_else(String::new, Checkpoint::shown)
            )
        })?;
    }
    for (n, &state) in present.iter().enumerate().skip(states.len()) {
        let (saved, converged) = solve(state, starts.start())?;
        states.push(saved);
        starts.solved(converged).map_err(|e| anyhow!(e))?;
        // 最後の状態の後は、呼ぶ側が出力を書いてチェックポイントを消す
        if let (Some(c), true) = (checkpoint, n + 1 < count) {
            c.save_states(&states, &starts)?;
        }
    }
    Ok(states)
}

#[cfg(test)]
mod tests;
