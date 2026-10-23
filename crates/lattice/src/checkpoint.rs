//! 長い計算（燃焼の点、分岐の状態）を途中で止めても、書き終えた点と状態から再開するためのチェックポイント（ADR-0196）。
//!
//! フィールド名付きの map 形式の MessagePack: `format` ("xlc-lattice-checkpoint"), `fingerprint`（`Fingerprint`）と、燃焼なら
//! `depletion: {history（`xlc-depletion-history` のバイト列）, next（次の点の組成 [燃焼の領域][チェーンの核種]）}`、分岐なら
//! `states: [{k, outputs（[出力][部分ノード] の `xlc-table` の状態 1 つのバイト列）, medium（細かい群の媒質か nil）}]`。

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, anyhow, bail};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use xlc_format::fewgroup::FewGroup;
use xlc_format::medium::MediumState;
use xlc_format::msgpack;
use xlc_format::table::Table;

use crate::depletion::{History, Partial};

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
    states: Option<Vec<StateFile>>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DepletionFile {
    history: rmpv::Value,
    next: Vec<Vec<f64>>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StateFile {
    k: f64,
    outputs: Vec<Vec<rmpv::Value>>,
    /// nil を書いても、キーは省けない（`deserialize_with` を付けると、serde は省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "xlc_format::msgpack::nil_or")]
    medium: Option<MediumFile>,
}

/// 媒質の 1 つの状態（`xlc-medium` の状態と同じ項目）
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MediumFile {
    total: Vec<f64>,
    scatter: Vec<Vec<Vec<f64>>>,
    fission: Vec<f64>,
    nu_fission: Vec<f64>,
    chi: Vec<f64>,
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
        let file: CheckpointFile =
            msgpack::from_value(msgpack::read_value(&bytes)?).with_context(|| self.shown())?;
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
        Ok(Some(file))
    }

    /// 一時ファイルに書いて `fsync` し、rename で置き換える（殺されても、ノードが止まっても、前の版か新しい版が残る）
    fn write(&self, file: &CheckpointFile) -> anyhow::Result<()> {
        let start = std::time::Instant::now();
        let bytes = msgpack::to_named(file)?;
        let temporary = self.temporary();
        let shown = || temporary.display().to_string();
        let mut out = File::create(&temporary).with_context(shown)?;
        out.write_all(&bytes).with_context(shown)?;
        out.sync_all().with_context(shown)?;
        std::fs::rename(&temporary, &self.path).with_context(|| self.shown())?;
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

    fn file(&self) -> CheckpointFile {
        CheckpointFile {
            format: FORMAT.into(),
            fingerprint: self.fingerprint.clone(),
            depletion: None,
            states: None,
        }
    }

    /// 燃焼の途中の履歴。なければ `None`
    pub fn depletion(&self) -> anyhow::Result<Option<Partial>> {
        let Some(file) = self.read()? else {
            return Ok(None);
        };
        let Some(DepletionFile { history, next }) = file.depletion else {
            bail!("{}: holds branch states, not a depletion", self.shown());
        };
        let rmpv::Value::Binary(bytes) = history else {
            bail!("{}: the history must be binary", self.shown());
        };
        let history = History::from_msgpack(&bytes).with_context(|| self.shown())?;
        Ok(Some(Partial { history, next }))
    }

    /// 燃焼の解き終えた点の履歴 `history` と、次の点の組成 `next` を書く
    pub fn save_depletion(&self, history: &History, next: &[Vec<f64>]) -> anyhow::Result<()> {
        self.write(&CheckpointFile {
            depletion: Some(DepletionFile {
                history: rmpv::Value::Binary(history.to_msgpack()?),
                next: next.to_vec(),
            }),
            ..self.file()
        })
    }

    /// 分岐の解き終えた状態（状態の順）。なければ空
    pub fn states(&self) -> anyhow::Result<Vec<SavedState>> {
        let Some(file) = self.read()? else {
            return Ok(Vec::new());
        };
        let Some(states) = file.states else {
            bail!("{}: holds a depletion, not branch states", self.shown());
        };
        states
            .into_iter()
            .enumerate()
            .map(|(i, s)| -> anyhow::Result<SavedState> {
                let outputs = s
                    .outputs
                    .into_iter()
                    .map(|subs| subs.into_iter().map(few_group).collect())
                    .collect::<anyhow::Result<Vec<Vec<FewGroup>>>>();
                let medium = s.medium.map(|m| MediumState {
                    total: m.total,
                    scatter: m.scatter,
                    fission: m.fission,
                    nu_fission: m.nu_fission,
                    chi: m.chi,
                });
                let outputs = outputs.with_context(|| format!("{}: state {i}", self.shown()))?;
                Ok(SavedState {
                    k: s.k,
                    outputs,
                    medium,
                })
            })
            .collect()
    }

    /// 分岐の解き終えた状態 `states`（状態 0 から）を書く
    pub fn save_states(&self, states: &[SavedState]) -> anyhow::Result<()> {
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
                let medium = s.medium.as_ref().map(|m| MediumFile {
                    total: m.total.clone(),
                    scatter: m.scatter.clone(),
                    fission: m.fission.clone(),
                    nu_fission: m.nu_fission.clone(),
                    chi: m.chi.clone(),
                });
                Ok(StateFile {
                    k: s.k,
                    outputs,
                    medium,
                })
            })
            .collect::<anyhow::Result<_>>()?;
        self.write(&CheckpointFile {
            states: Some(states),
            ..self.file()
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

/// 分岐の `count` 個の状態を順に解く。`checkpoint` があれば、書いてある状態を `restored` に渡して使い、続きの状態から `solve`
/// で解いて、状態を解き終えるたびに書く
pub fn solve_states(
    count: usize,
    checkpoint: Option<&Checkpoint>,
    mut restored: impl FnMut(usize, &SavedState),
    mut solve: impl FnMut(usize) -> anyhow::Result<SavedState>,
) -> anyhow::Result<Vec<SavedState>> {
    let mut states = match checkpoint {
        Some(c) => c.states()?,
        None => Vec::new(),
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
    for (state, saved) in states.iter().enumerate() {
        restored(state, saved);
    }
    for state in states.len()..count {
        states.push(solve(state)?);
        if let Some(c) = checkpoint {
            c.save_states(&states)?;
        }
    }
    Ok(states)
}

#[cfg(test)]
mod tests;
