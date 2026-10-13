//! 中性子のサブライブラリの評価ごとに、MF2 の共鳴の範囲の形式と、xs が読めるかを調べる。
//!
//! 使い方: survey_neutrons <中性子の ENDF のディレクトリ> <燃焼鎖の msgpack> [スレッド数]
//!
//! 評価ごとに 1 行（タブ区切り）を書く: ファイル、燃焼鎖の核種の名前（なければ `-`）、範囲（LRU/LRF/NRO/NAPS と el〜eh）、
//! xs の読み込みの結果。読み込みは、`process_nuclide` が多群の処理の前に読むもの（`read_evaluation`。0 K の再構成、MF1 MT452、
//! MF3〜6）、分離共鳴の上端とポテンシャル散乱、非分離共鳴の自己遮蔽（VERA の燃料の温度 565、600、900、1200 K と、
//! examples/vera/xs.toml の σ0 の 9 点）までで、散乱行列と群定数は作らない。
//! 最後に、最初の範囲の形式、NRO ≠ 0 の範囲、読めない理由ごとの数を書く。非分離共鳴の自己遮蔽が時間の大半で、557 評価に
//! 6 スレッドで約 1 時間かかる（2026-10-14、ほかの計算と重なった実行）ので、Slurm で流す。
//!
//! MF2 の範囲の見出しは、xs の `Mf2::parse`（読めない形式でエラーにして止まる）を通さずに、レコードを読み飛ばして
//! すべての範囲について集める。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use xlc_chain::endf::read_material;
use xlc_depletion::chain::Chain;
use xlc_xs::endf::{Reader, Tape};
use xlc_xs::pointwise::DEFAULT_TOLERANCE;
use xlc_xs::process::{read_evaluation, resolved, unresolved_shielding};

/// 非分離共鳴の自己遮蔽を確かめる温度 (K)。VERA の燃料の核種の温度（examples/vera/xs.toml）
const TEMPERATURES: [f64; 4] = [565.0, 600.0, 900.0, 1200.0];

/// 非分離共鳴の自己遮蔽を確かめる σ0 (barn)。examples/vera/xs.toml と同じ
const SIGMA0: [f64; 9] = [1e10, 1e5, 1e4, 1e3, 300.0, 100.0, 30.0, 10.0, 1.0];

/// 核種の (Z, A, 状態) と名前。
type Names = BTreeMap<(u32, u32, u32), String>;

/// MF2 の範囲の見出し。
#[derive(Debug, Clone)]
struct RangeHead {
    lru: i64,
    lrf: i64,
    nro: i64,
    naps: i64,
    el: f64,
    eh: f64,
}

impl RangeHead {
    fn label(&self) -> String {
        let kind = match (self.lru, self.lrf) {
            (0, _) => "radius".to_string(),
            (1, 1) => "SLBW".to_string(),
            (1, 2) => "MLBW".to_string(),
            (1, 3) => "RM".to_string(),
            (1, 4) => "AA".to_string(),
            (1, 7) => "RML".to_string(),
            (2, lrf) => format!("URR{lrf}"),
            (lru, lrf) => format!("LRU{lru}LRF{lrf}"),
        };
        format!(
            "{kind}(NRO={},NAPS={})[{:.4e},{:.4e}]",
            self.nro, self.naps, self.el, self.eh
        )
    }
}

struct Survey {
    file: String,
    /// 燃焼鎖の核種の名前（鎖になければ `None`）
    nuclide: Option<String>,
    ranges: Result<Vec<RangeHead>, String>,
    /// 読めない理由（読めれば空）
    failures: Vec<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let (dir, chain, threads) = match args.as_slice() {
        [_, dir, chain] => (dir, chain, 8),
        [_, dir, chain, threads] => (dir, chain, threads.parse()?),
        _ => {
            return Err(
                "usage: survey_neutrons <neutron ENDF dir> <chain msgpack> [threads]".into(),
            );
        }
    };
    let chain = Chain::from_msgpack(&std::fs::read(chain)?)?;
    let names: Names = chain
        .nuclides
        .iter()
        .map(|n| ((n.z, n.a, n.state), n.name.clone()))
        .collect();
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)?
        .map(|e| e.map(|e| e.path()))
        .collect::<Result<_, _>>()?;
    files.retain(|p| p.extension().is_some_and(|x| x == "endf"));
    files.sort();

    let next = AtomicUsize::new(0);
    let results = Mutex::new(Vec::with_capacity(files.len()));
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(path) = files.get(i) else { break };
                    let survey = survey(path, &names);
                    results
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .push((i, survey));
                }
            });
        }
    });
    let mut results = results.into_inner().unwrap_or_else(|e| e.into_inner());
    results.sort_by_key(|(i, _)| *i);

    println!("file\tnuclide\tranges\tstatus");
    for (_, s) in &results {
        let ranges = match &s.ranges {
            Ok(r) => r.iter().map(RangeHead::label).collect::<Vec<_>>().join(" "),
            Err(e) => format!("unreadable MF2 headers: {e}"),
        };
        let status = if s.failures.is_empty() {
            "ok".to_string()
        } else {
            s.failures.join("; ")
        };
        let nuclide = s.nuclide.as_deref().unwrap_or("-");
        println!("{}\t{nuclide}\t{ranges}\t{status}", s.file);
    }
    summarize(&results.into_iter().map(|(_, s)| s).collect::<Vec<_>>());
    Ok(())
}

fn survey(path: &Path, names: &Names) -> Survey {
    let file = path
        .file_stem()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let failed = |file: String, what: String| Survey {
        file,
        nuclide: None,
        ranges: Err(what.clone()),
        failures: vec![what],
    };
    let tape = match std::fs::read_to_string(path)
        .map_err(|e| e.to_string())
        .and_then(|t| Tape::parse(&t).map_err(|e| e.to_string()))
    {
        Ok(t) => t,
        Err(e) => return failed(file, format!("tape: {e}")),
    };
    let material = match read_material(&tape) {
        Ok(m) => m,
        Err(e) => return failed(file, format!("MF1 MT451: {e}")),
    };
    let nuclide = names
        .get(&(material.z, material.a, material.state))
        .cloned();
    Survey {
        file,
        nuclide,
        ranges: range_heads(&tape, material.mat),
        failures: read_failures(&tape, material.mat),
    }
}

/// MF2 MT151 のすべての範囲の見出し。本体はレコードの形だけを読み飛ばす。
fn range_heads(tape: &Tape, mat: i32) -> Result<Vec<RangeHead>, String> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(Vec::new());
    };
    let mut r = section.reader();
    let e = |e: xlc_xs::endf::Error| e.to_string();
    let head = r.cont().map_err(e)?;
    let mut heads = Vec::new();
    for _ in 0..head.n1 {
        let isotope = r.cont().map_err(e)?;
        let lfw = isotope.l2;
        for _ in 0..isotope.n1 {
            let range = r.cont().map_err(e)?;
            let h = RangeHead {
                lru: range.l1,
                lrf: range.l2,
                nro: range.n1,
                naps: range.n2,
                el: range.c1,
                eh: range.c2,
            };
            skip_range(&mut r, &h, lfw).map_err(|m| format!("{} after {:?}", m, h))?;
            heads.push(h);
        }
    }
    Ok(heads)
}

fn skip_range(r: &mut Reader<'_>, h: &RangeHead, lfw: i64) -> Result<(), String> {
    let e = |e: xlc_xs::endf::Error| e.to_string();
    if h.lru != 0 && h.nro != 0 {
        r.tab1().map_err(e)?;
    }
    match (h.lru, h.lrf) {
        (0, _) => {
            r.cont().map_err(e)?;
        }
        (1, 1..=3) => {
            let cont = r.cont().map_err(e)?;
            for _ in 0..cont.n1 {
                r.list().map_err(e)?;
            }
        }
        (1, 7) => {
            let cont = r.cont().map_err(e)?;
            r.list().map_err(e)?;
            for _ in 0..cont.n1 {
                let channels = r.list().map_err(e)?;
                r.list().map_err(e)?;
                // KBK、KPS（背景の R 行列、位相シフト）は読み飛ばし方を書いていない
                if channels.l1 != 0 || channels.l2 != 0 {
                    return Err(format!(
                        "RML spin group with KBK={} KPS={}",
                        channels.l1, channels.l2
                    ));
                }
            }
        }
        (2, 1) if lfw == 0 => {
            let cont = r.cont().map_err(e)?;
            for _ in 0..cont.n1 {
                r.list().map_err(e)?;
            }
        }
        (2, 1) => {
            let energies = r.list().map_err(e)?;
            for _ in 0..energies.n2 {
                let l = r.cont().map_err(e)?;
                for _ in 0..l.n1 {
                    r.list().map_err(e)?;
                }
            }
        }
        (2, 2) => {
            let cont = r.cont().map_err(e)?;
            for _ in 0..cont.n1 {
                let l = r.cont().map_err(e)?;
                for _ in 0..l.n1 {
                    r.list().map_err(e)?;
                }
            }
        }
        (lru, lrf) => return Err(format!("cannot skip LRU={lru} LRF={lrf}")),
    }
    Ok(())
}

/// xs の読み込みを順に試し、失敗したものの理由を返す。
fn read_failures(tape: &Tape, mat: i32) -> Vec<String> {
    let mut failures = Vec::new();
    if let Err(err) = read_evaluation(tape, DEFAULT_TOLERANCE) {
        failures.push(format!("evaluation: {err}"));
    }
    if let Err(err) = resolved(tape, mat) {
        failures.push(format!("resolved top: {err}"));
    }
    for t in TEMPERATURES {
        if let Err(err) = unresolved_shielding(tape, mat, t, &SIGMA0) {
            failures.push(format!("unresolved at {t} K: {err}"));
        }
    }
    failures
}

fn summarize(surveys: &[Survey]) {
    let count = |f: &dyn Fn(&Survey) -> bool| surveys.iter().filter(|s| f(s)).count();
    println!();
    println!("evaluations: {}", surveys.len());
    println!("in chain: {}", count(&|s| s.nuclide.is_some()));
    println!("readable: {}", count(&|s| s.failures.is_empty()));

    let mut first = BTreeMap::<String, usize>::new();
    let mut nro = BTreeMap::<String, Vec<&str>>::new();
    for s in surveys {
        let Ok(ranges) = &s.ranges else {
            *first.entry("unreadable headers".into()).or_default() += 1;
            continue;
        };
        let key = ranges.first().map_or("no MF2".to_string(), |h| {
            h.label().split('(').next().unwrap_or_default().to_string()
        });
        *first.entry(key).or_default() += 1;
        for h in ranges.iter().filter(|h| h.nro != 0) {
            let key = format!("LRU={} LRF={} NRO={} NAPS={}", h.lru, h.lrf, h.nro, h.naps);
            nro.entry(key).or_default().push(label(s));
        }
    }
    println!("first range:");
    for (k, n) in &first {
        println!("  {k}: {n}");
    }
    println!("ranges with NRO != 0:");
    for (k, names) in &nro {
        println!("  {k}: {} ({})", names.len(), names.join(" "));
    }

    // 理由の分類: エラーの文から行番号と数値を除いた形でまとめる
    let mut reasons = BTreeMap::<String, Vec<&str>>::new();
    for s in surveys {
        for f in &s.failures {
            reasons.entry(reason_class(f)).or_default().push(label(s));
        }
    }
    println!("failures by reason:");
    for (k, names) in &reasons {
        println!("  {k}: {} ({})", names.len(), names.join(" "));
    }
}

/// 評価の表示名（鎖の核種の名前、なければファイル）。
fn label(s: &Survey) -> &str {
    s.nuclide.as_deref().unwrap_or(&s.file)
}

/// 理由の分類: 行番号を除き、数（小数点と指数を含む）の並びを 1 つの `#` にまとめる。
fn reason_class(failure: &str) -> String {
    let without_line = failure
        .split_once("line ")
        .and_then(|(head, rest)| rest.split_once(": ").map(|(_, msg)| format!("{head}{msg}")))
        .unwrap_or_else(|| failure.to_string());
    let mut class = String::with_capacity(without_line.len());
    let mut in_number = false;
    for c in without_line.chars() {
        let numeric = c.is_ascii_digit() || (in_number && matches!(c, '.' | 'e' | '-' | '+'));
        if !numeric {
            class.push(c);
        } else if !in_number {
            class.push('#');
        }
        in_number = numeric;
    }
    class
}
