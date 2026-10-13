//! 中性子のサブライブラリの評価ごとに、MF2 の共鳴の範囲の形式と、xs が読めるかを調べる。
//!
//! 使い方: survey_neutrons <中性子の ENDF のディレクトリ> <燃焼鎖の msgpack> [スレッド数]
//!
//! 評価ごとに 1 行（タブ区切り）を書く: 核種、燃焼鎖にあるか、範囲（LRU/LRF/NRO/NAPS と el〜eh）、xs の読み込みの結果。
//! 読み込みは、0 K のポイントワイズの再構成（MF2 と MF3）、`process_nuclide` が読む MF1 MT452、MF4、MF5、MF6、
//! 分離共鳴の上端とポテンシャル散乱、非分離共鳴の自己遮蔽（VERA のライブラリの最も低い温度 565 K、無限希釈の 1 点）までで、多群の処理は走らせない。
//! 最後に、最初の範囲の形式、NRO ≠ 0 の範囲、読めない理由ごとの数を書く。
//!
//! MF2 の範囲の見出しは、xs の `Mf2::parse`（読めない形式でエラーにして止まる）を通さずに、レコードを読み飛ばして
//! すべての範囲について集める。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use xlc_xs::endf::{Mf4, Mf5, Mf6, Reader, Tape, parse_nubar};
use xlc_xs::multigroup::{Emission, emits_neutrons, uses_angular};
use xlc_xs::pointwise::{ColdParts, DEFAULT_TOLERANCE};
use xlc_xs::process::{resolved, unresolved_shielding};

/// 核種の (Z, A, 状態)。
type State = (i64, i64, i64);

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
    name: String,
    in_chain: bool,
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
    let chain = chain_states(Path::new(chain))?;
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
                    let survey = survey(path, &chain);
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

    println!("nuclide\tin_chain\tranges\tstatus");
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
        println!("{}\t{}\t{ranges}\t{status}", s.name, s.in_chain);
    }
    summarize(&results.into_iter().map(|(_, s)| s).collect::<Vec<_>>());
    Ok(())
}

/// 燃焼鎖の核種の (Z, A, 状態)。
fn chain_states(path: &Path) -> Result<Vec<State>, Box<dyn std::error::Error>> {
    let value = xlc_format::msgpack::read_value(&std::fs::read(path)?)?;
    let top = xlc_format::msgpack::map(&value)?;
    let field = |entries: &[(rmpv::Value, rmpv::Value)], key: &str| {
        entries
            .iter()
            .find(|(k, _)| k.as_str() == Some(key))
            .and_then(|(_, v)| v.as_i64())
            .ok_or_else(|| format!("chain nuclide without integer {key}"))
    };
    xlc_format::msgpack::array(top, "nuclides")
        .iter()
        .map(|n| {
            let entries = xlc_format::msgpack::map(n)?;
            Ok((
                field(entries, "z")?,
                field(entries, "a")?,
                field(entries, "state")?,
            ))
        })
        .collect()
}

fn survey(path: &Path, chain: &[State]) -> Survey {
    let name = path
        .file_stem()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let tape = match std::fs::read_to_string(path)
        .map_err(|e| e.to_string())
        .and_then(|t| Tape::parse(&t).map_err(|e| e.to_string()))
    {
        Ok(t) => t,
        Err(e) => {
            return Survey {
                name,
                in_chain: false,
                ranges: Err(e.clone()),
                failures: vec![format!("tape: {e}")],
            };
        }
    };
    let materials = tape.materials();
    let [mat] = materials[..] else {
        return Survey {
            name,
            in_chain: false,
            ranges: Err("not one material".into()),
            failures: vec![format!("materials {materials:?}")],
        };
    };
    let in_chain = identity(&tape, mat).is_some_and(|id| chain.contains(&id));
    let ranges = range_heads(&tape, mat);
    let failures = read_failures(&tape, mat);
    Survey {
        name,
        in_chain,
        ranges,
        failures,
    }
}

/// MF1 MT451 の (Z, A, LISO)。
fn identity(tape: &Tape, mat: i32) -> Option<State> {
    let section = tape.section(mat, 1, 451)?;
    let mut r = section.reader();
    let head = r.cont().ok()?;
    let second = r.cont().ok()?;
    let za = head.c1.round() as i64;
    Some((za / 1000, za % 1000, second.l2))
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
    if let Err(err) = ColdParts::read(tape, mat, DEFAULT_TOLERANCE) {
        failures.push(format!("reconstruct: {err}"));
    }
    if let Some(section) = tape.section(mat, 1, 452)
        && let Err(err) = parse_nubar(&section)
    {
        failures.push(format!("MF1 MT452: {err}"));
    }
    if let Some(section) = tape.section(mat, 5, 18)
        && let Err(err) = Mf5::parse(&section)
    {
        failures.push(format!("MF5 MT18: {err}"));
    }
    let mut emissions = BTreeMap::new();
    for (mt, section) in tape.sections(mat, 6) {
        if mt != 2 && emits_neutrons(mt) {
            match Mf6::parse(&section) {
                Ok(mf6) => {
                    emissions.insert(mt, Emission::Mf6(mf6));
                }
                Err(err) => failures.push(format!("MF6 MT{mt}: {err}")),
            }
        }
    }
    for (mt, section) in tape.sections(mat, 5) {
        if mt != 2 && emits_neutrons(mt) && !emissions.contains_key(&mt) {
            match Mf5::parse(&section) {
                Ok(mf5) => {
                    emissions.insert(mt, Emission::Mf5(mf5));
                }
                Err(err) => failures.push(format!("MF5 MT{mt}: {err}")),
            }
        }
    }
    for (mt, section) in tape.sections(mat, 4) {
        if uses_angular(mt, emissions.get(&mt))
            && let Err(err) = Mf4::parse(&section)
        {
            failures.push(format!("MF4 MT{mt}: {err}"));
        }
    }
    if let Err(err) = resolved(tape, mat) {
        failures.push(format!("resolved top: {err}"));
    }
    if let Err(err) = unresolved_shielding(tape, mat, 565.0, &[1e10]) {
        failures.push(format!("unresolved: {err}"));
    }
    failures
}

fn summarize(surveys: &[Survey]) {
    let count = |f: &dyn Fn(&Survey) -> bool| surveys.iter().filter(|s| f(s)).count();
    println!();
    println!("evaluations: {}", surveys.len());
    println!("in chain: {}", count(&|s| s.in_chain));
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
            nro.entry(key).or_default().push(&s.name);
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
            reasons.entry(reason_class(f)).or_default().push(&s.name);
        }
    }
    println!("failures by reason:");
    for (k, names) in &reasons {
        println!("  {k}: {} ({})", names.len(), names.join(" "));
    }
}

fn reason_class(failure: &str) -> String {
    let without_line = failure
        .split_once("line ")
        .and_then(|(head, rest)| rest.split_once(": ").map(|(_, msg)| format!("{head}{msg}")))
        .unwrap_or_else(|| failure.to_string());
    without_line
        .chars()
        .map(|c| {
            if c.is_ascii_digit() || c == '.' {
                '#'
            } else {
                c
            }
        })
        .collect::<String>()
        .split("##")
        .collect::<Vec<_>>()
        .join("#")
}
