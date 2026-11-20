//! 評価のファイルの TAB1 と TAB2 のレコードのうち、ln をとる補間則（INT=3、4、5）の区間に、ln をとる側の値が正でない
//! ものを数える（#0326）。
//!
//! 使い方: survey_log_interpolation <ENDF>...
//!
//! xs の読み込み（`Reader::tab1` の検査）を通さずに、行の形でレコードを見つける: ある行の N1、N2 を NR、NP（NZ）とし、続く
//! 行の整数の対が NBT の狭義の昇順で NP（NZ）に終わり、INT が 1〜6、11〜15、21〜25 なら TAB1 か TAB2 の見出しとする。続く
//! 行がすべて小数点のある実数で x が昇順なら TAB1、そうでなければ TAB2 とし、TAB2 は続く NZ 個の部分レコード（LIST、TAB1、
//! TAB2）の C2 を x（入射エネルギーなど）にする。レコードの種類を MF ごとの形から決めないので、ほかのレコードを TAB1 と
//! 取り違えうる。当たったものは、出力の行の位置で原文を見て確かめる。
//!
//! 当たった区間ごとに 1 行（タブ区切り）を書く: `HIT`、ファイル、MF、MT、見出しの行番号、種類（`TAB1` か `TAB2`）、補間則、
//! 分類（ln y の側の値が負なら `y<0`、0 なら `y=0`、ln x の側の x が正でなければ `x<=0`）、区間の点の番号（1 始まり）と
//! 値。TAB2 の ln y の側は、区間の両端の部分レコードの値（TAB1 は y、LIST は全部の値）を見る。最後に、見つけたレコードの
//! 数（`RECORDS`、MF、種類、ln をとる区間があるか、数）と、`SUMMARY`、MF、MT、種類、補間則、分類、当たったレコードの数、
//! ファイルの数を書く。読めないファイルは `PARSE` の 1 行にする。

use std::collections::{BTreeMap, BTreeSet};

/// 1 行の 6 つの欄（11 桁ずつ）。
fn fields(line: &str) -> [&str; 6] {
    std::array::from_fn(|i| line.get(i * 11..(i + 1) * 11).unwrap_or(""))
}

fn int(field: &str) -> Option<i64> {
    let s = field.trim();
    if s.is_empty() {
        return Some(0);
    }
    s.parse().ok()
}

/// 小数点のある実数の欄（整数の欄と区別するため、空欄と小数点のない欄は受け付けない）。
fn real(field: &str) -> Option<f64> {
    let s = field.trim();
    if !s.contains('.') {
        return None;
    }
    if s.contains(['e', 'E']) {
        return s.parse().ok();
    }
    let exponent = s
        .char_indices()
        .skip(1)
        .find(|&(_, c)| c == '+' || c == '-')
        .map(|(i, _)| i);
    match exponent {
        Some(i) => format!("{}e{}", &s[..i], &s[i..]).parse().ok(),
        None => s.parse().ok(),
    }
}

/// 見出しの C1、C2（空欄も 0 として受け付ける実数）。
fn head_real(field: &str) -> Option<f64> {
    if field.trim().is_empty() {
        Some(0.0)
    } else {
        real(field).or_else(|| int(field).map(|v| v as f64))
    }
}

/// 補間則の基本の番号（3、4、5 のどれかなら Some）。11〜15 と 21〜25 は下 1 桁。
fn log_law(code: i64) -> Option<i64> {
    let base = if code > 10 { code % 10 } else { code };
    matches!(base, 3..=5).then_some(base)
}

fn valid_code(code: i64) -> bool {
    matches!(code, 1..=6 | 11..=15 | 21..=25)
}

/// `start` の行から `count` 個の欄を読む（6 欄で改行）。
fn values<T>(
    lines: &[&str],
    start: usize,
    count: usize,
    parse: impl Fn(&str) -> Option<T>,
) -> Option<Vec<T>> {
    let rows = count.div_ceil(6);
    let mut out = Vec::with_capacity(count);
    for row in lines.get(start..start + rows)? {
        let f = fields(row);
        for field in f.iter().take(count - out.len()) {
            out.push(parse(field)?);
        }
    }
    Some(out)
}

/// 見出しの行 `at` の NR と、NBT・INT の対（NBT が狭義の昇順で `n` に終わる）。
fn ranges(lines: &[&str], at: usize) -> Option<(usize, Vec<(usize, i64)>)> {
    let f = fields(lines.get(at)?);
    if real(f[2]).is_some() || real(f[3]).is_some() {
        return None;
    }
    let nr = usize::try_from(int(f[4])?).ok()?;
    let n = usize::try_from(int(f[5])?).ok()?;
    if nr == 0 || n == 0 || nr > n {
        return None;
    }
    let pairs = values(lines, at + 1, 2 * nr, int)?;
    let mut out = Vec::with_capacity(nr);
    let mut previous = 0;
    for p in pairs.chunks(2) {
        let nbt = usize::try_from(p[0]).ok()?;
        if nbt <= previous || !valid_code(p[1]) {
            return None;
        }
        previous = nbt;
        out.push((nbt, p[1]));
    }
    (previous == n).then_some((n, out))
}

/// 区間 k（0 始まり）の補間則。
fn law_of(ranges: &[(usize, i64)], k: usize) -> i64 {
    ranges
        .iter()
        .find(|&&(nbt, _)| nbt >= k + 2)
        .map_or(2, |&(_, c)| c)
}

enum Record {
    Tab1 {
        ranges: Vec<(usize, i64)>,
        x: Vec<f64>,
        y: Vec<f64>,
    },
    Tab2 {
        ranges: Vec<(usize, i64)>,
        x: Vec<f64>,
        /// 部分レコードごとの値（TAB1 は y、LIST は全部、TAB2 は空）
        y: Vec<Vec<f64>>,
    },
}

/// 行 `at` から始まるレコードと、その行の数。
fn tab_record(lines: &[&str], at: usize) -> Option<(Record, usize)> {
    let (n, ranges) = ranges(lines, at)?;
    let after = at + 1 + (2 * ranges.len()).div_ceil(6);
    if let Some(points) = values(lines, after, 2 * n, real) {
        let x: Vec<f64> = points.iter().step_by(2).copied().collect();
        if x.windows(2).all(|w| w[0] <= w[1]) {
            let y = points.iter().skip(1).step_by(2).copied().collect();
            return Some((
                Record::Tab1 { ranges, x, y },
                after - at + (2 * n).div_ceil(6),
            ));
        }
    }
    // TAB2: NZ 個の部分レコード
    let mut j = after;
    let (mut x, mut y) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        let f = fields(lines.get(j)?);
        x.push(head_real(f[1])?);
        if let Some((record, len)) = tab_record(lines, j) {
            y.push(match record {
                Record::Tab1 { y, .. } => y,
                Record::Tab2 { .. } => Vec::new(),
            });
            j += len;
        } else {
            let npl = usize::try_from(int(f[4])?).ok()?;
            y.push(values(lines, j + 1, npl, real)?);
            j += 1 + npl.div_ceil(6);
        }
    }
    Some((Record::Tab2 { ranges, x, y }, after - at))
}

struct Hit {
    kind: &'static str,
    law: i64,
    class: &'static str,
    point: usize,
    detail: String,
}

fn hits(record: &Record) -> Vec<Hit> {
    let mut out = Vec::new();
    let (kind, ranges, x) = match record {
        Record::Tab1 { ranges, x, .. } => ("TAB1", ranges, x),
        Record::Tab2 { ranges, x, .. } => ("TAB2", ranges, x),
    };
    for k in 0..x.len().saturating_sub(1) {
        let code = law_of(ranges, k);
        let Some(law) = log_law(code) else {
            continue;
        };
        let (x0, x1) = (x[k], x[k + 1]);
        if matches!(law, 3 | 5) && x0 <= 0.0 {
            out.push(Hit {
                kind,
                law: code,
                class: "x<=0",
                point: k + 1,
                detail: format!("x = {x0:e}, {x1:e}"),
            });
        }
        if !matches!(law, 4 | 5) {
            continue;
        }
        let (negative, zero, detail) = match record {
            Record::Tab1 { y, .. } => {
                let (y0, y1) = (y[k], y[k + 1]);
                (
                    y0 < 0.0 || y1 < 0.0,
                    y0 == 0.0 || y1 == 0.0,
                    format!("(x, y) = ({x0:e}, {y0:e}), ({x1:e}, {y1:e})"),
                )
            }
            Record::Tab2 { y, .. } => {
                let ends = y[k].iter().chain(&y[k + 1]);
                let min = ends.clone().copied().fold(f64::INFINITY, f64::min);
                (
                    ends.clone().any(|&v| v < 0.0),
                    ends.clone().any(|&v| v == 0.0),
                    format!("x = {x0:e}, {x1:e}; smallest value {min:e}"),
                )
            }
        };
        let class = if negative {
            "y<0"
        } else if zero {
            "y=0"
        } else {
            continue;
        };
        out.push(Hit {
            kind,
            law: code,
            class,
            point: k + 1,
            detail,
        });
    }
    out
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        return Err("usage: survey_log_interpolation <ENDF>...".into());
    }
    // (MF, MT, 種類, 補間則, 分類) → (レコードの数, ファイル)
    type Key = (u32, u32, &'static str, i64, &'static str);
    let mut summary: BTreeMap<Key, (usize, BTreeSet<String>)> = BTreeMap::new();
    // (MF, 種類, ln をとる区間があるか) → 見つけたレコードの数（当たりの分母）
    let mut found: BTreeMap<(u32, &'static str, bool), usize> = BTreeMap::new();
    for path in &paths {
        let name: Vec<&str> = path.rsplit('/').take(2).collect();
        let name = format!("{}/{}", name[1], name[0]);
        let text = match std::fs::read_to_string(path) {
            Ok(text) => text,
            Err(error) => {
                println!("PARSE\t{name}\t{error}");
                continue;
            }
        };
        let lines: Vec<&str> = text.lines().collect();
        // (MF, MT) の連続した行。MT=0 の区切りの行は含めない
        let mut start = 0;
        while start < lines.len() {
            let tag = |l: &str| -> Option<(u32, u32)> {
                Some((
                    l.get(70..72)?.trim().parse().ok()?,
                    l.get(72..75)?.trim().parse().ok()?,
                ))
            };
            let Some((mf, mt)) = tag(lines[start]) else {
                start += 1;
                continue;
            };
            let end = (start..lines.len())
                .find(|&i| tag(lines[i]) != Some((mf, mt)))
                .unwrap_or(lines.len());
            if mt == 0 || mf == 1 && mt == 451 {
                start = end;
                continue;
            }
            let section = &lines[start..end];
            let mut at = 0;
            while at < section.len() {
                let Some((record, len)) = tab_record(section, at) else {
                    at += 1;
                    continue;
                };
                let (kind, ranges) = match &record {
                    Record::Tab1 { ranges, .. } => ("TAB1", ranges),
                    Record::Tab2 { ranges, .. } => ("TAB2", ranges),
                };
                let logarithmic = ranges.iter().any(|&(_, code)| log_law(code).is_some());
                *found.entry((mf, kind, logarithmic)).or_default() += 1;
                let mut counted = BTreeSet::new();
                for hit in hits(&record) {
                    println!(
                        "HIT\t{name}\t{mf}\t{mt}\t{}\t{}\t{}\t{}\t{}\t{}",
                        start + at + 1,
                        hit.kind,
                        hit.law,
                        hit.class,
                        hit.point,
                        hit.detail
                    );
                    let key = (mf, mt, hit.kind, hit.law, hit.class);
                    if counted.insert(key) {
                        let entry = summary.entry(key).or_default();
                        entry.0 += 1;
                        entry.1.insert(name.clone());
                    }
                }
                at += len;
            }
            start = end;
        }
    }
    for ((mf, kind, logarithmic), count) in &found {
        println!("RECORDS\t{mf}\t{kind}\t{logarithmic}\t{count}");
    }
    for ((mf, mt, kind, law, class), (records, files)) in &summary {
        println!(
            "SUMMARY\t{mf}\t{mt}\t{kind}\t{law}\t{class}\t{records}\t{}",
            files.len()
        );
    }
    Ok(())
}
