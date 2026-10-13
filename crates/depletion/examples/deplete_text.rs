//! 燃焼行列を標準入力の文字で受け取り、`xlc_depletion::deplete` の結果を 1 行に 1 核種で書く。
//! examples/depletion/compare_openmc.py が OpenMC の CRAM と同じ行列で比べるために呼ぶ。
//!
//! 入力（空白区切り、行の順は問わないが `n` を最初に置く）:
//!
//! ```text
//! n <核種の数>
//! order <16|48>
//! t <時間 s>
//! n0 <数密度 n 個>
//! entry <行> <列> <値>     # 何行でも。行が生まれる核種、列が親（DepletionMatrix::from_triplets）
//! ```

use std::io::Read;

use anyhow::{Context, Result, bail};
use xlc_depletion::{CramOrder, DepletionMatrix, deplete};

fn main() -> Result<()> {
    let mut text = String::new();
    std::io::stdin().read_to_string(&mut text)?;
    let mut n = None;
    let mut order = None;
    let mut t = None;
    let mut n0 = None;
    let mut entries = Vec::new();
    for line in text.lines() {
        let mut words = line.split_whitespace();
        let Some(key) = words.next() else { continue };
        let values: Vec<&str> = words.collect();
        match key {
            "n" => {
                n = Some(
                    values
                        .first()
                        .context("n needs a value")?
                        .parse::<usize>()?,
                )
            }
            "order" => {
                order = Some(match values.first().copied() {
                    Some("16") => CramOrder::Order16,
                    Some("48") => CramOrder::Order48,
                    other => bail!("order must be 16 or 48: {other:?}"),
                })
            }
            "t" => t = Some(values.first().context("t needs a value")?.parse::<f64>()?),
            "n0" => {
                n0 = Some(
                    values
                        .iter()
                        .map(|v| v.parse::<f64>())
                        .collect::<Result<Vec<_>, _>>()?,
                )
            }
            "entry" => {
                let [row, col, value] = values[..] else {
                    bail!("entry needs row, col and value: {line}")
                };
                entries.push((row.parse()?, col.parse()?, value.parse()?));
            }
            other => bail!("unknown key {other}"),
        }
    }
    let n = n.context("missing n")?;
    let a = DepletionMatrix::from_triplets(n, &entries)?;
    let result = deplete(
        &a,
        &n0.context("missing n0")?,
        t.context("missing t")?,
        order.context("missing order")?,
    )?;
    for v in result {
        println!("{v:.17e}");
    }
    Ok(())
}
