//! 燃焼鎖（`xlc_depletion::chain`）と 1 群の反応率から燃焼行列を組み立て、`deplete` の結果を鎖の核種の順に 1 行に
//! 1 核種（名前と数密度）で書く。examples/depletion/compare_chain_depletion.py が、OpenMC の鎖で同じ反応率を解いた結果と
//! 比べるために呼ぶ。
//!
//! 使い方: deplete_chain <鎖の msgpack>（反応率などは標準入力から）
//!
//! 入力（空白区切り、行の順は問わない。`order` と `t` は 1 回ずつ）:
//!
//! ```text
//! order <16|48>
//! t <時間 s>
//! n0 <核種> <数密度>
//! reaction <核種> <MT> <反応率 1/s>         （状態に分けない反応だけ）
//! fission <核種> <核分裂率 1/s> <収率のエネルギーの点ごとの重み…>
//! ```

use std::io::Read;

use anyhow::{Context, Result, bail};
use xlc_depletion::chain::Chain;
use xlc_depletion::{CramOrder, deplete};

/// 1 回だけ書ける項目に値を入れる
fn set<T>(slot: &mut Option<T>, key: &str, value: T) -> Result<()> {
    if slot.is_some() {
        bail!("{key} is given twice");
    }
    *slot = Some(value);
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [path] = &args[..] else {
        bail!("usage: deplete_chain <chain.msgpack>")
    };
    let chain = Chain::from_msgpack(&std::fs::read(path).with_context(|| path.clone())?)?;
    let find = |name: &str| {
        chain
            .nuclides
            .iter()
            .position(|n| n.name == name)
            .with_context(|| format!("{name} is not in the chain"))
    };
    let mut rates = chain.zero_rates();
    let mut n0 = vec![0.0; chain.nuclides.len()];
    let mut order = None;
    let mut t = None;
    // 同じ核種と反応を 2 度書いた行（打ち間違い）は止める
    let mut given = std::collections::HashSet::new();
    let mut once = |key: String| {
        if given.insert(key.clone()) {
            Ok(())
        } else {
            Err(anyhow::anyhow!("{key} is given twice"))
        }
    };
    let mut text = String::new();
    std::io::stdin().read_to_string(&mut text)?;
    for line in text.lines() {
        let words: Vec<&str> = line.split_whitespace().collect();
        match words[..] {
            [] => {}
            ["order", "16"] => set(&mut order, "order", CramOrder::Order16)?,
            ["order", "48"] => set(&mut order, "order", CramOrder::Order48)?,
            ["t", value] => set(&mut t, "t", value.parse::<f64>()?)?,
            ["n0", name, value] => {
                once(format!("n0 {name}"))?;
                n0[find(name)?] = value.parse()?;
            }
            ["reaction", name, mt, value] => {
                let i = find(name)?;
                let mt: u32 = mt.parse()?;
                once(format!("reaction {name} MT{mt}"))?;
                let k = chain.nuclides[i]
                    .reactions
                    .iter()
                    .position(|r| r.mt == mt)
                    .with_context(|| format!("{name} has no MT{mt}"))?;
                if rates[i].reactions[k].len() != 1 {
                    bail!("{name} MT{mt} has several channels: {line}");
                }
                rates[i].reactions[k][0] = value.parse()?;
            }
            ["fission", name, value, ref weights @ ..] => {
                once(format!("fission {name}"))?;
                let i = find(name)?;
                rates[i].fission = value.parse()?;
                rates[i].yield_weights = weights
                    .iter()
                    .map(|w| w.parse::<f64>())
                    .collect::<Result<_, _>>()?;
            }
            _ => bail!("cannot read the line: {line}"),
        }
    }
    let a = chain.matrix(&rates)?;
    let result = deplete(
        &a,
        &n0,
        t.context("missing t")?,
        order.context("missing order")?,
    )?;
    for (nuclide, v) in chain.nuclides.iter().zip(result) {
        println!("{} {v:.17e}", nuclide.name);
    }
    Ok(())
}
