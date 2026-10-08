//! 背景断面積 σ0 ごとの群平均断面積を、NJOY (GROUPR) の GENDF と比べる。
//!
//! 使い方: compare_njoy_shielding <ENDF> <GENDF> <MAT> <温度 K> [MT...]
//!
//! 群構造は WIMS 69群、重み関数は iwt=4 とする。GENDF は同じ温度に拡げた PENDF から、
//! 同じ群構造・重み関数で、σ0 の表を指定して作っておくこと（σ0 は GENDF の MF1 から読む）。
//! こちらはその温度に拡げたポイントワイズの断面積を使う。反応ごと・σ0 ごとに、全群での相対差の最大値を表示する。
//! 非分離共鳴があれば、その自己遮蔽（UNRESR に相当）も含める。GENDF も UNRESR を通して作っておくこと。
//! 環境変数 XLC_URR_ONLY を設定すると、非分離共鳴の範囲にかかる群だけを比べる。
//! 環境変数 XLC_SLOWING_DOWN を設定すると、分離共鳴の範囲の重みを超微細群の減速計算で作る
//! （GROUPR の flux calculator（重み関数の番号が負）と比べるとき）。GROUPR の sigpot は範囲より上の散乱源に
//! 使うので、こちらの分離共鳴のポテンシャル散乱 4πa²（a は s 波の散乱半径。Reich-Moore で l=0 の APL が
//! 0 でなければそれ）にそろえる。
//! GROUPR は自己遮蔽する反応（MT1, 2, 18, 19, 51, 102 など）だけを σ0 ごとに出力し、
//! それ以外は無限希釈の値だけを出すので、その反応は最初の σ0（無限希釈）だけを比べる。

use std::collections::BTreeMap;

use xlc_format::group::GroupStructure;
use xlc_xs::endf::{Section, Tape, parse_nubar};
use xlc_xs::multigroup::{Weight, shielded_averages};
use xlc_xs::pointwise::{DEFAULT_TOLERANCE, reconstruct_at};
use xlc_xs::process::{slowing_down_range, unresolved_shielding};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, gendf, mat, temperature, mts @ ..] = args.as_slice() else {
        return Err("usage: compare_njoy_shielding <ENDF> <GENDF> <MAT> <K> [MT...]".into());
    };
    let mat: i32 = mat.parse()?;
    let temperature: f64 = temperature.parse()?;
    let tape = Tape::parse(&std::fs::read_to_string(endf)?)?;
    let tables = reconstruct_at(&tape, mat, DEFAULT_TOLERANCE, temperature)?;
    let nubar = match tape.section(mat, 1, 452) {
        Some(section) => Some(parse_nubar(&section)?),
        None => None,
    };

    let reference = Tape::parse(&std::fs::read_to_string(gendf)?)?;
    let header = reference
        .section(mat, 1, 451)
        .ok_or("GENDF has no MF1 MT451")?;
    let sigma0 = read_sigma0(&header)?;
    let groups = GroupStructure::wims69();
    let ng = groups.len();
    let start = std::time::Instant::now();
    let shielding = unresolved_shielding(&tape, mat, temperature, &sigma0)?;
    println!(
        "unresolved shielding in {:.1} s",
        start.elapsed().as_secs_f64()
    );
    let slowing = match std::env::var_os("XLC_SLOWING_DOWN") {
        Some(_) => slowing_down_range(&tape, mat, &groups, &Weight::lwr())?,
        None => None,
    };
    let ours = shielded_averages(
        &groups,
        &Weight::lwr(),
        &tables,
        nubar.as_ref(),
        &sigma0,
        shielding.as_ref(),
        slowing.as_ref(),
    )?;

    let mts: Vec<u32> = if mts.is_empty() {
        reference.mts(mat, 3)
    } else {
        mts.iter().map(|m| m.parse()).collect::<Result<_, _>>()?
    };
    println!(
        "sigma0: {}",
        sigma0
            .iter()
            .map(|s| format!("{s:e}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    // 環境変数 XLC_URR_ONLY を設定すると、非分離共鳴の範囲にかかる群だけを比べる
    let urr_only = std::env::var_os("XLC_URR_ONLY").is_some();
    let in_scope = |g: usize| {
        let (lo, hi) = (groups.bounds()[g + 1], groups.bounds()[g]);
        !urr_only
            || shielding
                .as_ref()
                .is_some_and(|u| lo < u.range().1 && hi > u.range().0)
    };
    for mt in mts {
        let Some(section) = reference.section(mat, 3, mt) else {
            println!("MT{mt}: not in GENDF");
            continue;
        };
        let njoy = read_gendf_mf3(&section, ng, sigma0.len())?;
        let mut line = format!("MT{mt}:");
        let shielded = njoy[0].len();
        for (iz, avg) in ours.iter().enumerate().take(shielded) {
            let values: Vec<f64> = if mt == 452 {
                let (Some(nu_fission), Some(fission)) = (&avg.nu_fission, avg.sigma.get(&18))
                else {
                    return Err("MT452 needs fission".into());
                };
                nu_fission
                    .iter()
                    .zip(fission)
                    .map(|(nf, f)| if *f > 0.0 { nf / f } else { 0.0 })
                    .collect()
            } else {
                let Some(s) = avg.sigma.get(&mt) else {
                    line.push_str(" missing");
                    break;
                };
                s.clone()
            };
            let mut worst = (0.0_f64, 0);
            for g in (0..ng).filter(|&g| in_scope(g)) {
                let expected = njoy[g][iz];
                if expected == 0.0 {
                    continue;
                }
                let rel = (values[g] - expected).abs() / expected.abs();
                if rel > worst.0 {
                    worst = (rel, g + 1);
                }
            }
            line.push_str(&format!(" {:.1e}(g{})", worst.0, worst.1));
        }
        println!("{line}");
    }
    Ok(())
}

/// GENDF の MF1 MT451 から σ0 の表を読む。
fn read_sigma0(section: &Section<'_>) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    let (nz, ntw) = (usize::try_from(head.l2)?, usize::try_from(head.n2)?);
    let list = r.list()?;
    Ok(list.values[ntw..ntw + nz].to_vec())
}

/// GENDF の MF3 の1反応を、[群（高エネルギーから）][σ0] の P0 断面積として読む。
/// 自己遮蔽しない反応は、σ0 が無限希釈の1つだけになる。
fn read_gendf_mf3(
    section: &Section<'_>,
    ng: usize,
    nz_expected: usize,
) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    // 値は NL×NZ 個の中性子束のあとに断面積が並ぶ（ルジャンドル次数が速く変わる）
    let (nl, nz) = (usize::try_from(head.l1)?, usize::try_from(head.l2)?);
    if (nz != nz_expected && nz != 1) || usize::try_from(head.n2)? != ng {
        return Err(format!("unexpected GENDF layout: NZ {nz}, groups {}", head.n2).into());
    }
    let mut values = BTreeMap::new();
    while !r.at_end() {
        let list = r.list()?;
        // GENDF の群番号 IG は低エネルギーから 1 始まり
        let ig = usize::try_from(list.n2)?;
        let row: Vec<f64> = (0..nz).map(|iz| list.values[nl * nz + iz * nl]).collect();
        values.insert(ng - ig, row);
    }
    Ok((0..ng)
        .map(|g| values.get(&g).cloned().unwrap_or(vec![0.0; nz]))
        .collect())
}
