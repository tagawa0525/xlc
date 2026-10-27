//! UFG の S(α, β) の移行の表（ADR-0205）の、格子点の間の補間の誤差を測る。入射を格子点の間（レサジーの中点）に置き、表の累積の
//! 割合を、同じ S(α, β) の核の出射がそれ以下の割合と比べる。
//!
//! 使い方: sab_transfer_error <熱中性子散乱の ENDF> <温度 K> <UFG の下限 eV> <入射の間隔（レサジー）> <移行の間隔 eV>
//!
//! 出力: 入射ごとに、累積の割合の差の最大（出射は UFG の下限から上方散乱の届く幅までの 1000 点）と、WIMS 69 群の出射の群ごとの
//! 割合の相対差の最大（UFG の下限より上の、割合が 1e−4 以上の群）。最後に全体の最大と、表を作る時間。

use xlc_format::group::GroupStructure;
use xlc_xs::endf::Tape;
use xlc_xs::multigroup::{SabOutside, ThermalSab};
use xlc_xs::process::read_thermal_scattering;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, tsl, t, low, du, dt] = args.as_slice() else {
        return Err("usage: sab_transfer_error <TSL ENDF> <K> <low eV> <du> <dt eV>".into());
    };
    let (temperature, low, du, dt): (f64, f64, f64, f64) =
        (t.parse()?, low.parse()?, du.parse()?, dt.parse()?);
    let data = read_thermal_scattering(&Tape::parse(&std::fs::read_to_string(tsl)?)?)?;
    let sab = ThermalSab::new(
        &data,
        temperature,
        data.e_max,
        SabOutside::ShortCollisionTime,
    )?;
    let started = std::time::Instant::now();
    let table = sab
        .ufg_transfer_on(low, du, dt)?
        .ok_or("the S(α, β) limit is below the UFG reach")?;
    let seconds = started.elapsed().as_secs_f64();
    let reach = table.transfer[table.transfer.len() - 1];
    let groups = GroupStructure::wims69();
    let (mut worst_share, mut worst_group) = (0.0_f64, 0.0_f64);
    for w in table.incident.windows(2) {
        let e = (w[0] * w[1]).sqrt();
        if e < low {
            continue;
        }
        // 出射が UFG の下限から上方散乱の届く幅までの 1001 点を境にした群（最下位の群は 1e−5 eV から下限まで。群の外は端の
        // 群に入る）で核を積分し、下から足した割合を、表の累積と比べる
        let points: Vec<f64> = (0..=1000)
            .map(|k| ((low - e) + (reach - (low - e)) * k as f64 / 1000.0).min(reach))
            .collect();
        let bounds: Vec<f64> = points
            .iter()
            .rev()
            .map(|&t| e + t)
            .chain(std::iter::once(1.0e-5))
            .collect();
        let masses = sab.group_integrals(&GroupStructure::new(bounds)?, e);
        let total: f64 = masses.iter().sum();
        let mut below = masses[masses.len() - 1];
        let mut largest = 0.0_f64;
        for (k, &t) in points.iter().enumerate() {
            if k > 0 {
                below += masses[masses.len() - 1 - k];
            }
            let ours = table.cumulative_at(e, t).ok_or("outside the table")?;
            largest = largest.max((ours - below / total).abs());
        }
        let direct = sab.group_integrals(&groups, e);
        let total: f64 = direct.iter().sum();
        let mut group = 0.0_f64;
        for (g, d) in direct.iter().enumerate() {
            let (lo, hi) = groups.range(g);
            if lo < low || hi - e >= reach || *d < 1e-4 * total {
                continue;
            }
            let ours = table.cumulative_at(e, hi - e).ok_or("outside the table")?
                - table.cumulative_at(e, lo - e).ok_or("outside the table")?;
            group = group.max((ours / (d / total) - 1.0).abs());
        }
        println!("  {e:.4} eV: share {largest:.3e}, group {group:.3e}");
        worst_share = worst_share.max(largest);
        worst_group = worst_group.max(group);
    }
    println!(
        "{temperature} K, du {du}, dt {dt} eV: {} x {} points in {seconds:.1} s; max share {worst_share:.3e}, max group {worst_group:.3e}",
        table.incident.len(),
        table.transfer.len()
    );
    Ok(())
}
