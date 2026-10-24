//! S(α, β) の表による散乱の移行と、物質の温度の自由ガス（断面積を一定とした自由原子）の移行を、入射の群の中でレサジーについて
//! 平均して、出射の群ごとに比べる（ADR-0203: UFG の減速計算の H の核の誤差。H は質量比がほぼ 1 なので、S(α, β) の上限より下で
//! 短時間衝突近似（有効温度）の P0 の核は、物質の温度の自由ガスの核と倍精度で同じになる）。
//!
//! 使い方: sab_free_gas_rows <熱中性子散乱の ENDF> <温度 K> <入射の下端 eV> <入射の上端 eV> [入射の点の数]
//!
//! 群構造は WIMS 69群。表の外は THERMR と同じく短時間衝突近似で補い（ADR-0133）、上限は評価済みデータの B(4)。入射の範囲を
//! レサジーで等分した点の中点で、出射の群ごとの断面積 ∫_g' σ(E→E') dE'（barn）を求め、点の平均を出す。

use xlc_format::group::GroupStructure;
use xlc_xs::endf::Tape;
use xlc_xs::multigroup::{SabOutside, ThermalSab, free_gas_group_integrals};
use xlc_xs::physics::BK;
use xlc_xs::process::read_thermal_scattering;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let (tsl, temperature, low, high, points) = match args.as_slice() {
        [_, tsl, t, lo, hi] => (
            tsl,
            t.parse::<f64>()?,
            lo.parse::<f64>()?,
            hi.parse::<f64>()?,
            400,
        ),
        [_, tsl, t, lo, hi, n] => (
            tsl,
            t.parse::<f64>()?,
            lo.parse::<f64>()?,
            hi.parse::<f64>()?,
            n.parse::<usize>()?,
        ),
        _ => {
            return Err(
                "usage: sab_free_gas_rows <TSL ENDF> <K> <low eV> <high eV> [points]".into(),
            );
        }
    };
    let data = read_thermal_scattering(&Tape::parse(&std::fs::read_to_string(tsl)?)?)?;
    let sab = ThermalSab::new(
        &data,
        temperature,
        data.e_max,
        SabOutside::ShortCollisionTime,
    )?;
    let groups = GroupStructure::wims69();
    let kt = BK * temperature;
    let (mut table, mut free) = (vec![0.0; groups.len()], vec![0.0; groups.len()]);
    let width = (high / low).ln();
    for k in 0..points {
        let e = high * (-(k as f64 + 0.5) / points as f64 * width).exp();
        for (sum, v) in table.iter_mut().zip(sab.group_integrals(&groups, e)) {
            *sum += v / points as f64;
        }
        let row = free_gas_group_integrals(&groups, e, data.awr, kt);
        for (sum, v) in free.iter_mut().zip(row) {
            *sum += data.sigma_free * v / points as f64;
        }
    }
    let (sum_table, sum_free) = (table.iter().sum::<f64>(), free.iter().sum::<f64>());
    println!(
        "{temperature} K, incident {low}-{high} eV, B(4) {} eV, sigma_free {} b",
        data.e_max, data.sigma_free
    );
    println!(
        "row sum: S(a,b) {sum_table:.6} b, free gas {sum_free:.6} b, ratio - 1 {:+.4e}",
        sum_table / sum_free - 1.0
    );
    println!(
        "group (1-based)  bounds (eV)  S(a,b) (b)  free gas (b)  ratio - 1  share S(a,b)  share free gas"
    );
    for g in 0..groups.len() {
        if table[g].abs() < 1e-12 * sum_table && free[g].abs() < 1e-12 * sum_free {
            continue;
        }
        let (lo, hi) = groups.range(g);
        println!(
            "{:>3}  {lo:.4e}-{hi:.4e}  {:.6e}  {:.6e}  {:+.4e}  {:.6e}  {:.6e}",
            g + 1,
            table[g],
            free[g],
            table[g] / free[g] - 1.0,
            table[g] / sum_table,
            free[g] / sum_free
        );
    }
    Ok(())
}
