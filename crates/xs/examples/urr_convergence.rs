//! 非分離共鳴のサンプリングした自己遮蔽（ADR-0177）の統計誤差を測る。
//!
//! 使い方: urr_convergence <ENDF> <温度 (K)> <ラダーの数> <ラダーごとのエネルギーの数> <エネルギー (eV)>...
//!
//! エネルギーごとに、種 1〜16 の独立なサンプリングで、σ0 ごとの遮蔽の比（捕獲と全断面積の、実効断面積 / 同じサンプルの平均）の
//! 平均と、1 回のサンプリング（1 つの種）の相対標準偏差を書く。エネルギーごとに、捕獲と全断面積のサンプルの平均（種で平均）も書く。ライブラリは無限希釈の値にこの比（LSSF=1）か差（LSSF=0）を当てるので、比の誤差が
//! 群定数の統計誤差になる。

use xlc_xs::endf::{Mf2, RangeKind, Tape};
use xlc_xs::resonance::{Sampling, Sigma, UnresolvedXs};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, temperature, ladders, samples, energies @ ..] = args.as_slice() else {
        return Err(
            "usage: urr_convergence <ENDF> <temperature K> <ladders> <samples> <energy eV>..."
                .into(),
        );
    };
    let sampling = Sampling {
        ladders: ladders.parse()?,
        samples: samples.parse()?,
    };
    let temperature: f64 = temperature.parse()?;
    let tape = Tape::parse(&std::fs::read_to_string(endf)?)?;
    let mat = *tape.materials().first().ok_or("no material")?;
    let mf2 = Mf2::parse(&tape.section(mat, 2, 151).ok_or("no MF2")?)?;
    let urr = mf2
        .ranges
        .iter()
        .find_map(|r| match &r.kind {
            RangeKind::Unresolved(u) => Some(u),
            _ => None,
        })
        .ok_or("no unresolved range")?;
    let xs = UnresolvedXs::new(urr)?;
    let sigma0 = [1e10, 1e3, 100.0, 10.0, 1.0];
    println!("energy\tsigma0\tcapture_ratio\trel_std\ttotal_ratio\trel_std");
    for e in energies {
        let e: f64 = e.parse()?;
        let mut ratios = vec![Vec::new(); sigma0.len()];
        let mut means = [0.0; 2];
        for seed in 1..=16 {
            let b = xs.shielded(e, temperature, &sigma0, Sigma::default(), sampling, seed)?;
            means[0] += b.mean.capture / 16.0;
            means[1] += b.mean.total / 16.0;
            for (r, s) in ratios.iter_mut().zip(&b.shielded) {
                r.push([s.capture / b.mean.capture, s.total / b.mean.total]);
            }
        }
        println!(
            "# {e:e} eV: mean capture {:.6e} b, mean total {:.6e} b",
            means[0], means[1]
        );
        for (s0, r) in sigma0.iter().zip(&ratios) {
            let n = r.len() as f64;
            let stats = |k: usize| {
                let mean = r.iter().map(|v| v[k]).sum::<f64>() / n;
                let var = r.iter().map(|v| (v[k] - mean).powi(2)).sum::<f64>() / (n - 1.0);
                (mean, var.sqrt() / mean)
            };
            let ((c, ce), (t, te)) = (stats(0), stats(1));
            println!("{e:e}\t{s0:e}\t{c:.6}\t{ce:.2e}\t{t:.6}\t{te:.2e}");
        }
    }
    Ok(())
}
