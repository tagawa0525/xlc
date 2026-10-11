//! 格子の x を炉心の z とした軸方向の定数（軸方向の反射体の平板。ADR-0128）。

use xlc_format::fewgroup::FewGroup;

use crate::fewgroup::Node;

/// 平板の格子計算の定数 `few`（辺は（左、右、上、下、−z、+z））を、格子の x を炉心の z とした軸方向の定数にする。
///
/// 格子の（左、右）の辺の ADF、真空の比、f* の印を（−z、+z）に移し、径方向の4辺は ADF を 1、真空の比と f* の印を
/// なしにする。横（格子の y と炉心の径方向）に一様な平板では、径方向の面の不均質な中性子束も均質な解も横に一様で
/// ノード平均に等しく、GET の ADF は厳密に 1 なので、MOC の辺の値の離散化のゆらぎを持ち込まない。z の厚さは
/// 格子のノードの x の幅 `width`
pub fn axial(few: &FewGroup, width: f64) -> Result<FewGroup, String> {
    if few.pins.is_some() {
        return Err("axial constants have no pins (form functions are radial)".into());
    }
    few.validate().map_err(|e| e.to_string())?;
    let ng = few.groups();
    let mut out = few.clone();
    out.adf = vec![vec![1.0; ng]; 4];
    out.adf.extend([few.adf[0].clone(), few.adf[1].clone()]);
    out.vacuum_ratio = vec![None; 4];
    out.vacuum_ratio
        .extend([few.vacuum_ratio[0].clone(), few.vacuum_ratio[1].clone()]);
    out.normalized_adf = vec![false; 4];
    out.normalized_adf
        .extend([few.normalized_adf[0], few.normalized_adf[1]]);
    out.height = Some(width);
    out.validate().map_err(|e| e.to_string())?;
    Ok(out)
}

/// 平板の各ノードが格子の `rows` 行のすべてにわたることを確かめる（行の向きは炉心の径方向で、ノードは横に一様）
pub fn check_nodes(rows: usize, nodes: &[Node]) -> Result<(), String> {
    match nodes.iter().position(|node| node.rows != (0..rows)) {
        Some(n) => Err(format!(
            "axial node {n} must span every row of the lattice (rows 1 to {rows})"
        )),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests;
