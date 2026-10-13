//! 元素記号。

/// 原子番号 1〜118 の元素記号（添字は Z − 1）
const SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// ADR-0015 の形の核種名（`U-235`、`Am-242m`、`Hf-178m2`）。Z が範囲の外なら `None`。
pub fn nuclide_name(z: u32, a: u32, state: u32) -> Option<String> {
    let symbol = SYMBOLS.get(usize::try_from(z).ok()?.checked_sub(1)?)?;
    Some(match state {
        0 => format!("{symbol}-{a}"),
        1 => format!("{symbol}-{a}m"),
        s => format!("{symbol}-{a}m{s}"),
    })
}
