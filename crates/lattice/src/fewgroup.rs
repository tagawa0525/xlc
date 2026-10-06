//! 均質化した少数群定数（`lattice` から `core` に渡すファイル）。
//!
//! ```text
//! xlc-fewgroup 1
//! groups <G>
//! bounds <G+1 個の群境界 (eV)、降順>
//! k-inf <k∞>
//! total <G 個>
//! absorption <G 個>
//! nu-fission <G 個>
//! fission <G 個>
//! chi <G 個>
//! diffusion <G 個>
//! scatter
//! <G 行 × G 個: Σ(g→g')、行が入射の群>
//! ```

use xs::multigroup::GroupStructure;

use crate::material::Material;
use crate::moc::Solution;

/// セル全体で均質化した少数群定数。断面積は 1/cm、拡散係数は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroup {
    /// 群境界 (eV)、降順
    pub bounds: Vec<f64>,
    pub k_inf: f64,
    pub total: Vec<f64>,
    pub absorption: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    pub chi: Vec<f64>,
    pub diffusion: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
}

impl FewGroup {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    pub fn to_text(&self) -> String {
        let mut out = format!("{FORMAT} {VERSION}\n");
        out += &format!("groups {}\n", self.groups());
        out += &numbers_line("bounds", &self.bounds);
        out += &format!("k-inf {:e}\n", self.k_inf);
        out += &numbers_line("total", &self.total);
        out += &numbers_line("absorption", &self.absorption);
        out += &numbers_line("nu-fission", &self.nu_fission);
        out += &numbers_line("fission", &self.fission);
        out += &numbers_line("chi", &self.chi);
        out += &numbers_line("diffusion", &self.diffusion);
        out += "scatter\n";
        for row in &self.scatter {
            let values: Vec<String> = row.iter().map(|v| format!("{v:e}")).collect();
            out += &values.join(" ");
            out += "\n";
        }
        out
    }

    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut lines = Lines::new(text);
        let (line, header) = lines.next("format header")?;
        if header != [FORMAT, VERSION] {
            return Err(format!("line {line}: expected `{FORMAT} {VERSION}`"));
        }
        let (line, groups) = lines.keyword("groups", 1)?;
        let ng = groups[0]
            .parse::<usize>()
            .ok()
            .filter(|&n| n > 0 && n < usize::MAX)
            .ok_or_else(|| format!("line {line}: invalid number of groups `{}`", groups[0]))?;
        let bounds = lines.numbers_after("bounds", ng + 1)?;
        let k_inf = lines.numbers_after("k-inf", 1)?[0];
        let total = lines.numbers_after("total", ng)?;
        let absorption = lines.numbers_after("absorption", ng)?;
        let nu_fission = lines.numbers_after("nu-fission", ng)?;
        let fission = lines.numbers_after("fission", ng)?;
        let chi = lines.numbers_after("chi", ng)?;
        let diffusion = lines.numbers_after("diffusion", ng)?;
        lines.keyword("scatter", 0)?;
        let scatter = (0..ng)
            .map(|_| {
                let (line, tokens) = lines.next("a row of scattering cross sections")?;
                numbers(line, &tokens, ng)
            })
            .collect::<Result<_, _>>()?;
        if let Ok((line, _)) = lines.next("") {
            return Err(format!(
                "line {line}: unexpected data after the scattering matrix"
            ));
        }
        Ok(Self {
            bounds,
            k_inf,
            total,
            absorption,
            nu_fission,
            fission,
            chi,
            diffusion,
            scatter,
        })
    }
}

const FORMAT: &str = "xlc-fewgroup";
const VERSION: &str = "1";

fn numbers_line(keyword: &str, values: &[f64]) -> String {
    let mut out = keyword.to_string();
    for v in values {
        out += &format!(" {v:e}");
    }
    out + "\n"
}

fn numbers(line: usize, tokens: &[&str], n: usize) -> Result<Vec<f64>, String> {
    if tokens.len() != n {
        return Err(format!(
            "line {line}: expected {n} values, found {}",
            tokens.len()
        ));
    }
    tokens
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| format!("line {line}: invalid number `{t}`"))
        })
        .collect()
}

/// 空行を飛ばして、行番号付きでトークンに分けた行を順に返す。
struct Lines<'a> {
    lines: std::vec::IntoIter<(usize, Vec<&'a str>)>,
    end: usize,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        let lines: Vec<_> = text
            .lines()
            .enumerate()
            .map(|(i, l)| (i + 1, l.split_whitespace().collect::<Vec<_>>()))
            .filter(|(_, tokens)| !tokens.is_empty())
            .collect();
        let end = text.lines().count() + 1;
        Self {
            lines: lines.into_iter(),
            end,
        }
    }

    fn next(&mut self, expected: &str) -> Result<(usize, Vec<&'a str>), String> {
        self.lines.next().ok_or_else(|| {
            format!(
                "line {}: unexpected end of file, expected {expected}",
                self.end
            )
        })
    }

    /// `keyword` で始まり、続く語がちょうど n 個の行
    fn keyword(&mut self, keyword: &str, n: usize) -> Result<(usize, Vec<&'a str>), String> {
        let (line, tokens) = self.next(keyword)?;
        if tokens[0] != keyword || tokens.len() != n + 1 {
            return Err(format!("line {line}: expected `{keyword}` with {n} values"));
        }
        Ok((line, tokens[1..].to_vec()))
    }

    fn numbers_after(&mut self, keyword: &str, n: usize) -> Result<Vec<f64>, String> {
        let (line, tokens) = self.next(keyword)?;
        if tokens[0] != keyword {
            return Err(format!("line {line}: expected `{keyword}`"));
        }
        numbers(line, &tokens[1..], n)
    }
}

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は D = 1/(3Σt) とする（M0 の近似）。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<FewGroup, String> {
    let ng = fine.len();
    if materials.len() != solution.flux.len() {
        return Err(format!(
            "{} materials for {} regions",
            materials.len(),
            solution.flux.len()
        ));
    }
    // 少数群の境界を細かい群の境界の位置に対応させる
    let positions = coarse_bounds
        .iter()
        .map(|b| fine.bounds().iter().position(|f| f == b))
        .collect::<Option<Vec<usize>>>()
        .ok_or("few-group bounds must be bounds of the library groups")?;
    let covers = positions.first() == Some(&0) && positions.last() == Some(&ng);
    if !covers || positions.windows(2).any(|w| w[0] >= w[1]) {
        return Err("few-group bounds must be descending and span the library groups".into());
    }
    let coarse = positions.len() - 1;
    let mut group_of = vec![0; ng];
    for (big, w) in positions.windows(2).enumerate() {
        group_of[w[0]..w[1]].fill(big);
    }

    let mut flux = vec![0.0; coarse];
    let mut total = vec![0.0; coarse];
    let mut nu_fission = vec![0.0; coarse];
    let mut fission = vec![0.0; coarse];
    let mut scatter = vec![vec![0.0; coarse]; coarse];
    let mut chi = vec![0.0; coarse];
    for ((m, phi), volume) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
        for g in 0..ng {
            let w = volume * phi[g];
            let big = group_of[g];
            flux[big] += w;
            total[big] += w * m.total[g];
            nu_fission[big] += w * m.nu_fission[g];
            fission[big] += w * m.fission[g];
            for (to, s) in m.scatter[g].iter().enumerate() {
                scatter[big][group_of[to]] += w * s;
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (g, c) in term.chi.iter().enumerate() {
                chi[group_of[g]] += c * rate;
            }
        }
    }
    if flux.iter().any(|&f| f <= 0.0) {
        return Err("a few group has no flux".into());
    }
    let per_flux =
        |values: Vec<f64>| -> Vec<f64> { values.iter().zip(&flux).map(|(v, f)| v / f).collect() };
    let total = per_flux(total);
    let scatter: Vec<Vec<f64>> = scatter
        .into_iter()
        .zip(&flux)
        .map(|(row, f)| row.iter().map(|v| v / f).collect())
        .collect();
    let produced: f64 = chi.iter().sum();
    Ok(FewGroup {
        bounds: coarse_bounds.to_vec(),
        k_inf: solution.k,
        absorption: total
            .iter()
            .zip(&scatter)
            .map(|(t, row)| t - row.iter().sum::<f64>())
            .collect(),
        diffusion: total.iter().map(|t| 1.0 / (3.0 * t)).collect(),
        total,
        nu_fission: per_flux(nu_fission),
        fission: per_flux(fission),
        chi: if produced > 0.0 {
            chi.iter().map(|c| c / produced).collect()
        } else {
            chi
        },
        scatter,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    fn solved() -> (Solution, Vec<Material>, GroupStructure) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let moderator = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let materials = vec![fuel.clone(), moderator.clone(), moderator];
        let cell = PinCell::new(1.26, vec![0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 16, 0.05).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (
            solution,
            materials,
            GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]),
        )
    }

    #[test]
    fn one_group_collapse_preserves_k() {
        // 反応率を保存するので、1群の νΣf/Σa はセルの k∞ に一致する
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let one = homogenize(&solution, &refs, &fine, &[1.0e7, 1.0e-5]).unwrap();
        assert_eq!(one.groups(), 1);
        let k = one.nu_fission[0] / one.absorption[0];
        assert!(
            (k - solution.k).abs() < 1e-6 * solution.k,
            "{k} vs {}",
            solution.k
        );
        assert_eq!(one.k_inf, solution.k);
        assert!((one.chi[0] - 1.0).abs() < 1e-12);
        assert!((one.diffusion[0] - 1.0 / (3.0 * one.total[0])).abs() < 1e-15);
    }

    #[test]
    fn two_group_constants_reproduce_k_inf() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds()).unwrap();
        // 均質な2群の無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
        let a = [
            [fg.total[0] - fg.scatter[0][0], -fg.scatter[1][0]],
            [-fg.scatter[0][1], fg.total[1] - fg.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        let k = fg.nu_fission[0] * phi[0] + fg.nu_fission[1] * phi[1];
        assert!((k - solution.k).abs() < 1e-6 * solution.k);
        for g in 0..2 {
            let removal = fg.total[g] - fg.scatter[g].iter().sum::<f64>();
            assert!((fg.absorption[g] - removal).abs() < 1e-14);
        }
    }

    #[test]
    fn coarse_bounds_must_be_fine_bounds() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 1.0, 1.0e-5]).is_err());
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 0.625]).is_err());
    }

    #[test]
    fn text_round_trips() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds()).unwrap();
        let text = fg.to_text();
        assert!(text.starts_with("xlc-fewgroup 1\n"));
        assert_eq!(FewGroup::parse(&text).unwrap(), fg);
    }

    #[test]
    fn parse_errors_have_line_numbers() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let text = homogenize(&solution, &refs, &fine, fine.bounds())
            .unwrap()
            .to_text();
        let broken: String = text
            .lines()
            .map(|l| {
                if l.starts_with("diffusion") {
                    "diffusion 1.0"
                } else {
                    l
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let err = FewGroup::parse(&broken).unwrap_err();
        let line = text
            .lines()
            .position(|l| l.starts_with("diffusion"))
            .unwrap()
            + 1;
        assert!(err.starts_with(&format!("line {line}:")), "{err}");
        assert!(FewGroup::parse("xlc-fewgroup 9\n").is_err());
    }
}
