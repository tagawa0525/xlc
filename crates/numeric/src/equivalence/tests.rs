use super::*;
use crate::nem::{Edge, OneNode, OneNodeError, surface_fluxes};

const WIDTH: f64 = 21.5;

/// 群が結合しない2群の、源のない 1D の解 φ = a cosh κx + b sinh κx（x ∈ [0, L]）の
/// （ノード平均、−x の面の外向きの正味の中性子流、+x の面の外向きの正味の中性子流、−x の面の中性子束、+x の面の中性子束）
fn source_free(d: f64, removal: f64, b: f64) -> (f64, f64, f64, f64, f64) {
    let kappa = (removal / d).sqrt();
    let kl = kappa * WIDTH;
    let (ch, sh) = (kl.cosh(), kl.sinh());
    let a = 1.0;
    let left = -d * kappa * b;
    let right = -d * kappa * (a * sh + b * ch);
    let average = (a * sh + b * (ch - 1.0)) / kl;
    (average, -left, right, a, a * ch + b * sh)
}

fn node<'a>(
    diffusion: &'a [f64],
    total: &'a [f64],
    scatter: &'a [Vec<f64>],
    flux: &'a [f64],
    currents: [&'a [f64]; 4],
) -> OneNode<'a> {
    OneNode {
        diffusion,
        total,
        scatter,
        nu_fission: &[0.0, 0.0],
        chi: &[0.0, 0.0],
        k: 1.0,
        flux,
        widths: [WIDTH, WIDTH],
        currents,
        edges: [
            Edge::Reflective,
            Edge::Vacuum,
            Edge::Reflective,
            Edge::Reflective,
        ],
    }
}

#[test]
fn source_free_groups_take_the_closed_form() {
    // 源のない群では、両方の面の正味の中性子流 J と面の中性子束の比 t = φ+/φ− を、同じ不連続因子で再現する
    // D は閉じた形 D = Σr (L / acosh c)²、c = (1 + r t)/(r + t)、r = J+/J−（向きをそろえた正味の中性子流の比）
    let truth = [1.3, 0.35];
    let total = [0.6, 1.7];
    let scatter = vec![vec![0.57, 0.0], vec![0.0, 1.66]];
    let removal = [0.03, 0.04];
    let solutions = [
        source_free(truth[0], removal[0], -0.95),
        source_free(truth[1], removal[1], -0.97),
    ];
    let flux: Vec<f64> = solutions.iter().map(|s| s.0).collect();
    let minus: Vec<f64> = solutions.iter().map(|s| s.1).collect();
    let plus: Vec<f64> = solutions.iter().map(|s| s.2).collect();
    let zero = [0.0, 0.0];
    // 不均質な面の中性子束は、均質な面の中性子束に群ごとの同じ因子を掛けたもの
    let factors = [1.3, 0.7];
    let het_minus: Vec<f64> = (0..2).map(|g| factors[g] * solutions[g].3).collect();
    let het_plus: Vec<f64> = (0..2).map(|g| factors[g] * solutions[g].4).collect();
    let guess = [0.9, 0.25];
    let found = equivalent_diffusion(
        &node(
            &guess,
            &total,
            &scatter,
            &flux,
            [&minus, &plus, &zero, &zero],
        ),
        [&het_minus, &het_plus, &flux, &flux],
    )
    .unwrap();
    let Equivalence::Found(found) = found else {
        panic!("{found:?}");
    };
    for g in 0..2 {
        let r = plus[g] / -minus[g];
        let t = het_plus[g] / het_minus[g];
        let c = (1.0 + r * t) / (r + t);
        let closed = removal[g] * (WIDTH / c.acosh()).powi(2);
        assert!((closed - truth[g]).abs() < 1e-12 * truth[g], "{closed}");
        assert!((found[g] - truth[g]).abs() < 1e-10 * truth[g], "{found:?}");
    }
}

#[test]
fn coupled_groups_recover_the_generating_d() {
    // 下方散乱と上方散乱のある2群。拡散係数 `truth` の均質な解（ANM）の面の中性子束に群ごとの因子を掛けて
    // 不均質な面の中性子束とすると、等価な拡散係数は `truth` に戻る
    let truth = [1.45, 0.35];
    let total = [0.6215, 1.6711];
    let scatter = vec![vec![0.5956, 0.02371], vec![0.000806, 1.6323]];
    let width = WIDTH;
    // −x から入り +x へ出る（外向きの正味の中性子流）
    let minus = [-0.01136, -0.0005];
    let plus = [0.000927, 0.0002];
    let zero = [0.0, 0.0];
    // 1D のバランス (J+ + J−)/h + Σr φ̄_g − Σs(g'→g) φ̄_g' = 0 からノード平均
    let removal = [total[0] - scatter[0][0], total[1] - scatter[1][1]];
    let leak = [(plus[0] + minus[0]) / width, (plus[1] + minus[1]) / width];
    let det = removal[0] * removal[1] - scatter[1][0] * scatter[0][1];
    let flux = [
        (-leak[0] * removal[1] - scatter[1][0] * leak[1]) / det,
        (-leak[1] * removal[0] - scatter[0][1] * leak[0]) / det,
    ];
    let homogeneous = surface_fluxes(&node(
        &truth,
        &total,
        &scatter,
        &flux,
        [&minus, &plus, &zero, &zero],
    ))
    .unwrap();
    let factors = [0.92, 1.4];
    let het: Vec<Vec<f64>> = homogeneous
        .iter()
        .map(|face| (0..2).map(|g| factors[g] * face[g]).collect())
        .collect();
    let guess = [1.0, 0.3];
    let found = equivalent_diffusion(
        &node(
            &guess,
            &total,
            &scatter,
            &flux,
            [&minus, &plus, &zero, &zero],
        ),
        [&het[0], &het[1], &het[2], &het[3]],
    )
    .unwrap();
    let Equivalence::Found(found) = found else {
        panic!("{found:?}");
    };
    for g in 0..2 {
        assert!((found[g] - truth[g]).abs() < 1e-9 * truth[g], "{found:?}");
    }
}

#[test]
fn impossible_equivalences_are_missing_or_errors() {
    let total = [0.6, 1.7];
    let scatter = vec![vec![0.57, 0.0], vec![0.0, 1.66]];
    let flux = [0.02, 0.01];
    let zero = [0.0, 0.0];
    let guess = [1.0, 0.3];
    // 流れが −x から +x へ減るのに、面の中性子束が +x で大きい（t ≥ 1）: どの D でも同じ因子にならない
    let minus = [-0.01, -0.002];
    let plus = [0.001, 0.0005];
    let result = equivalent_diffusion(
        &node(
            &guess,
            &total,
            &scatter,
            &flux,
            [&minus, &plus, &zero, &zero],
        ),
        [&[0.01, 0.005], &[0.012, 0.004], &flux, &flux],
    );
    assert!(
        matches!(result, Ok(Equivalence::Missing { group: 0, .. })),
        "{result:?}"
    );
    // 正味の中性子流がない（対称な位置）: 向きが定まらない
    let result = equivalent_diffusion(
        &node(&guess, &total, &scatter, &flux, [&zero; 4]),
        [&flux, &flux, &flux, &flux],
    );
    assert!(
        matches!(result, Err(OneNodeError::NoNetCurrent)),
        "{result:?}"
    );
    // 核分裂のあるノードは対象でない
    let mut fissile = node(
        &guess,
        &total,
        &scatter,
        &flux,
        [&minus, &plus, &zero, &zero],
    );
    fissile.nu_fission = &[0.0, 0.1];
    fissile.chi = &[1.0, 0.0];
    let result = equivalent_diffusion(&fissile, [&flux, &flux, &flux, &flux]);
    assert!(
        matches!(result, Err(OneNodeError::Invalid(_))),
        "{result:?}"
    );
}

/// VERA 3A の軸方向の反射体の平板（ADR-0128）の格子計算のノードの1ノードの問題の入力（2群、x が炉心の z、y は
/// 鏡面反射で流れなし。examples/reflector の平板の失敗した実行のログの値、2026-10-11）。`diffusion` は中性子流で
/// 重みづけた D
struct Slab {
    diffusion: [f64; 2],
    total: [f64; 2],
    scatter: [[f64; 2]; 2],
    width: f64,
    flux: [f64; 2],
    heterogeneous: [[f64; 2]; 2],
    currents: [[f64; 2]; 2],
}

impl Slab {
    fn solve(
        &self,
        diffusion: &[f64],
        edges: [Edge<'_>; 2],
    ) -> Result<[Vec<f64>; 4], OneNodeError> {
        let scatter: Vec<Vec<f64>> = self.scatter.iter().map(|r| r.to_vec()).collect();
        let zero = [0.0, 0.0];
        surface_fluxes(&self.node(diffusion, &scatter, &zero, edges))
    }

    fn node<'a>(
        &'a self,
        diffusion: &'a [f64],
        scatter: &'a [Vec<f64>],
        zero: &'a [f64],
        edges: [Edge<'a>; 2],
    ) -> OneNode<'a> {
        OneNode {
            diffusion,
            total: &self.total,
            scatter,
            nu_fission: &[0.0, 0.0],
            chi: &[0.0, 0.0],
            k: 0.6342792405720961,
            flux: &self.flux,
            widths: [self.width, 1.0],
            currents: [&self.currents[0], &self.currents[1], zero, zero],
            edges: [edges[0], edges[1], Edge::Reflective, Edge::Reflective],
        }
    }

    fn choose(&self, edges: [Edge<'_>; 2]) -> Result<(Vec<f64>, Choice), OneNodeError> {
        let scatter: Vec<Vec<f64>> = self.scatter.iter().map(|r| r.to_vec()).collect();
        let zero = [0.0, 0.0];
        let node = self.node(&self.diffusion, &scatter, &zero, edges);
        let het = &self.heterogeneous;
        reflector_diffusion(&node, [&het[0], &het[1], &self.flux, &self.flux])
    }
}

#[test]
fn slab_gap_node_has_no_koebke_d_and_takes_the_current_weighted_d() {
    // 下のギャップと端栓（5.898 cm）: 熱群の源（減速）は燃料側で大きく、均質な面の比 φ+/φ− は正の面の範囲で
    // 1 より大きいが、不均質な比は t = 0.933。Koebke の D は正の面の範囲になく（ADR-0134）、中性子流で重みづけた D で
    // 両面の均質な中性子束が正
    let gap = Slab {
        diffusion: [1.6820777555877764, 0.763528784067662],
        total: [0.6440230468091337, 1.7422110599991552],
        scatter: [
            [0.6110549967237958, 0.031778774280522885],
            [0.0008769671192441397, 1.710908156843627],
        ],
        width: 5.898,
        flux: [0.011638191664828096, 0.006595898791517603],
        heterogeneous: [
            [0.0073158322245259885, 0.0054902840607720296],
            [0.0182116903170929, 0.005122180891888055],
        ],
        currents: [
            [0.001712212355002498, 0.00046632128751277076],
            [-0.0039410907280092, 0.0004972744001251937],
        ],
    };
    let zero = [0.0, 0.0];
    let edges = [
        Edge::Neighbor {
            width: 6.053,
            leakage: &zero,
        },
        Edge::Neighbor {
            width: 22.86,
            leakage: &zero,
        },
    ];
    let scatter: Vec<Vec<f64>> = gap.scatter.iter().map(|r| r.to_vec()).collect();
    let het = &gap.heterogeneous;
    let found = equivalent_diffusion(
        &gap.node(&gap.diffusion, &scatter, &zero, edges),
        [&het[0], &het[1], &gap.flux, &gap.flux],
    )
    .unwrap();
    let Equivalence::Missing {
        group,
        target,
        range,
    } = found
    else {
        panic!("{found:?}");
    };
    assert_eq!(group, 1);
    assert!((target - 0.005122180891888055 / 0.0054902840607720296).abs() < 1e-12);
    let [low, high] = range.unwrap();
    // 正の面の点の均質な比は、どれも t より大きい（D を大きくすると上から 1 に近づく）
    assert!(target < low && low <= high, "{low} {high}");
    let (diffusion, choice) = gap.choose(edges).unwrap();
    assert_eq!(diffusion, gap.diffusion.to_vec());
    assert!(
        matches!(choice, Choice::CurrentWeighted { group: 1, .. }),
        "{choice:?}"
    );
    let faces = gap.solve(&diffusion, edges).unwrap();
    assert!(faces.iter().flatten().all(|&f| f > 0.0), "{faces:?}");
}

#[test]
fn slab_water_node_at_vacuum_takes_the_koebke_d() {
    // 下の水（50 cm、−z が真空）: 中性子流で重みづけた D では真空の面の高速群の均質な中性子束が負になるが、Koebke の
    // D は正の面の範囲にあり、それを選ぶ
    let water = Slab {
        diffusion: [1.715312434162542, 0.28771112638663643],
        total: [0.7026021580140092, 1.948359353059257],
        scatter: [
            [0.6554684859649522, 0.045918789357451015],
            [0.0006212889915745182, 1.913432700752992],
        ],
        width: 50.0,
        flux: [0.00015580362145392823, 0.0001903211822100758],
        heterogeneous: [
            [1.247701469451303e-06, 4.3123280070015635e-07],
            [0.0012443350200998552, 0.0008074493486391072],
        ],
        currents: [
            [8.098077131024756e-07, 2.2708558646791588e-07],
            [-0.0003620774240241963, 2.512451034251237e-05],
        ],
    };
    let zero = [0.0, 0.0];
    let edges = [
        Edge::Vacuum,
        Edge::Neighbor {
            width: 5.0,
            leakage: &zero,
        },
    ];
    assert!(water.solve(&water.diffusion, edges).is_err());
    let (diffusion, choice) = water.choose(edges).unwrap();
    assert!(matches!(choice, Choice::Koebke), "{choice:?}");
    let faces = water.solve(&diffusion, edges).unwrap();
    assert!(faces.iter().flatten().all(|&f| f > 0.0), "{faces:?}");
    // 2つの面の ADF は等しい
    for (g, het) in water.heterogeneous[0].iter().enumerate() {
        let (a, b) = (het / faces[0][g], water.heterogeneous[1][g] / faces[1][g]);
        assert!((a - b).abs() < 1e-9 * a, "{a} {b}");
    }
}

#[test]
fn a_negative_face_root_is_missing_only_when_no_positive_ratio_reaches_the_target() {
    // 根で面の均質な中性子束が負のとき、正の面の点の均質な比の範囲に不均質な比 t が入っていなければ、正の面の範囲に
    // 根はない（Missing。ADR-0134）。入っていれば、正の面の根があるのに別の根に収束したので、黙って中性子流の D に
    // しないようエラーにする
    let missing = negative_root(1, 0.93, Some([1.0002, 11.98])).unwrap();
    assert!(
        matches!(missing, Equivalence::Missing { group: 1, .. }),
        "{missing:?}"
    );
    assert!(negative_root(1, 0.93, None).is_ok());
    assert!(matches!(
        negative_root(1, 2.0, Some([1.0002, 11.98])),
        Err(OneNodeError::Invalid(_))
    ));
}
