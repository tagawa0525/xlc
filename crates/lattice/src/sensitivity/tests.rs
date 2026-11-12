use super::*;

/// 1 群の、全断面積と自群の散乱 `scatter` だけを持つ物質（数密度 1 の核種の群定数の代わり）
fn one_group(total: f64, scatter: f64) -> Material {
    Material {
        total: vec![total],
        scatter: vec![vec![scatter]],
        anisotropic: Vec::new(),
        fission: vec![0.0],
        nu_fission: vec![0.0],
        fission_terms: Vec::new(),
    }
}

#[test]
fn sum_check_allows_rounding_of_cancelling_terms() {
    // 散乱の要素で核種の寄与が打ち消し合うと、和は項より桁が小さく、足す順の丸めの差が和に比べて大きい。許容は項の大きさで決める
    let list = vec![
        ("A".to_string(), 1.0, one_group(1.0, 0.1)),
        ("B".to_string(), 1.0, one_group(1.0, 0.2)),
        ("C".to_string(), 1.0, one_group(1.0, -0.3)),
    ];
    // 解いた物質は別の順に足したもの: 0.1 + (0.2 − 0.3) と (0.1 + 0.2) − 0.3 は相対で 1 ほど違う
    let solved = one_group(3.0, 0.1 + (0.2 - 0.3));
    assert_ne!(solved.scatter[0][0], (0.1 + 0.2) - 0.3);
    check_sum("fuel", &solved, &list).unwrap();
    // 項の大きさを超える食い違いは止める
    let wrong = one_group(3.0, 1e-6);
    assert!(check_sum("fuel", &wrong, &list).is_err());
}

#[test]
fn derivatives_tell_listed_nuclides_from_absent_ones() {
    // 数密度 0 で組成にある核種は、変化がなく感度 0（向きが定まらない）。どの物質の組成にもない核種は、名前の誤りの候補なので区別する
    let derivatives = Derivatives {
        materials: vec![
            Vec::new(),
            vec![
                ("U-235".to_string(), 7e-4, one_group(1.0, 0.1)),
                ("Xe-135".to_string(), 0.0, one_group(1.0, 0.1)),
            ],
        ],
    };
    assert_eq!(derivatives.changes("U-235").len(), 1);
    assert!(derivatives.changes("Xe-135").is_empty());
    assert!(derivatives.lists("Xe-135"));
    assert!(!derivatives.lists("Xe-136"));
}
