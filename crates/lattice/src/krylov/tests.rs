use super::*;

#[test]
fn gmres_stops_when_the_operator_maps_a_basis_vector_to_zero() {
    // 作用素が基底の値を 0 に移すと、Givens の回転と後退代入が 0 で割り、解が NaN になる
    let mut krylov = Krylov::new(vec![1.0, 0.0], 1e-12, 10);
    assert!(matches!(krylov.receive(vec![0.0, 0.0], 2), Err(Breakdown)));
}
