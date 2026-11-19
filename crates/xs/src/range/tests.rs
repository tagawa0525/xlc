use super::*;

#[test]
fn the_pointwise_range_file_round_trips() {
    for range in [Some((4.0, 9118.0)), Some((0.1 + 0.2, 1.0e5)), None] {
        let text = to_toml(range);
        assert!(text.contains("pointwise-range"), "{text}");
        assert_eq!(from_toml(&text).unwrap(), range);
    }
}

#[test]
fn rejects_a_bad_pointwise_range_file() {
    for text in [
        "",
        "pointwise-range = [4.0]",
        "pointwise-range = [4.0, 5.0, 6.0]",
        "pointwise-range = [5.0, 4.0]",
        "pointwise-range = [4.0, 4.0]",
        "pointwise-range = [0.0, 4.0]",
        "pointwise-range = [4.0, 1e5]\nextra = 1",
        "pointwise-range = \"4 1e5\"",
    ] {
        assert!(from_toml(text).is_err(), "{text}");
    }
}
