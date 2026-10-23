use std::path::PathBuf;

use xlc_format::medium::MediumState;

use super::*;
use crate::depletion::tests::{PIN, deplete_with, library};

/// テストごとの空のディレクトリ
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("xlc-checkpoint-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// 少数群定数と媒質を持つ状態（少数群定数は燃焼の点のもの）
fn states(count: usize) -> Vec<SavedState> {
    let text = format!("{PIN}few-group = {{ bounds = [1e7, 0.625, 1e-5] }}\n");
    let history = deplete_with(&text, &library());
    (0..count)
        .map(|i| {
            let point = &history.points[i % history.points.len()];
            let few = point.few_group.clone().unwrap();
            SavedState {
                k: point.k,
                outputs: vec![vec![few.clone()], vec![few.clone(), few]],
                medium: Some(MediumState {
                    total: vec![0.5 + i as f64, 1.0 / 3.0],
                    scatter: vec![vec![vec![0.25, 0.1], vec![0.0, 0.3]]],
                    fission: vec![0.01, 0.1],
                    nu_fission: vec![0.025, 0.24],
                    chi: vec![1.0, 0.0],
                }),
            }
        })
        .collect()
}

#[test]
fn fingerprints_tell_apart_where_the_bytes_belong() {
    let print = |parts: &[(&str, &[u8])]| {
        let mut f = Fingerprint::default();
        for (label, bytes) in parts {
            f.add(label, bytes);
        }
        f.finish()
    };
    let one = print(&[("input", b"ab"), ("library", b"c")]);
    assert_eq!(one, print(&[("input", b"ab"), ("library", b"c")]));
    assert_ne!(one, print(&[("input", b"a"), ("library", b"bc")]));
    assert_ne!(one, print(&[("input", b"ab"), ("chain", b"c")]));
    assert_eq!(one.len(), 64);
}

#[test]
fn without_a_checkpoint_the_run_starts_from_the_beginning() {
    let dir = scratch("missing");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    assert!(checkpoint.depletion().unwrap().is_none());
    assert!(checkpoint.states().unwrap().is_empty());
    // 消すものがなくてもよい（一気に流した燃焼は、点の数によっては 1 度も書かない）
    checkpoint.remove().unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn saved_states_read_back_bit_for_bit() {
    let dir = scratch("states");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    let saved = states(3);
    checkpoint.save_states(&saved).unwrap();
    let read = checkpoint.states().unwrap();
    assert_eq!(read, saved);
    for (r, s) in read.iter().zip(&saved) {
        assert_eq!(r.k.to_bits(), s.k.to_bits());
    }
    // 書き込みは一時ファイルを rename で置き換えるので、チェックポイントのほかに残らない
    let names: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    assert_eq!(names, vec![std::ffi::OsString::from("a.checkpoint")]);
    checkpoint.remove().unwrap();
    assert!(!checkpoint.path().exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_checkpoint_of_another_fingerprint_is_an_error() {
    let dir = scratch("fingerprint");
    let path = dir.join("a.checkpoint");
    Checkpoint::new(path.clone(), "old".into())
        .save_states(&states(1))
        .unwrap();
    let error = Checkpoint::new(path.clone(), "new".into())
        .states()
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("a.checkpoint") && error.contains("remove"),
        "{error}"
    );
    let history = deplete_with(PIN, &library());
    let next = history.points[1].densities.clone();
    Checkpoint::new(path.clone(), "old".into())
        .save_depletion(&history, &next)
        .unwrap();
    let error = Checkpoint::new(path, "new".into())
        .depletion()
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("a.checkpoint") && error.contains("remove"),
        "{error}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_checkpoint_of_the_other_kind_is_an_error() {
    let dir = scratch("kind");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint.save_states(&states(1)).unwrap();
    let error = checkpoint.depletion().unwrap_err().to_string();
    assert!(error.contains("states"), "{error}");
    let history = deplete_with(PIN, &library());
    let next = history.points[1].densities.clone();
    checkpoint.save_depletion(&history, &next).unwrap();
    let error = checkpoint.states().unwrap_err().to_string();
    assert!(error.contains("depletion"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn states_resume_after_the_saved_ones_and_match_the_straight_run() {
    let dir = scratch("solve");
    let all = states(4);
    let straight = solve_states(4, None, |_, _| {}, |i| Ok(all[i].clone())).unwrap();
    assert_eq!(straight, all);
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    // 状態 2 を解いている途中で止めた実行
    let error = solve_states(
        4,
        Some(&checkpoint),
        |_, _| {},
        |i| {
            if i == 2 {
                anyhow::bail!("killed");
            }
            Ok(all[i].clone())
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("killed"), "{error:#}");
    assert_eq!(checkpoint.states().unwrap().len(), 2);
    let mut restored = Vec::new();
    let mut solved = Vec::new();
    let resumed = solve_states(
        4,
        Some(&checkpoint),
        |i, s| restored.push((i, s.k)),
        |i| {
            solved.push(i);
            Ok(all[i].clone())
        },
    )
    .unwrap();
    assert_eq!(restored, vec![(0, all[0].k), (1, all[1].k)]);
    assert_eq!(solved, vec![2, 3]);
    assert_eq!(resumed, straight);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn more_saved_states_than_the_run_has_are_an_error() {
    let dir = scratch("count");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint.save_states(&states(3)).unwrap();
    let error = solve_states(
        2,
        Some(&checkpoint),
        |_, _| {},
        |_| Err(anyhow::anyhow!("no state is solved")),
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("3 states") && error.contains('2'), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}
