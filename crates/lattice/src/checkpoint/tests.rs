use std::path::PathBuf;

use xlc_format::medium::MediumState;

use super::*;
use crate::branches::Starts;
use crate::depletion::tests::{PIN, deplete_with, first_save, library};

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

/// 軸 `counts` の状態を `solved` 個解いた後の、反復を始める解（状態 i の解を i + 0.5 とする）
fn starts_after(counts: &[usize], solved: usize) -> Starts<f64> {
    let mut starts = Starts::new(counts);
    for i in 0..solved {
        starts.solved(i as f64 + 0.5);
    }
    starts
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
    let (states, starts) = checkpoint.states::<f64>(&[2, 3]).unwrap();
    assert!(states.is_empty());
    assert_eq!(starts.start(), None);
    // 消すものがなくてもよい（一気に流した燃焼は、点の数によっては 1 度も書かない）
    checkpoint.remove().unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn saved_states_and_starts_read_back_bit_for_bit() {
    let dir = scratch("states");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    let saved = states(3);
    let counts = [2, 3];
    checkpoint
        .save_states(&saved, &starts_after(&counts, 3))
        .unwrap();
    let (read, mut starts) = checkpoint.states::<f64>(&counts).unwrap();
    assert_eq!(read, saved);
    for (r, s) in read.iter().zip(&saved) {
        assert_eq!(r.k.to_bits(), s.k.to_bits());
    }
    // 状態 3（添字 [1, 0]）は状態 0 の解から、状態 4（[1, 1]）は状態 3 の解から始める
    let mut straight = starts_after(&counts, 3);
    assert_eq!(starts.start(), straight.start());
    assert_eq!(starts.start(), Some(&0.5));
    starts.solved(3.5);
    straight.solved(3.5);
    assert_eq!(starts.start(), straight.start());
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
        .save_states(&states(1), &starts_after(&[2], 1))
        .unwrap();
    let error = Checkpoint::new(path.clone(), "new".into())
        .states::<f64>(&[2])
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("a.checkpoint") && error.contains("remove"),
        "{error}"
    );
    let (history, next, start) = first_save(PIN, &library());
    Checkpoint::new(path.clone(), "old".into())
        .save_depletion(&history, &next, &start)
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
    checkpoint
        .save_states(&states(1), &starts_after(&[2], 1))
        .unwrap();
    let error = checkpoint.depletion().unwrap_err().to_string();
    assert!(error.contains("states"), "{error}");
    let (history, next, start) = first_save(PIN, &library());
    checkpoint.save_depletion(&history, &next, &start).unwrap();
    let error = checkpoint.states::<f64>(&[2]).unwrap_err().to_string();
    assert!(error.contains("depletion"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// 状態 `i` を、始める解 `start` から解いたふりをする。結果は始める解によるので、始める解を取り違えると変わる
fn solve_from(all: &[SavedState], i: usize, start: Option<&f64>) -> (SavedState, f64) {
    let value = start.map_or(0.0, |s| s * 0.5) + i as f64;
    let mut state = all[i].clone();
    state.k += value;
    (state, value)
}

#[test]
fn states_resume_after_the_saved_ones_and_match_the_straight_run() {
    let dir = scratch("solve");
    let all = states(6);
    let counts = [2, 3];
    let straight = solve_states(
        6,
        &counts,
        None,
        |_, _| Ok(()),
        |i, start: Option<&f64>| Ok(solve_from(&all, i, start)),
    )
    .unwrap();
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    // 状態 4 を解いている途中で止めた実行
    let error = solve_states(
        6,
        &counts,
        Some(&checkpoint),
        |_, _| Ok(()),
        |i, start: Option<&f64>| {
            if i == 4 {
                anyhow::bail!("killed");
            }
            Ok(solve_from(&all, i, start))
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("killed"), "{error:#}");
    assert_eq!(checkpoint.states::<f64>(&counts).unwrap().0.len(), 4);
    let mut restored = Vec::new();
    let mut solved = Vec::new();
    let resumed = solve_states(
        6,
        &counts,
        Some(&checkpoint),
        |i, s| {
            restored.push((i, s.k));
            Ok(())
        },
        |i, start: Option<&f64>| {
            solved.push(i);
            Ok(solve_from(&all, i, start))
        },
    )
    .unwrap();
    assert_eq!(
        restored,
        (0..4).map(|i| (i, straight[i].k)).collect::<Vec<_>>()
    );
    assert_eq!(solved, vec![4, 5]);
    assert_eq!(resumed, straight);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn more_saved_states_than_the_run_has_are_an_error() {
    let dir = scratch("count");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint
        .save_states(&states(3), &starts_after(&[2], 3))
        .unwrap();
    let error = solve_states(
        2,
        &[2],
        Some(&checkpoint),
        |_, _| Ok(()),
        |_, _: Option<&f64>| Err(anyhow::anyhow!("no state is solved")),
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("3 states") && error.contains('2'), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn saved_starts_must_follow_the_saved_states() {
    let dir = scratch("follow");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint
        .save_states(&states(2), &starts_after(&[4], 3))
        .unwrap();
    let error = format!("{:#}", checkpoint.states::<f64>(&[4]).unwrap_err());
    assert!(error.contains("2 states") && error.contains('3'), "{error}");
    // 軸の格子点の数が違う解も使わない
    checkpoint
        .save_states(&states(2), &starts_after(&[4], 2))
        .unwrap();
    let error = format!("{:#}", checkpoint.states::<f64>(&[2, 2]).unwrap_err());
    assert!(error.contains("axes"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn an_unreadable_checkpoint_names_its_file() {
    let dir = scratch("unreadable");
    let path = dir.join("a.checkpoint");
    std::fs::write(&path, b"").unwrap();
    let error = format!(
        "{:#}",
        Checkpoint::new(path, "f".into())
            .states::<f64>(&[2])
            .unwrap_err()
    );
    assert!(
        error.contains("a.checkpoint") && error.contains("remove"),
        "{error}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn restored_states_that_do_not_fit_the_run_stop_before_solving() {
    let dir = scratch("fit");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint
        .save_states(&states(2), &starts_after(&[4], 2))
        .unwrap();
    let mut solved = 0;
    let all = states(1);
    let error = solve_states(
        4,
        &[4],
        Some(&checkpoint),
        |i, _| Err(anyhow::anyhow!("state {i} has no medium")),
        |_, _: Option<&f64>| {
            solved += 1;
            Ok((all[0].clone(), 0.0))
        },
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("state 0"), "{error}");
    assert_eq!(solved, 0);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_checkpoint_needs_exactly_one_kind() {
    let dir = scratch("one-kind");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint.write(&checkpoint.file()).unwrap();
    let error = checkpoint.states::<f64>(&[2]).unwrap_err().to_string();
    assert!(error.contains("neither"), "{error}");
    let (history, next, start) = first_save(PIN, &library());
    let both = CheckpointFile {
        depletion: Some(DepletionFile {
            history: rmpv::Value::Binary(history.to_msgpack().unwrap()),
            next,
            start,
        }),
        states: Some(StatesFile {
            states: Vec::new(),
            starts: StartsFile {
                state: 0,
                solutions: Vec::new(),
                slots: Vec::new(),
            },
        }),
        ..checkpoint.file()
    };
    checkpoint.write(&both).unwrap();
    let error = checkpoint.depletion().unwrap_err().to_string();
    assert!(error.contains("both"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_last_state_is_not_saved() {
    // 最後の状態の後は出力を書いてチェックポイントを消すので、書くのはそれより前の状態だけ（燃焼の最後の点と同じ）
    let dir = scratch("last");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    let all = states(3);
    let mut saved_before = Vec::new();
    solve_states(
        3,
        &[3],
        Some(&checkpoint),
        |_, _| Ok(()),
        |i, start: Option<&f64>| {
            saved_before.push(checkpoint.states::<f64>(&[3]).unwrap().0.len());
            Ok(solve_from(&all, i, start))
        },
    )
    .unwrap();
    assert_eq!(saved_before, vec![0, 1, 2]);
    assert_eq!(checkpoint.states::<f64>(&[3]).unwrap().0.len(), 2);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_solution_shared_by_axes_is_written_once() {
    // 状態 0 の解は、すべての軸の次の状態が始める解になる。軸の数だけ書かない
    let dir = scratch("shared");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    let counts = [2, 2, 2];
    checkpoint
        .save_states(&states(1), &starts_after(&counts, 1))
        .unwrap();
    let bytes = std::fs::read(checkpoint.path()).unwrap();
    let file: CheckpointFile = msgpack::from_value(msgpack::read_value(&bytes).unwrap()).unwrap();
    let starts = file.states.unwrap().starts;
    assert_eq!(starts.solutions.len(), 1);
    assert_eq!(starts.slots, vec![Some(0), Some(0), Some(0)]);
    let (_, mut read) = checkpoint.states::<f64>(&counts).unwrap();
    let mut straight = starts_after(&counts, 1);
    for state in 1..8 {
        assert_eq!(read.start(), straight.start(), "state {state}");
        read.solved(state as f64 + 0.5);
        straight.solved(state as f64 + 0.5);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_missing_start_of_the_next_state_is_an_error() {
    // 状態 1 は状態 0 の解から始める。書いた解がなければ、平らな初期値から黙って始めずに止める
    let dir = scratch("missing-start");
    let checkpoint = Checkpoint::new(dir.join("a.checkpoint"), "f".into());
    checkpoint
        .save_states(&states(1), &starts_after(&[2], 1))
        .unwrap();
    let bytes = std::fs::read(checkpoint.path()).unwrap();
    let mut file: CheckpointFile =
        msgpack::from_value(msgpack::read_value(&bytes).unwrap()).unwrap();
    let starts = &mut file.states.as_mut().unwrap().starts;
    starts.slots = vec![None];
    starts.solutions.clear();
    checkpoint.write(&file).unwrap();
    let error = format!("{:#}", checkpoint.states::<f64>(&[2]).unwrap_err());
    assert!(error.contains("no start"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}
