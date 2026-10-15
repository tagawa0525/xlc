use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

use super::*;

/// 同時に走る `f` の数の最大を数える
struct Concurrency {
    now: AtomicUsize,
    max: AtomicUsize,
}

impl Concurrency {
    fn new() -> Self {
        Self {
            now: AtomicUsize::new(0),
            max: AtomicUsize::new(0),
        }
    }

    fn run(&self) {
        let now = self.now.fetch_add(1, Ordering::SeqCst) + 1;
        self.max.fetch_max(now, Ordering::SeqCst);
        std::thread::sleep(Duration::from_millis(2));
        self.now.fetch_sub(1, Ordering::SeqCst);
    }
}

#[test]
fn results_follow_the_index_order() {
    let values = in_parallel(100, |i| Ok::<_, ()>(i * i)).unwrap();
    assert_eq!(values, (0..100).map(|i| i * i).collect::<Vec<_>>());
}

#[test]
fn nested_calls_share_the_cores() {
    // 外の並列の中で呼んでも、内の計算が同時に走る数はコアの数を超えない（核種ごとの並列の中の格子点の並列）
    let cores = cores();
    let inner = Concurrency::new();
    in_parallel(cores, |_| {
        in_parallel(8, |_| {
            inner.run();
            Ok::<_, ()>(())
        })
    })
    .unwrap();
    let max = inner.max.load(Ordering::SeqCst);
    assert!(max <= cores, "{max} inner tasks at once on {cores} cores");
}
