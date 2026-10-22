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
fn nested_calls_stay_within_the_budget() {
    // 外の並列の中で呼んでも、内の計算が同時に走る数は、予算と呼んだスレッドの分を超えない（核種と温度の並列の中の格子点の並列）
    let budget = Budget::new(3);
    let inner = Concurrency::new();
    in_parallel_in(&budget, 4, |_| {
        in_parallel_in(&budget, 8, |_| {
            inner.run();
            Ok::<_, ()>(())
        })
    })
    .unwrap();
    let max = inner.max.load(Ordering::SeqCst);
    assert!(
        max <= 4,
        "{max} inner tasks at once on a budget of 3 and the caller"
    );
}

#[test]
#[ignore = "in_parallel takes its threads only when it starts"]
fn a_running_call_takes_threads_freed_later() {
    // 予算が空のときに始めた呼び出しも、途中で空いたスレッドを使う（ほかの核種や温度が先に終わったとき）
    let budget = Budget::new(3);
    let all = budget.take(3);
    let tasks = Concurrency::new();
    std::thread::scope(|scope| {
        let call = scope.spawn(|| {
            in_parallel_in(&budget, 100, |_| {
                tasks.run();
                Ok::<_, ()>(())
            })
        });
        std::thread::sleep(Duration::from_millis(20));
        drop(all);
        call.join().unwrap().unwrap();
    });
    let max = tasks.max.load(Ordering::SeqCst);
    assert!(max >= 2, "{max} tasks at once after the budget was freed");
}

#[test]
#[ignore = "in_parallel keeps its threads until every index is done"]
fn threads_out_of_indices_return_to_the_budget() {
    // 取る i がなくなったスレッドは予算に返り、まだ走っている重い i の中の並列が使える（最後に残った重い核種）
    let budget = Budget::new(2);
    let inner = Concurrency::new();
    in_parallel_in(&budget, 2, |i| {
        if i == 0 {
            return Ok(());
        }
        std::thread::sleep(Duration::from_millis(20));
        in_parallel_in(&budget, 100, |_| {
            inner.run();
            Ok::<_, ()>(())
        })
        .map(|_| ())
    })
    .unwrap();
    let max = inner.max.load(Ordering::SeqCst);
    assert!(
        max >= 2,
        "{max} inner tasks at once after the other index was done"
    );
}
