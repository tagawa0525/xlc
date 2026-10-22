//! 独立な計算の並列化。

use std::sync::LazyLock;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// `in_parallel` が使うスレッドの予算（プロセスで共有する）
static BUDGET: LazyLock<Budget> = LazyLock::new(|| Budget::new(cores()));

/// `f(i)`（i = 0..n）を、物理コアの数（ADR-0067。このプロセスに許された CPU の数で抑える）のスレッドで並列に求め、
/// i の順に並べて返す。
/// 1 つの i の重さは違いうるので、各スレッドは次の i を共有のカウンタから取る。
/// スレッドはプロセスで共有する予算（`BUDGET`）から取り、空きがなければ呼んだスレッドで順に求める。核種ごとのスレッドの中で
/// 格子点を並列にするときに、スレッドが核種の数 × コアの数に増えないようにする。
/// rayon を使わないのは、xs に依存を足さないため。
///
/// エラーが出たら、各スレッドは次に i を取る前にそれを見て止まる（見たあとに取った i は、スレッドごとに多くとも 1 つ
/// 余分に求める）。返すのは i の最も小さいエラーで、カウンタは増える一方なので、エラーの i より小さい i はすべて
/// 取られて最後まで求められており、返すエラーはスレッドの進み方によらない。
pub(crate) fn in_parallel<T: Send, E: Send>(
    n: usize,
    f: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    in_parallel_in(&BUDGET, n, f)
}

/// `in_parallel` を、予算 `budget` から取るスレッドで求める
fn in_parallel_in<T: Send, E: Send>(
    budget: &Budget,
    n: usize,
    f: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    let threads = budget.take(n);
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let work = || {
        let mut done = Vec::new();
        while !failed.load(Ordering::Relaxed) {
            let i = next.fetch_add(1, Ordering::Relaxed);
            if i >= n {
                break;
            }
            let value = f(i);
            if value.is_err() {
                failed.store(true, Ordering::Relaxed);
            }
            done.push((i, value));
        }
        done
    };
    let mut done: Vec<(usize, Result<T, E>)> = if threads.count == 0 {
        work()
    } else {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads.count).map(|_| scope.spawn(work)).collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        })
    };
    drop(threads);
    done.sort_unstable_by_key(|&(i, _)| i);
    done.into_iter().map(|(_, value)| value).collect()
}

/// 同時に使ってよいスレッドの数と、今使っている数
struct Budget {
    busy: AtomicUsize,
    limit: usize,
}

impl Budget {
    fn new(limit: usize) -> Self {
        Self {
            busy: AtomicUsize::new(0),
            limit,
        }
    }

    /// 空きのうち、多くとも `wanted` を取る（空きがなければ 0）
    fn take(&self, wanted: usize) -> Reservation<'_> {
        let mut busy = self.busy.load(Ordering::Relaxed);
        loop {
            let count = self.limit.saturating_sub(busy).min(wanted);
            match self.busy.compare_exchange_weak(
                busy,
                busy + count,
                Ordering::AcqRel,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    return Reservation {
                        budget: self,
                        count,
                    };
                }
                Err(now) => busy = now,
            }
        }
    }
}

/// `Budget` から取ったスレッドの数。落とすと（パニックのときも）返す
struct Reservation<'a> {
    budget: &'a Budget,
    count: usize,
}

impl Drop for Reservation<'_> {
    fn drop(&mut self) {
        self.budget.busy.fetch_sub(self.count, Ordering::AcqRel);
    }
}

/// 物理コアの数（ADR-0067）を、このプロセスに許された CPU の数で抑えたもの
fn cores() -> usize {
    let allowed = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    num_cpus::get_physical().min(allowed)
}

#[cfg(test)]
mod tests;
