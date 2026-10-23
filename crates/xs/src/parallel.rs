//! 独立な計算の並列化。

use std::any::Any;
use std::num::NonZeroUsize;
use std::panic::AssertUnwindSafe;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};
use std::thread::Scope;

/// `in_parallel` が呼んだスレッドのほかに起こすスレッドの予算（プロセスで共有する）。呼んだスレッドも計算するので、
/// 1 つの呼び出しが使うのは予算 + 1 で、`set_threads` の数か、なければ物理コアの数（ADR-0067。このプロセスに許された
/// CPU の数で抑える）にする
static BUDGET: OnceLock<Budget> = OnceLock::new();

fn budget() -> &'static Budget {
    BUDGET.get_or_init(|| Budget::new(cores().saturating_sub(1)))
}

/// `in_parallel` が同時に使うスレッドの数を `threads` にする。最初の `in_parallel` より前に 1 度だけ呼ぶ
pub fn set_threads(threads: NonZeroUsize) -> Result<(), ThreadsAlreadySet> {
    BUDGET
        .set(Budget::new(threads.get() - 1))
        .map_err(|_| ThreadsAlreadySet)
}

/// `set_threads` を、スレッドの数が決まった後に呼んだ
#[derive(Debug, thiserror::Error)]
#[error("the number of threads is already set")]
pub struct ThreadsAlreadySet;

/// `f(i)`（i = 0..n）を、呼んだスレッドと、プロセスで共有する予算（`BUDGET`）から取るスレッドで並列に求め、i の順に並べて
/// 返す。
/// 1 つの i の重さは違いうるので、各スレッドは次の i を共有のカウンタから取る。i を取るたびに、残りの i があり予算に空きが
/// あれば、スレッドを 1 つ足す。取る i がなくなった足したスレッドは、すぐに予算に返す（呼んだスレッドは、足したスレッドが
/// 終わるのを待つ間も、外側の呼び出しから取った予算を持ったまま）。こうして、入れ子の呼び出し（核種、
/// 温度、格子点）のどれかが空いたスレッドを、ほかの呼び出しが次に i を取るときに使う。始めるときに空きがなくても、
/// 呼んだスレッドで順に求めながら、空いたスレッドを足す（f(i) を求めている間に空いたスレッドは、その i が終わるまで
/// この呼び出しには足さない）。
/// rayon を使わないのは、xs に依存を足さないため。
///
/// エラーが出たら、各スレッドは次に i を取る前にそれを見て止まる（見たあとに取った i は、スレッドごとに多くとも 1 つ
/// 余分に求める）。返すのは i の最も小さいエラーで、カウンタは増える一方なので、エラーの i より小さい i はすべて
/// 取られて最後まで求められており、返すエラーはスレッドの進み方によらない。f(i) がパニックしたときは、パニックを投げ直す
/// （エラーとパニックが混ざると、どちらになるかはスレッドの進み方による）。
pub fn in_parallel<T: Send, E: Send>(
    n: usize,
    f: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    in_parallel_in(budget(), n, f)
}

/// `in_parallel` を、予算 `budget` から取るスレッドで求める
fn in_parallel_in<T: Send, E: Send, F: Fn(usize) -> Result<T, E> + Sync>(
    budget: &Budget,
    n: usize,
    f: F,
) -> Result<Vec<T>, E> {
    let call = Call {
        budget,
        n,
        f,
        next: AtomicUsize::new(0),
        failed: AtomicBool::new(false),
        done: Mutex::new(Vec::new()),
        panic: Mutex::new(None),
    };
    std::thread::scope(|scope| call.help(scope));
    if let Some(payload) = call
        .panic
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner)
    {
        std::panic::resume_unwind(payload);
    }
    let mut done = call
        .done
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    done.sort_unstable_by_key(|&(i, _)| i);
    done.into_iter().map(|(_, value)| value).collect()
}

/// 1 回の `in_parallel_in` の、スレッドが共有する状態
struct Call<'b, T, E, F> {
    budget: &'b Budget,
    n: usize,
    f: F,
    /// 次に取る i
    next: AtomicUsize,
    failed: AtomicBool,
    /// 求めた (i, f(i))
    done: Mutex<Vec<(usize, Result<T, E>)>>,
    /// パニックの最初のもの（呼んだスレッドで投げ直す）
    panic: Mutex<Option<Box<dyn Any + Send>>>,
}

impl<T: Send, E: Send, F: Fn(usize) -> Result<T, E> + Sync> Call<'_, T, E, F> {
    /// 取る i がなくなるか、エラーが出るまで求める。i を取るたびに、残りがあればスレッドを足す
    fn work<'scope>(&'scope self, scope: &'scope Scope<'scope, '_>) {
        let mut done = Vec::new();
        while !self.failed.load(Ordering::Relaxed) {
            let i = self.next.fetch_add(1, Ordering::Relaxed);
            if i >= self.n {
                break;
            }
            if i + 1 < self.n {
                let slot = self.budget.take(1);
                if slot.count == 1 {
                    // 足したスレッドでも、呼んだスレッドの span（核種の名前など）の中でログを出す
                    let span = tracing::Span::current();
                    scope.spawn(move || {
                        // スレッドが終わるまで予算を持つ
                        let _slot = slot;
                        let _span = span.enter();
                        self.help(scope);
                    });
                }
            }
            let value = (self.f)(i);
            if value.is_err() {
                self.failed.store(true, Ordering::Relaxed);
            }
            done.push((i, value));
        }
        self.done
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .extend(done);
    }

    /// `work` を実行し、パニックしたらほかのスレッドを止める（呼んだスレッドも足したスレッドも同じ）。パニックは、元の値のまま
    /// 呼んだスレッドで投げ直すために取っておく（`std::thread::scope` に任せると、値が「スレッドがパニックした」に置き換わる）
    fn help<'scope>(&'scope self, scope: &'scope Scope<'scope, '_>) {
        if let Err(payload) = std::panic::catch_unwind(AssertUnwindSafe(|| self.work(scope))) {
            self.failed.store(true, Ordering::Relaxed);
            self.panic
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .get_or_insert(payload);
        }
    }
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
