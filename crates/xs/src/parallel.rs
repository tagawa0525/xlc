//! 独立な計算の並列化。

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// `f(i)`（i = 0..n）を、物理コアの数（ADR-0067。このプロセスに許された CPU の数で抑える）のスレッドで並列に求め、
/// i の順に並べて返す。
/// 1 つの i の重さは違いうるので、各スレッドは次の i を共有のカウンタから取る。
/// rayon を使わないのは、xs が依存を持たないため（`OnceLock` の初期化の中で rayon のプールを待つと止まりうる。#0088）。
///
/// エラーが出たら、各スレッドは次に i を取る前にそれを見て止まる（見たあとに取った i は、スレッドごとに多くとも 1 つ
/// 余分に求める）。返すのは i の最も小さいエラーで、カウンタは増える一方なので、エラーの i より小さい i はすべて
/// 取られて最後まで求められており、返すエラーはスレッドの進み方によらない。
pub(crate) fn in_parallel<T: Send, E: Send>(
    n: usize,
    f: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    let allowed = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    let threads = num_cpus::get_physical().min(allowed).min(n);
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let mut done: Vec<(usize, Result<T, E>)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
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
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    done.sort_unstable_by_key(|&(i, _)| i);
    done.into_iter().map(|(_, value)| value).collect()
}
