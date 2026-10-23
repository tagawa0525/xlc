"""ADR-0189 の案 C（#0157）の試行: 核種を CPU 時間の長い順に、合計 16 スレッドの中で同時に実行する。

各核種のスレッドの数は ceil(CPU 時間 / T) を 1〜16 に抑えたもの（T は CPU 時間の和 / 16）。

Usage: python3 runner_trial.py <xs> <入力のディレクトリ> <核種と CPU 時間の一覧> <出力のディレクトリ>
"""

import math
import subprocess
import sys
import time

BUDGET = 16


def main() -> None:
    xs, directory, cpu_file, out = sys.argv[1:5]
    with open(cpu_file) as f:
        cpu = [(line.split()[0], float(line.split()[1])) for line in f]
    cpu.sort(key=lambda x: -x[1])
    unit = sum(c for _, c in cpu) / BUDGET
    queue = [(n, min(BUDGET, max(1, math.ceil(c / unit)))) for n, c in cpu]
    budget = BUDGET
    running = {}
    start = time.time()
    while queue or running:
        for item in list(queue):
            n, p = item
            if p <= budget:
                cmd = [
                    "/run/current-system/sw/bin/time",
                    "-f",
                    f"nuclide {n} threads={p} elapsed=%e user=%U sys=%S",
                    xs,
                    "--threads",
                    str(p),
                    "nuclide",
                    "xs.toml",
                    "../xs-range.toml",
                    n,
                    f"{out}/{n}.msgpack",
                ]
                # ログは量が読めないので、パイプでなくファイルに受ける（パイプが詰まると子が止まる）。
                # 子はファイルを開いたまま持つので、こちらはすぐ閉じてよい
                with open(f"{out}/{n}.log", "w") as log:
                    proc = subprocess.Popen(cmd, cwd=directory, stderr=log)
                running[proc] = (n, p)
                budget -= p
                queue.remove(item)
        time.sleep(0.2)
        for proc in list(running):
            if proc.poll() is not None:
                n, p = running.pop(proc)
                budget += p
                with open(f"{out}/{n}.log") as f:
                    last = [line for line in f if line.startswith("nuclide ")]
                last = [line.rstrip("\n") for line in last]
                print(f"{time.time() - start:.1f}", *last, flush=True)
                if proc.returncode != 0:
                    sys.exit(f"{n} failed")
    print(f"runner total elapsed={time.time() - start:.1f}")


if __name__ == "__main__":
    main()
