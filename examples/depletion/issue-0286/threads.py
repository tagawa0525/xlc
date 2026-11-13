"""#0286 のスレッドの数ごとの測定（`threads-*`、`ab-*`、`ab2-*`）を要約して標準出力に書く。

ログ（水準は `info,xlc_lattice::ufg=debug,xlc_lattice::shielding=debug`）ごとに、自己遮蔽の span の秒、1 ピンの段のスレッド秒
（`UFG slowing down` と `UFG transfers counted` を解いたピンの数で割ったもの）、自由ガスの散乱源の 1 回あたりの時間（`UFG free-gas
work` の秒を回数で割ったもの）、同じ名前の `.perf`（`perf stat -x,` の出力）のカウンタ、`.md5`（少数群定数の md5）と、1 状態の秒（`lattice state=` の行の `seconds`）を書く。

使い方: python3 threads.py <ディレクトリ> [<ログの glob>]（glob の既定は threads-*.log）
"""

import re
import sys
from pathlib import Path

D = Path(sys.argv[1])
F = re.compile(r"(\w+)=(\S+)")
UNIT = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "ns": 1e-9}


def secs(v):
    m = re.match(r"([\d.]+)(ms|s|µs|ns)", v)
    return float(m.group(1)) * UNIT[m.group(2)]


def fields(line):
    return {k: v for k, v in F.findall(line)}


rows = []
for log in sorted(D.glob(sys.argv[2] if len(sys.argv) > 2 else "threads-*.log")):
    name = log.stem
    st = {}
    pins = 0
    work = {}
    trans = 0.0
    shield = []
    total = None
    threads = None
    for line in log.open(errors="replace"):
        if "UFG slowing down" in line:
            pins += 1
            for k, v in fields(line).items():
                if k in (
                    "grid",
                    "table",
                    "evaluate",
                    "collision",
                    "slowing",
                    "free_gas",
                    "sab",
                ):
                    st[k] = st.get(k, 0.0) + float(v)
        elif "UFG free-gas work (thread seconds)" in line:
            for k, v in fields(line).items():
                try:
                    work[k] = work.get(k, 0.0) + float(v)
                except ValueError:
                    pass
        elif "UFG transfers counted" in line:
            trans += float(fields(line)["seconds"])
        elif " INFO shielding: xlc_lattice::shielding: close time.busy=" in line:
            shield.append(secs(line.split("time.busy=")[1].split()[0]))
        elif "thread pool threads=" in line:
            threads = fields(line)["threads"]
        elif "lattice state=" in line:
            total = fields(line).get("seconds")
    perf = {}
    p = log.with_suffix(".perf")
    if p.exists():
        for line in p.open():
            parts = line.strip().split(",")
            if len(parts) > 3 and parts[0] and not line.startswith("#"):
                try:
                    perf[parts[2]] = float(parts[0])
                except ValueError:
                    pass
    md5 = (
        log.with_suffix(".md5").read_text().split()[0][:8]
        if log.with_suffix(".md5").exists()
        else "-"
    )
    if pins == 0:
        continue
    per = {k: v / pins for k, v in st.items()}
    w = work
    ns = {
        "eval": 1e9 * w.get("quadratics", 0) / max(w.get("evaluations", 1), 1),
        "built": 1e9 * w.get("built", 0) / max(w.get("integrations_built", 1), 1),
        "held": 1e9 * w.get("held", 0) / max(w.get("integrations_held", 1), 1),
    }
    stage_sum = sum(
        st.get(k, 0) for k in ("grid", "table", "evaluate", "collision", "slowing")
    )
    print(
        f"== {name} threads={threads} md5={md5} pins={pins} shielding={shield} state={total}"
    )
    print(
        "  per pin (s): "
        + " ".join(f"{k}={v:.3f}" for k, v in per.items())
        + f" transfers={trans / pins:.3f} sum={(stage_sum + trans) / pins:.3f}"
    )
    print(
        "  ns/op: "
        + " ".join(f"{k}={v:.1f}" for k, v in ns.items())
        + f" fg_bounds={w.get('bounds', 0) / pins:.3f}"
    )
    if perf:
        tc = perf.get("task-clock", 0) / 1e3
        cu, ck = perf.get("cycles:u", 0), perf.get("cycles:k", 0)
        iu = perf.get("instructions:u", 0)
        print(
            f"  perf: task={tc:.0f}s cyc_u={cu:.3e} cyc_k={ck:.3e} ({ck / (cu + ck) * 100:.1f}%) ins_u={iu:.3e} IPC_u={iu / cu:.2f} "
            f"GHz={(cu + ck) / tc / 1e9:.2f} pf={perf.get('page-faults', 0):.3e} cs={perf.get('context-switches', 0):.3e} "
            f"mig={perf.get('cpu-migrations', 0):.3e} L1m={perf.get('L1-dcache-load-misses:u', 0):.3e} "
            f"L2={perf.get('ls_dmnd_fills_from_sys.local_l2:u', 0):.3e} ccx={perf.get('ls_dmnd_fills_from_sys.local_ccx:u', 0):.3e} "
            f"dram={perf.get('ls_dmnd_fills_from_sys.dram_io_all:u', 0):.3e}"
        )
