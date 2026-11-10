"""Compare the per-nuclide xs outputs before and after adding temperatures to a library.

usage: compare_library_temperatures.py <old dvc.lock> <new dvc.lock> <dvc cache dir> <stage dir>

The old outputs are read from the DVC cache by the md5 in the old lock, the new ones from the stage dir
(their md5 must match the new lock). Every field except the per-nuclide `temperatures` must be equal
after a msgpack re-pack, and every old temperature table must be equal to the new table at the same
temperature. Exits 1 on any mismatch.
"""

import hashlib
import os
import sys

import msgpack
import yaml


def load(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), strict_map_key=False)


def stages(lock):
    with open(lock) as f:
        return yaml.safe_load(f)["stages"]


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def main(old_lock, new_lock, cache, stage_dir):
    a = stages(old_lock)
    b = stages(new_lock)
    pk = msgpack.packb
    mismatch = []
    n_nuc = n_temp = 0
    new_temp_sets = {}
    for k, v in a.items():
        if not k.startswith("xs-nuclide@"):
            continue
        m = v["outs"][0]["md5"]
        rel = v["outs"][0]["path"]
        assert b[k]["outs"][0]["path"] == rel
        old = load(f"{cache}/files/md5/{m[:2]}/{m[2:]}")
        new_path = os.path.join(stage_dir, rel)
        assert md5(new_path) == b[k]["outs"][0]["md5"], k
        new = load(new_path)
        if set(old) != set(new):
            mismatch.append((k, "top keys", set(old) ^ set(new)))
            continue
        for key in old:
            if key in ("nuclides", "activation"):
                continue
            if pk(old[key]) != pk(new[key]):
                mismatch.append((k, key))
        for kind in ("nuclides", "activation"):
            if len(old[kind]) != len(new[kind]):
                mismatch.append((k, kind, "count"))
                continue
            for on, nn in zip(old[kind], new[kind]):
                n_nuc += 1
                if set(on) != set(nn):
                    mismatch.append((k, kind, on.get("name"), "keys"))
                    continue
                for f in on:
                    if f == "temperatures":
                        continue
                    if pk(on[f]) != pk(nn[f]):
                        mismatch.append((k, kind, on["name"], f))
                nt = {t["temperature"]: t for t in nn["temperatures"]}
                new_temp_sets.setdefault(tuple(sorted(nt)), []).append(on["name"])
                for t in on["temperatures"]:
                    n_temp += 1
                    T = t["temperature"]
                    if T not in nt:
                        mismatch.append((k, kind, on["name"], T, "missing"))
                        continue
                    for f in t:
                        if pk(t[f]) != pk(nt[T].get(f)):
                            mismatch.append((k, kind, on["name"], T, f))
                    if set(t) != set(nt[T]):
                        mismatch.append((k, kind, on["name"], T, "keys"))

    n_stages = sum(1 for k in a if k.startswith("xs-nuclide@"))
    print(f"stages={n_stages} nuclide entries={n_nuc} old temperature tables={n_temp}")
    for ts, names in new_temp_sets.items():
        listed = f" ({', '.join(names)})" if len(names) < 5 else ""
        print(f"new temperatures {list(ts)}: {len(names)} entries{listed}")
    print(f"mismatches={len(mismatch)}")
    for x in mismatch[:50]:
        print(x)
    return 1 if mismatch else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
