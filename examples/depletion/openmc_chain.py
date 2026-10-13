"""OpenMC の `Chain.from_endf` で、xlc の全核種の鎖（ADR-0164）と同じ ENDF/B-VIII.0 の decay、nfy、neutrons から鎖を作る。

使い方（references の devShell で。dvc.yaml の openmc-chain の段階が呼ぶ）:

    python3 openmc_chain.py <ENDF の根（ENDF-B-VIII.0_* の親）> <出力の XML>

反応は `from_endf` の既定（(n,2n)、(n,3n)、(n,4n)、(n,gamma)、(n,p)、(n,a) と核分裂）で、xlc の鎖と同じ。
"""

import sys
from pathlib import Path

from openmc.deplete import Chain


def endf_files(root: Path, part: str) -> list[str]:
    files = sorted(str(p) for p in (root / f"ENDF-B-VIII.0_{part}").glob("*.endf"))
    if not files:
        sys.exit(f"no *.endf in {root / f'ENDF-B-VIII.0_{part}'}")
    return files


def main():
    root, output = Path(sys.argv[1]), Path(sys.argv[2])
    chain = Chain.from_endf(
        endf_files(root, "decay"),
        endf_files(root, "nfy"),
        endf_files(root, "neutrons"),
        progress=False,
    )
    chain.export_to_xml(output)


if __name__ == "__main__":
    main()
