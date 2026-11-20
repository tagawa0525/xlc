"""Be-9 の MF6 MT16 の中性子（LAW=7）の入射エネルギーの補間則 INT=2 を、NJOY GROUPR が使う単位基底 INT=22 に替えたコピーを書く。

NJOY との差の原因を診断するためだけのもので、xs の本番の経路ではない（NJOY の選択を xs の機能として足さない。ADR-0271）。
"""

import sys
from pathlib import Path

src, dst = sys.argv[1], sys.argv[2]
lines = Path(src).read_text().splitlines(keepends=True)
start = next(i for i, l in enumerate(lines) if l[70:75] == " 6 16")
# 節の 2 行目が中性子の多重度の TAB1 の見出し（NR、NP）。その後に NBT/INT と点が続き、次が LAW=7 の外側の TAB2
nr, np_ = int(lines[start + 1][44:55]), int(lines[start + 1][55:66])
head = start + 2 + (2 * nr + 5) // 6 + (2 * np_ + 5) // 6
target = head + 1
assert lines[target][11:22].strip() == "2", lines[target]
lines[target] = lines[target][:11] + f"{22:>11}" + lines[target][22:]
Path(dst).write_text("".join(lines))
