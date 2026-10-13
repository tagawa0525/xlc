"""VERA の OpenMC の参照解に仕込む標準のタリーの組（ADR-0152）。

statepoint には実行する前に定義したタリーしか入らないので、xlc と細かく比べるのに要る量を、問題の種類ごとに決まった
組としてすべての参照解に入れる。組の中身（名前、フィルタ、スコア）を変えたら `VERSION` を上げる。タリーの名前は
`standard-<版>/<名前>` で、statepoint、summary.h5、段階の model.xml（`openmc-1B-model` の出力。ADR-0136）から版が分かる。

群は xlc のライブラリの群構造（WIMS 69 群）の下と上に、範囲の外の 2 つのビン（0〜1e-5 eV と 10 MeV〜1e10 eV）を足した
71 のビン（`XLC_BINS`。昇順）。細かい群は ECCO 1968 群（`FINE_BINS`。範囲の外は数えない）。

問題の種類ごとの組（`pincell`、`assembly`、`assembly_3d`、`core_2d`）:

- すべて: `material groups`（物質 × 71 ビンの反応率）、`material transfer`（物質 × 入射 × 出射 × Legendre P0〜P3 の
  nu-scatter）、`material fine`（物質 × ECCO 1968 群の中性子束、吸収、生成）、`nuclides`（物質 × 71 ビン × 追う核種の
  反応率）
- 集合体（2D と 3D）: `pin groups`（ピンのセル × 71 ビンの反応率。セルの中のすべての物質）、`pin nuclides` と
  `pin nuclides flux`（ピンのセル × 燃料の物質 × 71 ビンの、燃料の核種の反応率と中性子束）、`boundary`（集合体の外周の
  水のギャップの帯と角 × 71 ビンの中性子束）
- 3D の集合体: `layer groups`（軸方向の層 × 71 ビンの反応率）、`layer transfer`（層 × 入射 × 出射の nu-scatter）、
  `layer current`（層の境の面 × 71 ビンの正味の中性子流）。`boundary` は層ごと
- 2D の 1/4 炉心: `assembly groups`、`assembly transfer`、`assembly current`（集合体のメッシュで同じ）と、集合体の境の
  水のギャップの帯と角の `boundary`
"""

import openmc
import openmc.mgxs

VERSION = 1
PREFIX = f"standard-{VERSION}"

# WIMS 69 群の境界 (eV、昇順)。xlc の crates/format/src/group.rs の WIMS69_ASCENDING と同じ（openmc_tallies_check.py が
# 照らす）。xlc のライブラリ（vera.msgpack）から読まないのは、xs のコードを変えるたびに作り直されるライブラリを、
# OpenMC の段階の deps にしないため（ADR-0130 の帰結）
WIMS69 = (
    1.0e-5, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.042, 0.05, 0.058, 0.067, 0.08, 0.1,
    0.14, 0.18, 0.22, 0.25, 0.28, 0.3, 0.32, 0.35, 0.4, 0.5, 0.625, 0.78, 0.85, 0.91, 0.95, 0.972,
    0.996, 1.02, 1.045, 1.071, 1.097, 1.123, 1.15, 1.3, 1.5, 2.1, 2.6, 3.3, 4.0, 9.877, 15.968,
    27.7, 48.052, 75.501, 148.728, 367.262, 906.898, 1425.1, 2239.45, 3519.1, 5530.0, 9118.0,
    1.503e4, 2.478e4, 4.085e4, 6.734e4, 1.11e5, 1.83e5, 3.025e5, 5.0e5, 8.21e5, 1.353e6, 2.231e6,
    3.679e6, 6.0655e6, 1.0e7,
)  # fmt: skip
# 範囲の外のビン。上は核データの上限（150 MeV 程度）より十分に上
XLC_BINS = (0.0, *WIMS69, 1.0e10)
FINE_BINS = tuple(openmc.mgxs.GROUP_STRUCTURES["ECCO-1968"])

# 2 個以上の中性子を出す反応（ENDF の MT 11、16、17、24、25、30、37、41、42）。正味の吸収（吸収 − Σ (放出の数 − 1)
# × 率）を、xlc の全断面積 − 散乱行列の行の和と同じ量にするのに要る（openmc_balance.py と同じ）
MULTIPLYING = (
    "(n,2nd)",
    "(n,2n)",
    "(n,3n)",
    "(n,2na)",
    "(n,3na)",
    "(n,2n2a)",
    "(n,4n)",
    "(n,2np)",
    "(n,3np)",
)
# 領域（物質、ピン、層、集合体）ごとの反応率。total、absorption、fission、nu-fission は OpenMC が物質ごとに持つ巨視的
# 断面積で数えられ、MULTIPLYING のような反応は核種ごとに断面積を引くので、物質の組にだけ入れる
ZONE_SCORES = ("flux", "total", "scatter", "absorption", "fission", "nu-fission")
MATERIAL_SCORES = (*ZONE_SCORES, *MULTIPLYING)
FINE_SCORES = ("flux", "absorption", "nu-fission")
# 燃料と吸収体と減速材の、反応率を核種ごとに追う核種（モデルにあるものだけ入れる）
TRACKED = (
    "U234",
    "U235",
    "U236",
    "U238",
    "O16",
    "H1",
    "B10",
    "Ag107",
    "Ag109",
    "In113",
    "In115",
    "Cd113",
    "Gd155",
    "Gd157",
)
NUCLIDE_SCORES = (
    "total",
    "elastic",
    "absorption",
    "fission",
    "nu-fission",
    "(n,gamma)",
    "(n,2n)",
)
# 散乱行列の Legendre の次数。xlc の散乱は P3（ADR-0075）
LEGENDRE_ORDER = 3

# 入力のキー。`standard-tallies = false` は組を入れない（時間の比較のため）。書かなければ入れる（ADR-0152）
OPTIONAL = {"standard-tallies"}


def quiet(settings: openmc.Settings) -> None:
    """tallies.out（テキスト）を書かない。標準の組では数十 MB になり、statepoint と同じ値しか持たない。"""
    settings.output = {"tallies": False}


def enabled(run: dict) -> bool:
    """入力が標準のタリーの組を入れるか。"""
    value = run.get("standard-tallies", True)
    if not isinstance(value, bool):
        raise SystemExit(f"standard-tallies must be true or false, got {value!r}")
    return value


def name(kind: str) -> str:
    return f"{PREFIX}/{kind}"


def tally(
    kind: str,
    filters: list[openmc.Filter],
    scores: tuple[str, ...],
    nuclides: list[str] | None = None,
) -> openmc.Tally:
    t = openmc.Tally(name=name(kind))
    t.filters = filters
    t.scores = list(scores)
    if nuclides:
        t.nuclides = nuclides
    return t


def tracked(materials: list[openmc.Material]) -> list[str]:
    """`TRACKED` のうち、`materials` のどれかにある核種（`TRACKED` の順）。"""
    present = {n for m in materials for n in m.get_nuclides()}
    return [n for n in TRACKED if n in present]


class Filters:
    """1 つのモデルのタリーで共有するフィルタ。同じフィルタ（同じ ID）を使うタリーの間では、OpenMC が飛跡ごとの
    ビンの探索（メッシュの交差など）を 1 回で済ませる（src/tallies/filter.cpp の filter_matches）。"""

    def __init__(self, materials: list[openmc.Material]):
        self.materials = openmc.MaterialFilter(materials)
        self.groups = openmc.EnergyFilter(XLC_BINS)
        self.groups_out = openmc.EnergyoutFilter(XLC_BINS)
        self.fine = openmc.EnergyFilter(FINE_BINS)
        self.legendre = openmc.LegendreFilter(LEGENDRE_ORDER)
        self.meshes: dict[int, openmc.MeshFilter] = {}

    def mesh(self, mesh: openmc.MeshBase) -> openmc.MeshFilter:
        """`mesh` のフィルタ（メッシュごとに 1 つ。OpenMC は同じビンの別のフィルタを 1 つにまとめるが、ID を付け直す
        警告を出す）。"""
        return self.meshes.setdefault(mesh.id, openmc.MeshFilter(mesh))


def material_tallies(
    f: Filters, materials: list[openmc.Material]
) -> list[openmc.Tally]:
    """すべての問題に入れる、物質ごとの組。"""
    # 散乱の出来事ごとの推定（analog）が要る nu-scatter と EnergyoutFilter は別のタリーにする（同じタリーだと全体が
    # 出来事ごとの推定になる）
    return [
        tally("material groups", [f.materials, f.groups], MATERIAL_SCORES),
        tally(
            "material transfer",
            [f.materials, f.groups, f.groups_out, f.legendre],
            ("nu-scatter",),
        ),
        tally("material fine", [f.materials, f.fine], FINE_SCORES),
        tally("nuclides", [f.materials, f.groups], NUCLIDE_SCORES, tracked(materials)),
    ]


def zone_tallies(
    f: Filters, kind: str, mesh: openmc.MeshBase, transfer: bool, current: bool
) -> list[openmc.Tally]:
    """メッシュの領域（`kind` は名前の頭）ごとの反応率と、`transfer` なら P0 の nu-scatter の行列、`current` なら
    メッシュの面の正味の中性子流。"""
    zones = f.mesh(mesh)
    tallies = [tally(f"{kind} groups", [zones, f.groups], ZONE_SCORES)]
    if transfer:
        tallies.append(
            tally(f"{kind} transfer", [zones, f.groups, f.groups_out], ("nu-scatter",))
        )
    if current:
        tallies.append(
            tally(
                f"{kind} current",
                [openmc.MeshSurfaceFilter(mesh), f.groups],
                ("current",),
            )
        )
    return tallies


def pin_fuel_tallies(
    f: Filters, pins: openmc.MeshBase, fuels: list[openmc.Material]
) -> list[openmc.Tally]:
    """ピンのセル × 燃料の物質ごとの、燃料の核種の反応率と中性子束（ピンごとの実効断面積）。"""
    filters = [f.mesh(pins), openmc.MaterialFilter(fuels), f.groups]
    return [
        tally("pin nuclides", filters, NUCLIDE_SCORES, tracked(fuels)),
        tally("pin nuclides flux", filters, ("flux",)),
    ]


def strip_edges(low: float, faces: list[float], high: float, gap: float) -> list[float]:
    """`low`〜`high` を、面 `faces` の両側 `gap` の帯で分けるメッシュの境界（昇順）。範囲の外の帯は切る。"""
    edges = {low, high}
    for face in faces:
        edges |= {e for e in (face - gap, face + gap) if low < e < high}
    return sorted(edges)


def boundary_mesh(
    x: list[float], y: list[float], z: list[float] | None = None
) -> openmc.RectilinearMesh:
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = x
    mesh.y_grid = y
    # 2D の形状は z に一様で境界がないので、z は十分に広くとる
    mesh.z_grid = z if z is not None else [-1.0e6, 1.0e6]
    return mesh


def boundary_tally(f: Filters, mesh: openmc.MeshBase) -> openmc.Tally:
    """面の両側の水のギャップの帯と角の中性子束。鏡面反射の面では、帯の平均は面の中性子束と帯の幅の 2 乗の差しか
    違わない（面の法線の向きの勾配が 0）。集合体の間の面では、面を挟む対称な帯の平均が同じく 2 次で面の値になる。"""
    return tally("boundary", [f.mesh(mesh), f.groups], ("flux",))


def pincell(materials: list[openmc.Material]) -> list[openmc.Tally]:
    """ピンセル（1B〜1D）の組。ピンは 1 つなので、物質ごとの組だけ（`nuclides` が燃料の核種を含む）。"""
    return material_tallies(Filters(materials), materials)


def assembly(
    materials: list[openmc.Material],
    fuels: list[openmc.Material],
    pins: openmc.MeshBase,
    half: float,
    gap: float,
) -> list[openmc.Tally]:
    """外周が鏡面反射の 2D の集合体（2A〜2Q）の組。`pins` はピンのセルのメッシュ、`half` は外周の箱の半分の幅、
    `gap` は格子の外の水のギャップの幅。"""
    f = Filters(materials)
    edges = strip_edges(-half, [-half, half], half, gap)
    return [
        *material_tallies(f, materials),
        *zone_tallies(f, "pin", pins, transfer=False, current=False),
        *pin_fuel_tallies(f, pins, fuels),
        boundary_tally(f, boundary_mesh(edges, edges)),
    ]


def assembly_3d(
    materials: list[openmc.Material],
    fuels: list[openmc.Material],
    pins: openmc.MeshBase,
    half: float,
    gap: float,
    levels: list[float],
) -> list[openmc.Tally]:
    """外周が鏡面反射、上下が真空の 3D の集合体（3A、3B）の組。`pins` は燃料の全長で積分したピンのセルのメッシュ、
    `levels` は軸方向の層の境界（昇順。形状の下端から上端まで）。"""
    f = Filters(materials)
    layers = boundary_mesh([-half, half], [-half, half], levels)
    edges = strip_edges(-half, [-half, half], half, gap)
    return [
        *material_tallies(f, materials),
        *zone_tallies(f, "pin", pins, transfer=False, current=False),
        *pin_fuel_tallies(f, pins, fuels),
        *zone_tallies(f, "layer", layers, transfer=True, current=True),
        boundary_tally(f, boundary_mesh(edges, edges, levels)),
    ]


def core_2d(
    materials: list[openmc.Material],
    assemblies: openmc.RegularMesh,
    gap: float,
) -> list[openmc.Tally]:
    """2D の 1/4 炉心（5A〜5C-2D）の組。`assemblies` は集合体のメッシュ（x ≥ 0、y ≤ 0 の 1/4 で、中心の集合体の中心が
    原点）、`gap` は集合体の格子の外の水のギャップの幅（集合体の間の水は 2 × `gap`）。"""
    f = Filters(materials)
    n = assemblies.dimension[0]
    pitch = assemblies.width[0]
    faces = [(c + 0.5) * pitch for c in range(n)]
    x = strip_edges(0.0, faces, faces[-1] + gap, gap)
    y = [-e for e in reversed(x)]
    return [
        *material_tallies(f, materials),
        *zone_tallies(f, "assembly", assemblies, transfer=True, current=True),
        boundary_tally(f, boundary_mesh(x, y)),
    ]
