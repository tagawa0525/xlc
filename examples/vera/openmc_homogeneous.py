"""均質な無限媒質（`homogeneous_medium.py` の入力）を OpenMC の連続エネルギーモードで解き、k を書く。

使い方: python3 openmc_homogeneous.py <cross_sections.xml> <媒質 (TOML)> <作業ディレクトリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

形状は外周が鏡面反射の立方体 1 つ。H-1 があれば、軽水の H の熱中性子散乱（c_H_in_H2O）を使う（xlc の H-1 と
同じ扱い）。
"""

import sys
from pathlib import Path

import openmc
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import eigenvalue_settings, material


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    cross_sections = Path(sys.argv[1]).resolve()
    medium = tomllib.loads(Path(sys.argv[2]).read_text())
    work = Path(sys.argv[3]).resolve()
    particles = int(sys.argv[4]) if len(sys.argv) > 4 else 200_000
    batches = int(sys.argv[5]) if len(sys.argv) > 5 else 300
    inactive = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    work.mkdir(parents=True, exist_ok=True)

    densities = {k.replace("-", ""): v for k, v in medium["nuclides"].items()}
    mix = material("mix", densities, medium["temperature"])
    if "H1" in densities:
        mix.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([mix])
    materials.cross_sections = str(cross_sections)
    half = 1.0
    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    top = openmc.ZPlane(z0=half, boundary_type="reflective")
    bottom = openmc.ZPlane(z0=-half, boundary_type="reflective")
    geometry = openmc.Geometry([openmc.Cell(fill=mix, region=-box & -top & +bottom)])
    settings = eigenvalue_settings(particles, batches, inactive, seed, half)
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    model = openmc.Model(geometry, materials, settings)
    statepoint = model.run(cwd=work, output=False)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
    print(f"{sys.argv[2]} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
