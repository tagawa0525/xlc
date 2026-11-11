"""孤立した共鳴の実効共鳴積分の温度への依存のべきを、狭い共鳴の近似（NR）の J 関数で見る（#0274 の「燃料の温度への依存の
モデルと当てはめの式」）。

使い方: python3 doppler_exponents.py

J(ξ, β) = ∫ ψ(ξ, x) / (ψ(ξ, x) + β) dx（x = 2(E − E0)/Γ、ψ は Doppler で広がった形（Voigt）、ξ = Γ/Γ_D、
Γ_D = √(4 E0 kT / A)、β = σ0 / σ_peak）。実効共鳴積分は J に比例し、温度は ξ ∝ T^-1/2 だけを通して入る。共鳴（`RESONANCES`）と
β（`BETAS`）ごとに、`TEMPERATURES` の温度で、dJ/dT ∝ T^-n とみた局所のべき n = −d ln(dJ/dT) / d ln T と、J の 560 K から
1350 K までの相対変化を出す。n = 1/2 なら J は √T の 1 次、n = 0 なら T の 1 次、n = 1 なら ln T の 1 次。
"""

import math

from scipy.integrate import quad
from scipy.special import wofz

BOLTZMANN = 8.617333262e-5  # eV/K

# (名前, E0 (eV), Γ (eV), A（中性子の質量との比）)。ENDF/B-VIII.0 の MF2/MT151（Reich-Moore）の、中性子、捕獲、核分裂の幅の和
# （references/endf/ENDF-B-VIII.0_neutrons の n-092_U_238、n-094_Pu_239、n-094_Pu_240）
RESONANCES = [
    ("U-238 6.674 eV", 6.67428, 1.4923e-3 + 0.022711, 236.006),
    ("U-238 20.87 eV", 20.87177, 1.0076e-2 + 0.022753, 236.006),
    ("U-238 36.68 eV", 36.68283, 3.3592e-2 + 0.022264, 236.006),
    ("Pu-239 0.296 eV", 0.295624, 7.9e-5 + 0.039824 + 0.056197, 236.9986),
    ("Pu-240 1.056 eV", 1.05600, 2.454e-3 + 0.029148 + 8e-6, 237.9916),
]

BETAS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e3]
TEMPERATURES = [600.0, 900.0, 1200.0]


def psi(xi: float, x: float) -> float:
    """Doppler で広がった共鳴の形 ψ(ξ, x) = (√π ξ / 2) Re w(ξ (x + i) / 2)（Faddeeva の関数 w）"""
    return math.sqrt(math.pi) / 2 * xi * wofz(xi * (x + 1j) / 2).real


def j_function(xi: float, beta: float) -> float:
    # ψ は x について偶関数
    value, _ = quad(
        lambda x: psi(xi, x) / (psi(xi, x) + beta),
        0,
        math.inf,
        limit=500,
        epsabs=0,
        epsrel=1e-11,
    )
    return 2 * value


def xi_at(energy: float, width: float, mass: float, temperature: float) -> float:
    return width / math.sqrt(4 * energy * BOLTZMANN * temperature / mass)


def local_exponent(energy, width, mass, beta, temperature, step=0.02) -> float:
    """n = −d ln(dJ/dT) / d ln T（ln T の中心差分。dJ/dT も中心差分）"""

    def slope(t):
        h = t * step / 4
        return (
            j_function(xi_at(energy, width, mass, t + h), beta)
            - j_function(xi_at(energy, width, mass, t - h), beta)
        ) / (2 * h)

    up, down = temperature * math.exp(step), temperature * math.exp(-step)
    return -(math.log(abs(slope(up))) - math.log(abs(slope(down)))) / (2 * step)


def main() -> None:
    print(
        "共鳴、ξ(900 K)、β: n（"
        + "、".join(f"{t:g} K" for t in TEMPERATURES)
        + "）、J の 560→1350 K の相対変化"
    )
    for name, energy, width, mass in RESONANCES:
        for beta in BETAS:
            ns = [local_exponent(energy, width, mass, beta, t) for t in TEMPERATURES]
            low, high = (
                j_function(xi_at(energy, width, mass, t), beta) for t in (560.0, 1350.0)
            )
            print(
                f"{name}、{xi_at(energy, width, mass, 900.0):.3f}、{beta:g}: "
                + "、".join(f"{n:+.3f}" for n in ns)
                + f"、{(high - low) / low:+.3e}"
            )


if __name__ == "__main__":
    main()
