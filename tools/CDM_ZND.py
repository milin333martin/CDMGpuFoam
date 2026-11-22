#!/usr/bin/env python3
"""
CDM_ZND.py  --  part of CDMGpuFoam

CJ / von Neumann / ZND calculator for the Chemical-Diffusive Model (CDM)
of Kaplan et al. (2019) and Lu, Kaplan & Oran, Combust. Flame 230 (2021) 111417.

Model (single irreversible step R -> P, calorically perfect gas, same
gamma and Mw for reactant and product):

    dY/dt = -A rho Y exp(-Ta/T)          (Y = reactant mass fraction)
    heat of reaction q per unit mass of reactant consumed

What it does
    1. Converts the parameters as printed in the paper tables
       (A in cm^3/(g s), Ea/(R T0), q Mw/(R T0), kappa0 in g/(s cm K^0.7))
       into SI, and prints ready-to-paste OpenFOAM dictionaries.
    2. Computes D_CJ, the von Neumann state, the ZND structure,
       the half-reaction distance x_d, induction/heat-release lengths,
       the effective activation energy E = Ea/(R T_vn), and the cell size
       predicted by Eq. (11) of Lu et al. (2021).
    3. Optionally writes the ZND profile to a table that the
       setZNDProfile utility maps onto an OpenFOAM mesh.

Examples
    # Table 1 of Lu et al. (stoichiometric H2-air, x_d-calibrated)
    python3 CDM_ZND.py --paper --A 1.064e13 --EaByRT0 54.95 --qMwByRT0 46.36 \
        --kappa0 2.367e-5 --gamma 1.177 --Mw 24.36

    # Table 6, phi = 1 (lambda-calibrated), and write the profile
    python3 CDM_ZND.py --paper --A 1.01e12 --EaByRT0 39.98 --qMwByRT0 62.97 \
        --kappa0 2.00e-5 --gamma 1.125 --Mw 22.77 --write znd.dat
"""

import argparse
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

RU = 8314.462618  # J/(kmol K)  (OpenFOAM uses kmol-based molWeight)


# ---------------------------------------------------------------------------
# Thermodynamics of the CDM gas
# ---------------------------------------------------------------------------

class CDMGas:
    def __init__(self, A, Ta, q, kappa0, gamma, Mw, n=0.7, Pr=1.0, Le=1.0):
        self.A = A            # m^3/(kg s)
        self.Ta = Ta          # K
        self.q = q            # J/kg
        self.kappa0 = kappa0  # kg/(m s K^n)   (K = kappa0*cp*T^n)
        self.gamma = gamma
        self.Mw = Mw          # kg/kmol
        self.n = n
        self.Pr = Pr
        self.Le = Le
        self.R = RU / Mw
        self.cp = gamma / (gamma - 1.0) * self.R
        self.cv = self.cp / gamma

    # CJ detonation speed (exact for calorically perfect, same-gamma gas)
    def D_CJ(self, T0):
        g = self.gamma
        H = (g * g - 1.0) * self.q / (2.0 * g * self.R * T0)
        M = math.sqrt(1.0 + H) + math.sqrt(H)
        return M * math.sqrt(g * self.R * T0)

    # Post-shock (von Neumann) state for a shock of speed D into (p0, T0)
    def von_neumann(self, D, p0, T0):
        g = self.gamma
        c0 = math.sqrt(g * self.R * T0)
        M = D / c0
        rho0 = p0 / (self.R * T0)
        pr = 1.0 + 2.0 * g / (g + 1.0) * (M * M - 1.0)
        rr = (g + 1.0) * M * M / ((g - 1.0) * M * M + 2.0)
        p = p0 * pr
        rho = rho0 * rr
        T = p / (rho * self.R)
        u = D / rr  # speed relative to the shock
        return p, rho, T, u

    # State on the Rayleigh line / Hugoniot for progress lam = 1 - Y
    # (steady shock frame, gas entering at speed D)
    def state(self, lam, D, p0, T0):
        g = self.gamma
        rho0 = p0 / (self.R * T0)
        m = rho0 * D
        J = p0 + rho0 * D * D
        H = g / (g - 1.0) * p0 / rho0 + 0.5 * D * D + self.q * lam
        a = 0.5 - g / (g - 1.0)
        b = g / (g - 1.0) * J / m
        c = -H
        disc = max(b * b - 4.0 * a * c, 0.0)
        # a < 0 : the subsonic (shocked) branch is the smaller root
        u1 = (-b + math.sqrt(disc)) / (2.0 * a)
        u2 = (-b - math.sqrt(disc)) / (2.0 * a)
        u = min(u1, u2)
        rho = m / u
        p = J - m * u
        T = p / (rho * self.R)
        return p, rho, T, u

    def transport(self, T):
        k_over_cp = self.kappa0 * T ** self.n
        mu = self.Pr * k_over_cp
        K = k_over_cp * self.cp
        rhoD = k_over_cp / self.Le
        return mu, K, rhoD


# ---------------------------------------------------------------------------
# ZND integration
# ---------------------------------------------------------------------------

def znd(gas, p0, T0, D=None, lam_end=0.9999):
    DCJ = gas.D_CJ(T0)
    if D is None:
        D = DCJ * (1.0 + 1e-10)

    def rhs(x, y):
        lam = min(max(y[0], 0.0), 1.0)
        p, rho, T, u = gas.state(lam, D, p0, T0)
        return [gas.A * rho * (1.0 - lam) * math.exp(-gas.Ta / T) / u]

    def ev_half(x, y):
        return y[0] - 0.5
    ev_half.terminal = False

    def ev_end(x, y):
        return y[0] - lam_end
    ev_end.terminal = True
    ev_end.direction = 1

    # characteristic scale for the first step
    pvn, rvn, Tvn, uvn = gas.von_neumann(D, p0, T0)
    k0 = gas.A * rvn * math.exp(-gas.Ta / Tvn)
    L0 = uvn / k0

    sol = solve_ivp(
        rhs, [0.0, 1e4 * L0], [0.0],
        method="LSODA", events=[ev_half, ev_end],
        rtol=1e-10, atol=1e-13, dense_output=True,
        max_step=L0 / 20.0
    )
    if len(sol.t_events[0]) == 0:
        raise RuntimeError("half-reaction point not reached")

    xd = sol.t_events[0][0]
    x = sol.t
    lam = np.clip(sol.y[0], 0.0, 1.0)
    dlam = np.array([rhs(xi, [li])[0] for xi, li in zip(x, lam)])
    imax = int(np.argmax(dlam))
    Lind = x[imax]
    Lhr = 1.0 / dlam[imax]

    return dict(D=D, DCJ=DCJ, x=x, lam=lam, xd=xd, Lind=Lind, Lhr=Lhr,
                pvn=pvn, rhovn=rvn, Tvn=Tvn, uvn=uvn)


def lambda_over_xd(E):
    """Eq. (11) of Lu, Kaplan & Oran (2021)."""
    return (-5.34e-4 * E**4 + 3.63e-2 * E**3 - 9.41e-1 * E**2
            + 8.96 * E + 2.81)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="CJ/ZND calculator and unit converter for the CDM")
    ap.add_argument("--paper", action="store_true",
                    help="inputs are in the units of the Lu et al. tables: "
                         "A [cm^3/(g s)], EaByRT0, qMwByRT0, "
                         "kappa0 [g/(s cm K^0.7)]")
    ap.add_argument("--A", type=float, required=True)
    ap.add_argument("--Ta", type=float, help="activation temperature [K] (SI)")
    ap.add_argument("--EaByRT0", type=float, help="Ea/(R T0) (paper)")
    ap.add_argument("--q", type=float, help="heat release [J/kg] (SI)")
    ap.add_argument("--qMwByRT0", type=float, help="q Mw/(R T0) (paper)")
    ap.add_argument("--kappa0", type=float, required=True)
    ap.add_argument("--gamma", type=float, required=True)
    ap.add_argument("--Mw", type=float, required=True, help="kg/kmol (=g/mol)")
    ap.add_argument("--n", type=float, default=0.7)
    ap.add_argument("--Pr", type=float, default=1.0)
    ap.add_argument("--Le", type=float, default=1.0)
    ap.add_argument("--T0", type=float, default=298.0)
    ap.add_argument("--p0", type=float, default=101325.0)
    ap.add_argument("--nCellsPerXd", type=float, default=20.0)
    ap.add_argument("--write", type=str, default=None,
                    help="write ZND profile table (x_behind_shock p T u Y rho)")
    a = ap.parse_args()

    R = RU / a.Mw
    if a.paper:
        A = a.A * 1e-3
        kappa0 = a.kappa0 * 0.1
        if a.EaByRT0 is None or a.qMwByRT0 is None:
            sys.exit("--paper needs --EaByRT0 and --qMwByRT0")
        Ta = a.EaByRT0 * a.T0
        q = a.qMwByRT0 * R * a.T0
    else:
        A = a.A
        kappa0 = a.kappa0
        Ta = a.Ta if a.Ta is not None else a.EaByRT0 * a.T0
        q = a.q if a.q is not None else a.qMwByRT0 * R * a.T0

    gas = CDMGas(A, Ta, q, kappa0, a.gamma, a.Mw, a.n, a.Pr, a.Le)
    z = znd(gas, a.p0, a.T0)
    E = Ta / z["Tvn"]
    rho0 = a.p0 / (R * a.T0)

    print("\n=== CDM parameters (SI) ===")
    print(f"  A       = {A:.6e}  m^3/(kg s)")
    print(f"  Ta      = {Ta:.6g}  K          (Ea/R)")
    print(f"  q       = {q:.6e}  J/kg")
    print(f"  kappa0  = {kappa0:.6e}  kg/(m s K^{a.n})")
    print(f"  gamma   = {a.gamma},  Mw = {a.Mw} kg/kmol,  R = {R:.4f} J/(kg K)")
    print(f"  Cp      = {gas.cp:.6g} J/(kg K),  Cv = {gas.cv:.6g} J/(kg K)")
    print(f"  Pr = {a.Pr}, Le = {a.Le}")

    print("\n=== Initial state ===")
    print(f"  p0 = {a.p0} Pa, T0 = {a.T0} K, rho0 = {rho0:.6g} kg/m^3, "
          f"c0 = {math.sqrt(a.gamma*R*a.T0):.2f} m/s")

    print("\n=== CJ / von Neumann ===")
    print(f"  D_CJ    = {z['DCJ']:.2f} m/s")
    print(f"  p_vn    = {z['pvn']:.6g} Pa  ({z['pvn']/a.p0:.3f} p0)")
    print(f"  T_vn    = {z['Tvn']:.2f} K")
    print(f"  rho_vn  = {z['rhovn']:.6g} kg/m^3")
    pCJ, rCJ, TCJ, uCJ = gas.state(1.0, z["D"], a.p0, a.T0)
    print(f"  p_CJ    = {pCJ:.6g} Pa  ({pCJ/a.p0:.3f} p0),  T_CJ = {TCJ:.1f} K")
    print(f"  T_cv    = {a.T0 + q/gas.cv:.1f} K (const.-volume, complete reaction)")
    print(f"  T_b     = {a.T0 + q/gas.cp:.1f} K (const.-pressure, complete reaction)")

    print("\n=== ZND structure ===")
    print(f"  x_d     = {z['xd']*100:.5f} cm   (half-reaction distance)")
    print(f"  L_ind   = {z['Lind']*100:.5f} cm   (to max heat release)")
    print(f"  L_hr    = {z['Lhr']*100:.5f} cm   (1/max dlam/dx)")
    print(f"  E       = Ea/(R T_vn) = {E:.3f}")

    ratio = lambda_over_xd(E)
    print("\n=== Cell size estimate, Eq. (11) of Lu et al. 2021 ===")
    if 2.0 <= E <= 23.0:
        print(f"  lambda/x_d = {ratio:.2f}  ->  lambda ~ {ratio*z['xd']*100:.3f} cm")
    else:
        print(f"  E = {E:.2f} is outside the fitted range (2..~23)")

    dx = z["xd"] / a.nCellsPerXd
    print("\n=== Resolution (Choi et al.: xd/20, Lhr/5, Lind/8) ===")
    print(f"  dx <= xd/20   = {z['xd']/20*1e6:.3f} um")
    print(f"  dx <= Lhr/5   = {z['Lhr']/5*1e6:.3f} um")
    print(f"  dx <= Lind/8  = {z['Lind']/8*1e6:.3f} um")
    print(f"  chosen dx     = {dx*1e6:.3f} um  (xd/{a.nCellsPerXd:g})")

    print("\n=== OpenFOAM: constant/physicalProperties (mixture) ===")
    print(f"""    specie          {{ molWeight {a.Mw}; }}
    thermodynamics  {{ Cv {gas.cv:.8g}; hf 0; }}
    transport       {{ mu 0; Pr {a.Pr}; }}   // not used by CDMGpuFoam""")
    print("\n=== OpenFOAM: constant/CDMProperties ===")
    print(f"""    A        {A:.8g};
    Ta       {Ta:.8g};
    q        {q:.8g};
    kappa0   {kappa0:.8g};
    n        {a.n};
    Pr       {a.Pr};
    Le       {a.Le};""")

    if a.write:
        with open(a.write, "w") as f:
            f.write("# ZND profile, CDM, D = %.6f m/s (shock frame)\n" % z["D"])
            f.write("# p0 %.6g T0 %.6g D %.6f xd %.6e\n"
                    % (a.p0, a.T0, z["D"], z["xd"]))
            f.write("# columns: distance_behind_shock[m] p[Pa] T[K] "
                    "u_shockFrame[m/s] Y[-] rho[kg/m3]\n")
            # add the pre-shock point at negative distance for interpolation
            for xi, li in zip(z["x"], z["lam"]):
                p, rho, T, u = gas.state(li, z["D"], a.p0, a.T0)
                f.write(f"{xi:.10e} {p:.10e} {T:.10e} {u:.10e} "
                        f"{1.0-li:.10e} {rho:.10e}\n")
        print(f"\nWrote ZND profile ({len(z['x'])} points) to {a.write}")


if __name__ == "__main__":
    main()
