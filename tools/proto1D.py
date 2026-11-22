#!/usr/bin/env python3
"""
proto1D.py -- 1D numpy prototype of the CDMGpuFoam algorithm, used to
verify the numerical design (not the C++ code):

  MUSCL/vanLeer reconstruction of rho, U, e, p/rho, c, Y
  HLLC-P fluxes (same formulas as src/fluxSchemes/HLLCP)
  CDM transport (Pr = Le = 1), explicit
  SSP-RK2, Strang-split exponential sub-stepped CDM chemistry

Lab frame: ZND profile propagating into quiescent gas; measures D.
"""
import sys, math
import numpy as np
sys.path.insert(0, ".")
from CDM_ZND import CDMGas, znd

T0, p0 = 298.0, 101325.0
# Table 1 of Lu et al. 2021, converted to SI
Mw, gam = 24.36, 1.177
R = 8314.462618 / Mw
gas = CDMGas(A=1.064e10, Ta=54.95 * T0, q=46.36 * R * T0,
             kappa0=2.367e-6, gamma=gam, Mw=Mw)
if "--Evn" in sys.argv:   # stable test gas: set Ta/T_vn, keep xd similar
    Evn = float(sys.argv[sys.argv.index("--Evn") + 1])
    Tvn0 = gas.von_neumann(gas.D_CJ(T0), p0, T0)[2]
    gas.A *= math.exp(Evn - gas.Ta / Tvn0); gas.Ta = Evn * Tvn0
cv = gas.cv
z = znd(gas, p0, T0)
xd, D = z["xd"], z["D"]

viscous = "--inviscid" not in sys.argv
nPerXd = float(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 20
dx = xd / nPerXd
TRAVEL = float(sys.argv[sys.argv.index("--travel")+1]) if "--travel" in sys.argv else 45
L = (TRAVEL + 40) * xd
N = int(L / dx)
x = (np.arange(N) + 0.5) * dx
xs0 = 20 * xd

# ---- initial ZND profile, lab frame
prof_s, prof = z["x"], []
for li in z["lam"]:
    prof.append(gas.state(li, D, p0, T0))
prof = np.array(prof)  # p rho T u
dist = xs0 - x
rho = np.empty(N); u = np.empty(N); p = np.empty(N); Y = np.empty(N)
ahead = dist < 0
rho[ahead], u[ahead], p[ahead], Y[ahead] = p0 / (R * T0), 0, p0, 1
d = np.clip(dist[~ahead], 0, prof_s[-1])
lam = np.interp(d, prof_s, z["lam"])
for k, col in enumerate([0, 1, 3]):
    pass
pp = np.interp(d, prof_s, prof[:, 0]); rr = np.interp(d, prof_s, prof[:, 1])
uu = np.interp(d, prof_s, prof[:, 3])
p[~ahead], rho[~ahead], u[~ahead], Y[~ahead] = pp, rr, D - uu, 1 - lam

e_of_T = lambda T: cv * (T - T0)
def prim(q):
    r, m, E, rY = q
    u = m / r
    e = E / r - 0.5 * u * u
    T = T0 + e / cv
    Y = np.clip(rY / r, 0, 1)
    p = r * R * T
    return r, u, e, T, p, Y

T = p / (rho * R)
q = np.array([rho, rho * u, rho * (e_of_T(T) + 0.5 * u * u), rho * Y])

def vanleer(a):
    # ghost-padded array a (N+2); returns left/right face states at N+1 faces
    dm = a[1:-1] - a[:-2]; dp = a[2:] - a[1:-1]
    r = dm * dp
    slope = np.where(r > 0, 2 * r / (dm + dp + 1e-300), 0.0)
    aL = a[1:-1] + 0.5 * slope       # right face of each interior cell
    aR = a[1:-1] - 0.5 * slope       # left face
    # faces 0..N: left state = aL of cell i-1, right state = aR of cell i
    left = np.concatenate([[a[0]], aL])
    right = np.concatenate([aR, [a[-1]]])
    return left, right

def pad(a):  # transmissive boundaries
    return np.concatenate([[a[0]], a, [a[-1]]])

def rhs(q):
    r, u, e, T, p, Y = prim(q)
    c = np.sqrt(gam * p / r)
    rP = p / r
    rL, rR = vanleer(pad(r)); uL, uR = vanleer(pad(u))
    pL, pR = vanleer(pad(p)); YL, YR = vanleer(pad(Y))
    eL, eR = pL / rL / (gam - 1) - cv * T0, pR / rR / (gam - 1) - cv * T0  # e = p/(rho(g-1)) + eOffset, as in the C++
    cL, cR = np.sqrt(gam * pL / rL), np.sqrt(gam * pR / rR)
    YL, YR = np.clip(YL, 0, 1), np.clip(YR, 0, 1)
    EL, ER = eL + 0.5 * uL**2, eR + 0.5 * uR**2
    # pressure sensor
    pp = pad(p)
    rf = np.minimum(pp[:-1] / pp[1:], pp[1:] / pp[:-1])  # N+1 faces
    fc = np.minimum(rf[:-1], rf[1:]) ** 3                 # N cells
    fcp = pad(fc)
    f = 0.5 * (fcp[:-1] + fcp[1:])
    # HLLCP
    MaL, MaR = uL / cL, uR / cR
    wL = np.sqrt(rL) / (np.sqrt(rL) + np.sqrt(rR)); wR = 1 - wL
    MaT = wL * MaL + wR * MaR; uT = wL * uL + wR * uR; cT = wL * cL + wR * cR
    SL = np.minimum(np.minimum(uL - cL, uR - cR), uT - cT)
    SR = np.maximum(np.maximum(uR + cR, uL + cL), uT + cT)
    aL = rL * (SL - uL); aR = rR * (SR - uR)
    SS = (aR * uR - aL * uL + pL - pR) / (aR - aL)
    pS = (aR * pL - aL * pR - aL * aR * (uL - uR)) / (aR - aL)
    th = np.minimum(np.maximum(abs(MaL), abs(MaR)), 1)
    pSS = th * pS + (1 - th) * 0.5 * (pL + pR)
    pSSS = f * pSS + (1 - f) * pS
    phip = (f - 1) * SL * SR / (SR - SL) / (1 + abs(MaT)) * (pR - pL) / cT**2
    left = SS >= 0
    rK = np.where(left, rL, rR); uK = np.where(left, uL, uR)
    EK = np.where(left, EL, ER); pK = np.where(left, pL, pR)
    SK = np.where(left, SL, SR); YK = np.where(left, YL, YR)
    dS = SK - SS
    Fm = SS * rK * (SK - uK) / dS + phip
    Fu = (SS * (SK * rK * uK - (rK * uK * uK + pK)) + SK * pSSS) / dS + phip * uT
    FE = SS * (SK * rK * EK - (rK * EK + pK) * uK + SK * pS) / dS + 0.5 * phip * uT**2
    FY = Fm * np.where(Fm >= 0, YL, YR)
    supL, supR = SL >= 0, SR <= 0
    Fm = np.where(supL, rL * uL, np.where(supR, rR * uR, Fm))
    Fu = np.where(supL, rL * uL * uL + pL, np.where(supR, rR * uR * uR + pR, Fu))
    FE = np.where(supL, (rL * EL + pL) * uL, np.where(supR, (rR * ER + pR) * uR, FE))
    FY = np.where(supL, rL * uL * YL, np.where(supR, rR * uR * YR, FY))
    amax = np.max(np.maximum(abs(uL) + cL, abs(uR) + cR))
    if viscous:
        k = gas.kappa0 * T**0.7
        kf = pad(k); kf = 0.5 * (kf[:-1] + kf[1:])
        du = np.diff(pad(u)) / dx; dT = np.diff(pad(T)) / dx; dY = np.diff(pad(Y)) / dx
        uf = pad(u); uf = 0.5 * (uf[:-1] + uf[1:])
        tau = kf * (4.0 / 3.0) * du          # 1D: mu (du/dx + du/dx - 2/3 du/dx)
        Fu = Fu - tau
        FE = FE - tau * uf - kf * gas.cp * dT
        FY = FY - kf * dY
    F = np.array([Fm, Fu, FE, FY])
    return -(F[:, 1:] - F[:, :-1]) / dx, amax

def react(q, h):
    r, u, e, T, p, Y = prim(q)
    Y0 = Y.copy(); t = np.zeros_like(Y); Tl = T.copy(); Yl = Y.copy()
    for _ in range(1000):
        k = gas.A * r * np.exp(-gas.Ta / Tl)
        active = (t < h) & (Yl > 1e-15) & (k * (h - t) > 1e-15)
        if not active.any():
            break
        hs = np.where(active, np.minimum(h - t, 0.2 / k), 0)
        Yn = Yl * np.exp(-k * hs)
        Tl = Tl + gas.q * (Yl - Yn) / cv
        Yl, t = Yn, t + hs
    q = q.copy()
    q[2] += r * gas.q * (Y0 - Yl); q[3] = r * Yl
    return q

def xshock(q):
    r, u, e, T, p, Y = prim(q)
    return x[np.where(p > 5 * p0)[0].max()]

maxCo = 0.4
t, tend = 0.0, TRAVEL * xd / D
hist = []
nstep = 0
while t < tend:
    _, amax = rhs(q)
    dtc = maxCo * dx / amax
    if viscous:
        r_, u_, e_, T_, p_, Y_ = prim(q)
        nu = (gas.kappa0 * T_**0.7 / r_).max()
        dtc = min(dtc, maxCo * dx * dx / (2 * nu))
    dt = min(dtc, tend - t)
    q = react(q, 0.5 * dt)
    q0 = q.copy()
    L1, _ = rhs(q); q1 = q0 + dt * L1
    L2, _ = rhs(q1); q = 0.5 * q0 + 0.5 * (q1 + dt * L2)
    q = react(q, 0.5 * dt)
    t += dt; nstep += 1
    hist.append((t, xshock(q)))
    if np.any(prim(q)[0] <= 0) or not np.all(np.isfinite(q)):
        print("FAILED at step", nstep); sys.exit(1)

hist = np.array(hist)
n2 = len(hist) // 2
Dnum = np.polyfit(hist[n2:, 0], hist[n2:, 1], 1)[0]
r, u, e, T, p, Y = prim(q)
print(f"viscous={viscous} dx=xd/{nPerXd:g} cells={N} steps={nstep}")
print(f"D_CJ theory = {D:.1f} m/s,  D numerical (2nd half) = {Dnum:.1f} m/s,"
      f"  error = {100*(Dnum-D)/D:+.2f} %")
print(f"max p = {p.max()/p0:.2f} p0 (p_vn = {z['pvn']/p0:.2f} p0),"
      f"  min Y = {Y.min():.3e}, max T = {T.max():.0f} K")
# local speed history: fit over windows of ~5 xd of travel
w = max(int(len(hist) * 5 / TRAVEL), 5)
Dloc = [np.polyfit(hist[i:i+w,0], hist[i:i+w,1], 1)[0] for i in range(0, len(hist)-w, w)]
print("windowed D/D_CJ:", " ".join(f"{d/D:.2f}" for d in Dloc))
