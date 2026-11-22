/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: GPU back end
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

Description
    The CDM right-hand side on OpenFOAM's unstructured addressing, written as
    C++17 parallel algorithms so that one source serves every target:

      nvc++ -stdpar=gpu        -> CUDA kernels, managed memory
      nvc++ -stdpar=multicore  -> every CPU core
      g++   -ltbb -DCDM_STDPAR -> every CPU core, no NVIDIA SDK
      g++                      -> plain loops, for debugging and for machines
                                  with no GPU at all

    That portability is the point: the same libCDMGpu.so has to be as happy on
    a cluster card as on a workstation, and has to still build where there is
    no card.

    Gradients and divergences are gathers over a CSR cell-to-face map, not
    scatters over a face loop.  A scatter would need atomics, which are
    slower, and which stdpar does not express portably.

    This file includes no OpenFOAM header and no CUDA header.  Everything it
    knows about the mesh arrived through cdmGpuSetMesh.

\*---------------------------------------------------------------------------*/

#include "cdmGpu.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <vector>

#ifdef CDM_STDPAR
#include <execution>
#endif

// * * * * * * * * * * * * * * * * parallel  * * * * * * * * * * * * * * * * //

namespace
{

struct CountingIterator
{
    using value_type = int;
    using difference_type = std::ptrdiff_t;
    using pointer = const int*;
    using reference = const int&;
    using iterator_category = std::random_access_iterator_tag;

    int i;

    CountingIterator() : i(0) {}
    explicit CountingIterator(int v) : i(v) {}

    int operator*() const { return i; }
    int operator[](difference_type n) const { return i + int(n); }

    CountingIterator& operator++() { ++i; return *this; }
    CountingIterator& operator--() { --i; return *this; }
    CountingIterator operator++(int) { CountingIterator t = *this; ++i; return t; }
    CountingIterator operator--(int) { CountingIterator t = *this; --i; return t; }
    CountingIterator& operator+=(difference_type n) { i += int(n); return *this; }
    CountingIterator& operator-=(difference_type n) { i -= int(n); return *this; }

    friend CountingIterator operator+(CountingIterator a, difference_type n) { a += n; return a; }
    friend CountingIterator operator+(difference_type n, CountingIterator a) { a += n; return a; }
    friend CountingIterator operator-(CountingIterator a, difference_type n) { a -= n; return a; }
    friend difference_type operator-(CountingIterator a, CountingIterator b)
    { return difference_type(a.i) - difference_type(b.i); }

    friend bool operator==(CountingIterator a, CountingIterator b) { return a.i == b.i; }
    friend bool operator!=(CountingIterator a, CountingIterator b) { return a.i != b.i; }
    friend bool operator< (CountingIterator a, CountingIterator b) { return a.i <  b.i; }
    friend bool operator> (CountingIterator a, CountingIterator b) { return a.i >  b.i; }
    friend bool operator<=(CountingIterator a, CountingIterator b) { return a.i <= b.i; }
    friend bool operator>=(CountingIterator a, CountingIterator b) { return a.i >= b.i; }
};

template<class F>
inline void pfor(int n, F f)
{
#ifdef CDM_STDPAR
    std::for_each(std::execution::par_unseq, CountingIterator(0), CountingIterator(n), f);
#else
    for (int k = 0; k < n; ++k) { f(k); }
#endif
}

template<class T, class Op, class F>
inline T preduce(int n, T init, Op op, F f)
{
#ifdef CDM_STDPAR
    return std::transform_reduce
    (
        std::execution::par_unseq, CountingIterator(0), CountingIterator(n), init, op, f
    );
#else
    T a = init;
    for (int k = 0; k < n; ++k) { a = op(a, f(k)); }
    return a;
#endif
}

struct MaxOp { double operator()(double a, double b) const { return a > b ? a : b; } };
struct SumOp { double operator()(double a, double b) const { return a + b; } };


// * * * * * * * * * * * * * * * * limiters  * * * * * * * * * * * * * * * * //

inline double sq_(double x) { return x*x; }

//- L(r) in OpenFOAM's blending convention: 0 is upwind, 1 is central.
//  The ROUND limiters are those of Deng (2023), as vendored in
//  src/fluxSchemes/ROUND, reproduced here in the same convention.
inline double limiterL(int which, double r)
{
    if (!(r > 0.0)) { return 0.0; }

    switch (which)
    {
        case CDM_GPU_UPWIND:  return 0.0;
        case CDM_GPU_VANLEER: return 2.0*r/(1.0 + r);
        case CDM_GPU_MINMOD:  return r < 1.0 ? r : 1.0;

        case CDM_GPU_ROUNDF:
            return 2.0*r*(1.0 + 2.0*r)/(1.0 + 2.0*r + 3.0*r*r);

        case CDM_GPU_ROUNDA:
        case CDM_GPU_ROUNDAPLUS:
        {
            const double b1 = (which == CDM_GPU_ROUNDA) ? 500.0  : 1100.0;
            const double b2 = (which == CDM_GPU_ROUNDA) ? 1000.0 : 800.0;
            const double lm = (which == CDM_GPU_ROUNDA) ? 0.25   : 0.15;

            if (r > 1.0)
            {
                const double A = sq_(sq_(sq_(1.0 + r)));
                const double C = b2*sq_(sq_(1.0 - r)) + 32.0*sq_(sq_(1.0 + r));
                const double D = b2*sq_(sq_(1.0 - r)) + 16.0*sq_(sq_(1.0 + r));
                const double v =
                    (6.0*b2*(1.0 - lm)*sq_(sq_(r - 1.0))*C + 256.0*A*(2.0 + r))/(3.0*D*D);
                const double cap = 2.0 - 2.0*lm;
                return v < cap ? v : cap;
            }
            const double A = sq_(sq_(sq_(1.0 + r)));
            const double B = b1*sq_(sq_(0.5 - 0.5*r)) + sq_(sq_(1.0 + r));
            const double v = (A*(2.0 - 5.0*r) + 6.0*r*B*B)/(3.0*B*B);
            const double cap = 2.0*r;
            return v < cap ? v : cap;
        }

        case CDM_GPU_ROUNDL:
        {
            const double a  = r/(1.0 + r);
            const double a8 = sq_(sq_(sq_(1.0 + 5.0*sq_(a - 1.0))));
            const double a4 = sq_(sq_(1.0 + 12.0*sq_(a)));
            double ma = 18000.0*sq_(sq_(0.97 - a))*(0.97 - a)*sq_(a - 0.55)*(a - 0.55);
            double mb = 1100.0*sq_(0.47 - a)*(0.47 - a)*sq_(a - 0.05)*(a - 0.05);
            if (ma < 0.0) ma = 0.0;
            if (mb < 0.0) mb = 0.0;
            const double ms = 1.0/3.0 + 5.0*a/6.0 + ma + mb;
            return 2.0*(-a + (0.5 + 0.5*a)/a8
                 + (1.0 - 1.0/a8)*(3.0*a/2.0/a4 + (1.0 - 1.0/a4)*ms))*(1.0 + r);
        }
    }
    return 0.0;
}

inline double sign_(double x) { return x >= 0.0 ? 1.0 : -1.0; }

//- NVDTVD::r for a scalar: 2 (d.gradc)/(phiN - phiP) - 1, with OpenFOAM's
//  guard against a vanishing face jump.
inline double rScalar(double gradcf, double gradf)
{
    if (std::fabs(gradcf) >= 1000.0*std::fabs(gradf))
    {
        return 2.0*1000.0*sign_(gradcf)*sign_(gradf) - 1.0;
    }
    return 2.0*(gradcf/gradf) - 1.0;
}

} // namespace


// * * * * * * * * * * * * * * * * the case  * * * * * * * * * * * * * * * * //

struct cdmGpuCase
{
    // mesh
    int nCells = 0, nIF = 0, nBF = 0, nBC = 0;
    std::vector<int> bCell;

    //- Halo ghost cells, one per processor-patch face.  They occupy cell
    //  indices [nCells, nCells + nHalo) and exist only to be read as the
    //  neighbour of a processor face: they carry a primitive state and a
    //  gradient, which is the entire stencil a face needs, and they are
    //  never written by any pass because every pass loops over nCells.
    int nHalo = 0;
    std::vector<int> hCell;

    //- Contiguous staging for the boundary-cell exchange, so the host never
    //  touches a strided subset of a device-resident field
    std::vector<double> bcsRho, bcsU, bcsP, bcsT, bcsY;

    //- The same, for the cells whose state and gradient this rank sends to
    //  its neighbours
    std::vector<double> hsRho, hsU, hsP, hsT, hsY;
    std::vector<double> hsGRho, hsGU, hsGP, hsGY;

    std::vector<int> own, nei, bOwn, cfOff, cfFace;
    std::vector<signed char> cfSign;
    std::vector<double> Sf, magSf, wgt, dc, dlt, V, bSf, bMagSf, bDc;

    // model
    double gamma = 1.4, Cv = 718.0, eOffset = 0.0;
    double A = 0.0, Ta = 0.0, q = 0.0;
    double kappa0 = 0.0, nExp = 0.7, Pr = 1.0, Le = 1.0;
    double subFrac = 0.2;
    int maxSub = 100, flux = CDM_GPU_HLLCLM, lim = CDM_GPU_VANLEER, viscous = 1;

    // conserved, and the Runge-Kutta stage origin
    std::vector<double> rho, rhoU, rhoE, rhoY;
    std::vector<double> rho0, rhoU0, rhoE0, rhoY0;

    // primitives
    std::vector<double> U, p, T, Y, c, Qdot;

    // cell gradients: grad(rho), grad(U) (tensor), grad(p), grad(Y)
    std::vector<double> gRho, gU, gP, gY;

    // boundary primitives, evaluated by OpenFOAM on the host
    std::vector<double> bRho, bU, bP, bY;

    // face fluxes, area-integrated, and the wave speed
    std::vector<double> Fr, FrU, FrE, FrY, amaxSf;
    std::vector<double> bFr, bFrU, bFrE, bFrY, bAmaxSf;

    // cell right-hand sides
    std::vector<double> dRho, dRhoU, dRhoE, dRhoY;

    // per-face fallback flag, summed on the host
    std::vector<double> fallback;

    void allocate()
    {
        // nT covers the halo ghosts as well. Only the fields a face reads for its
        // neighbour are extended, the primitives and the gradients; the conserved
        // state, the volumes and the right-hand sides are written by passes that
        // loop over real cells, so a ghost has no use for them.
        const int nT = nCells + nHalo;

        rho.assign(nT, 0.0);      rhoU.assign(3*nCells, 0.0);
        rhoE.assign(nCells, 0.0); rhoY.assign(nCells, 0.0);
        rho0 = rho; rhoU0 = rhoU; rhoE0 = rhoE; rhoY0 = rhoY;

        U.assign(3*nT, 0.0);      p.assign(nT, 0.0);
        T.assign(nT, 0.0);        Y.assign(nT, 0.0);
        c.assign(nT, 0.0);        Qdot.assign(nCells, 0.0);

        gRho.assign(3*nT, 0.0); gU.assign(9*nT, 0.0);
        gP.assign(3*nT, 0.0);   gY.assign(3*nT, 0.0);

        bRho.assign(nBF, 0.0); bU.assign(3*nBF, 0.0);
        bP.assign(nBF, 0.0);   bY.assign(nBF, 0.0);

        Fr.assign(nIF, 0.0);  FrU.assign(3*nIF, 0.0);
        FrE.assign(nIF, 0.0); FrY.assign(nIF, 0.0);
        amaxSf.assign(nIF, 0.0); fallback.assign(nIF, 0.0);

        bFr.assign(nBF, 0.0);  bFrU.assign(3*nBF, 0.0);
        bFrE.assign(nBF, 0.0); bFrY.assign(nBF, 0.0);
        bAmaxSf.assign(nBF, 0.0);

        dRho.assign(nCells, 0.0);  dRhoU.assign(3*nCells, 0.0);
        dRhoE.assign(nCells, 0.0); dRhoY.assign(nCells, 0.0);
    }
};


// * * * * * * * * * * * * * * * * kernels * * * * * * * * * * * * * * * * * //

namespace
{

//- One face's inviscid flux, in the face-normal frame but returned as full
//  vectors.  Transcribed from src/fluxSchemes/HLLC and HLLCLM so that the two
//  back ends compute the same number for the same face states.
struct FaceFlux { double r, mU[3], E, Y, amax; };

inline FaceFlux riemann
(
    int scheme,
    double rl, const double* Ul, double pl, double Yl,
    double rr, const double* Ur, double pr, double Yr,
    const double* n, double gamma, double eOffset
)
{
    const double gm1 = gamma - 1.0;

    const double ul = Ul[0]*n[0] + Ul[1]*n[1] + Ul[2]*n[2];
    const double ur = Ur[0]*n[0] + Ur[1]*n[1] + Ur[2]*n[2];

    // The face internal energy carries OpenFOAM's energy datum, exactly as
    // fluxScheme.C's state(): e = p/(rho (gamma - 1)) + eOffset. Leaving it
    // out makes the energy flux wrong by rho u eOffset on every face, a 5%
    // error in the detonation speed for a datum of -575 kJ/kg.
    const double El = pl/(rl*gm1) + eOffset
                    + 0.5*(Ul[0]*Ul[0] + Ul[1]*Ul[1] + Ul[2]*Ul[2]);
    const double Er = pr/(rr*gm1) + eOffset
                    + 0.5*(Ur[0]*Ur[0] + Ur[1]*Ur[1] + Ur[2]*Ur[2]);

    const double cl = std::sqrt(gamma*pl/rl);
    const double cr = std::sqrt(gamma*pr/rr);

    FaceFlux F;
    F.amax = std::max(std::fabs(ul) + cl, std::fabs(ur) + cr);

    double SL, SR;

    if (scheme == CDM_GPU_HLLCLM)
    {
        const double sl = std::sqrt(rl), sr = std::sqrt(rr);
        const double rs = 1.0/(sl + sr);
        const double uHat = (ul*sl + ur*sr)*rs;
        const double cHat = std::sqrt
        (
            (cl*cl*sl + cr*cr*sr)*rs + 0.5*sl*sr*rs*rs*(ur - ul)*(ur - ul)
        );
        SL = std::min(ul - cl, uHat - cHat);
        SR = std::max(ur + cr, uHat + cHat);
    }
    else
    {
        const double sl = std::sqrt(rl), sr = std::sqrt(rr);
        const double wl = sl/(sl + sr), wr = 1.0 - wl;
        const double uT = wl*ul + wr*ur;
        const double cT = wl*cl + wr*cr;
        SL = std::min(std::min(ul - cl, ur - cr), uT - cT);
        SR = std::max(std::max(ur + cr, ul + cl), uT + cT);
    }

    // supersonic either way
    if (SL >= 0.0 || SR <= 0.0)
    {
        const bool L = (SL >= 0.0);
        const double rk = L ? rl : rr;
        const double uk = L ? ul : ur;
        const double pk = L ? pl : pr;
        const double Yk = L ? Yl : Yr;
        const double Ek = L ? El : Er;
        const double* Uk = L ? Ul : Ur;

        F.r = rk*uk;
        for (int d = 0; d < 3; ++d) { F.mU[d] = rk*Uk[d]*uk + pk*n[d]; }
        F.E = (rk*Ek + pk)*uk;
        F.Y = F.r*Yk;
        return F;
    }

    const double aL = rl*(SL - ul);
    const double aR = rr*(SR - ur);

    if (scheme == CDM_GPU_HLLCLM)
    {
        // Central form, Eqs. (19), (23)-(25) of Fleischmann et al.
        constexpr double MaLimit = 0.1;
        constexpr double piByTwo = 1.5707963267948966;

        const double SS = (pr - pl + aL*ul - aR*ur)/(aL - aR);
        const double kL = aL/(SL - SS), kR = aR/(SR - SS);

        const double ESL = kL*(El + (SS - ul)*(SS + pl/aL));
        const double ESR = kR*(Er + (SS - ur)*(SS + pr/aR));

        const double Ma = std::max(std::fabs(ul)/cl, std::fabs(ur)/cr);
        const double phi = std::sin(std::min(1.0, Ma/MaLimit)*piByTwo);
        const double sL = phi*SL, sR = phi*SR, sS = std::fabs(SS);

        const double FrL = rl*ul, FrR = rr*ur;
        F.r = 0.5*(FrL + FrR)
            + 0.5*(sL*(kL - rl) + sS*(kL - kR) + sR*(kR - rr));

        for (int d = 0; d < 3; ++d)
        {
            const double mSL = kL*(Ul[d] + (SS - ul)*n[d]);
            const double mSR = kR*(Ur[d] + (SS - ur)*n[d]);
            const double FL = rl*Ul[d]*ul + pl*n[d];
            const double FR = rr*Ur[d]*ur + pr*n[d];
            F.mU[d] = 0.5*(FL + FR)
                    + 0.5*(sL*(mSL - rl*Ul[d]) + sS*(mSL - mSR) + sR*(mSR - rr*Ur[d]));
        }

        const double FEL = (rl*El + pl)*ul, FER = (rr*Er + pr)*ur;
        F.E = 0.5*(FEL + FER)
            + 0.5*(sL*(ESL - rl*El) + sS*(ESL - ESR) + sR*(ESR - rr*Er));

        const double FYL = FrL*Yl, FYR = FrR*Yr;
        F.Y = 0.5*(FYL + FYR)
            + 0.5*(sL*(kL*Yl - rl*Yl) + sS*(kL*Yl - kR*Yr) + sR*(kR*Yr - rr*Yr));

        return F;
    }

    // HLLC: star-region flux on the side of the contact
    const double SS = (aR*ur - aL*ul + pl - pr)/(aR - aL);
    const double pS = (aR*pl - aL*pr - aL*aR*(ul - ur))/(aR - aL);

    const bool L = (SS >= 0.0);
    const double rk = L ? rl : rr;
    const double uk = L ? ul : ur;
    const double pk = L ? pl : pr;
    const double Yk = L ? Yl : Yr;
    const double Ek = L ? El : Er;
    const double Sk = L ? SL : SR;
    const double* Uk = L ? Ul : Ur;
    const double dS = Sk - SS;

    F.r = SS*rk*(Sk - uk)/dS;
    for (int d = 0; d < 3; ++d)
    {
        const double FmU = rk*Uk[d]*uk + pk*n[d];
        F.mU[d] = (SS*(Sk*rk*Uk[d] - FmU) + Sk*pS*n[d])/dS;
    }
    F.E = SS*(Sk*rk*Ek - (rk*Ek + pk)*uk + Sk*pS)/dS;
    F.Y = F.r*Yk;

    return F;
}

} // namespace


// * * * * * * * * * * * * * * * * C API * * * * * * * * * * * * * * * * * * //

extern "C" {

const char* cdmGpuBackend(void)
{
#if defined(CDM_STDPAR) && defined(__NVCOMPILER)
    return "stdpar (nvc++)";
#elif defined(CDM_STDPAR)
    return "stdpar (host threads)";
#else
    return "serial";
#endif
}

double cdmGpuWallTime(void)
{
    using clk = std::chrono::steady_clock;
    static const clk::time_point t0 = clk::now();
    return std::chrono::duration<double>(clk::now() - t0).count();
}

int cdmGpuAvailable(void)
{
#if defined(CDM_STDPAR) && defined(__NVCOMPILER)
    return 1;
#else
    return 0;
#endif
}

cdmGpuCase* cdmGpuCreate(void) { return new cdmGpuCase(); }
void cdmGpuDestroy(cdmGpuCase* c) { delete c; }


void cdmGpuSetMesh
(
    cdmGpuCase* c,
    int nCells, int nHalo, int nIF, int nBF,
    const int* owner, const int* neighbour,
    const double* Sf, const double* magSf, const double* weights,
    const double* deltaCoeffs, const double* delta, const double* V,
    const int* bOwner, const double* bSf, const double* bMagSf,
    const double* bDeltaCoeffs,
    const int* cellFaceOffset, const int* cellFace, const signed char* cellFaceSign
)
{
    c->nCells = nCells; c->nHalo = nHalo; c->nIF = nIF; c->nBF = nBF;

    c->own.assign(owner, owner + nIF);
    c->nei.assign(neighbour, neighbour + nIF);
    c->Sf.assign(Sf, Sf + 3*nIF);
    c->magSf.assign(magSf, magSf + nIF);
    c->wgt.assign(weights, weights + nIF);
    c->dc.assign(deltaCoeffs, deltaCoeffs + nIF);
    c->dlt.assign(delta, delta + 3*nIF);
    c->V.assign(V, V + nCells);

    c->bOwn.assign(bOwner, bOwner + nBF);
    c->bSf.assign(bSf, bSf + 3*nBF);
    c->bMagSf.assign(bMagSf, bMagSf + nBF);
    c->bDc.assign(bDeltaCoeffs, bDeltaCoeffs + nBF);

    c->cfOff.assign(cellFaceOffset, cellFaceOffset + nCells + 1);
    const int nnz = cellFaceOffset[nCells];
    c->cfFace.assign(cellFace, cellFace + nnz);
    c->cfSign.assign(cellFaceSign, cellFaceSign + nnz);

    c->allocate();
}


void cdmGpuSetModel
(
    cdmGpuCase* c,
    double gamma, double Cv, double eOffset,
    double A, double Ta, double q,
    double kappa0, double n, double Pr, double Le,
    double maxSubstepFraction, int maxSubsteps,
    int fluxScheme, int limiter, int viscous
)
{
    c->gamma = gamma; c->Cv = Cv; c->eOffset = eOffset;
    c->A = A; c->Ta = Ta; c->q = q;
    c->kappa0 = kappa0; c->nExp = n; c->Pr = Pr; c->Le = Le;
    c->subFrac = maxSubstepFraction; c->maxSub = maxSubsteps;
    c->flux = fluxScheme; c->lim = limiter; c->viscous = viscous;
}


void cdmGpuUploadState
(
    cdmGpuCase* c, const double* rho, const double* rhoU,
    const double* rhoE, const double* rhoY
)
{
    std::copy(rho,  rho  + c->nCells,   c->rho.begin());
    std::copy(rhoU, rhoU + 3*c->nCells, c->rhoU.begin());
    std::copy(rhoE, rhoE + c->nCells,   c->rhoE.begin());
    std::copy(rhoY, rhoY + c->nCells,   c->rhoY.begin());
}


void cdmGpuDownloadState
(
    cdmGpuCase* c, double* rho, double* rhoU, double* rhoE, double* rhoY
)
{
    // Counts, never end(). rho is nCells + nHalo long because a ghost needs a
    // density, and the host field it is copied into is nCells long; writing
    // the halo into it overruns the heap by nHalo doubles, invisible in
    // serial where nHalo is zero and fatal in parallel somewhere else.
    const int n = c->nCells;

    std::copy(c->rho.begin(),  c->rho.begin()  +   n, rho);
    std::copy(c->rhoU.begin(), c->rhoU.begin() + 3*n, rhoU);
    std::copy(c->rhoE.begin(), c->rhoE.begin() +   n, rhoE);
    std::copy(c->rhoY.begin(), c->rhoY.begin() +   n, rhoY);
}


void cdmGpuDownloadPrimitives
(
    cdmGpuCase* c, double* U, double* p, double* T, double* Y, double* Qdot
)
{
    // Counts, never end(): U, p, T and Y all carry nHalo ghosts past the
    // real cells, and the host fields do not.  See cdmGpuDownloadState.
    const int n = c->nCells;

    std::copy(c->U.begin(), c->U.begin() + 3*n, U);
    std::copy(c->p.begin(), c->p.begin() +   n, p);
    std::copy(c->T.begin(), c->T.begin() +   n, T);
    std::copy(c->Y.begin(), c->Y.begin() +   n, Y);
    std::copy(c->Qdot.begin(), c->Qdot.begin() + n, Qdot);
}


void cdmGpuDownloadQdot(cdmGpuCase* c, double* Qdot)
{
    std::copy(c->Qdot.begin(), c->Qdot.end(), Qdot);
}


void cdmGpuSetBoundaryCells(cdmGpuCase* c, int n, const int* cells)
{
    c->nBC = n;
    c->bCell.assign(cells, cells + n);

    c->bcsRho.assign(n, 0.0);
    c->bcsU.assign(3*n, 0.0);
    c->bcsP.assign(n, 0.0);
    c->bcsT.assign(n, 0.0);
    c->bcsY.assign(n, 0.0);
}


void cdmGpuDownloadBoundaryCells
(
    cdmGpuCase* c,
    double* rho, double* U, double* p, double* T, double* Y
)
{
    // Gathered by the device into contiguous staging buffers, and only those
    // buffers are read by the host. A host loop reading c->rho[bCell[i]] is
    // short in elements and ruinous in pages: these arrays are managed memory,
    // so a host read migrates whole 4 kB pages, and the boundary cells are the
    // mesh perimeter. On a 750 x 500 mesh that is some 36 MB migrated per
    // exchange instead of 20 kB, faulted straight back on the next kernel.
    const int n = c->nBC;
    if (n == 0) { return; }

    const int* bc = c->bCell.data();
    const double *rhoD = c->rho.data(), *UD = c->U.data(), *pD = c->p.data();
    const double *TD = c->T.data(), *YD = c->Y.data();

    double *sR = c->bcsRho.data(), *sU = c->bcsU.data(), *sP = c->bcsP.data();
    double *sT = c->bcsT.data(), *sY = c->bcsY.data();

    pfor(n, [=](int i)
    {
        const int ci = bc[i];
        sR[i] = rhoD[ci];
        sP[i] = pD[ci];
        sT[i] = TD[ci];
        sY[i] = YD[ci];
        sU[3*i + 0] = UD[3*ci + 0];
        sU[3*i + 1] = UD[3*ci + 1];
        sU[3*i + 2] = UD[3*ci + 2];
    });

    std::copy(c->bcsRho.begin(), c->bcsRho.end(), rho);
    std::copy(c->bcsU.begin(),   c->bcsU.end(),   U);
    std::copy(c->bcsP.begin(),   c->bcsP.end(),   p);
    std::copy(c->bcsT.begin(),   c->bcsT.end(),   T);
    std::copy(c->bcsY.begin(),   c->bcsY.end(),   Y);
}


void cdmGpuUploadBoundary
(
    cdmGpuCase* c, const double* bRho, const double* bU, const double* bP, const double* bY
)
{
    std::copy(bRho, bRho + c->nBF,   c->bRho.begin());
    std::copy(bU,   bU   + 3*c->nBF, c->bU.begin());
    std::copy(bP,   bP   + c->nBF,   c->bP.begin());
    std::copy(bY,   bY   + c->nBF,   c->bY.begin());
}


// * * * * * * * * * * * * * * * * halo exchange * * * * * * * * * * * * * * //
//
// A processor face is an internal face whose neighbour lives on another
// rank. Giving the device a ghost cell for that neighbour makes it an
// ordinary internal face and no kernel needs to know the difference, the
// same move the cyclic seam uses and for the same reason: one physical
// face must produce one flux, or mass leaks wherever the solution is not
// smooth across it.
//
// The reconstruction reads exactly two things per side, the cell value and
// the cell gradient, so a ghost needs those two and nothing else. Both
// ranks then evaluate the same expression from the same inputs and get a
// bit-identical flux, and a second layer of halo is never required.

void cdmGpuSetHaloCells(cdmGpuCase* c, int n, const int* cells)
{
    c->hCell.assign(cells, cells + n);

    c->hsRho.assign(n, 0.0);    c->hsU.assign(3*n, 0.0);
    c->hsP.assign(n, 0.0);      c->hsT.assign(n, 0.0);
    c->hsY.assign(n, 0.0);
    c->hsGRho.assign(3*n, 0.0); c->hsGU.assign(9*n, 0.0);
    c->hsGP.assign(3*n, 0.0);   c->hsGY.assign(3*n, 0.0);
}


void cdmGpuDownloadHaloCells
(
    cdmGpuCase* c,
    double* rho, double* U, double* p, double* T, double* Y,
    double* gRho, double* gU, double* gP, double* gY
)
{
    const int n = static_cast<int>(c->hCell.size());
    if (n == 0) { return; }

    const int* hc = c->hCell.data();

    const double *rhoD = c->rho.data(), *UD = c->U.data(), *pD = c->p.data();
    const double *TD = c->T.data(), *YD = c->Y.data();
    const double *gRD = c->gRho.data(), *gUD = c->gU.data();
    const double *gPD = c->gP.data(), *gYD = c->gY.data();

    double *sR = c->hsRho.data(), *sU = c->hsU.data(), *sP = c->hsP.data();
    double *sT = c->hsT.data(), *sY = c->hsY.data();
    double *sGR = c->hsGRho.data(), *sGU = c->hsGU.data();
    double *sGP = c->hsGP.data(), *sGY = c->hsGY.data();

    pfor(n, [=](int i)
    {
        const int ci = hc[i];
        sR[i] = rhoD[ci];
        sP[i] = pD[ci];
        sT[i] = TD[ci];
        sY[i] = YD[ci];
        for (int d = 0; d < 3; ++d)
        {
            sU[3*i + d]  = UD[3*ci + d];
            sGR[3*i + d] = gRD[3*ci + d];
            sGP[3*i + d] = gPD[3*ci + d];
            sGY[3*i + d] = gYD[3*ci + d];
        }
        for (int d = 0; d < 9; ++d) { sGU[9*i + d] = gUD[9*ci + d]; }
    });

    std::copy(c->hsRho.begin(),  c->hsRho.end(),  rho);
    std::copy(c->hsU.begin(),    c->hsU.end(),    U);
    std::copy(c->hsP.begin(),    c->hsP.end(),    p);
    std::copy(c->hsT.begin(),    c->hsT.end(),    T);
    std::copy(c->hsY.begin(),    c->hsY.end(),    Y);
    std::copy(c->hsGRho.begin(), c->hsGRho.end(), gRho);
    std::copy(c->hsGU.begin(),   c->hsGU.end(),   gU);
    std::copy(c->hsGP.begin(),   c->hsGP.end(),   gP);
    std::copy(c->hsGY.begin(),   c->hsGY.end(),   gY);
}


void cdmGpuUploadHalo
(
    cdmGpuCase* c,
    const double* rho, const double* U, const double* p, const double* T,
    const double* Y,
    const double* gRho, const double* gU, const double* gP, const double* gY
)
{
    const int nC = c->nCells, nH = c->nHalo;
    if (nH == 0) { return; }

    std::copy(rho,  rho  +   nH, c->rho.begin()  +   nC);
    std::copy(U,    U    + 3*nH, c->U.begin()    + 3*nC);
    std::copy(p,    p    +   nH, c->p.begin()    +   nC);
    std::copy(T,    T    +   nH, c->T.begin()    +   nC);
    std::copy(Y,    Y    +   nH, c->Y.begin()    +   nC);
    std::copy(gRho, gRho + 3*nH, c->gRho.begin() + 3*nC);
    std::copy(gU,   gU   + 9*nH, c->gU.begin()   + 9*nC);
    std::copy(gP,   gP   + 3*nH, c->gP.begin()   + 3*nC);
    std::copy(gY,   gY   + 3*nH, c->gY.begin()   + 3*nC);

    // The sound speed of a ghost follows from its own state, so it is derived
    // rather than sent, which also guarantees the two ranks agree exactly.
    const double gam = c->gamma;
    double* cs = c->c.data();
    const double* rhoD = c->rho.data();
    const double* pD = c->p.data();

    pfor(nH, [=](int h)
    {
        const int i = nC + h;
        cs[i] = std::sqrt(gam*pD[i]/rhoD[i]);
    });
}


/*- The cell-gradient pass on its own.  A parallel run calls this, exchanges
    the gradients of the cells behind each processor patch, uploads them into
    the ghosts, and only then forms the fluxes.

    cdmGpuComputeRHS recomputes these gradients, which costs a second pass
    but writes only indices below nCells - so the ghost gradients uploaded in
    between survive it untouched, and no kernel had to change. */
void cdmGpuComputeGradients(cdmGpuCase* c)
{
    const int nC = c->nCells, nIF = c->nIF;

    const int *own = c->own.data(), *nei = c->nei.data();
    const int *cfOff = c->cfOff.data(), *cfFace = c->cfFace.data();
    const signed char* cfSign = c->cfSign.data();
    const double *Sf = c->Sf.data(), *wgt = c->wgt.data();
    const double *V = c->V.data(), *bSf = c->bSf.data();

    const double *rho = c->rho.data(), *U = c->U.data(), *p = c->p.data();
    const double *Y = c->Y.data();
    const double *bRho = c->bRho.data(), *bU = c->bU.data();
    const double *bP = c->bP.data(), *bY = c->bY.data();

    double *gRho = c->gRho.data(), *gU = c->gU.data();
    double *gP = c->gP.data(), *gY = c->gY.data();

    pfor(nC, [=](int i)
    {
        double ar[3] = {0,0,0}, ap[3] = {0,0,0}, ay[3] = {0,0,0};
        double au[9] = {0,0,0,0,0,0,0,0,0};

        for (int k = cfOff[i]; k < cfOff[i+1]; ++k)
        {
            const int f = cfFace[k];
            const double s = double(cfSign[k]);

            double rf, pf, yf, uf[3], S[3];

            if (f < nIF)
            {
                const int P = own[f], N = nei[f];
                const double w = wgt[f], w1 = 1.0 - w;
                rf = w*rho[P] + w1*rho[N];
                pf = w*p[P]   + w1*p[N];
                yf = w*Y[P]   + w1*Y[N];
                for (int d = 0; d < 3; ++d)
                {
                    uf[d] = w*U[3*P+d] + w1*U[3*N+d];
                    S[d] = Sf[3*f+d];
                }
            }
            else
            {
                const int b = f - nIF;
                rf = bRho[b]; pf = bP[b]; yf = bY[b];
                for (int d = 0; d < 3; ++d)
                {
                    uf[d] = bU[3*b+d];
                    S[d] = bSf[3*b+d];
                }
            }

            for (int d = 0; d < 3; ++d)
            {
                const double sS = s*S[d];
                ar[d] += sS*rf;
                ap[d] += sS*pf;
                ay[d] += sS*yf;
                for (int e = 0; e < 3; ++e) { au[3*d+e] += sS*uf[e]; }
            }
        }

        const double iV = 1.0/V[i];
        for (int d = 0; d < 3; ++d)
        {
            gRho[3*i+d] = ar[d]*iV;
            gP[3*i+d]   = ap[d]*iV;
            gY[3*i+d]   = ay[d]*iV;
            for (int e = 0; e < 3; ++e) { gU[9*i+3*d+e] = au[3*d+e]*iV; }
        }
    });
}


void cdmGpuDownloadMassFlux(cdmGpuCase* c, double* phiI, double* phiB)
{
    if (phiI) { std::copy(c->Fr.begin(), c->Fr.end(), phiI); }
    if (phiB) { std::copy(c->bFr.begin(), c->bFr.end(), phiB); }
}


void cdmGpuDownloadRHS
(
    cdmGpuCase* c, double* dRho, double* dRhoU, double* dRhoE, double* dRhoY
)
{
    std::copy(c->dRho.begin(),  c->dRho.end(),  dRho);
    std::copy(c->dRhoU.begin(), c->dRhoU.end(), dRhoU);
    std::copy(c->dRhoE.begin(), c->dRhoE.end(), dRhoE);
    std::copy(c->dRhoY.begin(), c->dRhoY.end(), dRhoY);
}


void cdmGpuStoreStage(cdmGpuCase* c)
{
    double *r = c->rho.data(), *ru = c->rhoU.data(),
           *rE = c->rhoE.data(), *rY = c->rhoY.data();
    double *r0 = c->rho0.data(), *ru0 = c->rhoU0.data(),
           *rE0 = c->rhoE0.data(), *rY0 = c->rhoY0.data();
    const int n = c->nCells;

    pfor(n, [=](int i)
    {
        r0[i] = r[i]; rE0[i] = rE[i]; rY0[i] = rY[i];
        ru0[3*i] = ru[3*i]; ru0[3*i+1] = ru[3*i+1]; ru0[3*i+2] = ru[3*i+2];
    });
}


int cdmGpuDecode(cdmGpuCase* c)
{
    const double *r = c->rho.data(), *ru = c->rhoU.data(),
                 *rE = c->rhoE.data(), *rY = c->rhoY.data();
    double *U = c->U.data(), *p = c->p.data(), *T = c->T.data(),
           *Y = c->Y.data(), *cc = c->c.data();
    const double gam = c->gamma, gm1 = c->gamma - 1.0, Cv = c->Cv, eOff = c->eOffset;
    const int n = c->nCells;

    const double bad = preduce(n, 0.0, MaxOp(), [=](int i) -> double
    {
        const double ri = r[i];
        if (!(ri > 0.0)) { return double(i + 1); }

        const double ux = ru[3*i]/ri, uy = ru[3*i+1]/ri, uz = ru[3*i+2]/ri;
        const double e = rE[i]/ri - 0.5*(ux*ux + uy*uy + uz*uz);
        const double Ti = (e - eOff)/Cv;
        const double pi = ri*gm1*(e - eOff);

        U[3*i] = ux; U[3*i+1] = uy; U[3*i+2] = uz;
        T[i] = Ti;
        p[i] = pi;
        const double Yi = rY[i]/ri;
        Y[i] = Yi < 0.0 ? 0.0 : (Yi > 1.0 ? 1.0 : Yi);
        cc[i] = (pi > 0.0) ? std::sqrt(gam*pi/ri) : 0.0;

        return (!(Ti > 0.0) || !(pi > 0.0)) ? double(i + 1) : 0.0;
    });

    return int(bad);
}


void cdmGpuReact(cdmGpuCase* c, double dt)
{
    const double *r = c->rho.data();
    const double *Tc = c->T.data();
    double *rE = c->rhoE.data(), *rY = c->rhoY.data(), *Qd = c->Qdot.data();

    const double A = c->A, Ta = c->Ta, q = c->q, Cv = c->Cv;
    const double sub = c->subFrac;
    const int maxSub = c->maxSub, n = c->nCells;
    const double tiny = 1.0e-300;

    pfor(n, [=](int i)
    {
        const double ri = r[i];
        double Y0 = rY[i]/ri;
        Y0 = Y0 < 0.0 ? 0.0 : (Y0 > 1.0 ? 1.0 : Y0);

        // Leave Qdot alone on the early-out. CDMModel::react accumulates into
        // heatReleased across the two Strang half-steps and advance() zeroes it
        // once per step, so overwriting would discard the first half's release.
        if (Y0 <= 1.0e-12) { return; }

        double Y = Y0, Ti = Tc[i], t = 0.0;
        int ns = 0;

        while (t < dt && Y > 1.0e-12)
        {
            const double k = A*ri*std::exp(-Ta/Ti);
            if (k*(dt - t) < tiny) { break; }

            const double h = std::min(dt - t, sub/k);
            const double Yn = Y*std::exp(-k*h);
            Ti += q*(Y - Yn)/Cv;
            Y = Yn;
            t += h;

            if (++ns >= maxSub)
            {
                Y *= std::exp(-A*ri*std::exp(-Ta/Ti)*(dt - t));
                break;
            }
        }

        const double dQ = ri*q*(Y0 - Y);
        rE[i] += dQ;
        rY[i] = ri*Y;

        // Accumulated heat release, not a rate.  cdmGpuScaleQdot turns it
        // into a rate at the end of the step, as advance() does on the host.
        Qd[i] += dQ;
    });
}


void cdmGpuZeroQdot(cdmGpuCase* c)
{
    double* Qd = c->Qdot.data();
    pfor(c->nCells, [=](int i) { Qd[i] = 0.0; });
}


void cdmGpuScaleQdot(cdmGpuCase* c, double factor)
{
    double* Qd = c->Qdot.data();
    pfor(c->nCells, [=](int i) { Qd[i] *= factor; });
}


void cdmGpuRKStage(cdmGpuCase* c, double dt, double a, double b)
{
    double *r = c->rho.data(), *ru = c->rhoU.data(),
           *rE = c->rhoE.data(), *rY = c->rhoY.data();
    const double *r0 = c->rho0.data(), *ru0 = c->rhoU0.data(),
                 *rE0 = c->rhoE0.data(), *rY0 = c->rhoY0.data();
    const double *dr = c->dRho.data(), *dru = c->dRhoU.data(),
                 *drE = c->dRhoE.data(), *drY = c->dRhoY.data();
    const int n = c->nCells;

    pfor(n, [=](int i)
    {
        r[i]  = a*r0[i]  + b*(r[i]  + dt*dr[i]);
        rE[i] = a*rE0[i] + b*(rE[i] + dt*drE[i]);
        rY[i] = a*rY0[i] + b*(rY[i] + dt*drY[i]);
        for (int d = 0; d < 3; ++d)
        {
            ru[3*i+d] = a*ru0[3*i+d] + b*(ru[3*i+d] + dt*dru[3*i+d]);
        }
    });
}


double cdmGpuComputeRHS(cdmGpuCase* c, int* nFallback)
{
    const int nC = c->nCells, nIF = c->nIF, nBF = c->nBF;

    const int *own = c->own.data(), *nei = c->nei.data();
    const int *bOwn = c->bOwn.data();
    const int *cfOff = c->cfOff.data(), *cfFace = c->cfFace.data();
    const signed char* cfSign = c->cfSign.data();
    const double *Sf = c->Sf.data(), *magSf = c->magSf.data(), *wgt = c->wgt.data();
    const double *dcf = c->dc.data();
    const double *dlt = c->dlt.data(), *V = c->V.data();
    const double *bSf = c->bSf.data(), *bMagSf = c->bMagSf.data();
    const double *bDc = c->bDc.data();

    const double *rho = c->rho.data(), *U = c->U.data(), *p = c->p.data();
    const double *T = c->T.data(), *Y = c->Y.data();
    const double *bRho = c->bRho.data(), *bU = c->bU.data();
    const double *bP = c->bP.data(), *bY = c->bY.data();

    double *gRho = c->gRho.data(), *gU = c->gU.data();
    double *gP = c->gP.data(), *gY = c->gY.data();

    const double gam = c->gamma, Cp = c->Cv*c->gamma, eOff = c->eOffset;
    const double c_Cv = c->Cv;
    const double kappa0 = c->kappa0, nExp = c->nExp, Pr = c->Pr, Le = c->Le;
    const int scheme = c->flux, lim = c->lim, viscous = c->viscous;

    // ---- cell gradients, as a gather over the CSR cell-to-face map -------
    // Gauss linear: grad(phi)_P = (1/V) sum_f s_f Sf_f phi_f. A gather rather
    // than a face loop scattering into cells avoids atomics, which stdpar
    // cannot express portably.
    pfor(nC, [=](int i)
    {
        double ar[3] = {0,0,0}, ap[3] = {0,0,0}, ay[3] = {0,0,0};
        double au[9] = {0,0,0,0,0,0,0,0,0};

        for (int k = cfOff[i]; k < cfOff[i+1]; ++k)
        {
            const int f = cfFace[k];
            const double s = double(cfSign[k]);

            double rf, pf, yf, uf[3], S[3];

            if (f < nIF)
            {
                const int P = own[f], N = nei[f];
                const double w = wgt[f], w1 = 1.0 - w;
                rf = w*rho[P] + w1*rho[N];
                pf = w*p[P]   + w1*p[N];
                yf = w*Y[P]   + w1*Y[N];
                for (int d = 0; d < 3; ++d)
                {
                    uf[d] = w*U[3*P+d] + w1*U[3*N+d];
                    S[d] = Sf[3*f+d];
                }
            }
            else
            {
                const int b = f - nIF;
                rf = bRho[b]; pf = bP[b]; yf = bY[b];
                for (int d = 0; d < 3; ++d)
                {
                    uf[d] = bU[3*b+d];
                    S[d] = bSf[3*b+d];
                }
            }

            for (int d = 0; d < 3; ++d)
            {
                const double sS = s*S[d];
                ar[d] += sS*rf;
                ap[d] += sS*pf;
                ay[d] += sS*yf;
                for (int e = 0; e < 3; ++e) { au[3*d+e] += sS*uf[e]; }
            }
        }

        const double iV = 1.0/V[i];
        for (int d = 0; d < 3; ++d)
        {
            gRho[3*i+d] = ar[d]*iV;
            gP[3*i+d]   = ap[d]*iV;
            gY[3*i+d]   = ay[d]*iV;
            for (int e = 0; e < 3; ++e) { gU[9*i+3*d+e] = au[3*d+e]*iV; }
        }
    });

    // ---- internal faces --------------------------------------------------
    double *Fr = c->Fr.data(), *FrU = c->FrU.data();
    double *FrE = c->FrE.data(), *FrY = c->FrY.data();
    double *amaxSf = c->amaxSf.data(), *fallback = c->fallback.data();

    pfor(nIF, [=](int f)
    {
        const int P = own[f], N = nei[f];
        const double w = wgt[f], w1 = 1.0 - w;
        const double aSf = magSf[f];
        const double n[3] = {Sf[3*f]/aSf, Sf[3*f+1]/aSf, Sf[3*f+2]/aSf};
        const double d[3] = {dlt[3*f], dlt[3*f+1], dlt[3*f+2]};

        // --- MUSCL reconstruction, NVDTVD/NVDVTVDV as in OpenFOAM ---------
        double rL, rR, pL, pR, yL, yR, uL[3], uR[3];

        {
            const double gf = rho[N] - rho[P];
            const double gcP = d[0]*gRho[3*P] + d[1]*gRho[3*P+1] + d[2]*gRho[3*P+2];
            const double gcN = d[0]*gRho[3*N] + d[1]*gRho[3*N+1] + d[2]*gRho[3*N+2];
            rL = rho[P] + limiterL(lim, rScalar(gcP, gf))*w1*gf;
            rR = rho[N] - limiterL(lim, rScalar(gcN, gf))*w *gf;
        }
        {
            const double gf = p[N] - p[P];
            const double gcP = d[0]*gP[3*P] + d[1]*gP[3*P+1] + d[2]*gP[3*P+2];
            const double gcN = d[0]*gP[3*N] + d[1]*gP[3*N+1] + d[2]*gP[3*N+2];
            pL = p[P] + limiterL(lim, rScalar(gcP, gf))*w1*gf;
            pR = p[N] - limiterL(lim, rScalar(gcN, gf))*w *gf;
        }
        {
            const double gf = Y[N] - Y[P];
            const double gcP = d[0]*gY[3*P] + d[1]*gY[3*P+1] + d[2]*gY[3*P+2];
            const double gcN = d[0]*gY[3*N] + d[1]*gY[3*N+1] + d[2]*gY[3*N+2];
            yL = Y[P] + limiterL(lim, rScalar(gcP, gf))*w1*gf;
            yR = Y[N] - limiterL(lim, rScalar(gcN, gf))*w *gf;
        }
        {
            // NVDVTVDV: one limiter for the vector, from the projection of
            // the cell gradient onto the face jump
            double gfv[3], dgP[3] = {0,0,0}, dgN[3] = {0,0,0};
            for (int e = 0; e < 3; ++e) { gfv[e] = U[3*N+e] - U[3*P+e]; }
            for (int e = 0; e < 3; ++e)
            {
                for (int q = 0; q < 3; ++q)
                {
                    dgP[e] += d[q]*gU[9*P+3*q+e];
                    dgN[e] += d[q]*gU[9*N+3*q+e];
                }
            }
            const double gf = gfv[0]*gfv[0] + gfv[1]*gfv[1] + gfv[2]*gfv[2];
            const double gcP = gfv[0]*dgP[0] + gfv[1]*dgP[1] + gfv[2]*dgP[2];
            const double gcN = gfv[0]*dgN[0] + gfv[1]*dgN[1] + gfv[2]*dgN[2];
            const double LL = limiterL(lim, rScalar(gcP, gf));
            const double LR = limiterL(lim, rScalar(gcN, gf));
            for (int e = 0; e < 3; ++e)
            {
                uL[e] = U[3*P+e] + LL*w1*gfv[e];
                uR[e] = U[3*N+e] - LR*w *gfv[e];
            }
        }

        // --- first-order fallback, as cdmFoam's fluxScheme -----------
        const bool badL = !(rL > 0.0) || !(pL > 0.0);
        const bool badR = !(rR > 0.0) || !(pR > 0.0);
        fallback[f] = (badL || badR) ? 1.0 : 0.0;

        if (badL)
        {
            rL = rho[P]; pL = p[P]; yL = Y[P];
            for (int e = 0; e < 3; ++e) { uL[e] = U[3*P+e]; }
        }
        if (badR)
        {
            rR = rho[N]; pR = p[N]; yR = Y[N];
            for (int e = 0; e < 3; ++e) { uR[e] = U[3*N+e]; }
        }

        yL = yL < 0.0 ? 0.0 : (yL > 1.0 ? 1.0 : yL);
        yR = yR < 0.0 ? 0.0 : (yR > 1.0 ? 1.0 : yR);

        FaceFlux F = riemann(scheme, rL, uL, pL, yL, rR, uR, pR, yR, n, gam, eOff);

        double fr = F.r*aSf, fe = F.E*aSf, fy = F.Y*aSf;
        double fu[3] = {F.mU[0]*aSf, F.mU[1]*aSf, F.mU[2]*aSf};

        // --- CDM diffusive fluxes ----------------------------------------
        if (viscous)
        {
            const double dc = dcf[f];               // mesh.nonOrthDeltaCoeffs()

            const double kP = kappa0*std::pow(T[P], nExp);
            const double kN = kappa0*std::pow(T[N], nExp);
            const double kf = w*kP + w1*kN;
            const double mu = Pr*kf;

            double uf[3], gUf[9], tr = 0.0;
            for (int e = 0; e < 3; ++e) { uf[e] = w*U[3*P+e] + w1*U[3*N+e]; }
            for (int e = 0; e < 9; ++e) { gUf[e] = w*gU[9*P+e] + w1*gU[9*N+e]; }
            for (int e = 0; e < 3; ++e) { tr += gUf[3*e+e]; }

            // tauFlux = mu*(magSf*snGrad(U) + Sf & dev2(T(gradU)))
            double tau[3];
            for (int e = 0; e < 3; ++e)
            {
                double sg = 0.0;
                for (int qd = 0; qd < 3; ++qd)
                {
                    // dev2(transpose(gradU))[qd][e]
                    double M = gUf[3*e+qd];
                    if (qd == e) { M -= (2.0/3.0)*tr; }
                    sg += Sf[3*f+qd]*M;
                }
                tau[e] = mu*(aSf*dc*(U[3*N+e] - U[3*P+e]) + sg);
            }

            double tauU = 0.0;
            for (int e = 0; e < 3; ++e) { fu[e] -= tau[e]; tauU += tau[e]*uf[e]; }

            fe -= tauU + kf*Cp*dc*(T[N] - T[P])*aSf;
            fy -= (kf/Le)*dc*(Y[N] - Y[P])*aSf;
        }

        Fr[f] = fr; FrE[f] = fe; FrY[f] = fy;
        for (int e = 0; e < 3; ++e) { FrU[3*f+e] = fu[e]; }
        amaxSf[f] = F.amax*aSf;
    });

    // ---- boundary faces --------------------------------------------------
    // The patch values are the boundary state: OpenFOAM's own patch fields
    // have already evaluated every boundary condition in the package,
    // including the RDE injector and the conjugate temperature coupling.
    double *bFr = c->bFr.data(), *bFrU = c->bFrU.data();
    double *bFrE = c->bFrE.data(), *bFrY = c->bFrY.data();
    double *bAmax = c->bAmaxSf.data();

    pfor(nBF, [=](int b)
    {
        const double aSf = bMagSf[b];
        if (!(aSf > 0.0))        // empty patch
        {
            bFr[b] = bFrE[b] = bFrY[b] = bAmax[b] = 0.0;
            bFrU[3*b] = bFrU[3*b+1] = bFrU[3*b+2] = 0.0;
            return;
        }

        const double n[3] = {bSf[3*b]/aSf, bSf[3*b+1]/aSf, bSf[3*b+2]/aSf};
        const double ub[3] = {bU[3*b], bU[3*b+1], bU[3*b+2]};

        FaceFlux F = riemann
        (
            scheme, bRho[b], ub, bP[b], bY[b], bRho[b], ub, bP[b], bY[b],
            n, gam, eOff
        );

        double fr = F.r*aSf, fe = F.E*aSf, fy = F.Y*aSf;
        double fu[3] = {F.mU[0]*aSf, F.mU[1]*aSf, F.mU[2]*aSf};

        // --- boundary diffusive fluxes -----------------------------------
        // These are NOT zero on a wall or a symmetryPlane. A symmetryPlane on a
        // vector removes the normal component, so snGrad(U) there is
        // -n (n.U_cell) deltaCoeffs, which is exactly where a cellular
        // detonation's transverse waves reflect. Omitting the term leaves mass
        // and energy untouched but changes the peak pressure at the triple points.
        if (viscous)
        {
            const int P = bOwn[b];
            const double dcb = bDc[b];

            // kByCp at the patch, as linearInterpolate(kByCp) gives there.
            // T on the patch follows from p and rho: T = p/(rho R),
            // R = Cv (gamma - 1).
            const double R = c_Cv*(gam - 1.0);
            const double Tb = bP[b]/(bRho[b]*R);
            const double kb = kappa0*std::pow(Tb, nExp);
            const double mub = Pr*kb;

            // snGrad at the patch
            double snU[3];
            for (int e = 0; e < 3; ++e)
            {
                snU[e] = (ub[e] - U[3*P+e])*dcb;
            }
            const double snT = (Tb - T[P])*dcb;
            const double snY = (bY[b] - Y[P])*dcb;

            // gaussGrad::correctBoundaryConditions replaces the normal
            // component of the extrapolated cell gradient with snGrad:
            //   gB = gP + n (x) (snGrad - (n . gP))
            double gB[9];
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    gB[3*i+j] = gU[9*P+3*i+j];
                }
            }
            for (int j = 0; j < 3; ++j)
            {
                double ng = 0.0;
                for (int k = 0; k < 3; ++k) { ng += n[k]*gU[9*P+3*k+j]; }
                const double corr = snU[j] - ng;
                for (int i = 0; i < 3; ++i) { gB[3*i+j] += n[i]*corr; }
            }

            double trB = 0.0;
            for (int e = 0; e < 3; ++e) { trB += gB[3*e+e]; }

            double tau[3];
            for (int e = 0; e < 3; ++e)
            {
                double sg = 0.0;
                for (int qd = 0; qd < 3; ++qd)
                {
                    double M = gB[3*e+qd];
                    if (qd == e) { M -= (2.0/3.0)*trB; }
                    sg += bSf[3*b+qd]*M;
                }
                tau[e] = mub*(aSf*snU[e] + sg);
            }

            double tauU = 0.0;
            for (int e = 0; e < 3; ++e) { fu[e] -= tau[e]; tauU += tau[e]*ub[e]; }

            fe -= tauU + kb*Cp*snT*aSf;
            fy -= (kb/Le)*snY*aSf;
        }

        bFr[b] = fr; bFrE[b] = fe; bFrY[b] = fy;
        for (int e = 0; e < 3; ++e) { bFrU[3*b+e] = fu[e]; }
        bAmax[b] = F.amax*aSf;
    });

    // ---- divergence, the same CSR gather ---------------------------------
    double *dRho = c->dRho.data(), *dRhoU = c->dRhoU.data();
    double *dRhoE = c->dRhoE.data(), *dRhoY = c->dRhoY.data();

    const double *cFr = Fr, *cFrU = FrU, *cFrE = FrE, *cFrY = FrY;
    const double *cbFr = bFr, *cbFrU = bFrU, *cbFrE = bFrE, *cbFrY = bFrY;

    pfor(nC, [=](int i)
    {
        double ar = 0.0, ae = 0.0, ay = 0.0, au[3] = {0,0,0};

        for (int k = cfOff[i]; k < cfOff[i+1]; ++k)
        {
            const int f = cfFace[k];
            const double s = double(cfSign[k]);

            if (f < nIF)
            {
                ar += s*cFr[f]; ae += s*cFrE[f]; ay += s*cFrY[f];
                for (int e = 0; e < 3; ++e) { au[e] += s*cFrU[3*f+e]; }
            }
            else
            {
                const int b = f - nIF;
                ar += s*cbFr[b]; ae += s*cbFrE[b]; ay += s*cbFrY[b];
                for (int e = 0; e < 3; ++e) { au[e] += s*cbFrU[3*b+e]; }
            }
        }

        const double iV = -1.0/V[i];           // dU/dt = -div(F)
        dRho[i] = ar*iV; dRhoE[i] = ae*iV; dRhoY[i] = ay*iV;
        for (int e = 0; e < 3; ++e) { dRhoU[3*i+e] = au[e]*iV; }
    });

    // ---- reductions ------------------------------------------------------
    const double amaxI = preduce(nIF, 0.0, MaxOp(), [=](int f) -> double
        { return amaxSf[f]; });
    const double amaxB = preduce(nBF, 0.0, MaxOp(), [=](int b) -> double
        { return bAmax[b]; });

    if (nFallback)
    {
        *nFallback = int(preduce(nIF, 0.0, SumOp(), [=](int f) -> double
            { return fallback[f]; }));
    }

    return amaxI > amaxB ? amaxI : amaxB;
}

} // extern "C"
