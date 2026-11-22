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
    Validation of the GPU back end without OpenFOAM.

    Builds a one-dimensional mesh in OpenFOAM's own addressing - owner,
    neighbour, Sf, magSf, weights, delta, the CSR cell-to-face map - and runs
    the ZND1D_stable case through cdmGpuComputeRHS, cdmGpuReact and
    cdmGpuRKStage.  The detonation speed must come out at D_CJ = 1967.15 m/s,
    the same number the tutorial requires of CDMGpuFoam.

    This is the test that matters, because the unstructured face loop is the
    part of the back end that has no counterpart in a structured code and
    therefore the part most likely to be wrong.  It exercises the CSR
    gradient gather, the NVDTVD/NVDVTVDV reconstruction, the limiter, the
    Riemann solver, the diffusive terms, the divergence gather, the chemistry
    and the Runge-Kutta stage - everything except the OpenFOAM glue.

    Build and run with no GPU and no OpenFOAM:
        g++ -std=c++17 -O3 -I.. testZND.cpp ../cdmGpuKernels.cpp -o testZND
        ./testZND

\*---------------------------------------------------------------------------*/

#include "cdmGpu.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    // ---- the stable verification gas, as tutorials/ZND1D_stable ----------
    const double gamma = 1.177;
    const double Cv = 1928.3401;
    const double R = Cv*(gamma - 1.0);
    const double A = 336543.78;
    const double Ta = 4560.86;
    const double q = 4715378.9;
    const double kappa0 = 2.367e-6;
    const double nExp = 0.7;
    const double DCJ = 1967.15;

    int N = 8000;
    double L = 0.80;
    double endTime = 4.0e-4;
    const double CFL = 0.4;

    // OpenFOAM's sensible-energy datum, e = p/(rho(gamma-1)) + eOffset. It is
    // a pure datum shift, so D must not depend on it. Validating with a
    // non-zero value is what catches an energy flux that has dropped the
    // offset, which a run with eOffset = 0 cannot see.
    double eOffset = 0.0;

    for (int i = 1; i < argc - 1; ++i)
    {
        const std::string a(argv[i]);
        if      (a == "--N")       N = std::atoi(argv[++i]);
        else if (a == "--L")       L = std::atof(argv[++i]);
        else if (a == "--endTime") endTime = std::atof(argv[++i]);
        else if (a == "--eOffset") eOffset = std::atof(argv[++i]);
    }

    const double dx = L/N;
    const double Af = 1.0;              // unit face area
    const double Vc = Af*dx;

    // ---- mesh, in OpenFOAM addressing ------------------------------------
    const int nIF = N - 1;
    const int nBF = 2;

    std::vector<int> own(nIF), nei(nIF);
    std::vector<double> Sf(3*nIF, 0.0), magSf(nIF, Af), wgt(nIF, 0.5),
                        dcf(nIF, 1.0/dx), dlt(3*nIF, 0.0);
    for (int f = 0; f < nIF; ++f)
    {
        own[f] = f; nei[f] = f + 1;
        Sf[3*f] = Af;                   // +x, from owner to neighbour
        dlt[3*f] = dx;                  // C[nei] - C[own]
    }

    std::vector<double> V(N, Vc);

    // boundary 0 is the closed left end, boundary 1 the open right end
    std::vector<int> bOwn(nBF);
    std::vector<double> bSf(3*nBF, 0.0), bMagSf(nBF, Af);
    // cell centre to face centre is dx/2 on a boundary face
    std::vector<double> bDc(nBF, 2.0/dx);
    bOwn[0] = 0;      bSf[0] = -Af;     // outward normal points -x
    bOwn[1] = N - 1;  bSf[3] = +Af;

    // CSR cell-to-face map: +1 where the cell owns the face, -1 where it
    // neighbours it.  Face indices at or above nIF are boundary faces.
    std::vector<int> cfOff(N + 1, 0), cfFace;
    std::vector<signed char> cfSign;
    cfFace.reserve(2*N); cfSign.reserve(2*N);

    for (int i = 0; i < N; ++i)
    {
        cfOff[i] = int(cfFace.size());
        if (i == 0)      { cfFace.push_back(nIF + 0); cfSign.push_back(+1); }
        else             { cfFace.push_back(i - 1);   cfSign.push_back(-1); }
        if (i == N - 1)  { cfFace.push_back(nIF + 1); cfSign.push_back(+1); }
        else             { cfFace.push_back(i);       cfSign.push_back(+1); }
    }
    cfOff[N] = int(cfFace.size());

    // ---- back end --------------------------------------------------------
    cdmGpuCase* c = cdmGpuCreate();

    cdmGpuSetMesh
    (
        c, N, 0 /* nHalo: serial, no processor faces */, nIF, nBF,
        own.data(), nei.data(),
        Sf.data(), magSf.data(), wgt.data(), dcf.data(), dlt.data(), V.data(),
        bOwn.data(), bSf.data(), bMagSf.data(), bDc.data(),
        cfOff.data(), cfFace.data(), cfSign.data()
    );

    // The cells behind the two boundary faces. cdmGpuFoam exchanges only these
    // between Runge-Kutta stages, so the gather is cross-checked against the
    // full download on every step below.
    std::vector<int> bCell(2);
    bCell[0] = 0;
    bCell[1] = N - 1;
    cdmGpuSetBoundaryCells(c, 2, bCell.data());

    cdmGpuSetModel
    (
        c, gamma, Cv, eOffset, A, Ta, q, kappa0, nExp, 1.0, 1.0,
        0.2, 1000, CDM_GPU_HLLCLM, CDM_GPU_VANLEER, 1
    );

    // ---- initial condition ------------------------------------------------
    const double p0 = 101325.0, T0 = 298.0;
    const double xDriver = 0.02, pDriver = 50.0*p0, TDriver = 3000.0;

    std::vector<double> rho(N), rhoU(3*N, 0.0), rhoE(N), rhoY(N);
    for (int i = 0; i < N; ++i)
    {
        const double x = (i + 0.5)*dx;
        const bool drv = (x < xDriver);
        const double pp = drv ? pDriver : p0;
        const double TT = drv ? TDriver : T0;
        const double YY = drv ? 0.0 : 1.0;
        const double rr = pp/(R*TT);
        rho[i] = rr; rhoE[i] = rr*(Cv*TT + eOffset); rhoY[i] = rr*YY;
    }
    cdmGpuUploadState(c, rho.data(), rhoU.data(), rhoE.data(), rhoY.data());

    std::vector<double> U(3*N), p(N), T(N), Y(N), Qdot(N);
    std::vector<double> bRho(nBF), bU(3*nBF, 0.0), bP(nBF), bY(nBF);

    // Boundary conditions, exactly as OpenFOAM's patch fields would supply
    // them: a closed wall on the left, so the flux is pure pressure, and
    // zero-gradient on the right.
    std::vector<double> gRho(2), gU(6), gP(2), gT(2), gY(2);
    long nGatherBad = 0;

    auto setBCs = [&]()
    {
        cdmGpuDownloadPrimitives(c, U.data(), p.data(), T.data(), Y.data(), Qdot.data());
        cdmGpuDownloadState(c, rho.data(), rhoU.data(), rhoE.data(), rhoY.data());

        // The boundary-cell gather must agree exactly with the full download
        cdmGpuDownloadBoundaryCells
        (
            c, gRho.data(), gU.data(), gP.data(), gT.data(), gY.data()
        );

        // Compared against bOwn, not against bCell: checking the gather against
        // the list handed to it only proves it can read, and would pass even if
        // the list named the wrong cells. bOwn is the independent statement of
        // which cell sits behind each boundary face.
        for (int i = 0; i < 2; ++i)
        {
            const int ci = bOwn[i];
            if (gRho[i] != rho[ci] || gP[i] != p[ci] || gT[i] != T[ci]
             || gY[i] != Y[ci]
             || gU[3*i] != U[3*ci] || gU[3*i+1] != U[3*ci+1]
             || gU[3*i+2] != U[3*ci+2])
            {
                ++nGatherBad;
            }
        }

        bRho[0] = rho[0];     bP[0] = p[0];     bY[0] = Y[0];
        bU[0] = bU[1] = bU[2] = 0.0;

        bRho[1] = rho[N-1];   bP[1] = p[N-1];   bY[1] = Y[N-1];
        bU[3] = U[3*(N-1)];   bU[4] = 0.0;      bU[5] = 0.0;

        cdmGpuUploadBoundary(c, bRho.data(), bU.data(), bP.data(), bY.data());
    };

    std::printf("cdmGpu back end: %s\n", cdmGpuBackend());
    std::printf("mesh  %d cells, %d internal faces, %d boundary faces, dx = %.1f um\n",
                N, nIF, nBF, dx*1e6);
    std::printf("gas   gamma = %.3f  R = %.3f  Cv = %.3f   D_CJ = %.2f m/s\n",
                gamma, R, Cv, DCJ);
    std::printf("datum eOffset = %.1f J/kg\n\n", eOffset);

    // ---- time loop: Strang splitting around SSP-RK2 ----------------------
    double t = 0.0;
    long step = 0;
    int nFallbackTotal = 0;
    std::vector<double> hT, hX;

    if (cdmGpuDecode(c) != 0) { std::printf("bad initial state\n"); return 1; }
    setBCs();

    // Chemical energy audit. The heat release reported through Qdot must equal
    // the chemical energy consumed, q times the drop in reactant mass. With
    // Strang splitting the step releases heat in two half-steps, so this only
    // balances if Qdot accumulates across both.
    double rhoY0Total = 0.0;
    for (int i = 0; i < N; ++i) { rhoY0Total += rhoY[i]; }
    rhoY0Total *= Vc;

    double qReported = 0.0;

    while (t < endTime)
    {
        cdmGpuZeroQdot(c);

        const double amax = cdmGpuComputeRHS(c, nullptr);

        // Courant: dt = CFL*dx/(|u| + c), with amax = (|u| + c)*Af
        double dt = CFL*Vc/amax;

        // explicit diffusion limit
        double dmax = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const double nu = kappa0*std::pow(T[i], nExp)/rho[i];
            if (nu > dmax) { dmax = nu; }
        }
        const double dtd = 0.25*dx*dx/(dmax > 0.0 ? dmax : 1e-30);
        if (dtd < dt) { dt = dtd; }
        if (t + dt > endTime) { dt = endTime - t; }

        cdmGpuReact(c, 0.5*dt);
        if (cdmGpuDecode(c) != 0) { std::printf("bad state after react\n"); return 1; }
        setBCs();

        cdmGpuStoreStage(c);

        int nf = 0;
        cdmGpuComputeRHS(c, &nf); nFallbackTotal += nf;
        cdmGpuRKStage(c, dt, 0.0, 1.0);
        if (cdmGpuDecode(c) != 0) { std::printf("bad state, stage 1\n"); return 1; }
        setBCs();

        cdmGpuComputeRHS(c, &nf); nFallbackTotal += nf;
        cdmGpuRKStage(c, dt, 0.5, 0.5);
        if (cdmGpuDecode(c) != 0) { std::printf("bad state, stage 2\n"); return 1; }
        setBCs();

        cdmGpuReact(c, 0.5*dt);
        if (cdmGpuDecode(c) != 0) { std::printf("bad state after react 2\n"); return 1; }
        setBCs();

        // Qdot still holds the accumulated release for this step
        for (int i = 0; i < N; ++i) { qReported += Qdot[i]; }

        // Leave it as a rate, as advance() does with Qdot_ /= deltaT
        cdmGpuScaleQdot(c, 1.0/dt);

        t += dt;
        ++step;

        if (step % 20 == 0)
        {
            int lead = -1;
            for (int i = N - 1; i >= 0; --i) { if (p[i] > 5.0*p0) { lead = i; break; } }
            if (lead >= 0) { hT.push_back(t); hX.push_back((lead + 0.5)*dx); }
        }

        if (step % 4000 == 0)
        {
            double pmax = 0.0, Tmax = 0.0;
            for (int i = 0; i < N; ++i)
            {
                if (p[i] > pmax) pmax = p[i];
                if (T[i] > Tmax) Tmax = T[i];
            }
            std::printf("  step %6ld  t = %.4e  dt = %.3e  max p = %.4e  max T = %.1f\n",
                        step, t, dt, pmax, Tmax);
            std::fflush(stdout);
        }
    }

    // ---- detonation speed -------------------------------------------------
    auto fit = [&](size_t i0, size_t i1) -> double
    {
        if (i1 <= i0 + 2) { return 0.0; }
        const double n = double(i1 - i0);
        double st = 0, sx = 0, stt = 0, stx = 0;
        for (size_t i = i0; i < i1; ++i)
        {
            st += hT[i]; sx += hX[i]; stt += hT[i]*hT[i]; stx += hT[i]*hX[i];
        }
        return (n*stx - st*sx)/(n*stt - st*st);
    };

    const size_t H = hT.size();
    std::printf("\n  leading-shock speed by window\n");
    for (int w = 0; w < 5; ++w)
    {
        const size_t i0 = H*w/5, i1 = H*(w + 1)/5;
        std::printf("    t = %.3e .. %.3e s   D = %8.2f m/s\n",
                    hT[i0], hT[i1 - 1], fit(i0, i1));
    }

    const double D = fit(H*4/5, H);
    const double err = 100.0*(D - DCJ)/DCJ;

    // ---- chemical energy audit -------------------------------------------
    cdmGpuDownloadState(c, rho.data(), rhoU.data(), rhoE.data(), rhoY.data());

    double rhoYTotal = 0.0;
    for (int i = 0; i < N; ++i) { rhoYTotal += rhoY[i]; }
    rhoYTotal *= Vc;

    const double qConsumed = q*(rhoY0Total - rhoYTotal);
    qReported *= Vc;

    const double qErr = (qConsumed != 0.0)
        ? 100.0*(qReported - qConsumed)/qConsumed : 0.0;

    std::printf("\n  steps              %ld\n", step);
    std::printf("  first-order falls  %d\n", nFallbackTotal);
    std::printf("  measured D         %.2f m/s\n", D);
    std::printf("  D_CJ reference     %.2f m/s\n", DCJ);
    std::printf("  error              %+.3f %%\n", err);
    std::printf("\n  heat release reported through Qdot   %.6e J\n", qReported);
    std::printf("  chemical energy consumed, q d(rhoY)  %.6e J\n", qConsumed);
    std::printf("  Qdot accounting error                %+.3e %%\n", qErr);

    std::printf("  boundary-cell gather mismatches       %ld\n", nGatherBad);

    const bool speedOK = std::fabs(err) < 1.0;
    const bool qdotOK = std::fabs(qErr) < 1.0e-6;
    const bool gatherOK = (nGatherBad == 0);

    std::printf("  %s\n",
        (speedOK && qdotOK && gatherOK) ? "VALIDATION PASSED"
                                        : "VALIDATION FAILED");

    if (!qdotOK)
    {
        std::printf("  (Qdot does not account for the chemical energy; check"
                    " that cdmGpuReact accumulates)\n");
    }

    cdmGpuDestroy(c);
    return (speedOK && qdotOK && gatherOK) ? 0 : 1;
}
