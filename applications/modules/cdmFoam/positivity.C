/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam.
    GNU General Public License v3 or later; see COPYING.

Description
    Positivity-preserving cut-off flux limiter of Hu, Adams & Shu,
    "Positivity-preserving method for high-order conservative schemes solving
    compressible Euler equations", J. Comput. Phys. 242 (2013) 169-180,
    extended to unstructured finite volumes.

    The forward-Euler update of cell i (each SSP-RK stage),
        U_i' = U_i - dt/V_i sum_f s_f F_f,
    is written as the convex combination sum_f w_f U_i^f with
        w_f   = a_f |S_f|/W_i,   W_i = sum_f a_f |S_f|,
        U_i^f = U_i - dt W_i/(V_i a_f |S_f|) s_f F_f,
    a_f being the maximum wave speed of the face.  With the first-order
    Lax-Friedrichs flux F_LF every U_i^f has positive density and pressure
    for dt W_i/V_i <= 1, i.e. Co <= 0.5 with the Courant number of the
    solver.  Where the high-order flux F_HO would give rho(U_i^f) < rhoMin
    (or p < pMin), the face flux is replaced by
        F = F_LF + theta (F_HO - F_LF),   0 <= theta <= 1,
    with theta the largest value satisfying the constraint on both sides
    of the face (Hu et al., Section 2.2; the pressure is concave in U so the
    linear estimate of theta is conservative).  The same theta is applied to
    all the conserved quantities of the face, so conservation is exact.

    F_HO is the complete face flux (inviscid + viscous + diffusive); F_LF is
    the inviscid first-order Lax-Friedrichs flux of the cell states.

\*---------------------------------------------------------------------------*/

#include "cdmFoam.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::cdmFoam::limitPositivity
(
    const scalar deltaT,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi,
    surfaceScalarField& rhoYPhi
)
{
    const label nInternal = mesh.nInternalFaces();
    const label nBoundary = mesh.nFaces() - nInternal;

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& V = mesh.V().primitiveField();

    const volScalarField c(soundSpeed());

    // Conservative + primitive state of a cell or boundary face
    struct State
    {
        scalar rho;
        vector rhoU;
        scalar rhoE;
        scalar rhoY;
        vector U;
        scalar p;
        scalar c;
    };

    struct Flux
    {
        scalar rho;
        vector rhoU;
        scalar rhoE;
        scalar rhoY;
    };

    auto cellState = [&](const label i)
    {
        return State
        {
            rho_[i], rhoU_[i], rhoE_[i], rhoY_[i], U_[i], p_[i], c[i]
        };
    };

    // Area-weighted Euler flux of a state
    auto eulerFlux = [](const State& s, const vector& S)
    {
        const scalar un = s.U & S;
        return Flux
        {
            s.rho*un,
            s.rhoU*un + s.p*S,
            (s.rhoE + s.p)*un,
            s.rhoY*un
        };
    };

    // First-order Lax-Friedrichs flux and wave speed times area
    auto LF = [&](const State& o, const State& n, const vector& S,
        const scalar magS, scalar& aS)
    {
        const vector nf(S/magS);
        const scalar a = max(mag(o.U & nf) + o.c, mag(n.U & nf) + n.c);
        aS = a*magS;

        const Flux Fo(eulerFlux(o, S));
        const Flux Fn(eulerFlux(n, S));

        return Flux
        {
            0.5*(Fo.rho + Fn.rho) + 0.5*aS*(o.rho - n.rho),
            0.5*(Fo.rhoU + Fn.rhoU) + 0.5*aS*(o.rhoU - n.rhoU),
            0.5*(Fo.rhoE + Fn.rhoE) + 0.5*aS*(o.rhoE - n.rhoE),
            0.5*(Fo.rhoY + Fn.rhoY) + 0.5*aS*(o.rhoY - n.rhoY)
        };
    };

    // Pressure of a conservative state
    const scalar gm1 = gamma_ - 1;
    auto pressure = [&](const scalar r, const vector& rU, const scalar rE)
    {
        return gm1*(rE - 0.5*magSqr(rU)/r - r*eOffset_);
    };


    // Boundary states (coupled: neighbour cell values)
    List<State> bState(nBoundary);
    forAll(rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rho_.boundaryField()[patchi];
        const label start = prho.patch().start() - nInternal;

        if (prho.coupled())
        {
            const scalarField rN(prho.patchNeighbourField());
            const vectorField rUN
            (
                rhoU_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField rEN
            (
                rhoE_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField rYN
            (
                rhoY_.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField UN(U_.boundaryField()[patchi].patchNeighbourField());
            const scalarField pN(p_.boundaryField()[patchi].patchNeighbourField());
            const scalarField cN(c.boundaryField()[patchi].patchNeighbourField());

            forAll(prho, i)
            {
                bState[start + i] =
                    State{rN[i], rUN[i], rEN[i], rYN[i], UN[i], pN[i], cN[i]};
            }
        }
        else
        {
            forAll(prho, i)
            {
                bState[start + i] = State
                {
                    prho[i],
                    rhoU_.boundaryField()[patchi][i],
                    rhoE_.boundaryField()[patchi][i],
                    rhoY_.boundaryField()[patchi][i],
                    U_.boundaryField()[patchi][i],
                    p_.boundaryField()[patchi][i],
                    c.boundaryField()[patchi][i]
                };
            }
        }
    }


    // Lax-Friedrichs fluxes, face weights and cell sums W_i
    List<Flux> FLFi(nInternal);
    scalarField aSi(nInternal);
    List<Flux> FLFb(nBoundary);
    scalarField aSb(nBoundary, 0);
    scalarField W(mesh.nCells(), 0);

    for (label facei = 0; facei < nInternal; facei++)
    {
        FLFi[facei] = LF
        (
            cellState(own[facei]), cellState(nei[facei]),
            Sf[facei], magSf[facei], aSi[facei]
        );
        W[own[facei]] += aSi[facei];
        W[nei[facei]] += aSi[facei];
    }

    forAll(rho_.boundaryField(), patchi)
    {
        const fvPatch& pp = rho_.boundaryField()[patchi].patch();
        const labelUList& faceCells = pp.faceCells();
        const label start = pp.start() - nInternal;

        forAll(faceCells, i)
        {
            const label bi = start + i;
            FLFb[bi] = LF
            (
                cellState(faceCells[i]), bState[bi],
                Sf.boundaryField()[patchi][i],
                magSf.boundaryField()[patchi][i],
                aSb[bi]
            );
            W[faceCells[i]] += aSb[bi];
        }
    }


    // Limiting factor for the side of cell i with outward flux sign s
    // constraint: g(U_i - k s F) >= floor, g = rho or p
    auto thetaRho = [&](const label i, const scalar k, const scalar s,
        const Flux& FHO, const Flux& FLF)
    {
        const scalar rHO = rho_[i] - k*s*FHO.rho;
        if (rHO >= rhoFloor_) return scalar(1);
        const scalar rLF = rho_[i] - k*s*FLF.rho;
        if (rLF <= rhoFloor_) return scalar(0);
        return min(max((rLF - rhoFloor_)/(rLF - rHO), scalar(0)), scalar(1));
    };

    auto thetaP = [&](const label i, const scalar k, const scalar s,
        const Flux& FHO, const Flux& FLF)
    {
        const scalar pHO = pressure
        (
            rho_[i] - k*s*FHO.rho,
            rhoU_[i] - k*s*FHO.rhoU,
            rhoE_[i] - k*s*FHO.rhoE
        );
        if (pHO >= pFloor_) return scalar(1);
        const scalar pLF = pressure
        (
            rho_[i] - k*s*FLF.rho,
            rhoU_[i] - k*s*FLF.rhoU,
            rhoE_[i] - k*s*FLF.rhoE
        );
        if (pLF <= pFloor_) return scalar(0);
        return min(max((pLF - pFloor_)/(pLF - pHO), scalar(0)), scalar(1));
    };

    auto blend = [](const Flux& FLF, const Flux& FHO, const scalar theta)
    {
        return Flux
        {
            FLF.rho + theta*(FHO.rho - FLF.rho),
            FLF.rhoU + theta*(FHO.rhoU - FLF.rhoU),
            FLF.rhoE + theta*(FHO.rhoE - FLF.rhoE),
            FLF.rhoY + theta*(FHO.rhoY - FLF.rhoY)
        };
    };


    // Current high-order fluxes
    List<Flux> FHOi(nInternal);
    for (label facei = 0; facei < nInternal; facei++)
    {
        FHOi[facei] = Flux
        {
            rhoPhi[facei], rhoUPhi[facei], rhoEPhi[facei], rhoYPhi[facei]
        };
    }
    List<Flux> FHOb(nBoundary);
    forAll(rho_.boundaryField(), patchi)
    {
        const label start = rho_.boundaryField()[patchi].patch().start()
          - nInternal;
        forAll(rhoPhi.boundaryField()[patchi], i)
        {
            FHOb[start + i] = Flux
            {
                rhoPhi.boundaryField()[patchi][i],
                rhoUPhi.boundaryField()[patchi][i],
                rhoEPhi.boundaryField()[patchi][i],
                rhoYPhi.boundaryField()[patchi][i]
            };
        }
    }

    // Two passes: density, then pressure with the density-limited flux
    scalarField thetaI(nInternal, 1);
    scalarField thetaB(nBoundary, 1);

    for (label pass = 0; pass < 2; pass++)
    {
        auto theta = [&](const label i, const scalar k, const scalar sgn,
            const Flux& FHO, const Flux& FLF)
        {
            return pass == 0
              ? thetaRho(i, k, sgn, FHO, FLF)
              : thetaP(i, k, sgn, FHO, FLF);
        };

        scalarField tI(nInternal, 1);
        scalarField tB(nBoundary, 1);

        for (label facei = 0; facei < nInternal; facei++)
        {
            const label o = own[facei];
            const label n = nei[facei];
            const Flux F(blend(FLFi[facei], FHOi[facei], thetaI[facei]));
            const scalar ko = deltaT*W[o]/(V[o]*aSi[facei]);
            const scalar kn = deltaT*W[n]/(V[n]*aSi[facei]);

            tI[facei] = min
            (
                theta(o, ko, 1, F, FLFi[facei]),
                theta(n, kn, -1, F, FLFi[facei])
            );
        }

        forAll(rho_.boundaryField(), patchi)
        {
            const fvPatch& pp = rho_.boundaryField()[patchi].patch();
            const labelUList& faceCells = pp.faceCells();
            const label start = pp.start() - nInternal;

            forAll(faceCells, i)
            {
                const label bi = start + i;
                const label o = faceCells[i];
                const Flux F(blend(FLFb[bi], FHOb[bi], thetaB[bi]));
                const scalar ko = deltaT*W[o]/(V[o]*aSb[bi]);

                tB[bi] = theta(o, ko, 1, F, FLFb[bi]);
            }
        }

        // The other side of coupled faces is limited on the other processor
        scalarField tBN(tB);
        syncTools::swapBoundaryFaceList(mesh, tBN);

        forAll(rho_.boundaryField(), patchi)
        {
            const fvPatch& pp = rho_.boundaryField()[patchi].patch();
            if (pp.coupled())
            {
                const label start = pp.start() - nInternal;
                forAll(pp, i)
                {
                    tB[start + i] = min(tB[start + i], tBN[start + i]);
                }
            }
        }

        // theta_total = theta_rho*theta_p (F** - F_LF = theta_p(F* - F_LF))
        thetaI *= tI;
        thetaB *= tB;
    }


    // Apply
    label nLimited = 0;

    for (label facei = 0; facei < nInternal; facei++)
    {
        if (thetaI[facei] < 1)
        {
            const Flux F(blend(FLFi[facei], FHOi[facei], thetaI[facei]));
            rhoPhi[facei] = F.rho;
            rhoUPhi[facei] = F.rhoU;
            rhoEPhi[facei] = F.rhoE;
            rhoYPhi[facei] = F.rhoY;
            nLimited++;
        }
    }

    forAll(rho_.boundaryField(), patchi)
    {
        const label start = rho_.boundaryField()[patchi].patch().start()
          - nInternal;

        forAll(rhoPhi.boundaryField()[patchi], i)
        {
            const label bi = start + i;

            if (thetaB[bi] < 1)
            {
                const Flux F(blend(FLFb[bi], FHOb[bi], thetaB[bi]));
                rhoPhi.boundaryFieldRef()[patchi][i] = F.rho;
                rhoUPhi.boundaryFieldRef()[patchi][i] = F.rhoU;
                rhoEPhi.boundaryFieldRef()[patchi][i] = F.rhoE;
                rhoYPhi.boundaryFieldRef()[patchi][i] = F.rhoY;
                nLimited++;
            }
        }
    }

    nLimitedFaces_ += nLimited;
}


// ************************************************************************* //
