/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM.
    GNU General Public License v3 or later; see COPYING.

\*---------------------------------------------------------------------------*/

#include "cdmFoam.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::cdmFoam::waveSpeedSf() const
{
    tmp<surfaceScalarField> tamaxSf
    (
        surfaceScalarField::New
        (
            "amaxSf",
            mesh,
            dimensionedScalar(dimVolume/dimTime, 0)
        )
    );
    surfaceScalarField& amaxSf = tamaxSf.ref();

    const volScalarField c(soundSpeed());

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(own, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];

        amaxSf[facei] = max
        (
            mag(U_[o] & Sf[facei]) + c[o]*magSf[facei],
            mag(U_[n] & Sf[facei]) + c[n]*magSf[facei]
        );
    }

    forAll(amaxSf.boundaryField(), patchi)
    {
        const fvPatchVectorField& pU = U_.boundaryField()[patchi];
        const fvPatchScalarField& pc = c.boundaryField()[patchi];
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pmagSf = magSf.boundaryField()[patchi];
        fvsPatchScalarField& pa = amaxSf.boundaryFieldRef()[patchi];

        const vectorField UI(pU.patchInternalField());
        const scalarField cI(pc.patchInternalField());

        forAll(pa, i)
        {
            pa[i] = max
            (
                mag(UI[i] & pSf[i]) + cI[i]*pmagSf[i],
                mag(pU[i] & pSf[i]) + pc[i]*pmagSf[i]
            );
        }

        if (pU.coupled())
        {
            const vectorField UN(pU.patchNeighbourField());
            const scalarField cN(pc.patchNeighbourField());

            forAll(pa, i)
            {
                pa[i] = max
                (
                    pa[i],
                    mag(UN[i] & pSf[i]) + cN[i]*pmagSf[i]
                );
            }
        }
    }

    return tamaxSf;
}


void Foam::solvers::cdmFoam::correctCoNum
(
    const surfaceScalarField& amaxSf
)
{
    const scalarField& V = mesh.V().primitiveField();
    const scalar deltaT = runTime.deltaTValue();

    // Convective Courant number (as in shockFluid)
    const scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    const scalar convCo = 0.5*gMax(sumAmaxSf/V)*deltaT;

    // Explicit diffusion number: dt*nu*sum(|Sf| deltaCoeffs)/V
    scalar diffNo = 0;

    if (viscous_)
    {
        const scalarField sumDeltaSf
        (
            fvc::surfaceSum(mesh.magSf()*mesh.deltaCoeffs())()
           .primitiveField()
        );

        const scalarField nu
        (
            cdm_.maxDiffusivityFactor()
           *cdm_.kappaByCp(thermo_.T())().primitiveField()
           /rho_.primitiveField()
        );

        diffNo = gMax(nu*sumDeltaSf/V)*deltaT;
    }

    // Both limits are controlled by maxCo in controlDict
    CoNum_ = max(convCo, diffNo);

    Info<< "Courant Number max: " << convCo
        << "  Diffusion Number max: " << diffNo << endl;
}


void Foam::solvers::cdmFoam::updateDiagnostics()
{
    // Numerical smoke foil
    pMax_.primitiveFieldRef() =
        max(pMax_.primitiveField(), p_.primitiveField());
    pMax_.correctBoundaryConditions();

    // Refinement indicator: max over the faces of a cell of the relative
    // jumps of rho and p and the jump of Y; 1 in reacting gas
    {
        scalarField& ind = refinementIndicator_.primitiveFieldRef();
        ind = 0;

        const scalarField& Yi = Y_.primitiveField();
        const scalarField& Ti = thermo_.T().primitiveField();

        auto jump = [](const scalar a, const scalar b)
        {
            return mag(b - a)/max(min(a, b), vSmall);
        };

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        forAll(own, facei)
        {
            const label o = own[facei];
            const label n = nei[facei];

            const scalar j = max
            (
                max(jump(rho_[o], rho_[n]), jump(p_[o], p_[n])),
                mag(Yi[n] - Yi[o])
            );

            ind[o] = max(ind[o], j);
            ind[n] = max(ind[n], j);
        }

        forAll(rho_.boundaryField(), patchi)
        {
            const fvPatchScalarField& prho = rho_.boundaryField()[patchi];

            if (prho.coupled())
            {
                const fvPatchScalarField& pp = p_.boundaryField()[patchi];
                const fvPatchScalarField& pY = Y_.boundaryField()[patchi];

                const scalarField rI(prho.patchInternalField());
                const scalarField rN(prho.patchNeighbourField());
                const scalarField pI(pp.patchInternalField());
                const scalarField pN(pp.patchNeighbourField());
                const scalarField YI(pY.patchInternalField());
                const scalarField YN(pY.patchNeighbourField());

                const labelUList& faceCells = prho.patch().faceCells();

                forAll(faceCells, i)
                {
                    const scalar j = max
                    (
                        max(jump(rI[i], rN[i]), jump(pI[i], pN[i])),
                        mag(YN[i] - YI[i])
                    );

                    ind[faceCells[i]] = max(ind[faceCells[i]], j);
                }
            }
        }

        forAll(ind, celli)
        {
            if (Ti[celli] > TReacting_ && Yi[celli] > YReacting_)
            {
                ind[celli] = max(ind[celli], scalar(1));
            }
        }

        refinementIndicator_.correctBoundaryConditions();
    }

    // Leading-shock tracking
    if (shockTracking_)
    {
        const vectorField& C = mesh.C().primitiveField();
        const scalarField& pI = p_.primitiveField();

        scalar xShock = -great;

        forAll(pI, celli)
        {
            if (pI[celli] > shockPThreshold_)
            {
                xShock = max(xShock, C[celli] & shockDirection_);
            }
        }

        reduce(xShock, maxOp<scalar>());

        if (xShock > -great)
        {
            shockTimes_.append(runTime.value());
            shockPositions_.append(xShock);

            // Keep the last shockWindow_ samples
            if (shockTimes_.size() > shockWindow_)
            {
                const label nDrop = shockTimes_.size() - shockWindow_;

                for (label i = 0; i < shockWindow_; i++)
                {
                    shockTimes_[i] = shockTimes_[i + nDrop];
                    shockPositions_[i] = shockPositions_[i + nDrop];
                }

                shockTimes_.setSize(shockWindow_);
                shockPositions_.setSize(shockWindow_);
            }

            const label n = shockTimes_.size();

            if (n > 1)
            {
                const scalar D =
                    (shockPositions_[n - 1] - shockPositions_[0])
                   /(shockTimes_[n - 1] - shockTimes_[0] + vSmall);

                Info<< "Leading shock: x = " << xShock
                    << "  D(" << n << " steps) = " << D << " m/s" << endl;
            }
        }
    }
}


void Foam::solvers::cdmFoam::reportTimeSplit() const
{
    const double total =
        tRHS_ + tStage_ + tChem_ + tTail_ + tDiag_ + tSoot_ + tCourant_;

    if (total <= 0)
    {
        return;
    }

    const double pc = 100.0/total;

    // A flux evaluation that costs nothing did not happen. The way this goes
    // wrong is a derived solver declaring a member of the same name, which
    // shadows it: the derived advance() credits its own copy and the line
    // below reads zero.
    if (tRHS_ <= 0)
    {
        WarningInFunction
            << "computeRHS is not being timed - it reports zero while "
            << total << " s was tracked elsewhere." << nl
            << "    The split below is incomplete and does not add up to the"
            << " runtime.  Check that no derived solver shadows tRHS_."
            << endl;
    }

    Info<< nl << type() << " per-step time split" << nl
        << "    (tracked wall clock " << total
        << " s; percentages of that)" << nl
        << "    computeRHS          " << tRHS_ << " s  " << tRHS_*pc << " %"
        << nl
        << "    RK stage + decode   " << tStage_ << " s  " << tStage_*pc
        << " %" << nl
        << "    chemistry           " << tChem_ << " s  " << tChem_*pc << " %"
        << nl
        << "    step tail           " << tTail_ << " s  " << tTail_*pc << " %"
        << nl
        << "      (Qdot, hot spot, front tracking, min/max and totals)" << nl
        << "    updateDiagnostics   " << tDiag_ << " s  " << tDiag_*pc << " %"
        << nl
        << "      (pMax and the refinement indicator)" << nl
        << "    soot foil           " << tSoot_ << " s  " << tSoot_*pc << " %"
        << nl
        << "    Courant number      " << tCourant_ << " s  " << tCourant_*pc
        << " %" << nl;

    // decode() is called from the Runge-Kutta stages, from react() and from
    // preSolve(), so its segments are a breakdown across all of those rather
    // than a subset of any one line above.
    static const char* decodeNames[8] =
    {
        "U = rhoU/rho + BCs ",
        "Y bounded + BCs    ",
        "e = rhoE/rho - K   ",
        "checkState         ",
        "e BCs + constraints",
        "thermo_.correct()  ",
        "p = rho/psi + BCs  ",
        "encode()           "
    };

    double decodeTotal = 0;
    for (int i = 0; i < 8; i++)
    {
        decodeTotal += tDec_[i];
    }

    Info<< nl << "    decode() segments, summed over every call"
        << " (total " << decodeTotal << " s  " << decodeTotal*pc << " %)"
        << nl;

    for (int i = 0; i < 8; i++)
    {
        Info<< "      " << decodeNames[i] << "  " << tDec_[i] << " s  "
            << tDec_[i]*pc << " %" << nl;
    }

    Info<< endl;
}


// ************************************************************************* //
