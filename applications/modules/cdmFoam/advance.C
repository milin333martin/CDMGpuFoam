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
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "linear.H"
#include "globalMeshData.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::cdmFoam::computeRHS
(
    const scalar deltaT,
    scalarField& dRho,
    vectorField& dRhoU,
    scalarField& dRhoE,
    scalarField& dRhoY
)
{
    const volScalarField& Tfield = thermo_.T();

    tmp<surfaceScalarField> trhoPhi
    (
        surfaceScalarField::New
        (
            "rhoPhi",
            mesh,
            dimensionedScalar(dimMass/dimTime, 0)
        )
    );
    tmp<surfaceVectorField> trhoUPhi
    (
        surfaceVectorField::New
        (
            "rhoUPhi",
            mesh,
            dimensionedVector(dimForce, Zero)
        )
    );
    tmp<surfaceScalarField> trhoEPhi
    (
        surfaceScalarField::New
        (
            "rhoEPhi",
            mesh,
            dimensionedScalar(dimPower, 0)
        )
    );
    tmp<surfaceScalarField> trhoYPhi
    (
        surfaceScalarField::New
        (
            "rhoYPhi",
            mesh,
            dimensionedScalar(dimMass/dimTime, 0)
        )
    );
    tmp<surfaceScalarField> tamaxSf
    (
        surfaceScalarField::New
        (
            "amaxSf",
            mesh,
            dimensionedScalar(dimVolume/dimTime, 0)
        )
    );

    surfaceScalarField& rhoPhi = trhoPhi.ref();
    surfaceVectorField& rhoUPhi = trhoUPhi.ref();
    surfaceScalarField& rhoEPhi = trhoEPhi.ref();
    surfaceScalarField& rhoYPhi = trhoYPhi.ref();

    // Inviscid fluxes
    fluxScheme_->update
    (
        rho_, U_, p_, Y_, gamma_, eOffset_,
        rhoPhi, rhoUPhi, rhoEPhi, rhoYPhi, tamaxSf.ref()
    );

    // Viscous, conductive and diffusive fluxes (CDM transport)
    if (viscous_)
    {
        const volScalarField kByCp(cdm_.kappaByCp(Tfield));
        const surfaceScalarField kByCpf(linearInterpolate(kByCp));

        // mu = Pr K/Cp
        const surfaceScalarField muf("muf", cdm_.Pr()*kByCpf);

        const volTensorField gradU(fvc::grad(U_));

        // Sf & tau
        const surfaceVectorField tauFlux
        (
            "tauFlux",
            muf
           *(
                mesh.magSf()*fvc::snGrad(U_)
              + (mesh.Sf() & linearInterpolate(dev2(Foam::T(gradU))))
            )
        );

        rhoUPhi -= tauFlux;

        if (CHT_)
        {
            updateK(kByCp);
        }

        // Viscous work and heat conduction, K = kappa0 Cp T^n
        rhoEPhi -=
            (tauFlux & linearInterpolate(U_))
          + linearInterpolate(kByCp*thermo_.Cp())
           *fvc::snGrad(Tfield)*mesh.magSf();

        // Reactant diffusion, rho D = K/(Cp Le)
        rhoYPhi -= (kByCpf/cdm_.Le())*fvc::snGrad(Y_)*mesh.magSf();
    }

    if (positivityLimiter_)
    {
        limitPositivity(deltaT, rhoPhi, rhoUPhi, rhoEPhi, rhoYPhi);
    }

    phi_ = rhoPhi;

    dRho = -fvc::surfaceIntegrate(rhoPhi)().primitiveField();
    dRhoU = -fvc::surfaceIntegrate(rhoUPhi)().primitiveField();
    dRhoE = -fvc::surfaceIntegrate(rhoEPhi)().primitiveField();
    dRhoY = -fvc::surfaceIntegrate(rhoYPhi)().primitiveField();
}


void Foam::solvers::cdmFoam::react(const scalar deltaT)
{
    cdm_.react
    (
        deltaT,
        rho_.primitiveField(),
        thermo_.T().primitiveField(),
        thermo_.Cv().primitiveField(),
        rhoY_.primitiveFieldRef(),
        rhoE_.primitiveFieldRef(),
        Qdot_.primitiveFieldRef()
    );

    decode();
}


void Foam::solvers::cdmFoam::advance()
{
    const scalar deltaT = runTime.deltaTValue();

    // SSP Runge-Kutta coefficients (Shu & Osher 1988):
    //     q^(k) = a_k q^0 + b_k (q^(k-1) + dt L(q^(k-1)))
    static const scalar a[3][3] =
    {
        {0, 0, 0},
        {0, 0.5, 0},
        {0, 0.75, 1.0/3.0}
    };
    static const scalar b[3][3] =
    {
        {1, 0, 0},
        {1, 0.5, 0},
        {1, 0.25, 2.0/3.0}
    };
    const label s = nRKStages_ - 1;

    // Heat released during the step, converted to a rate at the end
    Qdot_.primitiveFieldRef() = 0;

    nLimitedFaces_ = 0;

    double tMark = CDMWallTime();

    if (chemistry_ && strang_)
    {
        react(0.5*deltaT);
    }

    tChem_ += CDMWallTime() - tMark;

    // Initial state of the flow step
    const scalarField rho0(rho_.primitiveField());
    const vectorField rhoU0(rhoU_.primitiveField());
    const scalarField rhoE0(rhoE_.primitiveField());
    const scalarField rhoY0(rhoY_.primitiveField());

    scalarField dRho(mesh.nCells());
    vectorField dRhoU(mesh.nCells());
    scalarField dRhoE(mesh.nCells());
    scalarField dRhoY(mesh.nCells());

    for (label stage = 0; stage < nRKStages_; stage++)
    {
        tMark = CDMWallTime();

        computeRHS(deltaT, dRho, dRhoU, dRhoE, dRhoY);

        const double tRHSDone = CDMWallTime();
        tRHS_ += tRHSDone - tMark;

        const scalar ak = a[s][stage];
        const scalar bk = b[s][stage];

        scalarField& rhoI = rho_.primitiveFieldRef();
        vectorField& rhoUI = rhoU_.primitiveFieldRef();
        scalarField& rhoEI = rhoE_.primitiveFieldRef();
        scalarField& rhoYI = rhoY_.primitiveFieldRef();

        rhoI = ak*rho0 + bk*(rhoI + deltaT*dRho);
        rhoUI = ak*rhoU0 + bk*(rhoUI + deltaT*dRhoU);
        rhoEI = ak*rhoE0 + bk*(rhoEI + deltaT*dRhoE);
        rhoYI = ak*rhoY0 + bk*(rhoYI + deltaT*dRhoY);

        decode();

        tStage_ += CDMWallTime() - tRHSDone;
    }

    tMark = CDMWallTime();

    if (chemistry_)
    {
        react(strang_ ? 0.5*deltaT : deltaT);
    }

    tChem_ += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    Qdot_.primitiveFieldRef() /= deltaT;
    Qdot_.correctBoundaryConditions();

    if (frontTracking_ && mesh.time().timeIndex() % frontInterval_ == 0)
    {
        trackFronts();
    }

    if (positivityLimiter_)
    {
        const label n = returnReduce(nLimitedFaces_, sumOp<label>());
        if (n > 0)
        {
            Info<< "Positivity limiter: " << n << " face fluxes limited"
                << endl;
        }
    }

    // Every per-step monitor in one pass. These were eight separate sweeps
    // over the cells plus four temporary fields of the mesh's own size,
    // costing 43 s of a 491 s run at 484k cells to print two lines per step.
    // The sums accumulate in gSum's order, so the values printed are identical
    // to the last bit.
    scalar pMin = great, pMax = -great;
    scalar TMin = great, TMax = -great;
    scalar QdotMax = -great;
    scalar mass = 0, energy = 0;
    label cellTMax = -1;

    {
        const scalarField& pI = p_.primitiveField();
        const scalarField& TI = thermo_.T().primitiveField();
        const scalarField& QI = Qdot_.primitiveField();
        const scalarField& rhoI = rho_.primitiveField();
        const scalarField& rhoEI = rhoE_.primitiveField();
        const scalarField& rhoYI = rhoY_.primitiveField();
        const scalarField& V = mesh.V().primitiveField();
        const scalar qCdm = cdm_.q();

        forAll(pI, celli)
        {
            const scalar pc = pI[celli];
            if (pc < pMin) { pMin = pc; }
            if (pc > pMax) { pMax = pc; }

            const scalar Tc = TI[celli];
            if (Tc < TMin) { TMin = Tc; }
            if (Tc > TMax) { TMax = Tc; cellTMax = celli; }

            if (QI[celli] > QdotMax) { QdotMax = QI[celli]; }

            const scalar Vc = V[celli];
            mass += rhoI[celli]*Vc;
            energy += (rhoEI[celli] + qCdm*rhoYI[celli])*Vc;
        }
    }

    const scalar TMaxGlobal = returnReduce(TMax, maxOp<scalar>());

    if (TMaxGlobal > TWarn_ && cellTMax >= 0 && TMax == TMaxGlobal)
    {
        Pout<< "Hot spot: T = " << TMax
            << " at " << mesh.C()[cellTMax]
            << " rho = " << rho_[cellTMax]
            << " p = " << p_[cellTMax]
            << " U = " << U_[cellTMax]
            << " V = " << mesh.V()[cellTMax] << endl;
    }

    Info<< "min/max(p) = " << returnReduce(pMin, minOp<scalar>()) << ", "
        << returnReduce(pMax, maxOp<scalar>())
        << "  min/max(T) = " << returnReduce(TMin, minOp<scalar>()) << ", "
        << TMaxGlobal
        << "  max(Qdot) = " << returnReduce(QdotMax, maxOp<scalar>()) << endl;

    // Conserved totals (sensible + kinetic + chemical energy)
    Info<< "Totals: mass = " << returnReduce(mass, sumOp<scalar>())
        << "  energy = " << returnReduce(energy, sumOp<scalar>())
        << "  cells = " << mesh.globalData().nTotalCells() << endl;

    tTail_ += CDMWallTime() - tMark;
}


// ************************************************************************* //
