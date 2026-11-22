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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::cdmFoam::encode()
{
    const volScalarField& e = thermo_.he();

    rhoU_ = rho_*U_;
    rhoE_ = rho_*(e + 0.5*magSqr(U_));
    rhoY_ = rho_*Y_;
}


void Foam::solvers::cdmFoam::decode(const bool evaluateBCs)
{
    const scalarField& rhoI = rho_.primitiveField();

    double tMark = CDMWallTime();

    // Velocity
    U_.primitiveFieldRef() = rhoU_.primitiveField()/rhoI;
    if (evaluateBCs)
    {
        U_.correctBoundaryConditions();
    }

    tDec_[0] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // Reactant mass fraction, bounded
    Y_.primitiveFieldRef() =
        min(max(rhoY_.primitiveField()/rhoI, scalar(0)), scalar(1));
    if (evaluateBCs)
    {
        Y_.correctBoundaryConditions();
    }

    tDec_[1] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // Sensible internal energy
    volScalarField& e = thermo_.he();
    e.primitiveFieldRef() =
        rhoE_.primitiveField()/rhoI - 0.5*magSqr(U_.primitiveField());

    tDec_[2] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // Stop with a description of the cell before the thermo fails
    checkState();

    tDec_[3] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // In a CHT run the coupled temperature condition needs the neighbouring
    // region, which does not exist while this solver is being constructed. The
    // fields have just been read, so their boundary values are already valid.
    if (evaluateBCs)
    {
        e.correctBoundaryConditions();
    }

    // Optional limits, e.g. limitTemperature in system/fvConstraints
    fvConstraints().constrain(e);

    tDec_[4] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // T, psi from e
    thermo_.correct();

    tDec_[5] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // Pressure from the equation of state
    const volScalarField& psi = thermo_.psi();
    p_.primitiveFieldRef() = rhoI/psi.primitiveField();
    if (evaluateBCs)
    {
        p_.correctBoundaryConditions();
    }
    rho_.boundaryFieldRef() == psi.boundaryField()*p_.boundaryField();

    tDec_[6] += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // Re-encode so that the conservative variables are consistent with
    // any clipping/constraints and with the boundary values
    encode();

    tDec_[7] += CDMWallTime() - tMark;
}


void Foam::solvers::cdmFoam::checkState() const
{
    const scalarField& rhoI = rho_.primitiveField();
    const scalarField& eI = thermo_.he().primitiveField();
    const scalar Cv = gAverage(thermo_.Cv().primitiveField());

    label nBad = 0;

    forAll(rhoI, celli)
    {
        const scalar Test = (eI[celli] - eOffset_)/Cv;

        if (rhoI[celli] <= 0 || Test <= 0 || Test > TFail_)
        {
            if (nBad < 5)
            {
                const cell& cFaces = mesh.cells()[celli];
                wordList patchNames;
                forAll(cFaces, i)
                {
                    if (!mesh.isInternalFace(cFaces[i]))
                    {
                        const label patchi =
                            mesh.boundaryMesh().whichPatch(cFaces[i]);
                        const word& name = mesh.boundaryMesh()[patchi].name();
                        if (findIndex(patchNames, name) == -1)
                        {
                            patchNames.append(name);
                        }
                    }
                }

                Pout<< "Bad cell " << celli
                    << " at " << mesh.C()[celli]
                    << " V = " << mesh.V()[celli]
                    << " rho = " << rhoI[celli]
                    << " U = " << U_[celli]
                    << " T(e) = " << Test
                    << " p(e) = "
                    << rhoI[celli]*(gamma_ - 1)*(eI[celli] - eOffset_)
                    << " Y = " << rhoY_[celli]/rhoI[celli]
                    << " boundary patches " << patchNames << endl;
            }
            nBad++;
        }
    }

    if (returnReduce(nBad, sumOp<label>()))
    {
        FatalErrorInFunction
            << returnReduce(nBad, sumOp<label>())
            << " cells with non-positive density or temperature, or"
            << " T > TFail = " << TFail_ << " at time " << runTime.name()
            << nl << "    See the 'Bad cell' lines above for the locations."
            << exit(FatalError);
    }
}


void Foam::solvers::cdmFoam::reportHotSpot() const
{
    const scalarField& TI = thermo_.T().primitiveField();

    scalar Tmax = -great;
    label cellMax = -1;
    forAll(TI, celli)
    {
        if (TI[celli] > Tmax)
        {
            Tmax = TI[celli];
            cellMax = celli;
        }
    }

    const scalar TmaxGlobal = returnReduce(Tmax, maxOp<scalar>());

    if (TmaxGlobal > TWarn_ && cellMax >= 0 && Tmax == TmaxGlobal)
    {
        Pout<< "Hot spot: T = " << Tmax
            << " at " << mesh.C()[cellMax]
            << " rho = " << rho_[cellMax]
            << " p = " << p_[cellMax]
            << " U = " << U_[cellMax]
            << " V = " << mesh.V()[cellMax] << endl;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::cdmFoam::soundSpeed() const
{
    return volScalarField::New
    (
        "c",
        sqrt(thermo_.Cp()/thermo_.Cv()/thermo_.psi())
    );
}


// ************************************************************************* //
