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
    Publishes the chemical-diffusive thermal conductivity

        K = kappa0 Cp T^n

    as a registered field, so that a CDMCoupledTemperature patch on an
    interface with a solid region uses exactly the conductivity the solver
    uses for its own heat conduction.  Only constructed when CHT is enabled.

\*---------------------------------------------------------------------------*/

#include "cdmFoam.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::cdmFoam::updateK(const volScalarField& kByCp)
{
    if (!KPtr_.valid())
    {
        KPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                kByCp*thermo_.Cp()
            )
        );
    }
    else
    {
        KPtr_() = kByCp*thermo_.Cp();
    }

    KPtr_().correctBoundaryConditions();
}


// ************************************************************************* //
