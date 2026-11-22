/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM.

    CDMGpuFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CDMGpuFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CDMGpuFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CDMModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "cpuLoad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
    defineTypeNameAndDebug(CDMModel, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::CDMModel::CDMModel(const fvMesh& mesh, const scalar Rspecific)
:
    IOdictionary
    (
        IOobject
        (
            "CDMProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    Rspecific_(Rspecific),
    cpuLoad_(false),
    A_(0),
    Ta_(0),
    q_(0),
    kappa0_(0),
    n_(0.7),
    Pr_(1),
    Le_(1),
    maxSubstepFraction_(0.2),
    maxSubsteps_(1000)
{
    read();
    printSummary();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::CDMModel::~CDMModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::CDM::CDMModel::read()
{
    const dictionary& dict = *this;

    // Reaction parameters: SI or the non-dimensional/CGS form of the paper
    if (dict.found("A"))
    {
        A_ = dict.lookup<scalar>("A");
        Ta_ = dict.lookup<scalar>("Ta");
        q_ = dict.lookup<scalar>("q");
        kappa0_ = dict.lookup<scalar>("kappa0");
    }
    else if (dict.found("A_cgs"))
    {
        const scalar T0 = dict.lookup<scalar>("T0");

        // cm^3/(g s) -> m^3/(kg s)
        A_ = 1e-3*dict.lookup<scalar>("A_cgs");

        Ta_ = dict.lookup<scalar>("EaByRT0")*T0;

        q_ = dict.lookup<scalar>("qMwByRT0")*Rspecific_*T0;

        // g/(s cm K^n) -> kg/(s m K^n)
        kappa0_ = 0.1*dict.lookup<scalar>("kappa0_cgs");
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Specify the CDM parameters either in SI units"
            << " (A, Ta, q, kappa0)" << nl
            << "or in the form of the tables of Lu et al. (2021)"
            << " (A_cgs, EaByRT0, qMwByRT0, kappa0_cgs, T0)"
            << exit(FatalIOError);
    }

    n_ = dict.lookupOrDefault<scalar>("n", 0.7);
    Pr_ = dict.lookupOrDefault<scalar>("Pr", 1);
    Le_ = dict.lookupOrDefault<scalar>("Le", 1);

    maxSubstepFraction_ =
        dict.lookupOrDefault<scalar>("maxSubstepFraction", 0.2);
    maxSubsteps_ = dict.lookupOrDefault<label>("maxSubsteps", 1000);

    cpuLoad_ = dict.lookupOrDefault<Switch>("cpuLoad", false);

    if (A_ < 0 || Ta_ < 0 || kappa0_ < 0 || Pr_ <= 0 || Le_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Unphysical CDM parameters" << exit(FatalIOError);
    }

    return true;
}


void Foam::CDM::CDMModel::printSummary() const
{
    Info<< nl << "Chemical-diffusive model (CDM)" << nl
        << "    A        = " << A_ << " m^3/(kg s)" << nl
        << "    Ta       = " << Ta_ << " K" << nl
        << "    q        = " << q_ << " J/kg" << nl
        << "    kappa0   = " << kappa0_ << " kg/(m s K^" << n_ << ")" << nl
        << "    Pr = " << Pr_ << ", Le = " << Le_ << nl
        << "    cpuLoad  = " << cpuLoad_ << nl
        << "    R (Cp-Cv)= " << Rspecific_ << " J/(kg K)" << nl << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::CDM::CDMModel::kappaByCp(const volScalarField& T) const
{
    return volScalarField::New
    (
        "kappaByCp",
        dimensionedScalar(dimDynamicViscosity, kappa0_)
       *pow(T/dimensionedScalar(dimTemperature, 1), n_)
    );
}


void Foam::CDM::CDMModel::react
(
    const scalar deltaT,
    const scalarField& rho,
    const scalarField& T,
    const scalarField& Cv,
    scalarField& rhoY,
    scalarField& rhoE,
    scalarField& heatReleased
) const
{
    optionalCpuLoad& chemistryCpuLoad
    (
        optionalCpuLoad::New("CDM:cpuLoad", mesh_, cpuLoad_)
    );

    chemistryCpuLoad.resetCpuTime();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Y0 = min(max(rhoY[celli]/rhoi, scalar(0)), scalar(1));

        if (Y0 <= small)
        {
            chemistryCpuLoad.cpuTimeIncrement(celli);
            continue;
        }

        scalar Y = Y0;
        scalar Ti = T[celli];
        scalar t = 0;
        label nSub = 0;

        while (t < deltaT && Y > small)
        {
            const scalar ki = k(rhoi, Ti);

            // Negligible reaction over the remaining interval
            if (ki*(deltaT - t) < small)
            {
                break;
            }

            const scalar h = min(deltaT - t, maxSubstepFraction_/ki);

            const scalar Ynew = Y*exp(-ki*h);

            // Constant-volume heat release raises the temperature
            Ti += q_*(Y - Ynew)/Cv[celli];

            Y = Ynew;
            t += h;

            if (++nSub >= maxSubsteps_)
            {
                // Finish the step with the frozen rate
                Y *= exp(-k(rhoi, Ti)*(deltaT - t));
                break;
            }
        }

        const scalar dQ = rhoi*q_*(Y0 - Y);

        rhoY[celli] = rhoi*Y;
        rhoE[celli] += dQ;
        heatReleased[celli] += dQ;

        chemistryCpuLoad.cpuTimeIncrement(celli);
    }
}


// ************************************************************************* //
