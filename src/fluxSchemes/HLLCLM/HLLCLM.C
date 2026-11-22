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

\*---------------------------------------------------------------------------*/

#include "HLLCLM.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLLCLM, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLCLM, fvMesh);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLCLM::HLLCLM(const fvMesh& mesh)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLCLM::~HLLCLM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::CDM::fluxScheme::faceFlux
Foam::CDM::fluxSchemes::HLLCLM::calculateFluxes
(
    const faceState& L,
    const faceState& R,
    const vector& n,
    const scalar magSf,
    const label facei,
    const label patchi
) const
{
    static const scalar MaLimit = 0.1;

    const scalar EL = L.e + 0.5*magSqr(L.U);
    const scalar ER = R.e + 0.5*magSqr(R.U);

    const scalar uL = L.U & n;
    const scalar uR = R.U & n;

    // Roe averages, Eqs. (7)-(8)
    const scalar sqrtRhoL = sqrt(L.rho);
    const scalar sqrtRhoR = sqrt(R.rho);
    const scalar rSum = 1/(sqrtRhoL + sqrtRhoR);
    const scalar uHat = (uL*sqrtRhoL + uR*sqrtRhoR)*rSum;
    const scalar cHat = sqrt
    (
        (sqr(L.c)*sqrtRhoL + sqr(R.c)*sqrtRhoR)*rSum
      + 0.5*sqrtRhoL*sqrtRhoR*sqr(rSum)*sqr(uR - uL)
    );

    // Einfeldt wave speeds, Eq. (6)
    const scalar SL = min(uL - L.c, uHat - cHat);
    const scalar SR = max(uR + R.c, uHat + cHat);

    // Fluxes of the left and right states
    const scalar FrhoL = L.rho*uL;
    const vector FrhoUL = L.rho*L.U*uL + L.p*n;
    const scalar FrhoEL = (L.rho*EL + L.p)*uL;
    const scalar FrhoYL = FrhoL*L.Y;

    const scalar FrhoR = R.rho*uR;
    const vector FrhoUR = R.rho*R.U*uR + R.p*n;
    const scalar FrhoER = (R.rho*ER + R.p)*uR;
    const scalar FrhoYR = FrhoR*R.Y;

    faceFlux F;

    if (SL >= 0)
    {
        F.rhoPhi = FrhoL;
        F.rhoUPhi = FrhoUL;
        F.rhoEPhi = FrhoEL;
        F.rhoYPhi = FrhoYL;
    }
    else if (SR <= 0)
    {
        F.rhoPhi = FrhoR;
        F.rhoUPhi = FrhoUR;
        F.rhoEPhi = FrhoER;
        F.rhoYPhi = FrhoYR;
    }
    else
    {
        // Contact speed, Eq. (9)
        const scalar aL = L.rho*(SL - uL);
        const scalar aR = R.rho*(SR - uR);
        const scalar SStar = (R.p - L.p + aL*uL - aR*uR)/(aL - aR);

        // Star states, Eq. (5)
        const scalar kL = aL/(SL - SStar);
        const scalar kR = aR/(SR - SStar);

        const scalar rhoSL = kL;
        const vector rhoUSL = kL*(L.U + (SStar - uL)*n);
        const scalar rhoESL =
            kL*(EL + (SStar - uL)*(SStar + L.p/aL));
        const scalar rhoYSL = kL*L.Y;

        const scalar rhoSR = kR;
        const vector rhoUSR = kR*(R.U + (SStar - uR)*n);
        const scalar rhoESR =
            kR*(ER + (SStar - uR)*(SStar + R.p/aR));
        const scalar rhoYSR = kR*R.Y;

        // Low Mach number correction of the acoustic dissipation,
        // Eqs. (23)-(25)
        const scalar MaLocal = max(mag(uL)/L.c, mag(uR)/R.c);
        const scalar phi = Foam::sin
        (
            min(scalar(1), MaLocal/MaLimit)
           *constant::mathematical::piByTwo
        );

        const scalar sL = phi*SL;
        const scalar sR = phi*SR;
        const scalar sS = mag(SStar);

        // Central form, Eq. (19)
        F.rhoPhi =
            0.5*(FrhoL + FrhoR)
          + 0.5*
            (
                sL*(rhoSL - L.rho)
              + sS*(rhoSL - rhoSR)
              + sR*(rhoSR - R.rho)
            );

        F.rhoUPhi =
            0.5*(FrhoUL + FrhoUR)
          + 0.5*
            (
                sL*(rhoUSL - L.rho*L.U)
              + sS*(rhoUSL - rhoUSR)
              + sR*(rhoUSR - R.rho*R.U)
            );

        F.rhoEPhi =
            0.5*(FrhoEL + FrhoER)
          + 0.5*
            (
                sL*(rhoESL - L.rho*EL)
              + sS*(rhoESL - rhoESR)
              + sR*(rhoESR - R.rho*ER)
            );

        F.rhoYPhi =
            0.5*(FrhoYL + FrhoYR)
          + 0.5*
            (
                sL*(rhoYSL - L.rho*L.Y)
              + sS*(rhoYSL - rhoYSR)
              + sR*(rhoYSR - R.rho*R.Y)
            );
    }

    F.rhoPhi *= magSf;
    F.rhoUPhi *= magSf;
    F.rhoEPhi *= magSf;
    F.rhoYPhi *= magSf;

    return F;
}


// ************************************************************************* //
