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

#include "HLLC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLC, fvMesh);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLC::HLLC(const fvMesh& mesh)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLC::~HLLC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::CDM::fluxScheme::faceFlux
Foam::CDM::fluxSchemes::HLLC::calculateFluxes
(
    const faceState& L,
    const faceState& R,
    const vector& n,
    const scalar magSf,
    const label facei,
    const label patchi
) const
{
    const scalar EL = L.e + 0.5*magSqr(L.U);
    const scalar ER = R.e + 0.5*magSqr(R.U);

    const scalar uL = L.U & n;
    const scalar uR = R.U & n;

    // Roe-type (sqrt(rho) weighted) averages for the wave-speed estimates
    const scalar sqrtRhoL = sqrt(L.rho);
    const scalar sqrtRhoR = sqrt(R.rho);
    const scalar wL = sqrtRhoL/(sqrtRhoL + sqrtRhoR);
    const scalar wR = 1 - wL;
    const scalar uT = wL*uL + wR*uR;
    const scalar cT = wL*L.c + wR*R.c;

    // Davis-Einfeldt bounds: include both one-sided estimates so that the
    // fastest waves are never underestimated (positivity of HLLC)
    const scalar SL = min(min(uL - L.c, uR - R.c), uT - cT);
    const scalar SR = max(max(uR + R.c, uL + L.c), uT + cT);

    const scalar aL = L.rho*(SL - uL);
    const scalar aR = R.rho*(SR - uR);

    const scalar SStar = (aR*uR - aL*uL + L.p - R.p)/(aR - aL);
    const scalar pStar = (aR*L.p - aL*R.p - aL*aR*(uL - uR))/(aR - aL);

    faceFlux F;

    if (SL >= 0)
    {
        F.rhoPhi = L.rho*uL;
        F.rhoUPhi = L.rho*L.U*uL + L.p*n;
        F.rhoEPhi = (L.rho*EL + L.p)*uL;
        F.rhoYPhi = F.rhoPhi*L.Y;
    }
    else if (SR <= 0)
    {
        F.rhoPhi = R.rho*uR;
        F.rhoUPhi = R.rho*R.U*uR + R.p*n;
        F.rhoEPhi = (R.rho*ER + R.p)*uR;
        F.rhoYPhi = F.rhoPhi*R.Y;
    }
    else
    {
        // Star-region flux on the side of the contact
        const bool left = SStar >= 0;
        const faceState& K = left ? L : R;
        const scalar SK = left ? SL : SR;
        const scalar uK = left ? uL : uR;
        const scalar EK = left ? EL : ER;
        const scalar dS = SK - SStar;

        const vector rhoUK(K.rho*K.U);
        const vector FrhoUK(rhoUK*uK + K.p*n);
        const scalar FrhoEK((K.rho*EK + K.p)*uK);

        F.rhoPhi = SStar*K.rho*(SK - uK)/dS;
        F.rhoUPhi = (SStar*(SK*rhoUK - FrhoUK) + SK*pStar*n)/dS;
        F.rhoEPhi = SStar*(SK*K.rho*EK - FrhoEK + SK*pStar)/dS;
        F.rhoYPhi = F.rhoPhi*K.Y;
    }

    F.rhoPhi *= magSf;
    F.rhoUPhi *= magSf;
    F.rhoEPhi *= magSf;
    F.rhoYPhi *= magSf;

    return F;
}


// ************************************************************************* //
