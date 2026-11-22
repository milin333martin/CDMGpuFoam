/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM and the HLLCP
    implementation of detonationFoam (Sun et al., CPC 2023).
    GNU General Public License v3 or later; see COPYING.

\*---------------------------------------------------------------------------*/

#include "HLLCP.H"
#include "surfaceInterpolate.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLLCP, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLCP, fvMesh);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLCP::HLLCP(const fvMesh& mesh)
:
    fluxScheme(mesh)
{
    WarningInFunction
        << "HLLCP is experimental in CDMGpuFoam: its pressure-control mass"
        << " flux (Xie et al. 2018, Eq. 16) drove cells to negative density"
        << " at shock collisions in 2D cellular detonations." << nl
        << "    Use HLLCLM (shock-stable) or HLLC." << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::fluxSchemes::HLLCP::~HLLCP()
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::CDM::fluxSchemes::HLLCP::preUpdate(const volScalarField& p)
{
    // Cell sensor: min over the faces of the cell of min(pP/pN, pN/pP)
    tmp<volScalarField> tfCells
    (
        volScalarField::New
        (
            "HLLCP:fCells",
            mesh_,
            dimensionedScalar(dimless, 1)
        )
    );
    volScalarField& fCells = tfCells.ref();
    scalarField& fc = fCells.primitiveFieldRef();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(own, facei)
    {
        const scalar pO = p[own[facei]];
        const scalar pN = p[nei[facei]];
        const scalar r = min(pO/pN, pN/pO);

        fc[own[facei]] = min(fc[own[facei]], r);
        fc[nei[facei]] = min(fc[nei[facei]], r);
    }

    forAll(p.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = p.boundaryField()[patchi];

        if (pp.coupled())
        {
            const scalarField pI(pp.patchInternalField());
            const scalarField pN(pp.patchNeighbourField());
            const labelUList& faceCells = pp.patch().faceCells();

            forAll(faceCells, i)
            {
                const scalar r = min(pI[i]/pN[i], pN[i]/pI[i]);
                fc[faceCells[i]] = min(fc[faceCells[i]], r);
            }
        }
    }

    fc = pow3(fc);

    // Boundary values: coupled patches are swapped by
    // correctBoundaryConditions, others take the adjacent cell value
    fCells.correctBoundaryConditions();

    forAll(fCells.boundaryField(), patchi)
    {
        fvPatchScalarField& pf = fCells.boundaryFieldRef()[patchi];

        if (!pf.coupled())
        {
            pf == pf.patchInternalField();
        }
    }

    f_ = linearInterpolate(fCells);
}


void Foam::CDM::fluxSchemes::HLLCP::postUpdate()
{
    f_.clear();
}


Foam::CDM::fluxScheme::faceFlux
Foam::CDM::fluxSchemes::HLLCP::calculateFluxes
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

    const scalar MaL = uL/L.c;
    const scalar MaR = uR/R.c;

    const scalar sqrtRhoL = sqrt(L.rho);
    const scalar sqrtRhoR = sqrt(R.rho);
    const scalar wL = sqrtRhoL/(sqrtRhoL + sqrtRhoR);
    const scalar wR = 1 - wL;

    const scalar MaT = wL*MaL + wR*MaR;
    const vector UT(wL*L.U + wR*R.U);
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

    // Low-Mach interface pressure
    const scalar theta = min(max(mag(MaL), mag(MaR)), scalar(1));
    const scalar pStarStar = theta*pStar + (1 - theta)*0.5*(L.p + R.p);

    // Pressure sensor
    const scalar f =
        patchi == -1
      ? f_()[facei]
      : f_().boundaryField()[patchi][facei];

    const scalar pStarStarStar = f*pStarStar + (1 - f)*pStar;

    // Pressure-difference driven mass-flux correction (zero where f = 1)
    const scalar phip =
        (f - 1)*SL*SR/(SR - SL)/(1 + mag(MaT))*(R.p - L.p)/sqr(cT);

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
        const bool left = SStar >= 0;
        const faceState& K = left ? L : R;
        const scalar SK = left ? SL : SR;
        const scalar uK = left ? uL : uR;
        const scalar EK = left ? EL : ER;
        const scalar dS = SK - SStar;

        const vector rhoUK(K.rho*K.U);
        const vector FrhoUK(rhoUK*uK + K.p*n);
        const scalar FrhoEK((K.rho*EK + K.p)*uK);

        F.rhoPhi = SStar*K.rho*(SK - uK)/dS + phip;
        F.rhoUPhi =
            (SStar*(SK*rhoUK - FrhoUK) + SK*pStarStarStar*n)/dS
          + phip*UT;
        F.rhoEPhi =
            SStar*(SK*K.rho*EK - FrhoEK + SK*pStar)/dS
          + 0.5*phip*magSqr(UT);

        // Reactant upwinded with the corrected mass flux
        F.rhoYPhi = F.rhoPhi*(F.rhoPhi >= 0 ? L.Y : R.Y);
    }

    F.rhoPhi *= magSf;
    F.rhoUPhi *= magSf;
    F.rhoEPhi *= magSf;
    F.rhoYPhi *= magSf;

    return F;
}


// ************************************************************************* //
