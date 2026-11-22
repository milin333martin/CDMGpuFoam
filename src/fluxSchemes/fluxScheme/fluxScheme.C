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

#include "fluxScheme.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
    defineTypeNameAndDebug(fluxScheme, 0);
    defineRunTimeSelectionTable(fluxScheme, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::fluxScheme::fluxScheme(const fvMesh& mesh)
:
    mesh_(mesh),
    pSchemeName_("reconstruct(p)")
{
    const dictionary& interpDict =
        mesh.schemes().dict().subDict("interpolationSchemes");

    if (!interpDict.found(pSchemeName_) && interpDict.found("reconstruct(T)"))
    {
        pSchemeName_ = "reconstruct(T)";
    }

    const wordList names({"rho", "U", "Y"});

    forAll(names, i)
    {
        const word schemeName("reconstruct(" + names[i] + ")");

        if (!interpDict.found(schemeName))
        {
            WarningInFunction
                << "No entry " << schemeName << " in interpolationSchemes."
                << nl << "    The default scheme will be used for the"
                << " MUSCL reconstruction of " << names[i]
                << "; use a limited scheme (e.g. vanLeer/vanLeerV)"
                << " to avoid oscillations." << endl;
        }
    }

    Info<< "    pressure reconstructed with " << pSchemeName_ << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::CDM::fluxScheme> Foam::CDM::fluxScheme::New
(
    const fvMesh& mesh
)
{
    const word schemeType
    (
        mesh.schemes().dict().lookupOrDefault<word>("fluxScheme", "HLLC")
    );

    Info<< "Selecting Riemann flux scheme " << schemeType << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(schemeType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluxScheme " << schemeType << nl << nl
            << "Valid fluxSchemes are : " << nl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::fluxScheme::~fluxScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CDM::fluxScheme::update
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& Y,
    const scalar gamma,
    const scalar eOffset,
    surfaceScalarField& rhoPhi,
    surfaceVectorField& rhoUPhi,
    surfaceScalarField& rhoEPhi,
    surfaceScalarField& rhoYPhi,
    surfaceScalarField& amaxSf
)
{
    // Direction fields for the upwind-biased (MUSCL) reconstruction.
    // Temporary, so nothing persists across topology changes.
    const tmp<surfaceScalarField> town
    (
        surfaceScalarField::New("own", mesh_, dimensionedScalar(dimless, 1))
    );
    const tmp<surfaceScalarField> tnei
    (
        surfaceScalarField::New("nei", mesh_, dimensionedScalar(dimless, -1))
    );
    const surfaceScalarField& own = town();
    const surfaceScalarField& nei = tnei();

    const surfaceScalarField rhoL(fvc::interpolate(rho, own, "reconstruct(rho)"));
    const surfaceScalarField rhoR(fvc::interpolate(rho, nei, "reconstruct(rho)"));

    const surfaceVectorField UL(fvc::interpolate(U, own, "reconstruct(U)"));
    const surfaceVectorField UR(fvc::interpolate(U, nei, "reconstruct(U)"));

    const surfaceScalarField pL(fvc::interpolate(p, own, pSchemeName_));
    const surfaceScalarField pR(fvc::interpolate(p, nei, pSchemeName_));

    const surfaceScalarField YL(fvc::interpolate(Y, own, "reconstruct(Y)"));
    const surfaceScalarField YR(fvc::interpolate(Y, nei, "reconstruct(Y)"));

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const scalar rGm1 = 1/(gamma - 1);

    // Thermodynamically consistent face state from (rho, U, p, Y)
    auto state = [&]
    (
        const scalar rhof,
        const vector& Uf,
        const scalar pf,
        const scalar Yf
    )
    {
        faceState s;
        s.rho = rhof;
        s.U = Uf;
        s.p = pf;
        s.e = pf/rhof*rGm1 + eOffset;
        s.c = sqrt(gamma*pf/rhof);
        s.Y = min(max(Yf, scalar(0)), scalar(1));
        return s;
    };

    preUpdate(p);

    // A reconstructed face value can lose positivity where the solution is very
    // steep. The cell value is positive by construction, so the reconstruction
    // drops to first order on that side of that face rather than taking the
    // square root of a negative number.
    label nFallback = 0;

    const labelUList& ownCell = mesh_.owner();
    const labelUList& neiCell = mesh_.neighbour();

    const scalarField& rhoI = rho.primitiveField();
    const vectorField& UI = U.primitiveField();
    const scalarField& pI = p.primitiveField();
    const scalarField& YI = Y.primitiveField();

    // Internal faces
    forAll(rhoL, facei)
    {
        const bool badL = rhoL[facei] <= 0 || pL[facei] <= 0;
        const bool badR = rhoR[facei] <= 0 || pR[facei] <= 0;

        if (badL || badR) nFallback++;

        const faceState L
        (
            badL
          ? state(rhoI[ownCell[facei]], UI[ownCell[facei]], pI[ownCell[facei]],
                  YI[ownCell[facei]])
          : state(rhoL[facei], UL[facei], pL[facei], YL[facei])
        );
        const faceState R
        (
            badR
          ? state(rhoI[neiCell[facei]], UI[neiCell[facei]], pI[neiCell[facei]],
                  YI[neiCell[facei]])
          : state(rhoR[facei], UR[facei], pR[facei], YR[facei])
        );

        const vector n(Sf[facei]/magSf[facei]);

        const faceFlux F
        (
            calculateFluxes(L, R, n, magSf[facei], facei, -1)
        );

        rhoPhi[facei] = F.rhoPhi;
        rhoUPhi[facei] = F.rhoUPhi;
        rhoEPhi[facei] = F.rhoEPhi;
        rhoYPhi[facei] = F.rhoYPhi;

        amaxSf[facei] =
            max(mag(L.U & n) + L.c, mag(R.U & n) + R.c)*magSf[facei];
    }

    // Boundary faces.  Empty patches have zero size and are skipped.
    forAll(rho.boundaryField(), patchi)
    {
        const fvsPatchScalarField& prhoL = rhoL.boundaryField()[patchi];
        const fvsPatchScalarField& prhoR = rhoR.boundaryField()[patchi];
        const fvsPatchVectorField& pUL = UL.boundaryField()[patchi];
        const fvsPatchVectorField& pUR = UR.boundaryField()[patchi];
        const fvsPatchScalarField& ppL = pL.boundaryField()[patchi];
        const fvsPatchScalarField& ppR = pR.boundaryField()[patchi];
        const fvsPatchScalarField& pYL = YL.boundaryField()[patchi];
        const fvsPatchScalarField& pYR = YR.boundaryField()[patchi];
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pmagSf = magSf.boundaryField()[patchi];

        fvsPatchScalarField& prhoPhi = rhoPhi.boundaryFieldRef()[patchi];
        fvsPatchVectorField& prhoUPhi = rhoUPhi.boundaryFieldRef()[patchi];
        fvsPatchScalarField& prhoEPhi = rhoEPhi.boundaryFieldRef()[patchi];
        fvsPatchScalarField& prhoYPhi = rhoYPhi.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pamaxSf = amaxSf.boundaryFieldRef()[patchi];

        const labelUList& faceCells =
            rho.boundaryField()[patchi].patch().faceCells();

        forAll(prhoL, facei)
        {
            const label celli = faceCells[facei];

            const bool badL = prhoL[facei] <= 0 || ppL[facei] <= 0;
            const bool badR = prhoR[facei] <= 0 || ppR[facei] <= 0;

            if (badL || badR) nFallback++;

            const faceState L
            (
                badL
              ? state(rhoI[celli], UI[celli], pI[celli], YI[celli])
              : state(prhoL[facei], pUL[facei], ppL[facei], pYL[facei])
            );
            const faceState R
            (
                badR
              ? state(rho.boundaryField()[patchi][facei],
                      U.boundaryField()[patchi][facei],
                      p.boundaryField()[patchi][facei],
                      Y.boundaryField()[patchi][facei])
              : state(prhoR[facei], pUR[facei], ppR[facei], pYR[facei])
            );

            const vector n(pSf[facei]/pmagSf[facei]);

            const faceFlux F
            (
                calculateFluxes(L, R, n, pmagSf[facei], facei, patchi)
            );

            prhoPhi[facei] = F.rhoPhi;
            prhoUPhi[facei] = F.rhoUPhi;
            prhoEPhi[facei] = F.rhoEPhi;
            prhoYPhi[facei] = F.rhoYPhi;

            pamaxSf[facei] =
                max(mag(L.U & n) + L.c, mag(R.U & n) + R.c)*pmagSf[facei];
        }
    }

    if (returnReduce(nFallback, sumOp<label>()) > 0)
    {
        Info<< "fluxScheme: " << returnReduce(nFallback, sumOp<label>())
            << " face reconstructions dropped to first order"
            << " (non-positive rho or p)" << endl;
    }

    postUpdate();
}


// ************************************************************************* //
