/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM and from the
    flux schemes of detonationFoam (Sun et al., Comput. Phys. Commun. 2023),
    themselves based on blastFoam (Synthetik Applied Technologies).

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

Class
    Foam::CDM::fluxScheme

Description
    Base class for approximate Riemann solvers computing the inviscid fluxes
    of mass, momentum, total energy and reactant mass.

    Left (owner) and right (neighbour) states are obtained by MUSCL
    reconstruction using the interpolation schemes

    \verbatim
        interpolationSchemes
        {
            default             linear;
            reconstruct(rho)    vanLeer;
            reconstruct(U)      vanLeerV;
            reconstruct(p)      vanLeer;
            reconstruct(Y)      vanLeer;
        }
    \endverbatim

    The scheme is selected by the fvSchemes entry

    \verbatim
        fluxScheme  HLLCLM;     // HLLC | HLLCLM | HLLCP (experimental)
    \endverbatim

    Primitive variables (rho, U, p, Y) are reconstructed; the internal
    energy and speed of sound at the faces follow from the calorically
    perfect equation of state of the CDM gas,
        e = p/(rho (gamma - 1)) + eOffset,   c = sqrt(gamma p/rho),
    so that the face states are thermodynamically consistent.  If
    reconstruct(p) is not given, reconstruct(T) is used for p.

    Any OpenFOAM limited interpolation scheme may be named, including the
    vendored ROUND schemes of Deng (2023), which reach beyond second order
    in the same compact stencil:

    \verbatim
        interpolationSchemes
        {
            default             linear;
            reconstruct(rho)    ROUNDAplus;
            reconstruct(U)      ROUNDAplusV;
            reconstruct(p)      ROUNDAplus;
            reconstruct(Y)      ROUNDAplus01;
        }
    \endverbatim

    See src/fluxSchemes/ROUND/README.md and docs/11-round.md.  Note that
    ROUNDL is compressive by design and reconstructs past the neighbour
    value for some gradient ratios; watch the first-order fallback count
    reported below if you select it.

    The object holds no persistent registered fields so that it is safe
    under mesh topology changes (refinement) and redistribution.

SourceFiles
    fluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CDM_fluxScheme_H
#define CDM_fluxScheme_H

#include "volFields.H"
#include "surfaceFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{

/*---------------------------------------------------------------------------*\
                          Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme
{
public:

    //- Reconstructed state on one side of a face
    struct faceState
    {
        scalar rho;
        vector U;
        scalar e;
        scalar p;
        scalar c;
        scalar Y;
    };

    //- Fluxes through a face (already multiplied by the face area)
    struct faceFlux
    {
        scalar rhoPhi;
        vector rhoUPhi;
        scalar rhoEPhi;
        scalar rhoYPhi;
    };


protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Interpolation scheme name used for p
        word pSchemeName_;


    // Protected Member Functions

        //- Hook called once per update before the face loop
        virtual void preUpdate(const volScalarField& p)
        {}

        //- Hook called once per update after the face loop
        virtual void postUpdate()
        {}

        //- Compute the fluxes through a face with unit normal n and area
        //  magSf.  patchi = -1 for internal faces.
        virtual faceFlux calculateFluxes
        (
            const faceState& L,
            const faceState& R,
            const vector& n,
            const scalar magSf,
            const label facei,
            const label patchi
        ) const = 0;


public:

    //- Runtime type information
    TypeName("fluxScheme");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            fvMesh,
            (const fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from mesh
        fluxScheme(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        fluxScheme(const fluxScheme&) = delete;


    // Selectors

        //- Select from the fvSchemes entry fluxScheme
        static autoPtr<fluxScheme> New(const fvMesh& mesh);


    //- Destructor
    virtual ~fluxScheme();


    // Member Functions

        //- Compute the inviscid face fluxes and the maximum wave speed
        //  times face area (for the Courant number).
        //  Virtual so that a back end can replace the whole pipeline -
        //  reconstruction included - and not merely the per-face arithmetic.
        //  See src/gpu.
        virtual void update
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
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluxScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
