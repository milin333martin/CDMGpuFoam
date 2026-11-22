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

Class
    Foam::CDM::fluxSchemes::HLLC

Description
    HLLC approximate Riemann solver (Toro, Spruce & Speares 1994) with
    Einfeldt/Roe-averaged wave-speed estimates.  The reactant is advected
    with the HLLC mass flux, upwinded on the side of the contact.

    Reference:
    \verbatim
        Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid
        Dynamics, 3rd ed., Springer, Ch. 10.
    \endverbatim

SourceFiles
    HLLC.C

\*---------------------------------------------------------------------------*/

#ifndef CDM_HLLC_H
#define CDM_HLLC_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                            Class HLLC Declaration
\*---------------------------------------------------------------------------*/

class HLLC
:
    public fluxScheme
{
protected:

    // Protected Member Functions

        //- Compute the HLLC fluxes
        virtual faceFlux calculateFluxes
        (
            const faceState& L,
            const faceState& R,
            const vector& n,
            const scalar magSf,
            const label facei,
            const label patchi
        ) const;


public:

    //- Runtime type information
    TypeName("HLLC");


    // Constructors

        //- Construct from mesh
        HLLC(const fvMesh& mesh);


    //- Destructor
    virtual ~HLLC();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxSchemes
} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
