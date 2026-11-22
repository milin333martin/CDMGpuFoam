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
    Foam::CDM::fluxSchemes::HLLCLM

Description
    HLLC-LM: shock-stable HLLC with low Mach number correction of the
    acoustic dissipation (Fleischmann, Adami & Adams, J. Comput. Phys. 423
    (2020) 109762).

    Central form of HLLC (their Eq. 19) for S_L < 0 < S_R:
    \verbatim
        F* = (F_L + F_R)/2
           + [phi S_L (U*_L - U_L) + |S*| (U*_L - U*_R) + phi S_R (U*_R - U_R)]/2
        phi = sin(min(1, Ma_local/Ma_limit) pi/2),
        Ma_local = max(|u_L|/c_L, |u_R|/c_R),   Ma_limit = 0.1
    \endverbatim
    where u is the face-normal velocity.  phi = 1 (classical HLLC) wherever
    Ma_local >= 0.1; the acoustic dissipation is reduced only on faces with a
    small normal velocity, i.e. the faces transverse to a grid-aligned shock
    where the carbuncle instability originates.  S*, U*_L and U*_R use the
    unmodified wave speeds.  Wave speeds: Einfeldt estimates with the
    Roe-averaged velocity and sound speed of their Eqs. (6)-(8).

    The reactant is carried as a passive scalar in the star states,
    U*_K(rhoY) = rho_K (S_K - u_K)/(S_K - S*) Y_K.

SourceFiles
    HLLCLM.C

\*---------------------------------------------------------------------------*/

#ifndef CDM_HLLCLM_H
#define CDM_HLLCLM_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                            Class HLLCLM Declaration
\*---------------------------------------------------------------------------*/

class HLLCLM
:
    public fluxScheme
{
protected:

    // Protected Member Functions

        //- Compute the HLLC-LM fluxes
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
    TypeName("HLLCLM");


    // Constructors

        //- Construct from mesh
        HLLCLM(const fvMesh& mesh);


    //- Destructor
    virtual ~HLLCLM();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxSchemes
} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
