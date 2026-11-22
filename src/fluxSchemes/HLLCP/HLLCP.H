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

Class
    Foam::CDM::fluxSchemes::HLLCP

Description
    HLLC-P: HLLC with a pressure-based shock-instability (carbuncle) fix and
    a low-Mach pressure correction, as in detonationFoam.

    A pressure sensor f = (min over the stencil of p_i/p_j)^3 is ~1 in smooth
    flow and <1 near strong pressure jumps.  Where f < 1 an anti-diffusive
    mass flux proportional to (1-f)(pR - pL) is added and the momentum
    interface pressure is blended, which suppresses the odd-even decoupling
    that pollutes cellular detonation structures.

    Differences to detonationFoam:
      - the sensor is computed face-wise and includes coupled (processor,
        cyclic) patches, so results do not depend on the decomposition;
      - on boundary patches the sensor is read from the patch field (the
        original indexed the internal field with the patch-local face index).

    Reference:
    \verbatim
        Xie, W., Zheng, R., Lai, J., Li, H. (2018).
        An accurate and robust HLLC-type Riemann solver for the compressible
        Euler system at various Mach numbers.
        Int. J. Numer. Methods Fluids, 89, 430-463.
    \endverbatim

SourceFiles
    HLLCP.C

\*---------------------------------------------------------------------------*/

#ifndef CDM_HLLCP_H
#define CDM_HLLCP_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
namespace fluxSchemes
{

/*---------------------------------------------------------------------------*\
                            Class HLLCP Declaration
\*---------------------------------------------------------------------------*/

class HLLCP
:
    public fluxScheme
{
    // Private Data

        //- Face pressure sensor, alive only during update()
        tmp<surfaceScalarField> f_;


protected:

    // Protected Member Functions

        //- Compute the pressure sensor
        virtual void preUpdate(const volScalarField& p);

        //- Clear the pressure sensor
        virtual void postUpdate();

        //- Compute the HLLC-P fluxes
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
    TypeName("HLLCP");


    // Constructors

        //- Construct from mesh
        HLLCP(const fvMesh& mesh);


    //- Destructor
    virtual ~HLLCP();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxSchemes
} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
