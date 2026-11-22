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

Class
    Foam::CDM::sootFoil

Description
    In-situ numerical soot foil on a uniform lab-frame grid, for 2D runs in
    the lab frame or in the shock-attached (detonation) frame.

    Following Ramachandran, Dushe & Yang, "Comparing shock-attached and lab
    frames of reference for numerical simulations of two-dimensional channel
    detonations", AIAA SciTech 2025 (AIAA 2025-1386), Eq. (13):

    \verbatim
        U_soot(x + V t, y) = max_{t' <= t} p(x, y, t')
    \endverbatim

    where x is the computational (frame) coordinate and V the velocity of
    the computational frame relative to the lab.  In the shock frame with
    fresh gas entering at -D_CJ (detonation moving along +x in the lab),
    V = (D_CJ 0 0).  With V = 0 the foil is a lab-frame foil that, unlike
    the pMax field, is not averaged when AMR coarsens cells.

    The paper reconstructs the foil a posteriori from plot files written at
    intervals of at most dx/D_CJ; here the foil is updated every time step,
    which always satisfies that sampling condition (CFL < 1).

    Every cell is rasterised onto the pixels whose centres lie inside its
    (shifted) bounding box; cells smaller than a pixel update the pixel
    containing their centre.  Each processor accumulates a local foil; the
    foils are combined with a max-reduction at write time only.

    Output, at every write time:
        postProcessing/sootFoil/<time>/sootFoil.vtk   (ParaView image)
        postProcessing/sootFoil/<time>/sootFoil.raw   (restart data)
    On start-up the .raw file of the start time is read if present.

Usage
    In system/fvSolution:
    \verbatim
    PIMPLE
    {
        sootFoil
        {
            frameVelocity   (1967.15 0 0);  // (0 0 0) for lab-frame runs
            min             (0 0);          // lab-frame box (x y) [m]
            max             (0.135 0.01);
            delta           2e-5;           // pixel size [m]
        }
    }
    \endverbatim

SourceFiles
    sootFoil.C

\*---------------------------------------------------------------------------*/

#ifndef CDM_sootFoil_H
#define CDM_sootFoil_H

#include "volFields.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{

/*---------------------------------------------------------------------------*\
                          Class sootFoil Declaration
\*---------------------------------------------------------------------------*/

class sootFoil
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Velocity of the computational frame relative to the lab
        vector frameVelocity_;

        //- Lab-frame box
        vector2D min_;
        vector2D max_;

        //- Pixel size
        scalar delta_;

        //- Steps between samples.  update() sweeps every cell and reads
        //  the mesh's cell-to-point addressing, which a topology change
        //  invalidates, so on an adaptive mesh it is one of the most
        //  expensive things in a time step: on the RDE case at 484k cells
        //  it was 53 s of a 560 s run, sampling every step.
        //
        //  The foil records a running maximum, so the sampling only has
        //  to be fast enough that the front does not cross a pixel
        //  between samples; beyond that, extra samples add nothing.  With
        //  the front advancing 0.14 pixels per step on that case, 1 is
        //  about seven times finer than it needs to be.
        label interval_;

        //- Steps since the last sample
        label nSinceSample_;

        //- Every cell is an axis-aligned box of side sqrt(V/thickness),
        //  so its extent follows from the centre and the volume and the
        //  cell-to-point addressing is not needed at all.
        //
        //  That addressing is what made this expensive on an adaptive
        //  mesh: a topology change invalidates it and the next sample
        //  rebuilds it, which measured 78 ms per change at 484k cells -
        //  29 s of a 535 s run, and an interval does not help because the
        //  first sample after each change pays it anyway.
        //
        //  Checked once against the real point bounds, for every cell,
        //  rather than assumed; -1 means not yet checked.
        label cartesian_;

        //- 1/thickness, for the cell half-width from the volume
        scalar rThickness_;

        //- Number of pixels
        label nx_;
        label ny_;

        //- Foil values (row-major, x fastest)
        scalarField foil_;


    // Private Member Functions

        //- Directory for the current time
        fileName outputDir() const;

        //- Read the restart data of the start time if present
        void readRestart();


public:

    //- Runtime type information
    ClassName("sootFoil");


    // Constructors

        //- Construct from mesh and dictionary
        sootFoil(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        sootFoil(const sootFoil&) = delete;


    //- Destructor
    ~sootFoil();


    // Member Functions

        //- Accumulate the field (normally p) at the current time
        void update(const volScalarField& p);

        //- Combine across processors and write
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sootFoil&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
