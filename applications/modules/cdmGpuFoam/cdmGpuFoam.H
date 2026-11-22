/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: a chemical-diffusive-model solver
   \\    /   O peration     | for gaseous detonations in OpenFOAM-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

Class
    Foam::solvers::cdmGpuFoam

Description
    cdmFoam with the right-hand side evaluated on a GPU.

    This is cdmFoam - the same class, subclassed - with one method
    replaced.  Everything else is inherited and unchanged: the same
    dictionaries, the same boundary conditions including the RDE injector and
    the conjugate temperature coupling, the same soot foil, front tracking,
    diagnostics, refinement, decomposition, function objects and output.  A
    case switches over by one line:

    \verbatim
        solver          cdmGpuFoam;
    \endverbatim

    and switches back the same way.  Nothing in a case directory has to
    change, so the two can be run against each other on identical input.

    \heading How it works

    computeRHS() uploads the conserved state and the boundary values, calls
    libCDMGpu to evaluate reconstruction, Riemann fluxes, the CDM diffusive
    terms and the divergence on the device, and downloads the cell
    right-hand sides.  The mesh is uploaded once and again after any topology
    change.

    Boundary conditions stay on the host.  OpenFOAM's own patch fields
    evaluate them exactly as they always did, and only the resulting patch
    values go to the device - which is why every boundary condition in the
    package keeps working without being ported.

    \heading What is different from cdmFoam

    \table
      Item                | CDMGpuFoam     | cdmGpuFoam
      Reconstruction      | host             | device
      Riemann solver      | host             | device
      Diffusive terms     | host             | device
      Divergence          | host             | device
      Boundary conditions | host             | host (unchanged)
      Chemistry           | host             | host (see below)
      Positivity limiter  | available        | NOT available
    \endtable

    The positivity-preserving flux limiter needs the face fluxes on the host,
    which the device form does not produce there.  The first-order
    reconstruction fallback is present and reported as usual, and is what
    guards admissibility.  The solver refuses to start if the limiter is
    switched on, rather than silently ignoring it.

    \heading Empty patches

    A two-dimensional OpenFOAM mesh is one cell thick with empty front and
    back patches: two boundary faces per cell, which for an 8.6M cell case is
    17M faces that take the device memory from 4.3 GB to 7.6 GB and do not
    belong in the discretisation at all.  They are excluded here, as OpenFOAM
    excludes them: emptyFvPatch has zero size, so building the boundary list
    from mesh.boundary() leaves them out, and the cell-to-face map skips any
    face that does not land in that list.

SourceFiles
    cdmGpuFoam.C

\*---------------------------------------------------------------------------*/

#ifndef cdmGpuFoam_H
#define cdmGpuFoam_H

#include "cdmFoam.H"

extern "C"
{
    #include "cdmGpu.h"
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                         Class cdmGpuFoam Declaration
\*---------------------------------------------------------------------------*/

class cdmGpuFoam
:
    public cdmFoam
{
    // Private Data

        //- The device-side case
        cdmGpuCase* gpu_;

        //- Mesh size the device currently holds, to detect a topology change
        label nCellsUploaded_;
        label nFacesUploaded_;

        //- Fingerprint of the uploaded mesh's face addressing.
        //  Counts alone do not identify a mesh: see cdmTopoHash.
        unsigned long long meshTopoHash_;

        //- Flattened boundary addressing, empty and coupled patches excluded.
        //  bPatch_/bFace_ give the patch and local face of each entry.
        labelList bPatch_;
        labelList bFace_;

        //- Cyclic faces, promoted to internal faces on the device and held
        //  here so the face mass flux can be written back to both sides of
        //  each pair.  cycPatch_/cycFace_ name the owner side; the device
        //  face index is mesh.nInternalFaces() + the position in these.
        labelList cycPatch_;
        labelList cycFace_;

        //- Processor faces, promoted to internal faces on the device in the
        //  same way, except that the neighbour is a ghost cell rather than
        //  a local one.  The device face index is mesh.nInternalFaces()
        //  + nCyc + the position here, and the ghost cell it reads is
        //  nCells + that position.
        //
        //  A processor face exists once physically and twice in the
        //  decomposition, so unlike a cyclic face only this rank's side is
        //  mapped back into phi_.  The two ranks still agree on the flux to
        //  the last bit: both evaluate the same Riemann problem from the
        //  same reconstructed states, and HLLC satisfies
        //  F(L,R,n) = -F(R,L,-n).
        labelList procPatch_;
        labelList procFace_;

        //- The local cell behind each processor face, in that same order:
        //  what this rank sends its neighbour every Runge-Kutta stage.
        labelList haloCell_;

        //- Where each processor patch's faces begin in the lists above,
        //  with a trailing total, and the rank on the other side of each.
        labelList procStart_;
        labelList procNbr_;

        //- Halo exchange buffers, sized once per topology change.  Twenty
        //  five doubles a face: rho, U, p, T, Y, and the gradients of rho,
        //  U, p and Y.  That is the entire stencil a face reads for its
        //  neighbour, which is why one layer of halo suffices where
        //  exchanging cell values alone would need two.
        scalarField hSend_;
        scalarField hRecv_;

        //- The same payload laid out patch by patch, which is what goes on
        //  the wire: a patch's slice is then contiguous and needs no
        //  per-exchange allocation.
        scalarField hPack_;
        scalarField hUnpack_;

        //- Number of processor faces, and so of ghost cells
        label nHalo_;

        //- Internal faces as the device counts them: the mesh's own plus
        //  one per cyclic face.  Larger than mesh.nInternalFaces(), which
        //  is why the mass flux cannot be downloaded straight into phi_.
        label nIFdev_;

        //- Receive buffer for the device's internal face mass flux
        scalarField phiScratch_;

        //- Cost of rebuilding and re-uploading the mesh after a topology
        //  change, and how often that happens.  With AMR this fires every
        //  refineInterval steps and rebuilds the CSR from scratch, which is
        //  the part a cellMap-based incremental update would remove.
        double tUploadMesh_;
        label nUploadMesh_;

        //- Total first-order fallbacks this step, for the log
        label nFallback_;

        //- Cv of the calorically perfect CDM gas.  Uniform, so one number;
        //  needed on the host to turn the device's T back into he.
        scalar Cv_;

        //- The cells behind a boundary face, de-duplicated and sorted.
        //  Within a time step these are the only cells that travel between
        //  host and device: OpenFOAM's patch fields need them to evaluate
        //  the boundary conditions, and everything else stays on the device.
        labelList bCell_;

        //- Receive buffers for those cells, sized once
        scalarField bcRho_;
        vectorField bcU_;
        scalarField bcP_;
        scalarField bcT_;
        scalarField bcY_;

        //- Wall-clock accumulators in seconds, for the time-split report
        //  written at the end of the run.  The cost is two clock reads per
        //  back-end call, so this is always on; tOutside_ is the gap between
        //  one computeRHS() returning and the next being entered, which is
        //  the host-side remainder: the Runge-Kutta updates, decode() with
        //  its boundary evaluation, the chemistry integration and the
        //  per-step diagnostics.
        double tUpState_;
        double tDecode_;
        double tUpBnd_;
        //- Named tBackRHS_, not tRHS_: cdmFoam has a tRHS_ of its
        //  own and a same-named member here shadows it, so advance()
        //  would silently credit the base class's counter to this one
        //  and the per-step split would report computeRHS as zero.
        double tBackRHS_;
        double tDownRHS_;
        double tDownPhi_;
        double tOutside_;
        double tLastExit_;
        label nProfCalls_;


    // Private Member Functions

        //- Build the flattened boundary list and the CSR cell-to-face map,
        //  and upload the mesh.  Called at construction and after any
        //  topology change.
        void uploadMesh();

        //- Upload the model constants from CDMProperties and fvSchemes
        void uploadModel();

        //- Gather the boundary primitive values and send them to the device
        void uploadBoundary();

        //- Bring the boundary-adjacent cells back from the device, let
        //  OpenFOAM's own patch fields evaluate every boundary condition
        //  from them, and send the resulting patch values back.  This is
        //  what allows the state to stay on the device across a whole step
        //  without porting a single boundary condition.
        //
        //  fetchFromDevice false skips the gather and only re-evaluates the
        //  patch fields from what the host already holds.  That is right at
        //  the top of a step, where the host's primitives came from its own
        //  decode() and the device's were just derived from the same
        //  conserved state - but the patch fields still have to be
        //  re-evaluated, because runTime has advanced and a time-dependent
        //  condition such as the RDE injector would otherwise stay a step
        //  behind.
        void exchangeBoundary(const bool fetchFromDevice = true);

        //- Fill the device's ghost cells from the neighbouring ranks.
        //
        //  Runs the cell-gradient pass, brings back the state and gradient
        //  of the cells behind each processor patch, swaps them with the
        //  rank on the other side, and uploads the result into the ghosts.
        //  Must happen before the fluxes are formed, which is why the
        //  gradient pass is a separate entry point: cdmGpuComputeRHS then
        //  recomputes those gradients but writes only real cells, so the
        //  ghosts it reads are the exchanged ones.
        //
        //  A no-op in serial, where there are no processor patches.
        //- One download / swap / upload of the ghost payload
        void haloRound();

        void exchangeHalo();

        //- Keep phi_ consistent with the device's face mass flux, for the
        //  flux-dependent boundary conditions
        void downloadMassFlux();

        //- Write the device/host time split
        void reportProfile() const;


public:

    //- Runtime type information
    TypeName("cdmGpuFoam");


    // Constructors

        //- Construct from a mesh
        cdmGpuFoam(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        cdmGpuFoam(const cdmGpuFoam&) = delete;


    //- Destructor
    virtual ~cdmGpuFoam();


    // Member Functions

        //- Evaluate the right-hand sides on the device
        virtual void computeRHS
        (
            const scalar deltaT,
            scalarField& dRho,
            vectorField& dRhoU,
            scalarField& dRhoE,
            scalarField& dRhoY
        );

        //- Run the whole time step on the device: the chemistry, the
        //  Runge-Kutta stages and the conservative-to-primitive decode,
        //  exchanging only the boundary-adjacent cells between stages.
        //
        //  The host owns the state at a step boundary and the device owns it
        //  within a step.  advance() uploads once at the top and runs the
        //  inherited decode() once at the bottom, so fvConstraints, the
        //  thermo and every diagnostic see exactly what they would under
        //  cdmFoam, and the device cannot drift from OpenFOAM's
        //  authoritative state across steps.
        virtual void advance();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cdmGpuFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
