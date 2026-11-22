/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: GPU back end
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

Description
    Plain-C interface to the GPU back end.

    Deliberately C, with no OpenFOAM types and no C++ in the ABI, so that this
    library can be compiled by a different compiler (nvc++) from the one that
    builds OpenFOAM (g++) without any risk of an ABI mismatch.  Nothing in
    libCDMGpu includes an OpenFOAM header; nothing in OpenFOAM includes a CUDA
    one.  The only thing crossing the boundary is this header.

    The back end holds the solution on the device.  OpenFOAM uploads the mesh
    once (and again after any topology change), uploads the state once at
    start-up, and from then on only exchanges boundary values each stage and
    downloads the fields at write times.  Moving the whole state across the
    bus every stage would cost more than the arithmetic saves, which is why
    this is a solver module rather than only a flux scheme.

    All arrays are the caller's, are read during the call only, and use
    OpenFOAM's own layout:

      - cell and face fields are contiguous, cell 0 first
      - vectors are interleaved xyz, so component c of entity i is at 3*i + c
      - tensors are interleaved row-major, 9 per entity
      - internal faces are indexed 0..nInternalFaces-1 in mesh.owner() order
      - boundary faces are flattened across patches into one 0..nBFaces-1 list

SourceFiles
    cdmGpuKernels.cpp

\*---------------------------------------------------------------------------*/

#ifndef CDM_GPU_H
#define CDM_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

/*- Opaque handle to one case's device state */
typedef struct cdmGpuCase cdmGpuCase;

/*- Flux schemes, matching the names in fvSchemes */
enum { CDM_GPU_HLLC = 0, CDM_GPU_HLLCLM = 1 };

/*- Limiters, matching the names in fvSchemes */
enum
{
    CDM_GPU_UPWIND = 0,
    CDM_GPU_VANLEER,
    CDM_GPU_MINMOD,
    CDM_GPU_ROUNDF,
    CDM_GPU_ROUNDA,
    CDM_GPU_ROUNDAPLUS,
    CDM_GPU_ROUNDL
};

/*- What this build can actually do.  Returns a static string such as
    "stdpar=gpu", "stdpar=multicore" or "serial", so that the solver can say
    so in its header rather than the user having to guess. */
const char* cdmGpuBackend(void);

/*- Non-zero if a usable device was found.  A build with no GPU still works;
    it runs the same kernels on the host. */
int cdmGpuAvailable(void);

/*- Monotonic wall-clock seconds since the first call.  Exposed so that the
    solver's time-split report and the back end read the same clock. */
double cdmGpuWallTime(void);

cdmGpuCase* cdmGpuCreate(void);
void cdmGpuDestroy(cdmGpuCase* c);

/*- Upload the mesh.  Call once, and again after any topology change.

    cellFaceOffset/cellFace/cellFaceSign are the cell-to-face map in CSR form:
    the faces of cell i are cellFace[cellFaceOffset[i] .. cellFaceOffset[i+1]),
    with sign +1 where the cell owns the face and -1 where it neighbours it.
    A face index below nInternalFaces is an internal face; at or above it, the
    boundary face (index - nInternalFaces).

    The CSR form exists so that gradients and divergences are pure gathers.
    The alternative - a face loop scattering into cells - needs atomics, which
    are both slower and less portable. */
void cdmGpuSetMesh
(
    cdmGpuCase* c,
    int nCells,
    int nHalo,                   /* ghost cells for processor faces; 0 in
                                    serial.  They occupy cell indices
                                    [nCells, nCells+nHalo) and are read only
                                    as the neighbour of a processor face. */
    int nInternalFaces,
    int nBoundaryFaces,
    const int* owner,            /* nInternalFaces */
    const int* neighbour,        /* nInternalFaces */
    const double* Sf,            /* 3*nInternalFaces */
    const double* magSf,         /* nInternalFaces */
    const double* weights,       /* nInternalFaces: mesh.weights(), the
                                    central-differencing weight on the owner */
    const double* deltaCoeffs,   /* nInternalFaces: mesh.nonOrthDeltaCoeffs(),
                                    so snGrad matches OpenFOAM on a skewed
                                    mesh and not only an orthogonal one */
    const double* delta,         /* 3*nInternalFaces: C[nei] - C[own] */
    const double* V,             /* nCells */
    const int* bOwner,           /* nBoundaryFaces */
    const double* bSf,           /* 3*nBoundaryFaces */
    const double* bMagSf,        /* nBoundaryFaces */
    const double* bDeltaCoeffs,  /* nBoundaryFaces: needed for the boundary
                                    viscous flux, which is NOT zero on a
                                    symmetryPlane or a wall */
    const int* cellFaceOffset,   /* nCells + 1 */
    const int* cellFace,         /* cellFaceOffset[nCells] */
    const signed char* cellFaceSign
);

/*- Model constants.  gamma, Cv and eOffset define the equation of state;
    the rest are the CDM parameters as CDMProperties gives them in SI. */
void cdmGpuSetModel
(
    cdmGpuCase* c,
    double gamma, double Cv, double eOffset,
    double A, double Ta, double q,
    double kappa0, double n, double Pr, double Le,
    double maxSubstepFraction, int maxSubsteps,
    int fluxScheme, int limiter, int viscous
);

/*- Conserved state, device-resident between calls */
void cdmGpuUploadState
(
    cdmGpuCase* c,
    const double* rho,   /* nCells */
    const double* rhoU,  /* 3*nCells */
    const double* rhoE,  /* nCells */
    const double* rhoY   /* nCells */
);

void cdmGpuDownloadState
(
    cdmGpuCase* c,
    double* rho, double* rhoU, double* rhoE, double* rhoY
);

/*- Primitive fields, for writing and for the host-side diagnostics */
void cdmGpuDownloadPrimitives
(
    cdmGpuCase* c,
    double* U, double* p, double* T, double* Y, double* Qdot
);

/*- Boundary values, evaluated by OpenFOAM's own patch fields on the host and
    uploaded each stage.  Boundary faces are a small fraction of the total, so
    this is cheap; and it means every boundary condition in the package keeps
    working unchanged, including the RDE injector and the conjugate patch. */
void cdmGpuUploadBoundary
(
    cdmGpuCase* c,
    const double* bRho,  /* nBoundaryFaces */
    const double* bU,    /* 3*nBoundaryFaces */
    const double* bP,    /* nBoundaryFaces */
    const double* bY     /* nBoundaryFaces */
);

/*- The cells that own a boundary face, de-duplicated, in the order the
    solver wants them back.  Uploaded with the mesh.

    This is what lets the state stay on the device across a whole time step.
    OpenFOAM's patch fields need the cell values behind each patch to
    evaluate a boundary condition; those cells are a thin shell - about
    15,000 of 8.6M for the RDE - so they come back each Runge-Kutta stage
    while the other 99.8% of the mesh never leaves the device. */
void cdmGpuSetBoundaryCells(cdmGpuCase* c, int nBoundaryCells, const int* cells);

/*- Primitives of those cells only, in cdmGpuSetBoundaryCells order */
void cdmGpuDownloadBoundaryCells
(
    cdmGpuCase* c,
    double* rho,  /* nBoundaryCells */
    double* U,    /* 3*nBoundaryCells */
    double* p,    /* nBoundaryCells */
    double* T,    /* nBoundaryCells */
    double* Y     /* nBoundaryCells */
);

/*- Heat release accumulates over a step, exactly as CDMModel::react does
    into heatReleased, so that Strang splitting's two half-steps both count.
    Zero it once per step and scale it to a rate at the end. */
void cdmGpuZeroQdot(cdmGpuCase* c);
void cdmGpuScaleQdot(cdmGpuCase* c, double factor);
void cdmGpuDownloadQdot(cdmGpuCase* c, double* Qdot);

/*- Save the conserved state as the Runge-Kutta stage origin */
void cdmGpuStoreStage(cdmGpuCase* c);

/*- The cells whose state and gradient this rank sends to its neighbours,
    de-duplicated, in the order the solver wants them back.  The mirror of
    cdmGpuSetBoundaryCells, for the halo rather than for patch evaluation. */
void cdmGpuSetHaloCells(cdmGpuCase* c, int nHaloCells, const int* cells);

/*- Primitives and gradients of those cells, in cdmGpuSetHaloCells order */
void cdmGpuDownloadHaloCells
(
    cdmGpuCase* c,
    double* rho, double* U, double* p, double* T, double* Y,
    double* gRho, double* gU, double* gP, double* gY
);

/*- Fill the ghosts from what the neighbouring ranks sent.  nHalo entries
    each, in the order the processor faces were given to cdmGpuSetMesh.  The
    sound speed is derived here rather than transferred, so both ranks agree
    on it exactly. */
void cdmGpuUploadHalo
(
    cdmGpuCase* c,
    const double* rho, const double* U, const double* p, const double* T,
    const double* Y,
    const double* gRho, const double* gU, const double* gP, const double* gY
);

/*- The cell-gradient pass alone.  A parallel step runs this, exchanges the
    halo, uploads it, then calls cdmGpuComputeRHS as usual: that recomputes
    these gradients but writes only indices below nCells, so the ghosts it
    then reads are the exchanged ones. */
void cdmGpuComputeGradients(cdmGpuCase* c);

/*- Evaluate the right-hand side from the current device state.  Returns the
    maximum of (|U.n| + c)*magSf over all faces, for the Courant number, and
    writes the number of face reconstructions that fell back to first order
    into nFallback. */
double cdmGpuComputeRHS(cdmGpuCase* c, int* nFallback);

/*- q = a*q0 + b*(q + dt*RHS), the SSP Runge-Kutta stage update */
void cdmGpuRKStage(cdmGpuCase* c, double dt, double a, double b);

/*- Integrate the reaction over dt at constant density and velocity */
void cdmGpuReact(cdmGpuCase* c, double dt);

/*- The face mass flux, area-integrated, so that OpenFOAM's phi stays
    consistent after a GPU right-hand side.  One double per face. */
void cdmGpuDownloadMassFlux(cdmGpuCase* c, double* phiInternal, double* phiBoundary);

/*- Cell right-hand sides from the last cdmGpuComputeRHS */
void cdmGpuDownloadRHS
(
    cdmGpuCase* c, double* dRho, double* dRhoU, double* dRhoE, double* dRhoY
);

/*- Decode conserved to primitive.  Returns 0 if every cell is admissible, or
    the index of the first cell that is not, plus one. */
int cdmGpuDecode(cdmGpuCase* c);

#ifdef __cplusplus
}
#endif

#endif

/* ************************************************************************* */
