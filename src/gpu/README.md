# libCDMGpu — the GPU back end

The CDM right-hand side evaluated on OpenFOAM's own unstructured addressing,
on a GPU, on every CPU core, or serially, from one source.

**Status: the kernels are written and validated. The OpenFOAM glue is not
written yet.** See *What is still missing* below. Nothing in CDMGpuFoam has
changed; this directory is additive and the existing solver is untouched.

## Why it is built separately

This is the one piece that may be compiled by `nvc++` rather than by the
compiler OpenFOAM itself was built with. The interface in `cdmGpu.h` is plain
C precisely so that the two can never disagree about an ABI: nothing in
`libCDMGpu` includes an OpenFOAM header, and nothing in OpenFOAM will include
a CUDA one.

It is written as C++17 parallel algorithms rather than CUDA so that one source
serves every target — a workstation card, a cluster card, a machine with no
card at all. `-gpu=ccall` emits code for every compute capability the SDK
supports plus PTX, so the same `.so` is portable across clusters.

```sh
make            # serial, system compiler.  Builds anywhere.
make parallel   # every CPU core, g++ + Intel TBB (libtbb-dev)
make multicore  # every CPU core, nvc++ -stdpar=multicore
make gpu        # one GPU,        nvc++ -stdpar=gpu
make test       # the ZND validation, on whichever build is present
make install    # -> $FOAM_USER_LIBBIN/libCDMGpu.so
```

`nvc++` comes with the **NVIDIA HPC SDK**, a separate free download that is
*not* part of the CUDA toolkit: <https://developer.nvidia.com/hpc-sdk>

## Validation

`make test` builds a one-dimensional mesh in OpenFOAM's own addressing —
`owner`, `neighbour`, `Sf`, `magSf`, `weights`, `delta`, and the CSR
cell-to-face map — and runs `ZND1D_stable` through the back end. The
detonation speed must come out at *D*<sub>CJ</sub> = 1967.15 m/s, the same
number `tutorials/ZND1D_stable` requires of CDMGpuFoam.

```
  leading-shock speed by window
    t = 3.813e-07 .. 8.684e-05 s   D = 1512.79 m/s
    t = 8.714e-05 .. 1.632e-04 s   D = 2036.91 m/s
    t = 1.635e-04 .. 2.416e-04 s   D = 1982.49 m/s
    t = 2.419e-04 .. 3.207e-04 s   D = 1973.70 m/s
    t = 3.210e-04 .. 4.000e-04 s   D = 1969.38 m/s

  steps              24260
  first-order falls  0
  measured D         1969.38 m/s     error +0.113 %      VALIDATION PASSED
```

Three independent implementations now agree on this number:

| implementation | addressing | *D* |
|---|---|---|
| 1D prototype | structured, 1D | 1969.23 m/s |
| 2D structured solver | structured, 2D | 1969.24 m/s |
| **this back end** | **OpenFOAM unstructured** | **1969.38 m/s** |

The third is the one that matters, because the unstructured face loop — the
CSR gradient gather, the NVDTVD/NVDVTVDV gradient-projection form of *r*, the
divergence gather — has no counterpart in a structured code and is therefore
where a mistake would hide.

## What the back end does and does not own

It owns the solution between calls: the conserved state lives in its arrays
and OpenFOAM only exchanges boundary values each stage and downloads fields at
write times. Moving the whole state across the bus every stage would cost more
than the arithmetic saves — which is why this has to be driven by a solver
module, not merely a flux scheme.

It does **not** own the boundary conditions. OpenFOAM's own patch fields
evaluate them on the host and upload the resulting patch values, so every
boundary condition in the package keeps working unchanged, including the RDE
injector and the conjugate temperature coupling.

## What is still missing

1. `fluxScheme::update()` must become `virtual` (`fluxScheme.H:206`). The
   expensive part of the RHS happens in that non-virtual base method, so
   nothing can replace it today.
2. `cdmFoam::computeRHS` and `react` must become `protected virtual`.
3. `applications/modules/cdmGpuFoam/` — the solver module that uploads the
   mesh, drives the back end, and falls back to the host path if no device is
   present.
4. The `Allwmake` step that builds this library and links it.

## Known limitations of this version

- **`snGrad` uses the orthogonal form** `1/|d|`, with no non-orthogonal
  correction. Exact on the hex meshes the RDE and cellular cases use; it would
  need the correction term for a skewed mesh.
- **No MPI.** One process, one GPU, the whole mesh — which is what an 8.64M
  cell case wants anyway. Processor patches would need either a device-to-host
  round trip each stage or GPU-aware MPI.
- **AMR re-uploads the mesh** on every topology change. It works, but with
  `refineInterval 2` the uploads would dominate; GPU runs should use a
  uniform mesh.
- **No positivity-preserving flux limiter.** The first-order reconstruction
  fallback is present and reported, as in cdmFoam; the Hu/Adams/Shu flux
  limiter is not.
