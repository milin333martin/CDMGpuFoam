# CDMGpuFoam

A density-based solver module for **OpenFOAM-13** (openfoam.org) for unsteady,
multidimensional, reactive, compressible Navier–Stokes flow with the
**chemical-diffusive model (CDM)** of Kaplan, Özgen & Oran (2019), as used by
Lu, Kaplan & Oran, *Combust. Flame* 230 (2021) 111417.

It is built for gaseous detonations: cellular structure, propagation through
obstacles and narrow channels, and flame arrestors, with 2D adaptive mesh
refinement, dynamic load balancing and optional conjugate heat transfer to the
solid walls.

| | |
|---|---|
| Physics | single-step CDM chemistry, power-law transport, calorically perfect gas |
| Numerics | MUSCL (van Leer) + HLLC / HLLC-LM / HLLC-P, SSP-RK2, Strang-split chemistry, positivity-preserving flux limiter |
| Meshes | 2D AMR (`refiner2D`/`hexRef2D`), run-time load balancing (zoltan, ptscotch or hierarchical) |
| Multi-region | conjugate heat transfer to solid walls through `foamMultiRun` |
| Diagnostics | in-situ soot foil, shock and flame tracking (global and per channel), conservation totals |
| Engines | plenum-fed slot injector for rotating detonation engines, cyclic circumference |

---

## 1. Requirements

- **OpenFOAM-13** from openfoam.org, with its environment sourced
  (`source /opt/openfoam13/etc/bashrc`, or your `of13` alias).
- A C++14 compiler (whatever OpenFOAM-13 was built with) and `python3`.
  **Build outside conda**: an active conda environment puts its own linker
  first, which cannot link against a system-built OpenFOAM
  (`undefined reference to hypot@GLIBC_2.xx`). `Allwmake` stops with that
  advice if it detects one. `conda deactivate` first; conda is fine for the
  Python tools afterwards.
- Optional: `numpy`, `scipy` for `tools/CDM_ZND.py` and the 1D prototypes;
  `matplotlib` for `tools/plotFronts.py`; ParaView for the fields and foils.
- Optional for parallel graph load balancing: `libzoltanDecomp.so` or
  `libptscotchDecomp.so` in your OpenFOAM installation (check with
  `ls $FOAM_LIBBIN/$FOAM_MPI | grep -E "zoltan|ptscotch"`). Without them, use
  `DISTRIBUTOR=hierarchical`.

## 2. Install

```sh
tar xzf CDMGpuFoam-1.4.0.tar.gz      # or unzip the .zip
cd CDMGpuFoam-1.4.0
source /opt/openfoam13/etc/bashrc      # your OpenFOAM-13
./Allwmake                             # add -j for a parallel build
```

This installs four libraries in `$FOAM_USER_LIBBIN` and three utilities in
`$FOAM_USER_APPBIN`:

```
libCDMModels  libCDMFluxSchemes  libCDMDynamicMesh2D  libcdmFoam
setZNDProfile  mapDetonation  decomposeRefinementHistory2D
```

`foamRun` loads the solver when a case's `controlDict` has
`solver cdmFoam;`. `./Allwclean` removes the build products.

Check the installation:

```sh
./Alltest            # ~5 minutes: a 1D detonation and a conjugate heat transfer slab
```

It reports the detonation speed against the CJ value and the conjugate
interface temperature against its analytic value, and prints PASS or FAIL for
each.

## 3. Quickstart

### A detonation in a channel (the usual starting point)

```sh
cd tutorials/cellular2D_T6_AMR
./Allrun                               # 2D H2-air cells, shock frame, AMR
```

Watch it with `tail -f log.foamRun`; the soot foil appears in
`postProcessing/sootFoil/<time>/sootFoil.vtk`.

### The whole study: a detonation, then a geometry

The detonation is grown in a frame that moves with the wave (small domain,
cheap), then transplanted into the geometry you care about, which must be in
the laboratory frame. One configuration file per case:

```sh
cd examples

../tools/runSource.sh source-T6.cfg          # 1. wave-fixed detonation
../tools/runSource.sh source-T6.cfg -s       #    progress and soot foil
../tools/runSource.sh source-T6.cfg -F -a arrestor-BR09.cfg
                                             # 2. measure the front and write
                                             #    the snapshot into the config
../tools/runArrestor.sh arrestor-BR09.cfg    # 3. the arrestor, lab frame
```

or all of it in one command, each stage waiting for the last:

```sh
../tools/runPipeline.sh source-T6.cfg arrestor-BR09.cfg arrestor-BR09-CHT.cfg
```

`docs/09-workflow.md` explains each stage and how to tell when the detonation
is quasi-steady enough to transplant.

`runArrestor.sh` does every step: mesh, region split and interface patches
(with `CHT=yes`), refinement zones, initial fields, the transplanted
detonation, decomposition, and the run itself. A copy of the configuration is
kept in the case as `case.cfg`.

Ready-made configurations in `examples/`:

| File | What it is |
|---|---|
| `source-T6.cfg` | the wave-fixed detonation that initialises the rest |
| `arrestor.cfg` | the arrestor template, every setting documented |
| `passage-single.cfg` | one 0.4 mm passage in a 4.5 cm channel |
| `arrestor-BR09.cfg` | 11 channels of 0.4 mm, blockage ratio 0.9 |
| `arrestor-BR09-CHT.cfg` | the same with conjugate heat transfer to steel plates |

### Looking at the results

```sh
./tools/plotFronts.py myCase --channels --arrestor 0.12 0.13   # shock, flame, separation
./tools/channels.sh myCase                                     # per-channel penetration
paraview myCase/postProcessing/sootFoil/<time>/sootFoil.vtk    # the soot foil
```

`docs/07-postprocessing.md` covers reconstruction, ParaView, movies, and the
traps involved.

## 4. Where things are

```
Allwmake Allwclean Alltest     build, clean, verify
package.sh clean-workspace.sh  make a release; tidy a working copy
src/CDMModels                  CDM parameters, transport, chemistry, soot foil,
                               CDMCoupledTemperature (the CHT patch field),
                               rdeInlet (the RDE injector patch fields)
src/fluxSchemes                HLLC, HLLC-LM, HLLC-P, positivity-safe face states
src/dynamicMesh2D              refiner2D, hexRef2D, refinementHistory2D
applications/modules/cdmFoam     the solver module
applications/utilities         setZNDProfile, mapDetonation,
                               decomposeRefinementHistory2D
tools/                         runSource.sh, runArrestor.sh, runPipeline.sh,
                               generators, plotting, CJ/ZND calculator
examples/                      ready-to-edit case configurations
tutorials/                     verification and demonstration cases
docs/                          the manual (PDF + LaTeX sources) and the
                               Markdown notes below
```

| Document | Contents |
|---|---|
| `docs/01-model.md` | the CDM equations, parameters, and how to set up a case by hand |
| `docs/02-numerics.md` | Riemann solvers, the positivity limiter, the soot foil |
| `docs/03-amr-parallel.md` | 2D AMR, refinement zones, load balancing, parallel running |
| `docs/04-transplant.md` | taking a quasi-steady detonation into another geometry |
| `docs/05-arrestor.md` | arrestor geometry, configurations, front tracking |
| `docs/06-cht.md` | conjugate heat transfer, including the arrestor case |
| `docs/07-postprocessing.md` | reconstruction, ParaView, movies, soot foils |
| `docs/08-verification.md` | what has been checked, and against what |
| `docs/09-workflow.md` | the end-to-end study: wave-fixed run, snapshot, geometry |
| `docs/10-rde.md` | rotating detonation engines: the slot injector and the unrolled annulus |
| `docs/11-round.md` | the ROUND high-resolution reconstruction schemes: which to use, and what to watch |
| `docs/detonation-transplant.html` | an illustrated, interactive walkthrough of the workflow |
| `docs/CDMGpuFoam-manual.pdf` | **the full manual**: theory, numerics, every dictionary key, all six tutorials, verification status |

## 5. Tutorials

| Case | What it checks |
|---|---|
| `ZND1D_stable` | 1D detonation on a stable gas: D within 0.1% of D_CJ |
| `AMR2D_blast` | inert blast with AMR and load balancing: conservation exact |
| `AMR2D_ZNDstrip` | detonation with the refined band following the front, serial and parallel |
| `cellular2D_H2air` | 2D cells, Table 1 of Lu et al. |
| `cellular2D_T6_AMR` | 2D cells, Table 6 (λ-calibrated): the usual source for a transplant |
| `quench2D_T6_AMR` | quenching in a single 0.4 mm passage |
| `RDE2D_T6` | rotating detonation engine, unrolled annulus with a slot injector |

## 6. References

- X. Lu, C. R. Kaplan, E. S. Oran, *Combust. Flame* **230** (2021) 111417.
- C. R. Kaplan, A. Özgen, E. S. Oran, *Combust. Theory Model.* **23** (2019) 67.
- S. Ramachandran, R. Dushe, S. Yang, AIAA 2025-1386 (shock-attached vs lab frames).
- Q. Xie et al., *Int. J. Numer. Methods Fluids* **89** (2018) 430 (HLLC-P).
- S. Fleischmann, S. Adami, N. A. Adams, *J. Comput. Phys.* **423** (2020) 109762 (HLLC-LM).
- X. Y. Hu, N. A. Adams, C.-W. Shu, *J. Comput. Phys.* **242** (2013) 169 (positivity limiter).
- H. Sun et al., detonationFoam, *Comput. Phys. Commun.* (2023).

License: GPL v3 (see `COPYING`). Derived from OpenFOAM and from detonationFoam.
