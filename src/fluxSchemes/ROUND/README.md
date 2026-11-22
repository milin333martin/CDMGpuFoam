# ROUND reconstruction schemes (vendored)

These four files are a verbatim copy of the limiter sources from
**libROUNDSchemes**, "An OpenFOAM Library of High-resolution
Structure-preserving Convection Schemes", by **Xi Deng** and **Bin Xie**,
redistributed under the GNU General Public License version 3 or later — the
same licence as CDMGpuFoam. The upstream licence text is in
`COPYING.upstream`.

Each file carries an added comment block at the top recording this
provenance. Nothing below that block has been altered: the limiter bodies are
byte-for-byte identical to the upstream sources, verified by hash at the time
of vendoring.

If you publish results obtained with these schemes, cite the upstream work:

> Deng, X. (2023) A unified framework for non-linear reconstruction schemes in
> a compact stencil. Part 1: beyond second order. *Journal of Computational
> Physics*, **481**, 112052.

> Deng, X. and Xie, B. An OpenFOAM library of high-resolution
> structure-preserving convection schemes.

## Why they live here

`fluxScheme` reconstructs the primitive variables through OpenFOAM's limited
interpolation machinery:

```cpp
fvc::interpolate(rho, own, "reconstruct(rho)")
fvc::interpolate(rho, nei, "reconstruct(rho)")
```

ROUND registers itself into exactly that machinery — it is a `LimiterFunc` used
with `LimitedScheme`, the same template `vanLeer` uses — so the schemes become
available with no change to `fluxScheme` at all.

They are compiled into `libCDMFluxSchemes.so` rather than built as a separate
library because the solver module already links that library. A separate
`libROUNDSchemes.so` would have to be named in a `libs` entry in every
`controlDict`, and a missing `libs` entry fails at run time with an unhelpful
"unknown scheme" message rather than at build time. This is the same failure
mode that makes a missing `"libzoltanDecomp.so"` report `Unknown distributor
zoltan`.

To use these schemes from a solver other than CDMGpuFoam, add

```
libs    ("libCDMFluxSchemes.so");
```

to that case's `controlDict`, or build the upstream library standalone.

## Scheme names registered

For each of `ROUNDF`, `ROUNDA`, `ROUNDAplus`, `ROUNDL`:

| name | type | use |
|---|---|---|
| `ROUNDx` | any | general scalar fields |
| `ROUNDxV` | vector | direction-aware limiter for `U` |
| `ROUNDx01` | scalar | fields bounded to [0, 1], such as `Y` |
| `limitedROUNDx` | scalar | explicit bounds, e.g. `limitedROUNDL 0 1;` |

See `docs/11-round.md` for which one to choose and what to watch.
