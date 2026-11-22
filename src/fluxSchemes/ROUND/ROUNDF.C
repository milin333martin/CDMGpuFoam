/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: vendored third-party scheme
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Provenance
    This file is part of libROUNDSchemes, "An OpenFOAM Library of
    High-resolution Structure-preserving Convection Schemes", by Xi Deng and
    Bin Xie, and is redistributed here under the GNU General Public License
    version 3 or later, the same licence as CDMGpuFoam itself.

    The limiter body below is a verbatim copy of the upstream source.  Only
    this comment block has been added.  Report numerical questions about the
    scheme to the upstream authors, not to CDMGpuFoam.

    If you use these schemes, cite:

      Deng, X. (2023) A unified framework for non-linear reconstruction
      schemes in a compact stencil. Part 1: beyond second order.
      Journal of Computational Physics, 481, 112052.

      Deng, X. and Xie, B. An OpenFOAM library of high-resolution
      structure-preserving convection schemes.

    Vendored into src/fluxSchemes/ROUND so that the schemes register when
    libCDMFluxSchemes.so loads, which the solver module already links.  No
    "libs" entry in controlDict is needed.  See docs/11-round.md.

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "ROUNDF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(ROUNDF, ROUNDFLimiter)
    makeLimitedVSurfaceInterpolationScheme(ROUNDFV, ROUNDFLimiter)
    
        makeLLimitedSurfaceInterpolationTypeScheme
    (
        limitedROUNDF,
        LimitedLimiter,
        ROUNDFLimiter,
        NVDTVD,
        magSqr,
        scalar
    )

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        ROUNDF01,
        Limited01Limiter,
        ROUNDFLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
