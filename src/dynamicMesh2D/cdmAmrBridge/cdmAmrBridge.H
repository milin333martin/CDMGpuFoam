/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: GPU adaptive mesh engine
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

Namespace
    cdmBridge

Description
    Turns a cdmAmr key set into OpenFOAM primitives, and builds the
    polyTopoChangeMap that carries fields across the change.

    This was developed inside applications/test/cdmAmrMesh, where it is
    validated against the real checkMesh and against the engine's own
    geometry element by element.  It lives here so that the solver's
    topoChanger and that test use one copy: the test is then a test of the
    code the solver runs, not of a copy of it.

    Everything here is host code.  The engine it reads is device code; the
    translation to points and faces is not, because a faceList is a list
    of lists and resetPrimitives wants it on the host anyway.

SourceFiles
    cdmAmrBridge.C

\*---------------------------------------------------------------------------*/

#ifndef cdmAmrBridge_H
#define cdmAmrBridge_H

#include "fvMesh.H"
// fvMesh.H only forward-declares surfaceScalarField and friends, so
// .primitiveField() on one is an incomplete type without these.
#include "volFields.H"
#include "surfaceFields.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "objectMap.H"
#include "HashSet.H"
#include "cyclicPolyPatch.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"

#include "cdmAmr.h"

// At file scope: including this inside a namespace would drop all of
// __gnu_parallel, and libstdc++'s internals with it, into that namespace.
#if defined(_OPENMP) && defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#define CDM_PSORT 1
#endif

#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <vector>

namespace cdmBridge
{

using namespace Foam;

// * * * * * * * * * * * * * * * * * Sorting * * * * * * * * * * * * * * * * //

// The two sorts in buildPrimitives, about 4 corners and 1 face per cell,
// are the bulk of its cost and the largest host stage left at half a
// million cells. libstdc++'s parallel sort needs a total order: see
// operator< on FaceSort.
template<class It>
inline void psort(It first, It last)
{
#ifdef CDM_PSORT
    __gnu_parallel::sort(first, last);
#else
    std::sort(first, last);
#endif
}

template<class It, class Compare>
inline void psort(It first, It last, Compare cmp)
{
#ifdef CDM_PSORT
    __gnu_parallel::sort(first, last, cmp);
#else
    std::sort(first, last, cmp);
#endif
}

// * * * * * * * * * * * * * * * * * Timing  * * * * * * * * * * * * * * * * //

// The stage timings come back through buildPrimitives' out-parameters, so
// the clock lives beside it rather than in the calling application.
using WClock = std::chrono::steady_clock;

inline double wall()
{
    return std::chrono::duration<double>
    (
        WClock::now().time_since_epoch()
    ).count();
}



// * * * * * * * * * * * * *  Integer cell geometry  * * * * * * * * * * * * //

// Every cell corner sits on a grid of spacing Lx/(nx*2^maxLevel) by
// Ly/(ny*2^maxLevel) whatever mixture of levels is present, so corners
// shared between levels are identical by construction and points
// deduplicate by sorting exact integers rather than comparing positions.
constexpr int cMaxLevel = cdmAmr::maxLevel;
constexpr std::uint64_t cFine = std::uint64_t(1) << cMaxLevel;   // per base cell

struct Box
{
    std::uint32_t ix0, iy0, s;      // corner and size, in finest units
};

inline Box cellBoxI(cdmAmr::Key k, int nx)
{
    const std::uint64_t base = cdmAmr::keyBase(k);
    const std::uint32_t level = cdmAmr::keyLevel(k);

    std::uint32_t i, j;
    cdmAmr::unMorton2D(cdmAmr::keyMorton(k), i, j);

    const std::uint32_t s =
        static_cast<std::uint32_t>(cFine >> level);

    Box b;
    b.ix0 = static_cast<std::uint32_t>((base % nx)*cFine) + i*s;
    b.iy0 = static_cast<std::uint32_t>((base / nx)*cFine) + j*s;
    b.s = s;
    return b;
}


// A point's name: its integer position, plus which of the two planes it is
// on.  Packed so that sorting the codes sorts the points.
inline std::uint64_t pointCode
(
    std::uint64_t ix, std::uint64_t iy, std::uint64_t iz, std::uint64_t strideX
)
{
    return ((iy*strideX + ix) << 1) | iz;
}


// * * * * * * * * * * * * * * *  Face records * * * * * * * * * * * * * * * //

//- A host face's identity, by geometry rather than by which cell made it.
//
//  The engine makes each face from the finer of its two cells, so the
//  creating cell is not a stable name: when one side refines, the creator
//  switches sides, and a face the neighbour used to make is now made by
//  this cell's children.  Naming faces that way leaves ~10% of them with
//  no counterpart across a refinement.  Geometry is stable: an x-normal
//  face sits at an integer x and spans an interval in y, whoever made it.
//
//  z-normal faces are the exception.  One belongs to exactly one cell, so
//  it follows the cell map instead of a geometric query; a coarse cell's
//  footprint covers children at two different y starts, which a single
//  ordered range query would only half find.
struct FaceKey
{
    int orient;             // 0 x-normal, 1 y-normal, 2 back, 3 front
    std::uint32_t pos;      // the coordinate along the normal
    std::uint32_t a0, a1;   // the span along the face
    std::uint32_t b0;       // second coordinate, for z faces only
    int side;               // 0..5, which side of its cell it sits on
    int sign;               // +1 if Sf points along the + axis
    int cell;               // owner, used for z faces
};


//- A face's identity packed into one integer, so sorting is a sort of
//  plain integers rather than an indirect comparison through FaceKey.
//  The ordering is (orient, pos, b0, a0), so the high bits are the "same
//  plane" prefix and the low 20 the position along it.
//
//  20 bits holds a coordinate up to 1048575 and the finest coordinate is
//  nx*2^maxLevel, so this is good for a base mesh up to 4095 cells wide.
//  Checked at run time rather than assumed.
struct FaceSort
{
    std::uint64_t key;
    label idx;
};

// A total order, not just on the key. Two faces should never share a key,
// but comparing only the key leaves equal elements in an order that
// depends on the sort, and the parallel sort below is not stable; idx
// costs nothing and makes the face numbering reproducible.
inline bool operator<(const FaceSort& a, const FaceSort& b)
{
    return a.key != b.key ? a.key < b.key : a.idx < b.idx;
}


inline std::uint64_t packFaceKey(const FaceKey& k)
{
    return (static_cast<std::uint64_t>(k.orient) << 60)
         | (static_cast<std::uint64_t>(k.pos)    << 40)
         | (static_cast<std::uint64_t>(k.b0)     << 20)
         |  static_cast<std::uint64_t>(k.a0);
}


struct FaceRec
{
    label patch;        // -1 internal, otherwise the patch index
    label own;
    label nei;          // -1 on a boundary face
    label skey;         // sort key within a boundary patch
    label eface;        // the engine face this came from, or -1
    FaceKey id;         // geometric identity, for the topology map
    label v[4];
};

// OpenFOAM wants internal faces first in upper-triangular order (by owner,
// then by neighbour), then boundary faces grouped by patch.
//
// skey is what makes the cyclic halves line up. A cyclic patch requires
// face i on one side to transform onto face i on the other, and the two
// sides have different owners, so owner order says nothing about the
// pairing. Sorting both halves by position along the seam does.
inline bool faceLess(const FaceRec& a, const FaceRec& b)
{
    if (a.patch != b.patch) { return a.patch < b.patch; }
    if (a.patch < 0)
    {
        if (a.own != b.own) { return a.own < b.own; }
        return a.nei < b.nei;
    }
    if (a.skey != b.skey)   { return a.skey < b.skey; }
    return a.own < b.own;
}

// * * * * * * * * * * * * * * * * Interface  * * * * * * * * * * * * * * * //

// Key set and its addressing -> points, faces, owner, neighbour and the
// patch slices, plus the two canonical-name arrays (faceKey, pointCodeOut)
// that the next call to buildTopoMaps joins against.
void buildPrimitives
(
    const cdmAmr::Mesh& m,
    const labelList& sidePatch,
    const label frontBack,
    const label cyclicMin,
    const label cyclicMax,
    const label nPatches,
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& faceEngine,
    std::vector<FaceKey>& faceKey,
    std::vector<std::uint64_t>& pointCodeOut,
    double& tPoints,
    double& tFaces,
    double& tSort
)
;

// The maps that carry cells, points and faces from one key set to the
// next. Every entity has a canonical name on both sides, so each map is
// a join rather than a geometric search.
void buildTopoMaps
(
    const cdmAmr::Mesh& m,
    const std::vector<cdmAmr::Key>& oldKeys,
    const std::vector<std::uint64_t>& oldCode,
    const std::vector<std::uint64_t>& newCode,
    const std::vector<FaceKey>& oldFace,
    const std::vector<FaceKey>& newFace,
    const std::vector<FaceSort>& oldSorted,
    std::vector<FaceSort>& newSorted,
    labelList& cellMap,
    List<objectMap>& cellsFromCells,
    labelList& pointMap,
    labelList& faceMap,
    List<objectMap>& facesFromFaces,
    labelList& revCell,
    labelList& revFace,
    labelList& revPoint,
    labelHashSet& flipFaceFlux
)
;

// Which patch each side of the base mesh belongs to, worked out from
// where the boundary faces actually are rather than from patch names.
void deducePatchLayout
(
    const fvMesh& mesh,
    labelList& sidePatch,       // [4], patch per in-plane side, or -1
    label& frontBack,
    label& cyclicMin,
    label& cyclicMax
)
;

} // End namespace cdmBridge

#endif

// ************************************************************************* //
