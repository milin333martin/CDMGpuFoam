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

Description
    A forest-of-quadtrees mesh whose topology is implicit in a sorted array
    of Morton keys.

    The motivation is measured.  In an RDE run with refineInterval 2, the
    refinement machinery was 74% of the runtime, and it divided as

        fvMesh::topoChange       41.5%   geometry rebuild and field mapping
        polyTopoChange           32.6%   renumbering points, faces, cells
        hexRef2D::setRefinement  22.7%   building the topology change
        everything else           3.2%

    Only the first of those is arithmetic.  The other 55% is bookkeeping:
    allocating entities, renumbering them and threading them onto lists.
    Porting that to a device would be a large amount of work for a data
    structure that does not suit one.

    The observation this file rests on is that cell-wise octree refinement
    does not need general topology change at all.  A refinement is always
    one cell into four, with known geometric relationships, so a cell can be
    named rather than stored:

        key = (base cell, level, position within that base cell)

    Refining replaces a key with its four children; coarsening replaces four
    siblings with their parent; and the neighbour of a cell is found by
    arithmetic on the key followed by a binary search. Nothing is allocated,
    renumbered or threaded. Every operation here is a map, a scan or a sort,
    which is what a device is for.

    The base mesh is assumed structured: nx by ny quadrilaterals, optionally
    cyclic in x, as blockMesh produces for the RDE and cellular cases.  The
    generalisation to an arbitrary base mesh needs a per-base-face table of
    neighbour and orientation; the engine is written so that is the only
    part that would change.

SourceFiles
    cdmAmr.cpp

\*---------------------------------------------------------------------------*/

#ifndef CDM_AMR_H
#define CDM_AMR_H

#include <cstdint>
#include <string>
#include <vector>

namespace cdmAmr
{

//- A cell's name.  Bits 0..3 are the level, bits 4..19 the Morton code of
//  the cell's (i, j) position within its base cell, left-shifted to a
//  common finest level, and bits 20 and up the base cell.
//
//  Putting z above level is what makes this a linear octree rather than a
//  list sorted by level.  A cell at level L with Morton code m owns the z
//  interval [m << 2(maxLevel-L), that + 2^2(maxLevel-L)), and its four
//  children tile that same interval.  So the children sort exactly where
//  their parent was: refining is an in-place replacement that leaves the
//  array sorted, with no sort to pay for, and the leaf covering a point is
//  the last one whose z does not exceed it - one binary search instead of
//  a walk up the levels.
using Key = std::uint64_t;

constexpr int maxLevel = 8;

inline Key makeKey(std::uint64_t base, std::uint32_t level, std::uint32_t morton)
{
    const std::uint64_t z =
        static_cast<std::uint64_t>(morton) << (2*(maxLevel - level));
    return (base << 20) | (z << 4) | level;
}

inline std::uint64_t keyBase(Key k)  { return k >> 20; }
inline std::uint32_t keyLevel(Key k) { return static_cast<std::uint32_t>(k & 0xFu); }
inline std::uint32_t keyZ(Key k)     { return static_cast<std::uint32_t>((k >> 4) & 0xFFFFu); }

inline std::uint32_t keyMorton(Key k)
{
    return keyZ(k) >> (2*(maxLevel - keyLevel(k)));
}

//- The width of a cell's z interval
inline std::uint32_t keyZSpan(Key k)
{
    return 1u << (2*(maxLevel - keyLevel(k)));
}

//- Interleave the low bits of i and j, j in the odd positions
std::uint32_t morton2D(std::uint32_t i, std::uint32_t j);

//- The inverse
void unMorton2D(std::uint32_t m, std::uint32_t& i, std::uint32_t& j);


/*---------------------------------------------------------------------------*\
                            Class Mesh Declaration
\*---------------------------------------------------------------------------*/

//- A structured base mesh of nx by ny quadrilaterals over [0,Lx] x [0,Ly],
//  with the refinement state held as a sorted key array.
class Mesh
{
public:

    // Base mesh

        int nx = 0;
        int ny = 0;
        double Lx = 0;
        double Ly = 0;
        double thickness = 1;

        //- Cyclic in x, as an unrolled annulus is
        bool cyclicX = false;

    // Refinement state

        //- Sorted, unique, and covering the base mesh exactly once
        std::vector<Key> keys;

    // Addressing, rebuilt from the keys

        std::vector<int> owner;
        std::vector<int> neighbour;

        //- The cell that made each face and the side of that cell it sits
        //  on, 0 = -x, 1 = +x, 2 = -y, 3 = +y.  A face is always made by
        //  the finer of the two cells, so this fixes the face's extent
        //  exactly: it is one whole side of that cell.  owner/neighbour
        //  cannot say this, because enforcing owner < neighbour swaps the
        //  pair and loses which side the face came from.  Nothing can
        //  build the face's vertices without it.
        std::vector<int> faceCell;
        std::vector<int> faceSide;

        //- Outward from the owner, times the face area
        std::vector<double> Sf;      // 2 per face, the mesh being planar
        std::vector<double> magSf;

        //- Interpolation weight on the owner
        std::vector<double> weights;

        //- 1/|d| with d from owner to neighbour centre
        std::vector<double> deltaCoeffs;

        std::vector<double> V;
        std::vector<double> C;       // 2 per cell

        //- Faces on the outer boundary, by the side they sit on
        std::vector<int> bCell;
        std::vector<int> bSide;      // 0 = -x, 1 = +x, 2 = -y, 3 = +y


    // Construction

        //- Declared out of line, and defined in cdmAmr.cpp, so that the
        //  member vectors are allocated and freed by the same binary.
        //
        //  This matters because libCDMAmr.so is linked into OpenFOAM
        //  alongside libCDMGpu.so, which nvc++ builds with -stdpar=gpu
        //  and its own allocator.  With the destructor inline, the
        //  vectors were allocated by code in libCDMAmr and freed by code
        //  in whatever library held the Mesh.  It does not remove the
        //  one that mattered, which releaseWithoutFreeing describes.
        Mesh();
        ~Mesh();

        //- Hand every buffer to a holder that is never deleted, leaving
        //  this Mesh empty.
        //
        //  For process teardown only.  nvc++ -stdpar=gpu places the
        //  process heap under CUDA managed memory, and in a process that
        //  has loaded libCDMGpu the *last* free of these buffers is
        //  rejected: glibc aborts with "free(): invalid pointer", or
        //  "munmap_chunk(): invalid pointer", after the run has already
        //  printed "End".
        //
        //  Every free during the run succeeds.  A production run made 369
        //  mesh changes, reallocating these arrays as the cell count grew
        //  to 484,245, and reached the end of the time loop cleanly; only
        //  the teardown aborted.  So the fault is in the final release and
        //  not in the allocator pairing, which was the first four guesses.
        //  Emptying these thirteen vectors in the debugger made that same
        //  run exit normally, and the four large std::vectors keyRefiner
        //  frees moments later were never affected - which places the
        //  whole of the problem in this object, and nowhere else.
        //
        //  The memory is not recovered, deliberately.  The only caller is
        //  the destructor of an object whose process is on its way out,
        //  and the kernel reclaims the pages regardless.  Calling this at
        //  any other point would be a real leak.
        void releaseWithoutFreeing();

        //- Every base cell unrefined
        void initialise(int nx, int ny, double Lx, double Ly, bool cyclicX);


    // Topology

        //- Replace each flagged key with its four children.  Returns the
        //  number of cells refined.
        //
        //  If newIndex is given it is filled with, for each old cell, the
        //  position in the new array of that cell or of the first of its
        //  four children.  A caller that holds a second per-cell flag -
        //  a coarsening flag, say - needs this to carry that flag across,
        //  because refine() shifts every index after the first refined
        //  cell.  Deriving the shift from the flag array instead means
        //  duplicating the level guard above, and the two copies then
        //  have to be kept in step; this keeps it in one place.
        int refine
        (
            const std::vector<char>& flag,
            std::vector<int>* newIndex = nullptr
        );

        //- Replace complete flagged sibling groups with their parent.  A
        //  group only coarsens if all four siblings exist and all are
        //  flagged, which is the same rule hexRef2D applies.  Returns the
        //  number of groups coarsened.
        int unrefine(const std::vector<char>& flag);

        //- Refine until no two neighbours differ by more than one level.
        //  Returns the number of cells added.
        int balance();

        //- Rebuild owner/neighbour, the boundary list and all geometry.
        void buildAddressing();


    // Field transfer

        //- For each cell of the new mesh, one cell of the old mesh it
        //  overlaps: its ancestor if it was refined, its lowest-Morton
        //  child if it was coarsened, itself otherwise.
        //
        //  This is a nearest-value map and it is NOT conservative on
        //  coarsening: it keeps one of the four children and discards the
        //  other three, and always the same corner one, so the error is
        //  biased rather than cancelling.  OpenFOAM's refiner maps fields
        //  volume-weighted on unrefinement, so this does not reproduce it.
        //  Use mapRange/mapField for anything that has to be conserved,
        //  which for this solver is everything the flux scheme touches.
        //  This is kept only for quantities where nearest-value is wanted
        //  deliberately, such as an integer marker.
        static void parentMap
        (
            const std::vector<Key>& oldKeys,
            const std::vector<Key>& newKeys,
            std::vector<int>& map
        );

        //- For each cell of the new mesh, the range [first, first+count)
        //  of old cells that overlap it.
        //
        //  The range is always contiguous, and that is a property of the
        //  key layout rather than an assumption: keys sort by (base, z),
        //  and a cell's children tile its own z interval, so the leaves
        //  overlapping any z interval are adjacent in the array.  No index
        //  list is needed, and a device reads each cell's contributions as
        //  one coalesced run with nothing to serialise.
        static void mapRange
        (
            const std::vector<Key>& oldKeys,
            const std::vector<Key>& newKeys,
            std::vector<int>& first,
            std::vector<int>& count
        );

        //- Map a cell field conservatively:
        //
        //      phiNew[c] = sum_i |V_old,i ^ V_new,c| phiOld[i] / |V_new,c|
        //
        //  One expression covers all three cases.  Under refinement there
        //  is a single overlapping old cell and the overlap is the whole
        //  new cell, so the child inherits the parent's value exactly.
        //  Under coarsening the four children lie wholly inside, so the
        //  parent gets their volume-weighted mean.  Unchanged cells are
        //  the identity.  The integral of phi over the domain is therefore
        //  invariant, which is the property the solver needs and the one
        //  parentMap does not have.
        static void mapField
        (
            const std::vector<Key>& oldKeys,
            const std::vector<Key>& newKeys,
            const std::vector<int>& first,
            const std::vector<int>& count,
            const double* phiOld,
            double* phiNew
        );


    // Queries

        int nCells() const { return static_cast<int>(keys.size()); }
        int nFaces() const { return static_cast<int>(owner.size()); }
        int nBoundaryFaces() const { return static_cast<int>(bCell.size()); }

        //- Index of a key, or -1
        int find(Key k) const;

        //- The cell's extent, for tests and for geometry
        void cellBox(int i, double& x0, double& y0, double& x1, double& y1) const;
};


//- Whether the key set is a valid leaf set: sorted, unique, and covering
//  every base cell exactly once.  Returns an empty string if so.
std::string checkKeys(const Mesh& m);

//- Whether the addressing is a valid mesh: every face referenced by two
//  distinct cells, every cell closed to rounding, volumes summing to the
//  domain, and no two neighbours more than one level apart.
std::string checkMesh(const Mesh& m);

} // End namespace cdmAmr

#endif

// ************************************************************************* //
