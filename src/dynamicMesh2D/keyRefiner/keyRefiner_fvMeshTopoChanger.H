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

Class
    Foam::fvMeshTopoChangers::keyRefiner

Description
    2D adaptive refinement driven by a Morton key set rather than by
    polyTopoChange.

    Drop-in replacement for refiner2D: same dictionary entries, same
    criterion, same hysteresis, same buffer layers.  What changes is
    underneath.  refiner2D asks hexRef2D to build a topology change,
    polyTopoChange to renumber every entity through it, and fvMesh to map
    the fields; it does that twice per check, once to refine and once to
    unrefine.  Here the mesh is named rather than stored - a cell is
    (base cell, level, Morton position) - so refining is replacing a key
    with four keys, coarsening is replacing four with one, and both happen
    in one pass.  The whole mesh is then regenerated from the keys and
    handed to polyMesh::resetPrimitives, and the fields are carried across
    by a polyTopoChangeMap built by joining the two key sets.

    Measured on the RDE base mesh: 1.4 us per cell per check against
    refiner2D's 4.9, with every face's geometry agreeing with the engine
    to roundoff and the integral of a mapped field invariant to 1e-11.

Restrictions
    - Serial only.  The key engine's base mesh is one structured nx by ny
      block, which a decomposed case does not have.  Unlike keyRefiner,
      cdmGpuFoam supports multiple MPI ranks and GPUs with the refiner2D
      backend.  keyRefiner does not support distributed refinement or
      redistribution.

    - The run must start from the unrefined base mesh.  A key set cannot
      be recovered from an already-refined mesh: the refinement history is
      not in the mesh, and two different key sets can produce the same
      cells in a different order.  So the keys are written alongside the
      mesh at every write, and a restart reads them back.  Starting from a
      refined time directory with no keys file is refused rather than
      guessed at.

    - A side of the base mesh whose boundary faces span more than one
      patch is refused.  Handling it needs a per-base-face patch table,
      which is a small extension and not written yet.

Usage
    In constant/dynamicMeshDict, replacing refiner2D:
    \verbatim
    topoChanger
    {
        type            keyRefiner;
        libs            ("libCDMDynamicMesh2D.so");

        refineInterval  2;

        field           refinementIndicator;
        lowerRefineLevel 0.05;
        upperRefineLevel 1e6;
        unrefineLevel   0.015;

        nBufferLayers   2;
        maxRefinement   2;
        maxCells        4000000;

        // phi is recomputed by the solver every stage
        correctFluxes   ((phi none));

        dumpLevel       true;

        // Optional; deduced from the mesh if absent, and checked against it if given.
        // nx            750;
        // ny            180;
    }
    \endverbatim

SourceFiles
    keyRefiner_fvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef keyRefiner_fvMeshTopoChanger_H
#define keyRefiner_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "Switch.H"
#include "cdmAmrBridge.H"

#include <cstdint>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                 Class fvMeshTopoChangers::keyRefiner Declaration
\*---------------------------------------------------------------------------*/

class keyRefiner
:
    public fvMeshTopoChanger
{
    // Private Member Data

        //- Refinement control dictionary
        dictionary dict_;

        label refineInterval_;
        label maxCells_;
        label nBufferLayers_;

        //- Dump cellLevel as a volScalarField
        Switch dumpLevel_;

        //- Fluxes to map, as refiner2D reads them.  "none" means the
        //  field needs no correction because the solver recomputes it.
        HashTable<word> correctFluxes_;

        //- The mesh, as a sorted key array.  This is the authority: the
        //  fvMesh is regenerated from it, never the other way round.
        cdmAmr::Mesh amr_;

        //- Which patch each side of the base mesh belongs to, settled
        //  once from the mesh's own geometry
        labelList sidePatch_;
        label frontBack_;
        label cyclicMin_;
        label cyclicMax_;
        label nPatches_;

        //- The current mesh's canonical names, carried from one change to
        //  the next.  A topology map is a statement about how the new
        //  mesh relates to the old one, so the old one's names have to
        //  survive the change that replaces it.
        std::vector<cdmBridge::FaceKey> faceKey_;
        std::vector<std::uint64_t> pointCode_;
        std::vector<cdmBridge::FaceSort> faceSorted_;

        //- The key set the current mesh was built from.  Distinct from
        //  amr_.keys, which update() has already moved on by the time the
        //  map is built: the map needs both sides.
        std::vector<cdmAmr::Key> prevKeys_;

        //- The patch slices of the current mesh, which the map needs as
        //  the old mesh's
        labelList patchSizes_;
        labelList patchStarts_;

        label nRefinementIterations_;
        mutable bool changedSinceWrite_;
        label timeIndex_;


    // Wall-clock accumulators, in seconds, reported at the end of the run.
    // Deliberately the same breakdown refiner2D reports, so the two can be
    // put side by side without reinterpreting anything.

        double tCriterion_;
        double tBuffer_;
        double tKeys_;
        double tPrimitives_;
        double tMaps_;
        double tReset_;

        //- fvMesh::topoChange: field mapping, the geometry and addressing
        //  caches, and every meshObject and functionObject hook
        double tTopoChange_;

        double tFluxes_;
        double tGeometry_;

        double geomTouch_;
        label nUpdateCalls_;
        label nChanges_;
        label nRefined_;
        label nCoarsened_;


    // Private Member Functions

        //- Read the parameters from the dictionary
        void readDict();

        //- Write the stage timings
        void reportSplit() const;

        //- Work out nx, ny and the extent from the mesh, and check that
        //  the engine's base mesh is the same mesh cell for cell.  This
        //  is the one assumption the whole scheme rests on - that base
        //  cell j*nx + i is the fvMesh's cell j*nx + i - and it is cheap
        //  to verify, so it is verified rather than assumed.
        void fitBaseMesh();

        //- nx and ny from the key file's header, without reading the keys
        bool peekKeysHeader(label& nx, label& ny) const;

        //- Check the engine's mesh is the case's mesh cell for cell
        void verifyCentres() const;

        //- Seed the canonical names from the mesh as it stands
        void seedNames();

        //- The file the keys are written to and read from
        fileName keysPath(const word& instance) const;

        //- Write the key array beside the mesh
        void writeKeys() const;

        //- Read a key array written by a previous run.  Returns false if
        //  there is no file to read.
        bool readKeys();

        //- Cells the criterion wants refined, by the same route refiner2D
        //  takes: error band on a field, optionally per cellZone
        void selectRefineCandidates(PackedBoolList&, label& maxRefinement) const;

        //- Cells protected from coarsening: those marked for refinement,
        //  plus the unrefineLevel hysteresis band
        void selectKeepCells(const PackedBoolList& refineCells, PackedBoolList&) const;

        //- Grow a marked set by one cell-face-cell layer.  Done on the
        //  engine's own owner/neighbour arrays, which are already built,
        //  rather than on mesh().cells(), which is not.  The cyclic seam
        //  is internal to the engine, so the layer crosses it without a
        //  patch synchronisation.
        void extendMarked(std::vector<char>& marked) const;

        //- Rebuild the fvMesh from the keys and carry every field across
        void rebuild();

        //- Rescale mapped surface fluxes by the area ratio, for any flux
        //  field not declared "none"
        void correctFluxes(const labelList& faceMap, const scalarField& oldMagSf);

        //- Write cellLevel as a volScalarField, as refiner2D's dumpLevel
        //  does
        bool writeLevel() const;


public:

    //- Runtime type information
    TypeName("keyRefiner");


    // Constructors

        //- Construct from fvMesh and dictionary
        keyRefiner(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        keyRefiner(const keyRefiner&) = delete;


    //- Destructor
    virtual ~keyRefiner();


    // Member Functions

        //- Refinement level of every cell, from the keys
        labelList cellLevel() const;

        //- Update the mesh for topology change
        virtual bool update();

        //- Update corresponding to the given map.  Nothing to do: the
        //  keys are the authority and were updated before the map was
        //  built from them.
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Writing

        //- Write the key set beside the mesh
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const keyRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
