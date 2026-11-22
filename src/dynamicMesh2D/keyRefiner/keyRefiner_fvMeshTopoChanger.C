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

\*---------------------------------------------------------------------------*/

#include "keyRefiner_fvMeshTopoChanger.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"
#include "polyTopoChangeMap.H"
#include "objectMap.H"
#include "PackedBoolList.H"
#include "Pair.H"
#include "cellZone.H"
#include "Time.H"
#include "OFstream.H"
#include "IFstream.H"
#include "OSspecific.H"

#include <algorithm>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(keyRefiner, 0);
    addToRunTimeSelectionTable
    (
        fvMeshTopoChanger,
        keyRefiner,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::keyRefiner::readDict()
{
    refineInterval_ = dict_.lookupOrDefault<label>("refineInterval", 1);
    maxCells_ = dict_.lookupOrDefault<label>("maxCells", labelMax);
    nBufferLayers_ = dict_.lookupOrDefault<label>("nBufferLayers", 0);
    dumpLevel_ = dict_.lookupOrDefault("dumpLevel", Switch(false));

    if (refineInterval_ < 0)
    {
        FatalErrorInFunction
            << "refineInterval should be >= 0, is " << refineInterval_
            << exit(FatalError);
    }

    if (dict_.found("correctFluxes"))
    {
        const List<Pair<word>> fluxVelocities
        (
            dict_.lookup<List<Pair<word>>>("correctFluxes")
        );

        forAll(fluxVelocities, i)
        {
            correctFluxes_.insert
            (
                fluxVelocities[i][0],
                fluxVelocities[i][1]
            );
        }
    }

    if (dict_.isDict("refinementRegions"))
    {
        FatalErrorInFunction
            << "refinementRegions is not supported by keyRefiner yet."
            << "  Per-zone maximum levels need a per-cell level limit"
            << " threaded through the key refinement, which is a small"
            << " extension and not written." << nl
            << "  Use refiner2D for a case that needs it."
            << exit(FatalError);
    }
}


void Foam::fvMeshTopoChangers::keyRefiner::fitBaseMesh()
{
    const fvMesh& m = mesh();

    cdmBridge::deducePatchLayout
    (
        m, sidePatch_, frontBack_, cyclicMin_, cyclicMax_
    );

    // The polyBoundaryMesh, not the fvBoundaryMesh: resetPrimitives resets
    // every polyPatch, including the zero-size internal one the
    // finite-volume view omits.
    nPatches_ = m.boundaryMesh().size();

    const boundBox& bb = m.bounds();
    const scalar Lx = bb.span().x();
    const scalar Ly = bb.span().y();
    const scalar Lz = bb.span().z();
    const scalar tol = 1e-9*bb.mag();

    // The engine's cells are at absolute positions: base cell (i, j) spans
    // [i*dx, (i+1)*dx). An offset base mesh would need the engine to carry
    // an origin, so it is refused rather than built a domain width away.
    if (mag(bb.min().x()) > tol || mag(bb.min().y()) > tol)
    {
        FatalErrorInFunction
            << "the mesh starts at " << bb.min()
            << " but the key engine places base cell (i, j) at"
            << " [i*dx, (i+1)*dx) x [j*dy, (j+1)*dy), so the base mesh"
            << " must start at x = 0, y = 0." << nl
            << "  Translate the mesh, or extend cdmAmr::Mesh with an"
            << " origin." << exit(FatalError);
    }

    // nx and ny, from the dictionary if given, otherwise deduced from the
    // first cell's centre. Either way they are checked against every cell,
    // so a wrong value is caught here rather than permuting every field.
    label nx = dict_.lookupOrDefault<label>("nx", 0);
    label ny = dict_.lookupOrDefault<label>("ny", 0);

    // A refined mesh cannot be counted back to its base, so on a restart
    // nx and ny come from the key file's own header.
    const bool restarting = isFile(keysPath(mesh().time().name()));

    if ((nx <= 0 || ny <= 0) && restarting)
    {
        peekKeysHeader(nx, ny);
    }

    if (nx <= 0 || ny <= 0)
    {
        if (m.nCells() == 0)
        {
            FatalErrorInFunction << "empty mesh" << exit(FatalError);
        }

        const point& c0 = m.C()[0];

        if (c0.x() <= 0 || c0.y() <= 0)
        {
            FatalErrorInFunction
                << "cannot deduce the base mesh: cell 0 is centred at "
                << c0 << ", which is not in the first corner of the"
                << " domain.  Give nx and ny in the topoChanger dict."
                << exit(FatalError);
        }

        nx = label(Lx/(2*c0.x()) + 0.5);
        ny = label(Ly/(2*c0.y()) + 0.5);
    }

    if (nx <= 0 || ny <= 0 || (!restarting && nx*ny != m.nCells()))
    {
        FatalErrorInFunction
            << "the base mesh does not fit: nx = " << nx << ", ny = " << ny
            << ", nx*ny = " << nx*ny << ", but the mesh has "
            << m.nCells() << " cells." << nl
            << "  The run must start from the unrefined base mesh: a key"
            << " set cannot be recovered from an already-refined mesh,"
            << " because the refinement history is not in the mesh." << nl
            << "  Restart from time 0, or from a time directory whose"
            << " polyMesh carries the cdmAmrKeys file this changer"
            << " writes." << exit(FatalError);
    }

    amr_.thickness = Lz;
    amr_.initialise(nx, ny, Lx, Ly, cyclicMin_ >= 0);
    amr_.buildAddressing();

    // The assumption the whole scheme rests on: base cell j*nx + i is the
    // fvMesh's cell j*nx + i. If it were false every field would be permuted
    // on the first rebuild while the mesh stayed perfectly valid, so it is
    // checked here against every cell. On a restart the case mesh is refined,
    // and verifyCentres() makes the stronger comparison once the keys are read.
    if (restarting)
    {
        Info<< "keyRefiner: base mesh " << nx << " x " << ny
            << ", restarting from a key file" << endl;
        return;
    }

    scalar worst = 0;
    label worstCell = -1;
    const scalar h = min(Lx/nx, Ly/ny);

    forAll(m.C(), celli)
    {
        const scalar d =
            mag(m.C()[celli].x() - amr_.C[2*celli])
          + mag(m.C()[celli].y() - amr_.C[2*celli + 1]);

        if (d > worst)
        {
            worst = d;
            worstCell = celli;
        }
    }

    if (worst > 1e-8*h)
    {
        FatalErrorInFunction
            << "the key engine's base mesh is not the case's base mesh"
            << " cell for cell." << nl
            << "  Worst disagreement at cell " << worstCell << ": the mesh"
            << " has it at " << m.C()[worstCell] << ", the engine at ("
            << amr_.C[2*worstCell] << " " << amr_.C[2*worstCell + 1]
            << "), off by " << worst << " against a cell size of " << h
            << "." << nl
            << "  The engine numbers base cell (i, j) as j*nx + i over a"
            << " uniform grid.  A graded or multi-block base mesh needs a"
            << " per-base-cell geometry table, which is not written."
            << exit(FatalError);
    }

    Info<< "keyRefiner: base mesh " << nx << " x " << ny << " over "
        << Lx << " x " << Ly << " x " << Lz
        << ", cyclic in x: " << Switch(cyclicMin_ >= 0) << nl
        << "    cell centres agree with the engine to " << worst
        << " (cell size " << h << ")" << endl;
}


void Foam::fvMeshTopoChangers::keyRefiner::seedNames()
{
    // The names of the mesh as it stands. The primitives that come back are
    // discarded, but the canonical names are what the first topology map
    // joins against; without them it would call every face new.
    pointField points;
    faceList faces;
    labelList owner, neighbour, faceEngine;
    double tP = 0, tF = 0, tS = 0;

    patchSizes_.setSize(nPatches_);
    patchStarts_.setSize(nPatches_);

    cdmBridge::buildPrimitives
    (
        amr_, sidePatch_, frontBack_, cyclicMin_, cyclicMax_, nPatches_,
        points, faces, owner, neighbour, patchSizes_, patchStarts_,
        faceEngine, faceKey_, pointCode_, tP, tF, tS
    );

    if (faces.size() != mesh().nFaces())
    {
        FatalErrorInFunction
            << "the bridge builds " << faces.size() << " faces from the"
            << " keys but the case mesh has " << mesh().nFaces()
            << ".  The base mesh or the patch layout has been misread."
            << exit(FatalError);
    }

    faceSorted_.clear();
    for (std::size_t f = 0; f < faceKey_.size(); ++f)
    {
        if (faceKey_[f].orient >= 2) { continue; }
        faceSorted_.push_back
        (
            cdmBridge::FaceSort
            {
                cdmBridge::packFaceKey(faceKey_[f]),
                static_cast<label>(f)
            }
        );
    }
    psort(faceSorted_.begin(), faceSorted_.end());
}


// * * * * * * * * * * * * * * *  Key file I/O  * * * * * * * * * * * * * * //

Foam::fileName Foam::fvMeshTopoChangers::keyRefiner::keysPath
(
    const word& instance
) const
{
    return
        mesh().time().path()/instance/mesh().dbDir()
       /polyMesh::meshSubDir/"cdmAmrKeys";
}


void Foam::fvMeshTopoChangers::keyRefiner::writeKeys() const
{
    const fileName f(keysPath(mesh().time().name()));
    mkDir(f.path());

    OFstream os(f);

    // Written as three components rather than the packed 64-bit key, so every
    // number fits a label and the file survives a build with a different
    // label size. makeKey reassembles them exactly.
    os  << "// cdmAmrKeys 1" << nl
        << amr_.nx << token::SPACE << amr_.ny << token::SPACE
        << label(amr_.keys.size()) << nl;

    for (std::size_t i = 0; i < amr_.keys.size(); ++i)
    {
        const cdmAmr::Key k = amr_.keys[i];
        os  << label(cdmAmr::keyBase(k)) << token::SPACE
            << label(cdmAmr::keyLevel(k)) << token::SPACE
            << label(cdmAmr::keyMorton(k)) << nl;
    }
}


bool Foam::fvMeshTopoChangers::keyRefiner::readKeys()
{
    const fileName f(keysPath(mesh().time().name()));

    if (!isFile(f)) { return false; }

    IFstream is(f);

    label nx = 0, ny = 0, n = 0;
    is >> nx >> ny >> n;

    if (nx != amr_.nx || ny != amr_.ny)
    {
        FatalErrorInFunction
            << f << " is for a " << nx << " x " << ny << " base mesh but"
            << " this one is " << amr_.nx << " x " << amr_.ny
            << exit(FatalError);
    }

    std::vector<cdmAmr::Key> keys;
    keys.reserve(n);

    for (label i = 0; i < n; ++i)
    {
        label base = 0, level = 0, morton = 0;
        is >> base >> level >> morton;
        keys.push_back
        (
            cdmAmr::makeKey
            (
                static_cast<std::uint64_t>(base),
                static_cast<std::uint32_t>(level),
                static_cast<std::uint32_t>(morton)
            )
        );
    }

    amr_.keys.swap(keys);

    const std::string err = cdmAmr::checkKeys(amr_);
    if (!err.empty())
    {
        FatalErrorInFunction
            << "the key set in " << f << " is not a valid leaf set: "
            << err.c_str() << exit(FatalError);
    }

    amr_.buildAddressing();

    Info<< "keyRefiner: read " << n << " keys from " << f.name()
        << ", " << amr_.nCells() << " cells" << endl;

    return true;
}


// * * * * * * * * * * * * * * *  Cell selection * * * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::keyRefiner::selectRefineCandidates
(
    PackedBoolList& candidateCells,
    label& maxRefinement
) const
{
    const word fieldName(dict_.lookup("field"));
    const volScalarField& vFld =
        mesh().lookupObject<volScalarField>(fieldName);

    const scalar lowerRefineLevel =
        dict_.lookup<scalar>("lowerRefineLevel");
    const scalar upperRefineLevel =
        dict_.lookup<scalar>("upperRefineLevel");

    maxRefinement = dict_.lookup<label>("maxRefinement");

    if (maxRefinement <= 0)
    {
        FatalErrorInFunction
            << "Illegal maximum refinement level " << maxRefinement
            << exit(FatalError);
    }

    if (maxRefinement >= cdmAmr::maxLevel)
    {
        FatalErrorInFunction
            << "maxRefinement " << maxRefinement << " is at or above the"
            << " engine's level limit of " << cdmAmr::maxLevel
            << ".  A key holds the level in four bits and the Morton code"
            << " in sixteen, so " << cdmAmr::maxLevel - 1
            << " is the most this build supports."
            << exit(FatalError);
    }

    // The same error band refiner2D uses: in range is a candidate, out of
    // range is not.
    forAll(vFld, celli)
    {
        const scalar err =
            min(vFld[celli] - lowerRefineLevel, upperRefineLevel - vFld[celli]);

        if (err > 0)
        {
            candidateCells.set(celli, 1);
        }
    }

    if (dict_.found("cellZone"))
    {
        const word zoneName(dict_.lookup("cellZone"));
        const label zoneID = mesh().cellZones().findIndex(zoneName);

        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "cannot find cellZone " << zoneName << exit(FatalError);
        }

        PackedBoolList inZone(mesh().nCells());
        const labelList& zoneCells = mesh().cellZones()[zoneID];
        forAll(zoneCells, i) { inZone.set(zoneCells[i], 1); }

        forAll(vFld, celli)
        {
            if (!inZone.get(celli)) { candidateCells.unset(celli); }
        }
    }
}


void Foam::fvMeshTopoChangers::keyRefiner::selectKeepCells
(
    const PackedBoolList& refineCells,
    PackedBoolList& keepCells
) const
{
    forAll(keepCells, celli)
    {
        if (refineCells.get(celli)) { keepCells.set(celli, 1); }
    }

    if (!dict_.found("unrefineLevel")) { return; }

    // Hysteresis, as refiner2D does it: a cell still above unrefineLevel is
    // kept though no longer a refinement candidate. Without it a cell whose
    // indicator depends on its own size oscillates forever.
    const scalar unrefineLevel = dict_.lookup<scalar>("unrefineLevel");

    const volScalarField& vFld =
        mesh().lookupObject<volScalarField>(dict_.lookup<word>("field"));

    forAll(vFld, celli)
    {
        if (vFld[celli] >= unrefineLevel) { keepCells.set(celli, 1); }
    }
}


void Foam::fvMeshTopoChangers::keyRefiner::extendMarked
(
    std::vector<char>& marked
) const
{
    // One pass over the engine's internal faces. owner and neighbour are
    // already in hand and a cyclic direction is internal to the engine, so
    // the layer crosses the seam with nothing to synchronise.
    std::vector<char> grown(marked);
    const int nF = amr_.nFaces();

    for (int f = 0; f < nF; ++f)
    {
        const int o = amr_.owner[f];
        const int n = amr_.neighbour[f];

        if (marked[o] || marked[n])
        {
            grown[o] = 1;
            grown[n] = 1;
        }
    }

    marked.swap(grown);
}


// * * * * * * * * * * * * * * * *  The rebuild  * * * * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::keyRefiner::rebuild()
{
    fvMesh& m = mesh();
    double t0 = 0;

    // ---- the state the map is built against, before the mesh moves ----
    const label oldNPoints = m.nPoints();
    const label oldNFaces = m.nFaces();
    const label oldNCells = m.nCells();

    autoPtr<scalarField> oldVPtr(new scalarField(m.V()));

    // magSf of the old mesh, by face, for the flux rescaling. Empty patches
    // carry no finite-volume faces, so their entries stay zero and the
    // rescaling skips them.
    scalarField oldMagSf(oldNFaces, 0);

    // Only if something actually wants rescaling. Every entry being none is
    // the normal case, since the solver recomputes its flux.
    bool needFluxes = false;
    forAllConstIter(HashTable<word>, correctFluxes_, iter)
    {
        if (iter() != "none") { needFluxes = true; break; }
    }

    if (needFluxes)
    {
        const surfaceScalarField& ms = m.magSf();
        const scalarField& msi = ms.primitiveField();
        forAll(msi, f) { oldMagSf[f] = msi[f]; }

        forAll(m.boundary(), patchi)
        {
            const label st = m.boundary()[patchi].start();
            const fvsPatchScalarField& pf = ms.boundaryField()[patchi];
            forAll(pf, i) { oldMagSf[st + i] = pf[i]; }
        }
    }

    // ---- keys -> primitives ----
    t0 = cdmBridge::wall();

    pointField points;
    faceList faces;
    labelList owner, neighbour, faceEngine;
    labelList newPatchSizes(nPatches_), newPatchStarts(nPatches_);
    std::vector<cdmBridge::FaceKey> newFaceKey;
    std::vector<std::uint64_t> newPointCode;
    double tP = 0, tF = 0, tS = 0;

    cdmBridge::buildPrimitives
    (
        amr_, sidePatch_, frontBack_, cyclicMin_, cyclicMax_, nPatches_,
        points, faces, owner, neighbour, newPatchSizes, newPatchStarts,
        faceEngine, newFaceKey, newPointCode, tP, tF, tS
    );

    tPrimitives_ += cdmBridge::wall() - t0;

    // ---- the topology map, built while both meshes still have names ----
    t0 = cdmBridge::wall();

    labelList cellMap, pointMap, faceMap, revCell, revFace, revPoint;
    List<objectMap> cellsFromCells, facesFromFaces, pointsFromPoints;
    labelHashSet flipFaceFlux;
    std::vector<cdmBridge::FaceSort> newSorted;

    cdmBridge::buildTopoMaps
    (
        amr_, prevKeys_, pointCode_, newPointCode, faceKey_, newFaceKey,
        faceSorted_, newSorted,
        cellMap, cellsFromCells, pointMap,
        faceMap, facesFromFaces,
        revCell, revFace, revPoint, flipFaceFlux
    );

    tMaps_ += cdmBridge::wall() - t0;

    // ---- hand the primitives over ----
    t0 = cdmBridge::wall();

    m.resetPrimitives
    (
        std::move(points),
        std::move(faces),
        std::move(owner),
        std::move(neighbour),
        newPatchSizes,
        newPatchStarts,
        true
    );

    tReset_ += cdmBridge::wall() - t0;

    // fvMesh::topoChange can be called on a mesh resetPrimitives has already
    // replaced, which is not obvious and is the whole reason this works. Its
    // first act is polyMesh::topoChange, pure bookkeeping that reorders
    // nothing, so a map whose new side is already installed is what it wants.
    // Everything after that is needed and none of it reachable from outside:
    // clearFvGeomNotOldVol drops fvMesh's geometry while keeping V0, V00 and
    // the mesh flux a plain clearOut would destroy, then mapFields,
    // surfaceInterpolation::clearOut, clearAddressing, the fvMesh and lduMesh
    // meshObject hooks, the function objects, and the stitcher, distributor
    // and mover.
    // Doing resetPrimitives, clearOut, mapFields and the two hooks by hand
    // builds a valid mesh and maps every field correctly, then segfaults a
    // step later: a mesh object registered on the polyMesh survives clearOut
    // as updateable and is never told the topology changed.
    t0 = cdmBridge::wall();
    {
        labelListList patchPointMap(nPatches_);
        labelList oldPatchNMeshPoints(nPatches_, 0);

        polyTopoChangeMap tcMap
        (
            m,
            oldNPoints, oldNFaces, oldNCells,
            std::move(pointMap), std::move(pointsFromPoints),
            std::move(faceMap), std::move(facesFromFaces),
            std::move(cellMap), std::move(cellsFromCells),
            std::move(revPoint), std::move(revFace), std::move(revCell),
            std::move(flipFaceFlux),
            std::move(patchPointMap),
            std::move(patchSizes_),
            std::move(patchStarts_),
            std::move(oldPatchNMeshPoints),
            std::move(oldVPtr)
        );

        m.topoChange(tcMap);
        tTopoChange_ += cdmBridge::wall() - t0;

        if (needFluxes)
        {
            t0 = cdmBridge::wall();
            correctFluxes(tcMap.faceMap(), oldMagSf);
            tFluxes_ += cdmBridge::wall() - t0;
        }
    }

    // Force the geometry back, as the solver's first access would. This
    // attributes the cost rather than changing it; left alone, fvMesh
    // rebuilds V and Sf lazily and the bill lands on the first reader.
    t0 = cdmBridge::wall();
    geomTouch_ += mag(m.V()[0]) + mag(m.magSf().primitiveField()[0]);
    tGeometry_ += cdmBridge::wall() - t0;

    // ---- the state the next map will be built against ----
    prevKeys_ = amr_.keys;
    faceKey_.swap(newFaceKey);
    pointCode_.swap(newPointCode);
    faceSorted_.swap(newSorted);
    patchSizes_ = newPatchSizes;
    patchStarts_ = newPatchStarts;

    ++nChanges_;
}


void Foam::fvMeshTopoChangers::keyRefiner::correctFluxes
(
    const labelList& faceMap,
    const scalarField& oldMagSf
)
{
    // A flux is extensive in area: mapping gives a new face its ancestor's
    // value, wrong by the area ratio for one that split or merged. Scaling by
    // that ratio is exact when the velocity is uniform across the ancestor.
    // A field declared none is left alone, which is what the solver's own
    // cases ask for: phi is recomputed from the Riemann solver every stage,
    // so correcting it here would be overwritten before anything read it.
    const fvMesh& m = mesh();
    const surfaceScalarField& ms = m.magSf();

    forAllConstIter(HashTable<word>, correctFluxes_, iter)
    {
        const word& fluxName = iter.key();

        if (iter() == "none") { continue; }

        if (!m.foundObject<surfaceScalarField>(fluxName)) { continue; }

        surfaceScalarField& phi =
            m.lookupObjectRef<surfaceScalarField>(fluxName);

        scalarField& phii = phi.primitiveFieldRef();
        const scalarField& msi = ms.primitiveField();

        forAll(phii, facei)
        {
            const label oldFacei = faceMap[facei];

            if (oldFacei >= 0 && oldMagSf[oldFacei] > vSmall)
            {
                phii[facei] *= msi[facei]/oldMagSf[oldFacei];
            }
        }

        forAll(phi.boundaryField(), patchi)
        {
            fvsPatchScalarField& pf = phi.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& pms = ms.boundaryField()[patchi];
            const label st = m.boundary()[patchi].start();

            forAll(pf, i)
            {
                const label oldFacei = faceMap[st + i];

                if (oldFacei >= 0 && oldMagSf[oldFacei] > vSmall)
                {
                    pf[i] *= pms[i]/oldMagSf[oldFacei];
                }
            }
        }
    }
}


Foam::labelList Foam::fvMeshTopoChangers::keyRefiner::cellLevel() const
{
    labelList level(amr_.nCells());

    forAll(level, celli)
    {
        level[celli] = label(cdmAmr::keyLevel(amr_.keys[celli]));
    }

    return level;
}


bool Foam::fvMeshTopoChangers::keyRefiner::writeLevel() const
{
    // Built here and discarded rather than registered on the mesh: a
    // registered field would be mapped on every rebuild and then immediately
    // overwritten from the keys.
    volScalarField scalarCellLevel
    (
        IOobject
        (
            "cellLevel",
            mesh().time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh(),
        dimensionedScalar(dimless, 0)
    );

    forAll(scalarCellLevel, celli)
    {
        scalarCellLevel[celli] = cdmAmr::keyLevel(amr_.keys[celli]);
    }

    return scalarCellLevel.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::keyRefiner::keyRefiner
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    refineInterval_(1),
    maxCells_(labelMax),
    nBufferLayers_(0),
    dumpLevel_(false),
    sidePatch_(4, -1),
    frontBack_(-1),
    cyclicMin_(-1),
    cyclicMax_(-1),
    nPatches_(0),
    nRefinementIterations_(0),
    changedSinceWrite_(false),
    timeIndex_(-1),
    tCriterion_(0),
    tBuffer_(0),
    tKeys_(0),
    tPrimitives_(0),
    tMaps_(0),
    tReset_(0),
    tTopoChange_(0),
    tFluxes_(0),
    tGeometry_(0),
    geomTouch_(0),
    nUpdateCalls_(0),
    nChanges_(0),
    nRefined_(0),
    nCoarsened_(0)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "keyRefiner runs in serial only." << nl
            << "  The key engine's base mesh is one structured nx by ny"
            << " block, which a decomposed case does not have, and a"
            << " redistribution would move cells between key sets that"
            << " know nothing of each other." << nl
            << "  For a parallel run use refiner2D." << exit(FatalError);
    }

    readDict();
    fitBaseMesh();

    // A restart from a refined mesh is only possible if the keys were written
    // with it; fitBaseMesh has already refused a refined mesh without them.
    readKeys();

    // After the keys, so that a restart compares the engine's refined
    // mesh with the case's refined mesh.
    verifyCentres();

    seedNames();
    prevKeys_ = amr_.keys;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::keyRefiner::~keyRefiner()
{
    if (nUpdateCalls_ > 0)
    {
        reportSplit();
    }

    // The engine's arrays cannot be returned to the heap at this point in a
    // GPU run: see cdmAmr::Mesh::releaseWithoutFreeing for what goes wrong.
    amr_.releaseWithoutFreeing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::keyRefiner::update()
{
    // Only once per time step
    if (timeIndex_ != mesh().time().timeIndex())
    {
        timeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        return false;
    }

    if (refineInterval_ == 0) { return false; }

    if
    (
        mesh().time().timeIndex() <= 0
     || mesh().time().timeIndex() % refineInterval_ != 0
    )
    {
        return false;
    }

    ++nUpdateCalls_;

    // ---- what the criterion wants ----
    double t0 = cdmBridge::wall();

    PackedBoolList refineCells(mesh().nCells());
    label maxRefinement = 0;
    selectRefineCandidates(refineCells, maxRefinement);

    tCriterion_ += cdmBridge::wall() - t0;

    // ---- buffer layers, and the cells kept from coarsening ----
    t0 = cdmBridge::wall();

    const label nCells = amr_.nCells();

    if (nCells != mesh().nCells())
    {
        FatalErrorInFunction
            << "the engine has " << nCells << " cells and the fvMesh "
            << mesh().nCells() << ".  Something changed the mesh behind"
            << " this changer's back; the criterion's per-cell values"
            << " would be read against the wrong cells."
            << exit(FatalError);
    }

    std::vector<char> toRefine(nCells, 0), keep(nCells, 0);

    {
        PackedBoolList keepCells(nCells);
        selectKeepCells(refineCells, keepCells);

        for (label c = 0; c < nCells; ++c)
        {
            toRefine[c] = refineCells.get(c) ? 1 : 0;
            keep[c] = keepCells.get(c) ? 1 : 0;
        }
    }

    for (label i = 0; i < nBufferLayers_; ++i)
    {
        extendMarked(toRefine);
        extendMarked(keep);
    }

    // A cell at the level limit is not a candidate, and nothing refines once
    // the cell count is spent. Coarsening is always allowed, which is what
    // keeps the budget available where it is needed.
    const bool roomToRefine = (amr_.nCells() < maxCells_);

    for (label c = 0; c < nCells; ++c)
    {
        const int level = int(cdmAmr::keyLevel(amr_.keys[c]));

        if (!roomToRefine || level >= maxRefinement) { toRefine[c] = 0; }
        if (toRefine[c]) { keep[c] = 1; }
    }

    tBuffer_ += cdmBridge::wall() - t0;

    // ---- apply it to the keys, in one pass each way ----
    t0 = cdmBridge::wall();

    std::vector<char> toCoarsen(nCells, 0);
    for (label c = 0; c < nCells; ++c)
    {
        toCoarsen[c] =
            (!keep[c] && cdmAmr::keyLevel(amr_.keys[c]) > 0) ? 1 : 0;
    }

    std::vector<int> newIndex;
    const int nRef = amr_.refine(toRefine, &newIndex);

    // Carry the coarsening flag onto the refined array: refine() shifts every
    // index after the first refined cell. A just-refined cell is not a
    // coarsening candidate and its four children inherit no flag.
    std::vector<char> crs(amr_.nCells(), 0);
    for (std::size_t c = 0; c < newIndex.size(); ++c)
    {
        const int w = newIndex[c];

        if (w >= 0 && w < amr_.nCells() && amr_.keys[w] == prevKeys_[c])
        {
            crs[w] = toCoarsen[c];
        }
    }

    const int nCrs = amr_.unrefine(crs);
    const int nBal = amr_.balance();

    if (nRef == 0 && nCrs == 0)
    {
        tKeys_ += cdmBridge::wall() - t0;
        return false;
    }

    amr_.buildAddressing();
    tKeys_ += cdmBridge::wall() - t0;

    nRefined_ += nRef + nBal;
    nCoarsened_ += nCrs;

    // ---- regenerate the mesh and carry the fields ----
    rebuild();

    ++nRefinementIterations_;
    changedSinceWrite_ = true;

    Info<< "keyRefiner: refined " << nRef << " cells";
    if (nBal) { Info<< " (+" << nBal << " to keep 2:1)"; }
    Info<< ", coarsened " << nCrs << " groups, now "
        << amr_.nCells() << " cells" << endl;

    return true;
}


void Foam::fvMeshTopoChangers::keyRefiner::topoChange
(
    const polyTopoChangeMap&
)
{
    // Nothing to do. The keys are the authority and were updated before the
    // map was built from them, so there is no numbering of ours to correct.
}


void Foam::fvMeshTopoChangers::keyRefiner::mapMesh(const polyMeshMap&)
{
    FatalErrorInFunction
        << "keyRefiner cannot follow a change to an entirely new mesh:"
        << " the key set describes the mesh it built, and a mesh that"
        << " arrived from elsewhere has no key set."
        << exit(FatalError);
}


void Foam::fvMeshTopoChangers::keyRefiner::distribute
(
    const polyDistributionMap&
)
{
    FatalErrorInFunction
        << "keyRefiner cannot be redistributed; it is serial only."
        << exit(FatalError);
}


bool Foam::fvMeshTopoChangers::keyRefiner::write(const bool) const
{
    // The keys go beside the mesh every time it is written: they are the only
    // record of the refinement, so a restart from a written time has no other
    // way back to a key set.
    if (!changedSinceWrite_) { return true; }

    writeKeys();

    bool ok = true;
    if (dumpLevel_) { ok = writeLevel(); }

    changedSinceWrite_ = false;
    return ok;
}


void Foam::fvMeshTopoChangers::keyRefiner::reportSplit() const
{
    const double total =
        tCriterion_ + tBuffer_ + tKeys_ + tPrimitives_ + tMaps_
      + tReset_ + tTopoChange_ + tFluxes_ + tGeometry_;

    Info<< nl
        << "keyRefiner stage split over " << nUpdateCalls_
        << " refinement checks (" << nChanges_ << " changed the mesh)"
        << nl
        << "    (tracked wall clock " << total << " s)" << nl;

    #define REPORT(t, what)                                               \
        Info<< "    " << what << "  " << t << " s  "                      \
            << (total > 0 ? 100.0*t/total : 0.0) << " %" << nl;

    REPORT(tCriterion_,   "criterion        ")
    REPORT(tBuffer_,      "buffer layers    ")
    REPORT(tKeys_,        "keys + addressing")
    REPORT(tPrimitives_,  "keys -> faces    ")
    REPORT(tMaps_,        "topology maps    ")
    REPORT(tReset_,       "resetPrimitives  ")
    REPORT(tTopoChange_,  "fvMesh topoChange")
    REPORT(tFluxes_,      "flux rescaling   ")
    REPORT(tGeometry_,    "fvMesh geometry  ")
    #undef REPORT

    Info<< "    ----" << nl
        << "    " << nRefined_ << " cells refined, " << nCoarsened_
        << " groups coarsened" << nl;

    if (nChanges_ > 0)
    {
        Info<< "    " << 1000.0*total/nChanges_ << " ms per change, "
            << 1.0e6*total/(scalar(nChanges_)*max(amr_.nCells(), 1))
            << " us per cell" << nl;
    }

    Info<< endl;
}


// ************************************************************************* //


bool Foam::fvMeshTopoChangers::keyRefiner::peekKeysHeader
(
    label& nx,
    label& ny
) const
{
    const fileName f(keysPath(mesh().time().name()));

    if (!isFile(f)) { return false; }

    IFstream is(f);
    label n = 0;
    is >> nx >> ny >> n;

    return nx > 0 && ny > 0;
}


void Foam::fvMeshTopoChangers::keyRefiner::verifyCentres() const
{
    // The assumption the whole scheme rests on: the engine's cell i is the
    // fvMesh's cell i. If it were false every field would be permuted on the
    // first rebuild while the mesh stayed valid, and nothing could detect it.
    const fvMesh& m = mesh();

    if (amr_.nCells() != m.nCells())
    {
        FatalErrorInFunction
            << "the engine has " << amr_.nCells() << " cells but the mesh"
            << " has " << m.nCells() << "." << nl
            << "  On a restart that means the key file and the mesh in"
            << " the same time directory are not the same mesh."
            << exit(FatalError);
    }

    const scalar h = min(amr_.Lx/amr_.nx, amr_.Ly/amr_.ny);

    scalar worst = 0;
    label worstCell = -1;

    forAll(m.C(), celli)
    {
        const scalar d =
            mag(m.C()[celli].x() - amr_.C[2*celli])
          + mag(m.C()[celli].y() - amr_.C[2*celli + 1]);

        if (d > worst)
        {
            worst = d;
            worstCell = celli;
        }
    }

    if (worst > 1e-8*h)
    {
        FatalErrorInFunction
            << "the key engine's mesh is not the case's mesh cell for"
            << " cell." << nl
            << "  Worst at cell " << worstCell << ": the mesh has it at "
            << m.C()[worstCell] << ", the engine at ("
            << amr_.C[2*worstCell] << " " << amr_.C[2*worstCell + 1]
            << "), off by " << worst << " against a base cell size of "
            << h << "." << exit(FatalError);
    }

    Info<< "    cell centres agree with the engine to " << worst
        << " over " << m.nCells() << " cells (base cell " << h << ")"
        << endl;
}


// ************************************************************************* //
