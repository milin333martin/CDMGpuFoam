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

Application
    cdmAmrMesh

Description
    Builds an OpenFOAM polyMesh from a cdmAmr key set, and measures what
    that costs.

    The engine in src/amr names cells with Morton keys and derives its own
    addressing from them.  On the device that is 67 ms per refinement check
    at a million cells, against 1.44 s for hexRef2D plus polyTopoChange
    plus fvMesh::topoChange.  But the solver's boundary conditions, its
    output and everything that is not the GPU flux path need a real
    fvMesh, so the keys have to be turned back into points, faces, owner
    and neighbour and handed to polyMesh::resetPrimitives.

    This application measures exactly that, and nothing else:

        keys -> primitives      how long to build points/faces/owner/nei
        resetPrimitives         how long OpenFOAM takes to accept them

    The second number is the one that decides the design.  Geometry was
    measured at 20 ms per topology change and the faceList at 32 ms, but
    initMesh's reconstruction of cells from owner/neighbour has never been
    measured, and if it is large then pushing the device side further is
    pointless.

    Correctness is not checked here by assertions of mine.  The mesh is
    written out, and the real checkMesh utility is run on it.  That checks
    closure, face ordering, pyramid volumes, non-orthogonality and skewness
    - far more than I would have thought to assert, and it is the code the
    rest of OpenFOAM trusts.

    A note on the cyclic seam.  The engine treats a cyclic direction as
    internal: a face across the seam is an ordinary internal face with an
    owner and a neighbour, not a pair of coupled patch faces.  That is
    deliberate.  Matching cyclic patch faces one-to-one would force equal
    refinement on both sides of the seam, which defeats the point of
    adaptive refinement there.  Keeping the seam internal costs nothing,
    removes an interpolation, and leaves every other patch - and so every
    other boundary condition - exactly as it was.  The host mesh built here
    therefore has no cyclic patch at all when cyclicX is set.

Usage
    cdmAmrMesh [-case dir]

    Reads system/amrMeshDict.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OSspecific.H"
#include "zeroGradientFvPatchFields.H"
#include "wallPolyPatch.H"

#include "cdmAmrBridge.H"

#include <algorithm>
#include <chrono>
#include <utility>
#include <vector>

using namespace Foam;
using namespace cdmBridge;



// * * * * * * * * * * * * * * *  Refinement * * * * * * * * * * * * * * * * //

// A band of refinement centred on xc, as a detonation front would be. The
// point is to exercise refinement and coarsening together at a realistic
// rate, not to be physical.
static void flagBand
(
    const cdmAmr::Mesh& m,
    double xc,
    double halfWidth,
    int maxRefine,
    std::vector<char>& toRefine,
    std::vector<char>& toCoarsen
)
{
    const int n = m.nCells();
    toRefine.assign(n, 0);
    toCoarsen.assign(n, 0);

    for (int c = 0; c < n; ++c)
    {
        const double x = m.C[2*c];
        const bool in = (std::fabs(x - xc) < halfWidth);
        const int level = static_cast<int>(cdmAmr::keyLevel(m.keys[c]));

        if (in && level < maxRefine) { toRefine[c] = 1; }
        if (!in && level > 0)        { toCoarsen[c] = 1; }
    }
}


// * * * * * * * * * * * * * * * * * *  main * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary dict
    (
        IOobject
        (
            "amrMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ
        )
    );

    const label nx = dict.lookup<label>("nx");
    const label ny = dict.lookup<label>("ny");
    const scalar Lx = dict.lookup<scalar>("Lx");
    const scalar Ly = dict.lookup<scalar>("Ly");
    const scalar thickness = dict.lookup<scalar>("thickness");
    const Switch cyclicX(dict.lookup("cyclicX"));
    const label nSteps = dict.lookup<label>("nSteps");
    const label maxRefine = dict.lookup<label>("maxRefinement");
    const Switch writeMesh(dict.lookup("writeMesh"));

    // How many registered fields to carry. mapFields costs roughly its
    // per-field cost times this, and the solver carries about fifteen, so
    // measuring with one understates it badly.
    const label nTestFields = dict.lookupOrDefault<label>("nTestFields", 1);

    Info<< "Base mesh " << nx << " x " << ny
        << " over " << Lx << " x " << Ly
        << ", cyclic in x: " << cyclicX
        << ", refining to level " << maxRefine << nl << endl;

    // readCase reads the case's own mesh, with its own patches, and works the
    // base mesh and patch layout out from it: the production path, where the
    // solver's mesh exists already. Without it the application builds its own
    // uniform mesh with invented patch names, which exercises the bridge but
    // not the fitting.
    const Switch readCase(dict.lookupOrDefault("readCase", Switch(false)));

    autoPtr<fvMesh> meshPtr;

    labelList sidePatch(4, -1);
    label frontBack = -1, cyclicMin = -1, cyclicMax = -1, nPatches = 0;

    // the effective geometry; taken from the case when reading one
    scalar LxUse = Lx, LyUse = Ly, thickUse = thickness;
    bool cyclicXUse = cyclicX;

    // only used when building our own mesh
    DynamicList<word> patchName;
    DynamicList<label> patchKind;

    if (readCase)
    {
        meshPtr.set
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        const fvMesh& cm = meshPtr();

        deducePatchLayout(cm, sidePatch, frontBack, cyclicMin, cyclicMax);
        // the polyBoundaryMesh: resetPrimitives resets every polyPatch,
        // including the zero-size internal one the fvBoundaryMesh omits
        nPatches = cm.boundaryMesh().size();

        const boundBox& bb = cm.bounds();
        LxUse = bb.span().x();
        LyUse = bb.span().y();
        thickUse = bb.span().z();
        cyclicXUse = (cyclicMin >= 0);

        if (cm.nCells() != nx*ny)
        {
            FatalErrorInFunction
                << "the case has " << cm.nCells() << " cells but nx*ny is "
                << nx*ny << ".  The engine must start from the unrefined"
                << " base mesh: a key set cannot be recovered from an"
                << " already-refined mesh, because the refinement history"
                << " is not in the mesh." << abort(FatalError);
        }

        Info<< "Read the case mesh: " << cm.nCells() << " cells, "
            << cm.nFaces() << " faces, " << nPatches << " patches" << nl
            << "Base " << nx << " x " << ny << " over " << LxUse
            << " x " << LyUse << " x " << thickUse
            << ", cyclic in x: " << Switch(cyclicXUse) << nl << endl;
    }
    else
    {
        // A cyclic direction gets no patch of its own: the engine makes those
        // faces internal and the bridge turns each into a matched cyclic pair.
        // Concrete constructors rather than polyPatch::New; the cyclic pair goes
        // through the dictionary constructor, the path a boundary file takes.
        // kind: 0 generic, 1 wall, 2 empty, 3 cyclic
        if (cyclicX)
        {
            cyclicMin = patchName.size();
            patchName.append("xMin"); patchKind.append(3);
            cyclicMax = patchName.size();
            patchName.append("xMax"); patchKind.append(3);
        }
        else
        {
            sidePatch[0] = patchName.size();
            patchName.append("xMin"); patchKind.append(0);
            sidePatch[1] = patchName.size();
            patchName.append("xMax"); patchKind.append(0);
        }

        sidePatch[2] = patchName.size();
        patchName.append("yMin"); patchKind.append(1);
        sidePatch[3] = patchName.size();
        patchName.append("yMax"); patchKind.append(1);

        frontBack = patchName.size();
        patchName.append("frontAndBack");
        patchKind.append(2);

        nPatches = patchName.size();
    }

    // Reading a case, it starts unrefined to match the mesh on disk and the
    // loop below does the refining. Building our own, it is refined first so
    // the run does not begin from a trivial mesh.
    cdmAmr::Mesh m;
    m.initialise(nx, ny, LxUse, LyUse, cyclicXUse);
    m.thickness = thickUse;
    m.buildAddressing();

    std::vector<char> toRefine, toCoarsen;

    if (!readCase)
    {
        for (int pass = 0; pass < maxRefine; ++pass)
        {
            flagBand(m, 0.5*LxUse, 0.08*LxUse, maxRefine, toRefine, toCoarsen);
            m.refine(toRefine);
            m.balance();
            m.buildAddressing();
        }
    }

    Info<< "Starting mesh " << m.nCells() << " cells, "
        << m.nFaces() << " internal faces" << nl << endl;


    // ---- build the primitives and construct the mesh once ----
    pointField points;
    faceList faces;
    labelList owner, neighbour, patchSizes, patchStarts, faceEngine;
    std::vector<FaceKey> faceKeyNew;
    std::vector<std::uint64_t> pointCodeNew;
    double tP = 0, tF = 0, tS = 0;

    buildPrimitives
    (
        m, sidePatch, frontBack, cyclicMin, cyclicMax, nPatches,
        points, faces, owner, neighbour, patchSizes, patchStarts,
        faceEngine, faceKeyNew, pointCodeNew, tP, tF, tS
    );

    Info<< "First build: " << points.size() << " points, "
        << faces.size() << " faces" << nl;

    double t0 = wall();

    // An fvMesh, not a polyMesh. polyMesh::resetPrimitives is the hook, but
    // fvMesh caches surface interpolation weights, deltaCoeffs and V0 on top
    // of polyMesh's addressing, and clearAddressing(true) does not reach them.
    if (!readCase)
    {
    meshPtr.set
    (
        new fvMesh
        (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points),
        std::move(faces),
        std::move(owner),
        std::move(neighbour)
        )
    );

    fvMesh& newMesh = meshPtr();
    const double tConstruct = wall() - t0;

    // The patches must exist before resetPrimitives, which resets their
    // size and start but never adds or removes one.
    {
        List<polyPatch*> p(nPatches);

        forAll(p, i)
        {
            switch (patchKind[i])
            {
                case 3:
                {
                    // Built from a dictionary holding exactly what a boundary file holds for
                    // a translational cyclic, so this goes down the same code path as reading
                    // one. The patch at x = 0 separates by -Lx.
                    dictionary pd;
                    pd.add("type", cyclicPolyPatch::typeName);
                    pd.add("nFaces", patchSizes[i]);
                    pd.add("startFace", patchStarts[i]);
                    pd.add("matchTolerance", scalar(1e-4));
                    pd.add
                    (
                        "neighbourPatch",
                        (i == cyclicMin) ? patchName[cyclicMax]
                                         : patchName[cyclicMin]
                    );
                    pd.add("transformType", word("translational"));
                    pd.add
                    (
                        "separation",
                        vector((i == cyclicMin) ? -LxUse : LxUse, 0, 0)
                    );

                    p[i] = new cyclicPolyPatch
                    (
                        patchName[i], pd, i, newMesh.boundaryMesh(),
                        cyclicPolyPatch::typeName
                    );
                    break;
                }

                case 1:
                    p[i] = new wallPolyPatch
                    (
                        patchName[i], patchSizes[i], patchStarts[i], i,
                        newMesh.boundaryMesh(), wallPolyPatch::typeName
                    );
                    break;

                case 2:
                    p[i] = new emptyPolyPatch
                    (
                        patchName[i], patchSizes[i], patchStarts[i], i,
                        newMesh.boundaryMesh(), emptyPolyPatch::typeName
                    );
                    break;

                default:
                    p[i] = new polyPatch
                    (
                        patchName[i], patchSizes[i], patchStarts[i], i,
                        newMesh.boundaryMesh(), polyPatch::typeName
                    );
                    break;
            }
        }

        // addFvPatches, not addPatches: on an fvMesh the fvBoundaryMesh has to be
        // built alongside the polyBoundaryMesh, or weights() reads a boundary
        // that does not exist.
        t0 = wall();
        newMesh.addFvPatches(p);
        Info<< "addFvPatches " << 1e3*(wall() - t0) << " ms" << nl << endl;
    }

    Info<< "construct from primitives " << 1e3*tConstruct << " ms" << nl
        << endl;
    }   // end of !readCase

    fvMesh& mesh = meshPtr();

    // When reading a case the primitives built above describe the base mesh
    // the case already has, so they are not installed; they only seed the
    // state the first topology map is built against.
    if (readCase && mesh.nFaces() != faces.size())
    {
        WarningInFunction
            << "the case mesh has " << mesh.nFaces() << " faces but the"
            << " bridge builds " << faces.size() << " for the same cells."
            << "  The face ordering will differ from the case's, which is"
            << " harmless here - the first resetPrimitives replaces it -"
            << " but it means the first topology map relates the bridge's"
            << " own ordering, not the case's." << endl;
    }

    // A registered field, so fvMesh::mapFields has something real to map and
    // the result can be checked: the integral of a cell field over the domain
    // must not move when the mesh changes under it.
    volScalarField phiTest
    (
        IOobject
        (
            "phiTest",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    {
        scalarField& pf = phiTest.primitiveFieldRef();
        forAll(pf, c)
        {
            const scalar x = mesh.C()[c].x();
            pf[c] = 1.0 + 5.0*Foam::exp(-sqr((x - 0.5*LxUse)/(0.03*LxUse)));
        }
        phiTest.correctBoundaryConditions();
    }

    // The rest of the registered fields. Never read; they exist so mapFields
    // does the work it would do in the solver, which carries roughly fifteen.
    // Each is registered on the mesh, so mapFields finds it without being told.
    PtrList<volScalarField> ballast
    (
        (nTestFields > 1) ? nTestFields - 1 : 0
    );

    forAll(ballast, i)
    {
        ballast.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "ballast" + Foam::name(i),
                    runTime.name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phiTest
            )
        );
    }

    Info<< "Carrying " << nTestFields << " registered volScalarFields"
        << nl << endl;

    // ---- the measurement: refine, rebuild, resetPrimitives ----
    Info<< "Timing " << nSteps << " refinement checks" << nl << endl;

    double sKeys = 0, sPoints = 0, sFaces = 0, sSort = 0, sReset = 0;
    double sFv = 0, sClear = 0, sMaps = 0, sMapFields = 0;

    // The previous state, which each topology map is built against. The map is
    // a statement about how the new mesh relates to the old one, so the old
    // one's names have to survive the change.
    std::vector<cdmAmr::Key> prevKeys = m.keys;
    std::vector<FaceKey> prevFaceKey = faceKeyNew;
    std::vector<std::uint64_t> prevCode = pointCodeNew;
    labelList prevPatchSizes(patchSizes), prevPatchStarts(patchStarts);
    label prevNPoints = mesh.nPoints();
    label prevNFaces = mesh.nFaces();
    label prevNCells = mesh.nCells();
    scalarField prevV(mesh.V());

    scalar worstMapDrift = 0;
    label nFromFaces = 0, nFromCells = 0, nFlip = 0, nInsertedFace = 0;

    // The old mesh's in-plane faces, sorted. Carried across steps rather than
    // rebuilt: the old list is the previous step's new list, already in order,
    // so one of the two sorts per step disappears.
    std::vector<FaceSort> prevSorted, newSorted;

    // Seeded from the mesh as first built, or the first map would find no
    // old faces at all and call every face new.
    for (std::size_t f = 0; f < prevFaceKey.size(); ++f)
    {
        if (prevFaceKey[f].orient >= 2) { continue; }
        prevSorted.push_back
        (
            FaceSort{packFaceKey(prevFaceKey[f]), static_cast<label>(f)}
        );
    }
    std::sort(prevSorted.begin(), prevSorted.end());

    // Worst relative disagreement between what fvMesh reports after
    // resetPrimitives and what the engine computed from the keys. Agreement to
    // roundoff means fvMesh really did rebuild; a stale cache from the
    // previous mesh would show here and nowhere else.
    scalar worstV = 0, worstMagSf = 0, worstSf = 0, worstW = 0, worstDC = 0;
    label nFvChecked = 0;

    for (label step = 0; step < nSteps; ++step)
    {
        const label nDen = (nSteps > 1) ? (nSteps - 1) : 1;
        const double xc = LxUse*(0.2 + 0.6*scalar(step)/scalar(nDen));

        // the engine's own work
        t0 = wall();
        flagBand(m, xc, 0.08*LxUse, maxRefine, toRefine, toCoarsen);

        // Refine first, then carry the coarsening flag across. refine() shifts
        // every index after the first refined cell, so a flag built before the
        // refinement coarsens the wrong cells: in bounds, so nothing warns, and
        // the only symptom is a mesh that refines about half as much as asked.
        {
            std::vector<int> newIndex;
            m.refine(toRefine, &newIndex);

            std::vector<char> crs(m.nCells(), 0);
            for (std::size_t c = 0; c < newIndex.size(); ++c)
            {
                const int w = newIndex[c];
                if (w >= 0 && w < m.nCells() && m.keys[w] == prevKeys[c])
                {
                    crs[w] = toCoarsen[c];
                }
            }
            m.unrefine(crs);
        }
        m.balance();
        m.buildAddressing();
        sKeys += wall() - t0;

        buildPrimitives
        (
            m, sidePatch, frontBack, cyclicMin, cyclicMax, nPatches,
            points, faces, owner, neighbour, patchSizes, patchStarts,
            faceEngine, faceKeyNew, pointCodeNew, tP, tF, tS
        );
        sPoints += tP;
        sFaces += tF;
        sSort += tS;

        // ---- the topology change map, built before the mesh moves ----
        t0 = wall();

        labelList cellMap, pointMap, faceMap, revCell, revFace, revPoint;
        List<objectMap> cellsFromCells, facesFromFaces, pointsFromPoints;
        labelHashSet flipFaceFlux;

        buildTopoMaps
        (
            m, prevKeys, prevCode, pointCodeNew, prevFaceKey, faceKeyNew,
            prevSorted, newSorted,
            cellMap, cellsFromCells, pointMap,
            faceMap, facesFromFaces,
            revCell, revFace, revPoint, flipFaceFlux
        );

        sMaps += wall() - t0;
        nFromCells += cellsFromCells.size();
        nFromFaces += facesFromFaces.size();
        nFlip += flipFaceFlux.size();
        forAll(faceMap, f) { if (faceMap[f] < 0) { ++nInsertedFace; } }

        const scalar Ibefore =
            sum(prevV*phiTest.primitiveField());

        t0 = wall();
        mesh.resetPrimitives
        (
            std::move(points),
            std::move(faces),
            std::move(owner),
            std::move(neighbour),
            patchSizes,
            patchStarts,
            true
        );

        // Force the cells back, as the solver's first access would
        const label nc = mesh.cells().size();
        sReset += wall() - t0;

        // polyMesh::resetPrimitives calls polyMesh::clearAddressing, which does
        // not reach fvMesh's own caches: V_, the surface interpolation weights
        // and deltaCoeffs survive and describe the mesh as it was before. Nothing
        // warns and checkMesh still passes, so they are dropped explicitly.
        // Blunt: clearOut also discards V0, which an Euler ddt scheme needs across
        // a topology change and which OpenFOAM's own topoChange maps rather than
        // throws away. That matters for the solver, not for the geometry measured
        // here.
        t0 = wall();
        mesh.clearOut();
        sClear += wall() - t0;

        // fvMesh::mapFields is public, so the mapping can be driven from a map
        // built here rather than one polyTopoChange produces. oldCellVolumes is
        // what makes coarsening volume weighted; without it cellMapper falls back
        // to a uniform average.
        t0 = wall();
        {
            autoPtr<scalarField> oldVPtr(new scalarField(prevV));
            labelListList patchPointMap(nPatches);
            labelList oldPatchNMeshPoints(nPatches, 0);

            polyTopoChangeMap tcMap
            (
                mesh,
                prevNPoints, prevNFaces, prevNCells,
                std::move(pointMap), std::move(pointsFromPoints),
                std::move(faceMap), std::move(facesFromFaces),
                std::move(cellMap), std::move(cellsFromCells),
                std::move(revPoint), std::move(revFace), std::move(revCell),
                std::move(flipFaceFlux),
                std::move(patchPointMap),
                std::move(prevPatchSizes),
                std::move(prevPatchStarts),
                std::move(oldPatchNMeshPoints),
                std::move(oldVPtr)
            );

            mesh.mapFields(tcMap);
        }
        sMapFields += wall() - t0;

        // The integral must not have moved.
        {
            const scalarField& Vn = mesh.V();
            const scalar Iafter = sum(Vn*phiTest.primitiveField());
            worstMapDrift = max
            (
                worstMapDrift, mag(Iafter - Ibefore)/max(mag(Ibefore), vSmall)
            );
        }

        // the state the next map is built against
        prevKeys = m.keys;
        prevFaceKey = faceKeyNew;
        prevSorted.swap(newSorted);
        prevCode = pointCodeNew;
        prevPatchSizes = patchSizes;
        prevPatchStarts = patchStarts;
        prevNPoints = mesh.nPoints();
        prevNFaces = mesh.nFaces();
        prevNCells = mesh.nCells();
        prevV = mesh.V();

        // ---- does fvMesh's own geometry follow? ----
        t0 = wall();

        // primitiveField(), because a GeometricField has no operator[]: its
        // internal field and its boundary field are separate objects. Each
        // accessor returns a cached reference, so taking the internal field of it
        // does not dangle.
        const scalarField& Vfv = mesh.V();
        const scalarField& magSffv = mesh.magSf().primitiveField();
        const vectorField& Sffv = mesh.Sf().primitiveField();
        const scalarField& wfv = mesh.weights().primitiveField();
        const scalarField& dcfv = mesh.deltaCoeffs().primitiveField();

        sFv += wall() - t0;

        if (Vfv.size() != m.nCells())
        {
            FatalErrorInFunction
                << "fvMesh reports " << Vfv.size() << " cells, engine has "
                << m.nCells() << abort(FatalError);
        }

        forAll(Vfv, c)
        {
            worstV = max(worstV, mag(Vfv[c] - m.V[c])/m.V[c]);
        }

        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            const label e = faceEngine[f];
            if (e < 0) { continue; }
            ++nFvChecked;

            const scalar a = m.magSf[e];
            worstMagSf = max(worstMagSf, mag(magSffv[f] - a)/a);
            worstW = max(worstW, mag(wfv[f] - m.weights[e]));
            worstDC =
                max(worstDC, mag(dcfv[f] - m.deltaCoeffs[e])/m.deltaCoeffs[e]);

            const vector se(m.Sf[2*e], m.Sf[2*e + 1], 0);
            worstSf = max(worstSf, mag(Sffv[f] - se)/a);
        }

        if (step == nSteps - 1)
        {
            Info<< "Final mesh " << nc << " cells, "
                << mesh.nFaces() << " faces ("
                << mesh.nInternalFaces() << " internal)" << nl << endl;
        }
    }

    const scalar perStep = 1e3/scalar(nSteps);

    Info<< "Per refinement check, averaged over " << nSteps << ":" << nl
        << "    keys + addressing   " << sKeys*perStep << " ms" << nl
        << "    points              " << sPoints*perStep << " ms" << nl
        << "    faces               " << sFaces*perStep << " ms" << nl
        << "    order and pack      " << sSort*perStep << " ms" << nl
        << "    resetPrimitives     " << sReset*perStep << " ms" << nl
        << "    fvMesh clearOut     " << sClear*perStep << " ms" << nl
        << "    topology maps       " << sMaps*perStep << " ms" << nl
        << "    mapFields           " << sMapFields*perStep << " ms" << nl
        << "    fvMesh geometry     " << sFv*perStep << " ms" << nl
        << "    ----" << nl
        << "    host total          "
        << (sPoints + sFaces + sSort + sReset + sClear + sFv
          + sMaps + sMapFields)*perStep << " ms" << nl
        << "    everything          "
        << (sKeys + sPoints + sFaces + sSort + sReset + sClear + sFv
          + sMaps + sMapFields)*perStep << " ms"
        << nl << endl;

    Info<< "field mapping through OpenFOAM, " << nTestFields
        << " fields, over " << nSteps << " rebuilds:" << nl
        << "    cells from several cells  " << nFromCells << nl
        << "    faces from several faces  " << nFromFaces << nl
        << "    faces with no predecessor " << nInsertedFace << nl
        << "    faces with flipped flux   " << nFlip << nl
        << "    worst drift in the integral of a mapped field  "
        << worstMapDrift << nl << endl;

    Info<< "fvMesh against the engine, over " << nFvChecked
        << " internal faces and " << nSteps << " rebuilds:" << nl
        << "    V            " << worstV << nl
        << "    magSf        " << worstMagSf << nl
        << "    Sf           " << worstSf << nl
        << "    weights      " << worstW << nl
        << "    deltaCoeffs  " << worstDC << nl
        << "  (deltaCoeffs may differ by construction: the engine uses"
        << " 1/|d|," << nl
        << "   OpenFOAM may project d onto the face normal, which differs"
        << " at a" << nl
        << "   hanging node.  V, magSf, Sf and weights must agree.)"
        << nl << endl;

    if (writeMesh)
    {
        Info<< "Writing the mesh; run checkMesh on it" << nl << endl;
        mesh.setInstance(runTime.constant());
        mesh.write();
    }

    Info<< "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
