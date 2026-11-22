/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam.
    GNU General Public License v3 or later; see COPYING.

Application
    mapDetonation

    Use -region <name> for the fluid region of a multi-region (CHT) case;
    the source case is read from the region of the same name.

Description
    Initialises a (lab-frame) case with a quasi-steady detonation taken from
    another run, typically a shock-frame channel run: in the spirit of the
    RESIDENT workflow (run the detonation in a channel until it is
    quasi-steady, then use it as the initial condition of the configuration
    of interest).

    Steps
      1. Optional pre-refinement of a box of the (2D) target mesh with
         hexRef2D, so that the cellular structure is not lost on a coarse
         base mesh.  The refinement data (cellLevel, pointLevel, level0Edge,
         refinementHistory2D) are written so that refiner2D continues from
         this state (and can unrefine it).
      2. For every target cell centre c the source point is
             x_s = c - shift,
         optionally mirror-tiled in y inside [yMin, yMax] of the source
         (valid for sources bounded by symmetry planes), and the source cell
         containing it is found.
             x_s beyond the source in +x (ahead of the front): fresh gas
             x_s before the source in -x (behind): clamped to the source
                 inlet of the burnt side, i.e. the burnt state extended
      3. U_target = U_source + frameVelocity (source frame -> lab frame);
         the transverse velocity changes sign in mirrored copies.

    The source case must be reconstructed (reconstructPar -time <t>).

    system/mapDetonationDict:
    \verbatim
        sourceCase      "../run_T6";
        sourceTime      6e-05;
        frameVelocity   (1970.06 0 0);  // shock-frame run: (D_CJ 0 0)
        shift           (0.02 0 0);     // target = source + shift

        mirrorY         yes;            // tile source in y (symmetry walls)
        mapMethod       nearest;        // nearest | containing
        sourceRegion    region0;        // region of the source case
        behindX         0.02;           // x below which a uniform burnt
                                        // state is used (optional)

        fresh                           // state ahead of the front
        {
            p           101325;
            T           298;
            U           (0 0 0);
            Y           1;
        }

        preRefine                       // optional
        {
            min         (0.04 0);       // box in the target (x y)
            max         (0.10 0.09);
            levels      2;
        }
    \endverbatim

    Fields written: p, T, U, Y in the start time directory of the target.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshSearch.H"
#include "hexRef2D.H"
#include "polyTopoChange.H"
#include "polyTopoChangeMap.H"
#include "vector2D.H"

#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSpecifiedMeshNoChangers.H"

    IOdictionary dict
    (
        IOobject
        (
            "mapDetonationDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName sourceCase(dict.lookup<fileName>("sourceCase"));
    sourceCase.expand();
    if (!sourceCase.isAbsolute())
    {
        sourceCase = runTime.path()/sourceCase;
    }
    sourceCase.clean();

    const scalar sourceTime = dict.lookup<scalar>("sourceTime");
    const vector frameVelocity(dict.lookup<vector>("frameVelocity"));
    const vector shift(dict.lookup<vector>("shift"));
    const bool mirrorY = dict.lookupOrDefault<bool>("mirrorY", false);

    // The source is normally a single-region case, even when the target is
    // the fluid region of a multi-region (CHT) case
    const word sourceRegion
    (
        dict.lookupOrDefault<word>("sourceRegion", polyMesh::defaultRegion)
    );

    // containing: the source cell that contains the point (exact, but slow
    // on AMR meshes, where containment needs a tet decomposition)
    // nearest:    the source cell whose centre is nearest (much faster and,
    //             for cell-averaged fields of similar cell size, equivalent)
    const word mapMethod(dict.lookupOrDefault<word>("mapMethod", "nearest"));

    // Target cells with x < behindX are given a uniform burnt state, the
    // y-average of the source at its burnt end, instead of the mapped fields:
    // the products far behind the front carry cellular structure the indicator
    // would resolve at needless cost. Default: map everything.
    const scalar behindX = dict.lookupOrDefault<scalar>("behindX", -great);

    if (mapMethod != "nearest" && mapMethod != "containing")
    {
        FatalIOErrorInFunction(dict)
            << "mapMethod must be nearest or containing"
            << exit(FatalIOError);
    }

    const dictionary& freshDict = dict.subDict("fresh");
    const scalar pFresh = freshDict.lookup<scalar>("p");
    const scalar TFresh = freshDict.lookup<scalar>("T");
    const vector UFresh(freshDict.lookup<vector>("U"));
    const scalar YFresh = freshDict.lookup<scalar>("Y");


    // 1. Pre-refinement of the target mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (dict.isDict("preRefine"))
    {
        const dictionary& rDict = dict.subDict("preRefine");
        const vector2D bbMin(rDict.lookup<vector2D>("min"));
        const vector2D bbMax(rDict.lookup<vector2D>("max"));
        const label levels = rDict.lookup<label>("levels");

        const fileName oldInstance = mesh.pointsInstance();

        hexRef2D meshCutter(mesh);

        for (label level = 1; level <= levels; level++)
        {
            const vectorField& C = mesh.C().primitiveField();
            const labelList& cellLevel = meshCutter.cellLevel();

            DynamicList<label> candidates;
            forAll(C, celli)
            {
                if
                (
                    C[celli].x() >= bbMin.x() && C[celli].x() <= bbMax.x()
                 && C[celli].y() >= bbMin.y() && C[celli].y() <= bbMax.y()
                 && cellLevel[celli] < level
                )
                {
                    candidates.append(celli);
                }
            }

            // Maintain 2:1
            const labelList cellsToRefine
            (
                meshCutter.consistentRefinement(labelList(candidates), true)
            );

            Info<< "Pre-refinement to level " << level << ": refining "
                << returnReduce(cellsToRefine.size(), sumOp<label>())
                << " cells" << endl;

            polyTopoChange meshMod(mesh);
            meshCutter.setRefinement(cellsToRefine, meshMod);

            autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh);
            mesh.topoChange(map);
            meshCutter.topoChange(map);

            Info<< "    cells: " << mesh.globalData().nTotalCells() << endl;
        }

        // Overwrite the mesh in constant together with the refinement data
        mesh.setInstance(oldInstance);
        meshCutter.setInstance(oldInstance);

        Info<< nl << "Writing refined mesh and hexRef2D data to "
            << oldInstance << nl << endl;

        mesh.write();
        meshCutter.write();
    }


    // 2. Source case
    // ~~~~~~~~~~~~~~

    Info<< "Source case " << sourceCase << endl;

    Time runTimeSource
    (
        Time::controlDictName,
        sourceCase.path(),
        sourceCase.name()
    );

    const instantList sourceTimes = runTimeSource.times();
    const label timei =
        Time::findClosestTimeIndex(sourceTimes, sourceTime);

    if (timei < 0)
    {
        FatalErrorInFunction
            << "No time directories in " << sourceCase
            << " (reconstructPar first)" << exit(FatalError);
    }

    runTimeSource.setTime(sourceTimes[timei], timei);

    Info<< "    time " << runTimeSource.name() << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.name(),
            runTimeSource,
            IOobject::MUST_READ
        )
    );

    Info<< "    " << meshSource.nCells() << " cells" << endl;

    volScalarField pS
    (
        IOobject("p", runTimeSource.name(), meshSource, IOobject::MUST_READ),
        meshSource
    );
    volScalarField TS
    (
        IOobject("T", runTimeSource.name(), meshSource, IOobject::MUST_READ),
        meshSource
    );
    volVectorField US
    (
        IOobject("U", runTimeSource.name(), meshSource, IOobject::MUST_READ),
        meshSource
    );
    volScalarField YS
    (
        IOobject("Y", runTimeSource.name(), meshSource, IOobject::MUST_READ),
        meshSource
    );

    const boundBox sourceBb(meshSource.bounds());
    const scalar xMin = sourceBb.min().x();
    const scalar xMax = sourceBb.max().x();
    const scalar yMin = sourceBb.min().y();
    const scalar yMax = sourceBb.max().y();
    const scalar zMid = 0.5*(sourceBb.min().z() + sourceBb.max().z());
    const scalar H = yMax - yMin;
    const scalar tol = 1e-9*max(xMax - xMin, H);

    Info<< "    bounds " << sourceBb << nl << endl;

    // meshSearch's interface differs between OpenFOAM-13 patch releases in
    // opposite directions: TAMU HPRC's module has no meshSearch::New but
    // leaves the plain constructor public, while the 13-18870c24d21c source
    // build has New and makes that constructor protected. A derived type may
    // invoke either, so this one form builds against both with no version
    // test and no #if.
    struct cdmMeshSearch : public meshSearch
    {
        explicit cdmMeshSearch(const polyMesh& m) : meshSearch(m) {}
    };
    const cdmMeshSearch search(meshSource);

    // Uniform burnt state: average of the source cells within one cell
    // width of its burnt (low x) end
    scalar pBurnt = 0, TBurnt = 0, YBurnt = 0;
    vector UBurnt(Zero);
    if (behindX > -great)
    {
        const scalar dx = 2*(xMax - xMin)/max(label(1), label(Foam::cbrt(scalar(meshSource.nCells()))));
        const vectorField& Cs = meshSource.C().primitiveField();
        label n = 0;
        forAll(Cs, celli)
        {
            if (Cs[celli].x() < xMin + dx)
            {
                pBurnt += pS[celli]; TBurnt += TS[celli];
                YBurnt += YS[celli]; UBurnt += US[celli];
                n++;
            }
        }
        if (!n)
        {
            FatalErrorInFunction
                << "No source cells found at the burnt end" << exit(FatalError);
        }
        pBurnt /= n; TBurnt /= n; YBurnt /= n; UBurnt /= n;
        UBurnt += frameVelocity;

        Info<< "    uniform burnt state for x < " << behindX << ": p = "
            << pBurnt << ", T = " << TBurnt << ", Y = " << YBurnt
            << ", U = " << UBurnt << " (average of " << n << " cells)"
            << nl << endl;
    }


    // 3. Map
    // ~~~~~~

    volScalarField p
    (
        IOobject("p", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField T
    (
        IOobject("T", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volVectorField U
    (
        IOobject("U", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField Y
    (
        IOobject("Y", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );

    const vectorField& C = mesh.C().primitiveField();

    label nFresh = 0, nMapped = 0, nExtended = 0, nNearest = 0, nBurnt = 0;

    forAll(C, celli)
    {
        scalar xs = C[celli].x() - shift.x();
        scalar ys = C[celli].y() - shift.y();

        if (C[celli].x() < behindX)
        {
            p[celli] = pBurnt;
            T[celli] = TBurnt;
            U[celli] = UBurnt;
            Y[celli] = YBurnt;
            nBurnt++;
            continue;
        }

        if (xs > xMax)
        {
            p[celli] = pFresh;
            T[celli] = TFresh;
            U[celli] = UFresh;
            Y[celli] = YFresh;
            nFresh++;
            continue;
        }

        if (xs < xMin)
        {
            xs = xMin + tol;
            nExtended++;
        }

        // Mirror tiling in y (source bounded by symmetry planes)
        bool mirrored = false;
        if (mirrorY)
        {
            scalar yy = std::fmod(ys - yMin, 2*H);
            if (yy < 0) yy += 2*H;
            if (yy > H)
            {
                yy = 2*H - yy;
                mirrored = true;
            }
            ys = yMin + yy;
        }

        ys = min(max(ys, yMin + tol), yMax - tol);

        const point ps(xs, ys, zMid);

        label cs = -1;
        if (mapMethod == "containing")
        {
            cs = search.findCell(ps);
        }
        if (cs < 0)
        {
            cs = search.findNearestCell(ps);
            nNearest++;
        }

        p[celli] = pS[cs];
        T[celli] = TS[cs];
        Y[celli] = YS[cs];

        vector Ui(US[cs]);
        if (mirrored)
        {
            Ui.y() = -Ui.y();
        }
        U[celli] = Ui + frameVelocity;

        nMapped++;
    }

    Info<< "Method " << mapMethod << ": mapped " << nMapped << " cells ("
        << nExtended
        << " behind the source window, " << nNearest
        << " by nearest cell), " << nFresh << " fresh, "
        << nBurnt << " uniform burnt" << nl << endl;

    p.write();
    T.write();
    U.write();
    Y.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
