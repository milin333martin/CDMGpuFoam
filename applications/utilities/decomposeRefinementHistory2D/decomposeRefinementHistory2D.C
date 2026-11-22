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
    decomposeRefinementHistory2D

Description
    Decomposes constant/polyMesh/refinementHistory2D (written by
    mapDetonation's pre-refinement) into the processor directories.

    decomposePar decomposes cellLevel, pointLevel and level0Edge but not
    refinementHistory2D, which it does not know.  Without it, the cells that
    were refined before decomposition have no parents on the processors and
    can never be unrefined.  Run after decomposePar:

    \verbatim
        decomposePar
        decomposeRefinementHistory2D
    \endverbatim

    For each processor the history is the subset of the serial history for
    the cells listed in processorN/constant/polyMesh/cellProcAddressing.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "labelIOList.H"
#include "refinementHistory2D.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const word regionName
    (
        args.optionFound("region")
      ? args.optionRead<word>("region")
      : polyMesh::defaultRegion
    );

    const fileName meshDir
    (
        regionName == polyMesh::defaultRegion
      ? fileName(polyMesh::meshSubDir)
      : fileName(regionName)/polyMesh::meshSubDir
    );

    const word historyName("refinementHistory2D");

    IOobject historyIO
    (
        historyName,
        runTime.constant(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!historyIO.headerOk())
    {
        Info<< "No " << historyIO.relativeObjectPath()
            << ": nothing to decompose" << nl << endl;
        return 0;
    }

    const refinementHistory2D history(historyIO);

    Info<< "Read " << historyIO.relativeObjectPath() << ": "
        << history.visibleCells().size() << " cells, "
        << history.splitCells().size() << " split entries" << nl << endl;

    label proci = 0;

    while (isDir(args.path()/("processor" + Foam::name(proci))))
    {
        const fileName procCase
        (
            args.caseName()/("processor" + Foam::name(proci))
        );

        Time procTime(Time::controlDictName, args.rootPath(), procCase);

        const labelIOList cellMap
        (
            IOobject
            (
                "cellProcAddressing",
                procTime.constant(),
                meshDir,
                procTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        autoPtr<refinementHistory2D> procHistory
        (
            history.clone
            (
                IOobject
                (
                    historyName,
                    procTime.constant(),
                    meshDir,
                    procTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                cellMap
            )
        );

        procHistory->write();

        Info<< "    processor" << proci << ": " << cellMap.size()
            << " cells, " << procHistory->splitCells().size()
            << " split entries" << endl;

        proci++;
    }

    if (proci == 0)
    {
        FatalErrorInFunction
            << "No processor directories found: run decomposePar first"
            << exit(FatalError);
    }

    Info<< nl << "Decomposed " << historyName << " into " << proci
        << " processors" << nl << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
