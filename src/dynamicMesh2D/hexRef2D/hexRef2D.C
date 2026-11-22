/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "hexRef2D.H"
#include "emptyPolyPatch.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "polyDistributionMap.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef2D, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<label value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hexRef2D::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorInFunction << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


Foam::label Foam::hexRef2D::getPatchIndex(const label facei) const
{
    label patchID = -1;

    if (!mesh_.isInternalFace(facei))
    {
        patchID = mesh_.boundaryMesh().whichPatch(facei);
    }

    return patchID;
}


Foam::label Foam::hexRef2D::addFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    const label patchID = getPatchIndex(facei);

    label newFacei = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFacei = meshMod.addFace
        (
            newFace,                    // face
            own,                        // owner
            nei,                        // neighbour
            facei,                      // master face for addition
            false,                      // flux flip
            patchID                     // patch for face
        );
    }
    else
    {
        // Reverse owner/neighbour
        newFacei = meshMod.addFace
        (
            newFace.reverseFace(),      // face
            nei,                        // owner
            own,                        // neighbour
            facei,                      // master face for addition
            false,                      // flux flip
            patchID                     // patch for face
        );
    }
    return newFacei;
}


Foam::label Foam::hexRef2D::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFacei,
    const label meshPointi,
    const face& newFace,
    const label own,
    const label nei
) const
{
    if (mesh_.isInternalFace(meshFacei))
    {
        return meshMod.addFace
        (
            newFace,                    // face
            own,                        // owner
            nei,                        // neighbour
            meshFacei,                  // master face for addition
            false,                      // flux flip
            -1                          // patch for face
        );
    }
    else
    {
        // For now create out of nothing
        //   problem: field does not get mapped.
        return meshMod.addFace
        (
            newFace,                    // face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master face for addition
            false,                      // flux flip
            -1                          // patch for face
        );
    }
}


void Foam::hexRef2D::modifyFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
) const
{
    const label patchID = getPatchIndex(facei);

    if
    (
        (own != mesh_.faceOwner()[facei])
     || (
            mesh_.isInternalFace(facei)
         && (nei != mesh_.faceNeighbour()[facei])
        )
     || (newFace != mesh_.faces()[facei])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.modifyFace
            (
                newFace,            // modified face
                facei,              // label of face being modified
                own,                // owner
                nei,                // neighbour
                false,              // face flip
                patchID             // patch for face
            );
        }
        else
        {
            meshMod.modifyFace
            (
                newFace.reverseFace(),  // modified face
                facei,                  // label of face being modified
                nei,                    // owner
                own,                    // neighbour
                false,                  // face flip
                patchID                 // patch for face
            );
        }
    }
}


Foam::scalar Foam::hexRef2D::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar great2 = sqr(great);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, great2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, celli)
        {
            const label cLevel = cellLevel_[celli];

            const labelList& cEdges = mesh_.cellEdges(celli);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // hexRef2D: only edges in the plane of the mesh carry the level
        {
            boolList isEmptyFace;
            boolList isInPlaneEdge;
            calcEmptyFacesAndInPlaneEdges(isEmptyFace, isInPlaneEdge);

            forAll(edgeLevel, edgeI)
            {
                if (!isInPlaneEdge[edgeI])
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef2D::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -great2);

    forAll(cellLevel_, celli)
    {
        const label cLevel = cellLevel_[celli];

        const labelList& cEdges = mesh_.cellEdges(celli);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef2D::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == great2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexRef2D::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < great2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexRef2D::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorInFunction
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


Foam::label Foam::hexRef2D::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label celli,
    const label facei,
    const label pointi
) const
{
    if (cellAnchorPoints[celli].size())
    {
        label index = findIndex(cellAnchorPoints[celli], pointi);

        if (index != -1)
        {
            return cellAddedCells[celli][index];
        }


        // pointi is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[facei];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[celli], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[celli][index];
            }
        }

        // Problem.
        dumpCell(celli);
        Perr<< "cell:" << celli << " anchorPoints:" << cellAnchorPoints[celli]
            << endl;

        FatalErrorInFunction
            << "Could not find point " << pointi
            << " in the anchorPoints for cell " << celli << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return celli;
    }
}


void Foam::hexRef2D::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label facei,
    const label pointi,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[facei],
        facei,
        pointi
    );

    if (mesh_.isInternalFace(facei))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[facei],
            facei,
            pointi
        );
    }
    else
    {
        nei = -1;
    }
}


Foam::label Foam::hexRef2D::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


Foam::label Foam::hexRef2D::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexRef2D::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexRef2D::dumpCell(const label celli) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(celli) + ".obj");
    Pout<< "hexRef2D : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[celli];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointi = f[fp];
            label nexPointi = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointi]+1
                << ' ' << pointToObjVert[nexPointi]+1 << nl;
        }
    }
}


Foam::label Foam::hexRef2D::findLevel
(
    const label facei,
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointi = f[fp];

        if (pointLevel_[pointi] < wantedLevel)
        {
            dumpCell(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                dumpCell(mesh_.faceNeighbour()[facei]);
            }

            FatalErrorInFunction
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointi] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    dumpCell(mesh_.faceOwner()[facei]);
    if (mesh_.isInternalFace(facei))
    {
        dumpCell(mesh_.faceNeighbour()[facei]);
    }

    FatalErrorInFunction
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


Foam::label Foam::hexRef2D::faceLevel(const label facei) const
{
    const face& f = mesh_.faces()[facei];

    if (f.size() <= 4)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[facei]];

        if (countAnchors(f, ownLevel) == 4)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 4)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::hexRef2D::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    const face compactFace(identityMap(newFace.size()));
    const pointField compactPoints(meshMod.points(), newFace);

    const vector a(compactFace.area(compactPoints));

    const vector dir(neiPt - ownPt);

    if ((dir & a) < 0)
    {
        FatalErrorInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    const scalar s = (fcToOwn & a) / (dir & a);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexRef2D::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label celli,
    const label facei,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    const face compactFace(identityMap(newFace.size()));
    const pointField compactPoints(meshMod.points(), newFace);

    const vector a(compactFace.area(compactPoints));

    const vector dir(boundaryPt - ownPt);

    if ((dir & a) < 0)
    {
        FatalErrorInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    const vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    const scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningInFunction
            << "cell:" << celli << " old face:" << facei
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
    }
}


void Foam::hexRef2D::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


Foam::label Foam::hexRef2D::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label celli,
    const label facei,
    const bool faceOrder,
    const label edgeMidPointi,
    const label anchorPointi,
    const label faceMidPointi,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointi);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointi, edge(anchorPointi, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointi != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointi;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointi);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointi, edge(faceMidPointi, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointi != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointi;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointi];
        const edge& faceMids = midPointToFaceMids[edgeMidPointi];

        label otherFaceMidPointi = faceMids.otherVertex(faceMidPointi);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cellMid and face mids.

        DynamicList<label> newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[facei] == celli))
        {
            newFaceVerts.append(faceMidPointi);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointi,  // edge between faceMid and
                edgeMidPointi,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                otherFaceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointi);
            newFaceVerts.append(cellMidPoint[celli]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointi,
                edgeMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointi);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointi,
                faceMidPointi,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointi);
            newFaceVerts.append(cellMidPoint[celli]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchorPointi
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            celli,
            facei,
            anchors.otherVertex(anchorPointi)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointi];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            neiPt = mesh_.points()[anchorPointi];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointi];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointi)];
                neiPt = mesh_.points()[anchorPointi];
            }

            checkInternalOrientation
            (
                meshMod,
                celli,
                facei,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            facei,
            anchorPointi,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


void Foam::hexRef2D::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label celli,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[celli];
    const label cLevel = cellLevel_[celli];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(24);
    // From edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        const face& f = mesh_.faces()[facei];
        const labelList& fEdges = mesh_.faceEdges(facei, storage);

        // We are on the celli side of face f. The face will have 1 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointi = -1;

        label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel+1 point
            label edgeMid = findLevel
            (
                facei,
                f,
                f.fcIndex(anchorFp),
                true,
                cLevel+1
            );
            label faceMid = findLevel
            (
                facei,
                f,
                f.fcIndex(edgeMid),
                true,
                cLevel+1
            );

            faceMidPointi = f[faceMid];
        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting.

            faceMidPointi = faceMidPoint[facei];
        }
        else
        {
            dumpCell(mesh_.faceOwner()[facei]);
            if (mesh_.isInternalFace(facei))
            {
                dumpCell(mesh_.faceNeighbour()[facei]);
            }

            FatalErrorInFunction
                << "nAnchors:" << nAnchors
                << " facei:" << facei
                << abort(FatalError);
        }



        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all.

        forAll(f, fp0)
        {
            label point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // Anchor.

                // Walk forward
                // ~~~~~~~~~~~~
                // to cLevel+1 or edgeMidPoint of this level.


                label edgeMidPointi = -1;

                label fp1 = f.fcIndex(fp0);

                if (pointLevel_[f[fp1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fp0];

                    edgeMidPointi = edgeMidPoint[edgeI];

                    if (edgeMidPointi == -1)
                    {
                        dumpCell(celli);

                        const labelList& cPoints = mesh_.cellPoints(celli);

                        FatalErrorInFunction
                            << "cell:" << celli << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << facei
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[facei]
                            << " faceMidPoint:" << faceMidPoint[facei]
                            << " faceMidPointi:" << faceMidPointi
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel+1
                    label edgeMid = findLevel(facei, f, fp1, true, cLevel+1);

                    edgeMidPointi = f[edgeMid];
                }

                label newFacei = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    celli,
                    facei,
                    true,                   // mid point after anchor
                    edgeMidPointi,          // edgemid
                    point0,                 // anchor
                    faceMidPointi,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFacei != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }



                // Walk backward
                // ~~~~~~~~~~~~~

                label fpMin1 = f.rcIndex(fp0);

                if (pointLevel_[f[fpMin1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fpMin1];

                    edgeMidPointi = edgeMidPoint[edgeI];

                    if (edgeMidPointi == -1)
                    {
                        dumpCell(celli);

                        const labelList& cPoints = mesh_.cellPoints(celli);

                        FatalErrorInFunction
                            << "cell:" << celli << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << facei
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[facei]
                            << " faceMidPoint:" << faceMidPoint[facei]
                            << " faceMidPointi:" << faceMidPointi
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back to clevel+1
                    label edgeMid = findLevel
                    (
                        facei,
                        f,
                        fpMin1,
                        false,
                        cLevel+1
                    );

                    edgeMidPointi = f[edgeMid];
                }

                newFacei = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    celli,
                    facei,
                    false,                  // mid point before anchor
                    edgeMidPointi,          // edgemid
                    point0,                 // anchor
                    faceMidPointi,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFacei != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }
            }   // done anchor
        }   // done face

        if (nFacesAdded == 12)
        {
            break;
        }
    }
}


void Foam::hexRef2D::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


void Foam::hexRef2D::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label facei,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges(facei);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


Foam::label Foam::hexRef2D::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


void Foam::hexRef2D::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[facei];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorInFunction
                << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label facei = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);

            dumpCell(own);
            FatalErrorInFunction
                << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << facei
                << " on patch " << patchi << " "
                << mesh_.boundaryMesh()[patchi].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


void Foam::hexRef2D::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexRef2D::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


void Foam::hexRef2D::collectLevelPoints
(
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= level)
        {
            points.append(f[fp]);
        }
    }
}


void Foam::hexRef2D::collectLevelPoints
(
    const labelList& meshPoints,
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        label pointi = meshPoints[f[fp]];
        if (pointLevel_[pointi] <= level)
        {
            points.append(pointi);
        }
    }
}


bool Foam::hexRef2D::matchHexShape
(
    const label celli,
    const label cellLevel,
    DynamicList<face>& quads
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    // Work arrays
    DynamicList<label> verts(4);
    quads.clear();


    // 1. pick up any faces with four cellLevel points

    forAll(cFaces, i)
    {
        label facei = cFaces[i];
        const face& f = mesh_.faces()[facei];

        verts.clear();
        collectLevelPoints(f, cellLevel, verts);
        if (verts.size() == 4)
        {
            if (mesh_.faceOwner()[facei] != celli)
            {
                reverse(verts);
            }
            quads.append(face(0));
            labelList& quadVerts = quads.last();
            quadVerts.transfer(verts);
        }
    }


    if (quads.size() < 6)
    {
        Map<labelList> pointFaces(2*cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];
            const face& f = mesh_.faces()[facei];

            // Pick up any faces with only one level point.
            // See if there are four of these where the common point
            // is a level+1 point. This common point is then the mid of
            // a split face.

            verts.clear();
            collectLevelPoints(f, cellLevel, verts);
            if (verts.size() == 1)
            {
                // Add to pointFaces for any level+1 point (this might be
                // a midpoint of a split face)
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (pointLevel_[pointi] == cellLevel+1)
                    {
                        Map<labelList>::iterator iter =
                            pointFaces.find(pointi);
                        if (iter != pointFaces.end())
                        {
                            labelList& pFaces = iter();
                            if (findIndex(pFaces, facei) == -1)
                            {
                                pFaces.append(facei);
                            }
                        }
                        else
                        {
                            pointFaces.insert
                            (
                                pointi,
                                labelList(1, facei)
                            );
                        }
                    }
                }
            }
        }

        // 2. Check if we've collected any midPoints.
        forAllConstIter(Map<labelList>, pointFaces, iter)
        {
            const labelList& pFaces = iter();

            if (pFaces.size() == 4)
            {
                // Collect and orient.
                faceList fourFaces(pFaces.size());
                forAll(pFaces, pFacei)
                {
                    label facei = pFaces[pFacei];
                    const face& f = mesh_.faces()[facei];
                    if (mesh_.faceOwner()[facei] == celli)
                    {
                        fourFaces[pFacei] = f;
                    }
                    else
                    {
                        fourFaces[pFacei] = f.reverseFace();
                    }
                }

                primitivePatch bigFace
                (
                    SubList<face>(fourFaces, fourFaces.size()),
                    mesh_.points()
                );
                const labelListList& edgeLoops = bigFace.edgeLoops();

                if (edgeLoops.size() == 1)
                {
                    // Collect the 4 cellLevel points
                    verts.clear();
                    collectLevelPoints
                    (
                        bigFace.meshPoints(),
                        bigFace.edgeLoops()[0],
                        cellLevel,
                        verts
                    );

                    if (verts.size() == 4)
                    {
                        quads.append(face(0));
                        labelList& quadVerts = quads.last();
                        quadVerts.transfer(verts);
                    }
                }
            }
        }
    }

    return (quads.size() == 6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRef2D::hexRef2D(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementHistory2D",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // All cells visible if not read or readHistory = false
        (readHistory ? mesh_.nCells() : 0)
    ),
    faceRemover_(mesh_, great),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells "
            << history_.visibleCells().size() << " in "
            << history_.objectPath()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "cellLevel file " << cellLevel_.objectPath() << endl
            << "pointLevel file " << pointLevel_.objectPath() << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    // if (debug)
    {
        checkMesh();
    }
}


Foam::hexRef2D::hexRef2D
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory2D& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar
        (
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory2D",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        history
    ),
    faceRemover_(mesh_, great),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    // if (debug)
    {
        checkMesh();
    }
}


Foam::hexRef2D::hexRef2D
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar
        (
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory2D",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<refinementHistory2D::splitCell8>(0),
        labelList(0),
        false
    ),
    faceRemover_(mesh_, great),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    // if (debug)
    {
        checkMesh();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hexRef2D::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : deselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef2D::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell.get(celli))
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef2D::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialise it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, celli)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[celli] = refinementData
        (
            maxFaceDiff*(cellLevel_[celli]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[celli]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].count() = allCellInfo[celli].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[facei]];

        if (mesh_.isInternalFace(facei))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[facei]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(facei);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[facei] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces in between different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];
            label nei = faceNeighbour[facei];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(facei);
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                facei,
                own,
                allCellInfo[own],
                FaceCellWave<refinementData, int>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(facei);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexRef2D::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointi)
        {
            label& pLevel = maxPointCount[pointi];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                refinementData& cellInfo = allCellInfo[celli];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointi]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            facei,
                            celli,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[facei].count())
                        {
                            changedFacesInfo.insert(facei, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorInFunction
                    << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << facei << " faceRefData:" << allFaceInfo[facei]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << facei
                    << " refData:" << allFaceInfo[facei]
                    << " on patch " << patchi << " "
                    << mesh_.boundaryMesh()[patchi].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].isRefined())
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef2D::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef2D::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorInFunction
            << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[celli],
            cellLevel_[celli]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, celli)
    {
        if (!allCellInfo[celli].valid(dummyTrackData))
        {
            allCellInfo[celli] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[celli],
                cellLevel_[celli]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label facei = facesToCheck[i];

        if (allFaceInfo[facei].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorInFunction
                << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << facei << " occurs at positions "
                << findIndices(facesToCheck, facei)
                << abort(FatalError);
        }

        label own = faceOwner[facei];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(facei))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[facei];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[facei].updateFace
            (
                mesh_,
                facei,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(facei);
        seedFacesInfo.append(allFaceInfo[facei]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces in between different refinement levels for now
    forAll(faceNeighbour, facei)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            label own = faceOwner[facei];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[facei];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(facei);
                allFaceInfo[facei].updateFace
                (
                    mesh_,
                    facei,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[facei]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    // if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().name(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::NO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar(dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, celli)
    //    {
    //        wantedLevel[celli] = allCellInfo[celli].wantedLevel(cc[celli]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        allCellInfo[celli].originLevel() = sizeof(label)*8-2;
        allCellInfo[celli].origin() = cc[celli];
    }

    // 2. Extend to 2:1. I don't understand yet why this is not done
    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, celli)
    {
        label wanted = allCellInfo[celli].wantedLevel(cc[celli]);

        if (wanted > cellLevel_[celli]+1)
        {
            refineCell.set(celli);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef2D::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell[celli])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, celli)
    {
        if (refineCell[celli])
        {
            newCellsToRefine[nRefined++] = celli;
        }
    }

    if (debug)
    {
        Pout<< "hexRef2D::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexRef2D::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexRef2D::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli))
                {
                    cellsOut2.insert(celli);
                }
            }
            Pout<< "hexRef2D::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, celli)
            {
                if (refineCell.get(celli) && !savedRefineCell.get(celli))
                {
                    dumpCell(celli);
                    FatalErrorInFunction
                        << "Cell:" << celli << " cc:"
                        << mesh_.cellCentres()[celli]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}




void Foam::hexRef2D::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointi = pointsToStore[i];
        savedPointLevel_.insert(pointi, pointLevel_[pointi]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label celli = cellsToStore[i];
        savedCellLevel_.insert(celli, cellLevel_[celli]);
    }
}


void Foam::hexRef2D::topoChange(const polyTopoChangeMap& map)
{
    Map<label> dummyMap(0);

    topoChange(map, dummyMap, dummyMap, dummyMap);
}


void Foam::hexRef2D::topoChange
(
    const polyTopoChangeMap& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef2D::topoChange :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef2D::topoChange :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef2D that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCelli)
            {
                label oldCelli = cellMap[newCelli];

                if (oldCelli == -1)
                {
                    newCellLevel[newCelli] = -1;
                }
                else
                {
                    newCellLevel[newCelli] = cellLevel_[oldCelli];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCelli = iter.key();
            label storedCelli = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCelli);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new cell "
                    << newCelli << " but cannot find old cell " << storedCelli
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCelli] = fnd();
        }
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef2D that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];

                if (oldPointi == -1)
                {
                    newPointLevel[newPointi] = -1;
                }
                else
                {
                    newPointLevel[newPointi] = pointLevel_[oldPointi];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointi = iter.key();
            label storedPointi = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointi);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorInFunction
                    << "Problem : trying to restore old value for new point "
                    << newPointi << " but cannot find old point "
                    << storedPointi
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointi] = fnd();
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.topoChange(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.topoChange(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef2D::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef2D::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningInFunction
            << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCelli)
        {
            newCellLevel[newCelli] = cellLevel_[cellMap[newCelli]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointi)
        {
            newPointLevel[newPointi] = pointLevel_[pointMap[newPointi]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorInFunction
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    // faceRemover_.subset(pointMap, faceMap, cellMap);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef2D::distribute(const polyDistributionMap& map)
{
    if (debug)
    {
        Pout<< "hexRef2D::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);

    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef2D::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexRef2D::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<>>
                    cellToFace(2*pp.size());

                label facei = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[facei];
                    label bFacei = facei-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFacei]), facei))
                    {
                        dumpCell(own);
                        FatalErrorInFunction
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << facei
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFacei]
                            << " has two faces connected to it:"
                            << facei << " and "
                            << cellToFace[labelPair(own, nei[bFacei])]
                            << abort(FatalError);
                    }

                    facei++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mesh_.magFaceAreas()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const scalar magArea = mesh_.magFaceAreas()[facei];

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[facei];
                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                dumpCell(mesh_.faceOwner()[facei]);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label facei = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[facei];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label facei = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[facei];
            const face& f = mesh_.faces()[facei];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[facei]);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << facei
                    << " on patch " << patchi
                    << " " << mesh_.boundaryMesh()[patchi].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef2D::checkMesh : Returning" << endl;
    }
}


void Foam::hexRef2D::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexRef2D::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorInFunction
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << facei << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[facei];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchi = mesh_.boundaryMesh().whichPatch(facei);

                FatalErrorInFunction
                    << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << facei
                    << " on patch " << patchi << " "
                    << mesh_.boundaryMesh()[patchi].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronised
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointi)
        {
            if (pointLevel_[pointi] != syncPointLevel[pointi])
            {
                FatalErrorInFunction
                    << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointi << " coord:" << mesh_.points()[pointi]
                    << " has level " << pointLevel_[pointi]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointi]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointi)
        {
            const labelList& pCells = mesh_.pointCells(pointi);

            label& pLevel = maxPointLevel[pointi];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointi = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointi);

            forAll(pCells, i)
            {
                label celli = pCells[i];

                if
                (
                    mag(cellLevel_[celli]-maxPointLevel[pointi])
                  > maxPointDiff
                )
                {
                    dumpCell(celli);

                    FatalErrorInFunction
                        << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << celli
                        << " cellLevel:" << cellLevel_[celli]
                        << " uses point:" << pointi
                        << " coord:" << mesh_.points()[pointi]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointi]
                        << abort(FatalError);
                }
            }
        }
    }


    //- Gives problems after first splitting off inside mesher.
    //// Hanging points
    //{
    //    // Any patches with points having only two edges.
    //
    //    boolList isHangingPoint(mesh_.nPoints(), false);
    //
    //    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    //
    //    forAll(patches, patchi)
    //    {
    //        const polyPatch& pp = patches[patchi];
    //
    //        const labelList& meshPoints = pp.meshPoints();
    //
    //        forAll(meshPoints, i)
    //        {
    //            label pointi = meshPoints[i];
    //
    //            const labelList& pEdges = mesh_.pointEdges()[pointi];
    //
    //            if (pEdges.size() == 2)
    //            {
    //                isHangingPoint[pointi] = true;
    //            }
    //        }
    //    }
    //
    //    syncTools::syncPointList
    //    (
    //        mesh_,
    //        isHangingPoint,
    //        andEqOp<bool>(),        // only if all decide it is hanging point
    //        true,                   // null
    //        false                   // no separation
    //    );
    //
    //    // OFstream str(mesh_.time().path()/"hangingPoints.obj");
    //
    //    label nHanging = 0;
    //
    //    forAll(isHangingPoint, pointi)
    //    {
    //        if (isHangingPoint[pointi])
    //        {
    //            nHanging++;
    //
    //            Pout<< "Hanging boundary point " << pointi
    //                << " at " << mesh_.points()[pointi]
    //                << endl;
    //            // meshTools::writeOBJ(str, mesh_.points()[pointi]);
    //        }
    //    }
    //
    //    if (returnReduce(nHanging, sumOp<label>()) > 0)
    //    {
    //        FatalErrorInFunction
    //            << "Detected a point used by two edges only (hanging point)"
    //            << nl << "This is not allowed"
    //            << abort(FatalError);
    //    }
    //}
}


const Foam::cellShapeList& Foam::hexRef2D::cellShapes() const
{
    if (cellShapesPtr_.empty())
    {
        if (debug)
        {
            Pout<< "hexRef2D::cellShapes() : calculating splitHex cellShapes."
                << " cellLevel:" << cellLevel_.size()
                << " pointLevel:" << pointLevel_.size()
                << endl;
        }

        const cellShapeList& meshShapes = mesh_.cellShapes();
        cellShapesPtr_.reset(new cellShapeList(meshShapes));

        label nSplitHex = 0;
        label nUnrecognised = 0;

        forAll(cellLevel_, celli)
        {
            if (meshShapes[celli].model().index() == 0)
            {
                label level = cellLevel_[celli];

                // Return true if we've found 6 quads
                DynamicList<face> quads;
                bool haveQuads = matchHexShape
                (
                    celli,
                    level,
                    quads
                );

                if (haveQuads)
                {
                    faceList faces(move(quads));
                    cellShapesPtr_()[celli] = degenerateMatcher::match(faces);
                    nSplitHex++;
                }
                else
                {
                    nUnrecognised++;
                }
            }
        }
        if (debug)
        {
            Pout<< "hexRef2D::cellShapes() :"
                << " nCells:" << mesh_.nCells() << " of which" << nl
                << "    primitive:" << (mesh_.nCells()-nSplitHex-nUnrecognised)
                << nl
                << "    split-hex:" << nSplitHex << nl
                << "    poly     :" << nUnrecognised << nl
                << endl;
        }
    }
    return cellShapesPtr_();
}




Foam::labelList Foam::hexRef2D::consistentUnrefinement
(
    const labelList& pointsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef2D::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorInFunction
            << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // Loop, modifying pointsToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : deselect points to refine
    // maxSet = true: select points to refine

    // Maintain boolList for pointsToUnrefine and cellsToUnrefine
    PackedBoolList unrefinePoint(mesh_.nPoints());

    forAll(pointsToUnrefine, i)
    {
        label pointi = pointsToUnrefine[i];

        unrefinePoint.set(pointi);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefinePoint, pointi)
        {
            if (unrefinePoint.get(pointi))
            {
                const labelList& pCells = mesh_.pointCells(pointi);

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[facei];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorInFunction
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label facei = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[facei];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorInFunction
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef2D::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any point whose cell neighbour cannot be unrefined.
        forAll(unrefinePoint, pointi)
        {
            if (unrefinePoint.get(pointi))
            {
                const labelList& pCells = mesh_.pointCells(pointi);

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefinePoint.unset(pointi);
                        break;
                    }
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefinePoint, pointi)
    {
        if (unrefinePoint.get(pointi))
        {
            nSet++;
        }
    }

    labelList newPointsToUnrefine(nSet);
    nSet = 0;

    forAll(unrefinePoint, pointi)
    {
        if (unrefinePoint.get(pointi))
        {
            newPointsToUnrefine[nSet++] = pointi;
        }
    }

    return newPointsToUnrefine;
}




bool Foam::hexRef2D::write(const bool write) const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Size of cellLevel:" << cellLevel_.size()
            << " does not equal number of cells in mesh:" << mesh_.nCells()
            << abort(FatalError);
    }

    if (pointLevel_.size() != mesh_.nPoints())
    {
        FatalErrorInFunction
            << "Size of pointLevel:" << pointLevel_.size()
            << " does not equal number of points in mesh:" << mesh_.nPoints()
            << abort(FatalError);
    }

    bool writeOk =
        cellLevel_.write(write)
     && pointLevel_.write(write)
     && level0Edge_.write(write);

    if (history_.active())
    {
        writeOk = writeOk && history_.write(write);
    }

    return writeOk;
}


// * * * * * * * * * * * * 2D helpers (CDMGpuFoam)  * * * * * * * * * * * //

void Foam::hexRef2D::check2D() const
{
    if (mesh_.nGeometricD() != 2)
    {
        FatalErrorInFunction
            << "hexRef2D requires a 2D mesh: a single layer of hexahedra"
            << " between empty patches." << nl
            << "    The mesh has " << mesh_.nGeometricD()
            << " geometric directions." << exit(FatalError);
    }
}


void Foam::hexRef2D::calcEmptyFacesAndInPlaneEdges
(
    boolList& isEmptyFace,
    boolList& isInPlaneEdge
) const
{
    isEmptyFace.setSize(mesh_.nFaces());
    isEmptyFace = false;

    isInPlaneEdge.setSize(mesh_.nEdges());
    isInPlaneEdge = false;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        if (isA<emptyPolyPatch>(patches[patchi]))
        {
            const polyPatch& pp = patches[patchi];

            forAll(pp, i)
            {
                const label facei = pp.start() + i;

                isEmptyFace[facei] = true;

                const labelList& fEdges = mesh_.faceEdges(facei);

                forAll(fEdges, j)
                {
                    isInPlaneEdge[fEdges[j]] = true;
                }
            }
        }
    }
}


Foam::label Foam::hexRef2D::normalPartner
(
    const label pointi,
    const boolList& isInPlaneEdge
) const
{
    const labelList& pEdges = mesh_.pointEdges()[pointi];

    forAll(pEdges, i)
    {
        if (!isInPlaneEdge[pEdges[i]])
        {
            return mesh_.edges()[pEdges[i]].otherVertex(pointi);
        }
    }

    FatalErrorInFunction
        << "Point " << pointi << " at " << mesh_.points()[pointi]
        << " has no edge normal to the plane of the mesh." << nl
        << "    hexRef2D requires a single layer of cells between"
        << " empty patches." << abort(FatalError);

    return -1;
}


Foam::label Foam::hexRef2D::faceEdge
(
    const label facei,
    const label p,
    const label q
) const
{
    const labelList& fEdges = mesh_.faceEdges(facei);
    const edgeList& edges = mesh_.edges();

    forAll(fEdges, i)
    {
        const edge& e = edges[fEdges[i]];

        if ((e[0] == p && e[1] == q) || (e[0] == q && e[1] == p))
        {
            return fEdges[i];
        }
    }

    return -1;
}


void Foam::hexRef2D::extendedFace
(
    const label facei,
    const labelList& edgeMidPoint,
    const boolList& isInPlaneEdge,
    DynamicList<label>& pts,
    DynamicList<bool>& isMid
) const
{
    const face& f = mesh_.faces()[facei];

    pts.clear();
    isMid.clear();

    forAll(f, fp)
    {
        const label p = f[fp];
        const label prev = f[f.rcIndex(fp)];
        const label next = f[f.fcIndex(fp)];

        const label ePrev = faceEdge(facei, prev, p);
        const label eNext = faceEdge(facei, p, next);

        // An existing point is a mid-point of an in-plane edge if both of
        // its edges in this face lie in the plane (meaningful for faces
        // normal to the plane only)
        pts.append(p);
        isMid.append
        (
            ePrev >= 0 && eNext >= 0
         && isInPlaneEdge[ePrev] && isInPlaneEdge[eNext]
        );

        if (eNext >= 0 && edgeMidPoint[eNext] >= 0)
        {
            pts.append(edgeMidPoint[eNext]);
            isMid.append(true);
        }
    }
}


// * * * * * * * * * * * * 2D refinement (CDMGpuFoam) * * * * * * * * * * //

Foam::labelListList Foam::hexRef2D::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    check2D();

    if (debug)
    {
        Pout<< "hexRef2D::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();

    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, celli)
    {
        newCellLevel.append(cellLevel_[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointi)
    {
        newPointLevel.append(pointLevel_[pointi]);
    }

    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const edgeList& edges = mesh_.edges();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const pointField& points = mesh_.points();
    const label nOldPoints = mesh_.nPoints();

    boolList isEmptyFace;
    boolList isInPlaneEdge;
    calcEmptyFacesAndInPlaneEdges(isEmptyFace, isInPlaneEdge);

    boolList refineCell(mesh_.nCells(), false);
    forAll(cellLabels, i)
    {
        refineCell[cellLabels[i]] = true;
    }

    // Positions of the added points (for orienting new faces)
    Map<point> addedPoints(8*cellLabels.size());

    auto pointPos = [&](const label pointi)
    {
        return pointi < nOldPoints ? points[pointi] : addedPoints[pointi];
    };


    // 1. Split the in-plane edges of the cells to refine
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const labelList& cEdges = mesh_.cellEdges(celli);

        forAll(cEdges, j)
        {
            const label edgei = cEdges[j];
            const edge& e = edges[edgei];

            if
            (
                isInPlaneEdge[edgei]
             && pointLevel_[e[0]] <= cellLevel_[celli]
             && pointLevel_[e[1]] <= cellLevel_[celli]
            )
            {
                edgeMidPoint[edgei] = 12345;    // mark need for splitting
            }
        }
    }

    // Any split takes precedence across coupled patches
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );

    {
        pointField edgeMids(mesh_.nEdges(), point(-great, -great, -great));

        forAll(edgeMidPoint, edgei)
        {
            if (edgeMidPoint[edgei] >= 0)
            {
                edgeMids[edgei] = edges[edgei].centre(points);
            }
        }

        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-great, -great, -great)
        );

        forAll(edgeMidPoint, edgei)
        {
            if (edgeMidPoint[edgei] >= 0)
            {
                const edge& e = edges[edgei];

                edgeMidPoint[edgei] = meshMod.addPoint
                (
                    edgeMids[edgei],    // point
                    e[0],               // master point
                    true                // supports a cell
                );

                addedPoints.insert(edgeMidPoint[edgei], edgeMids[edgei]);

                newPointLevel(edgeMidPoint[edgei]) =
                    max(pointLevel_[e[0]], pointLevel_[e[1]]) + 1;
            }
        }
    }


    // 2. Mid-points of the two empty faces of every cell to refine
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList faceMidPoint(mesh_.nFaces(), -1);
    labelListList cellEmptyFaces(mesh_.nCells());

    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const cell& cFaces = cells[celli];

        labelList& ef = cellEmptyFaces[celli];

        forAll(cFaces, j)
        {
            if (isEmptyFace[cFaces[j]])
            {
                ef.append(cFaces[j]);
            }
        }

        if (ef.size() != 2)
        {
            FatalErrorInFunction
                << "Cell " << celli << " at " << mesh_.cellCentres()[celli]
                << " has " << ef.size() << " faces on empty patches"
                << " instead of 2." << nl
                << "    hexRef2D requires a single layer of cells between"
                << " empty patches." << abort(FatalError);
        }

        forAll(ef, j)
        {
            const label facei = ef[j];
            const point fc = mesh_.faceCentres()[facei];

            faceMidPoint[facei] = meshMod.addPoint
            (
                fc,                     // point
                faces[facei][0],        // master point
                true                    // supports a cell
            );

            addedPoints.insert(faceMidPoint[facei], fc);

            newPointLevel(faceMidPoint[facei]) = cellLevel_[celli] + 1;
        }
    }


    // 3. Split every cell into 4: faces on the empty patches and the 4 new
    //    internal faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per refined cell the 4 cells (original at 0)
    labelListList cellAddedCells(mesh_.nCells());

    // Per refined cell the boundary points on its empty faces that belong
    // to a single child (corners and the finer points on the half-edges,
    // not the edge mid-points) and the index of that child
    labelListList cellChildPoints(mesh_.nCells());
    labelListList cellChildIndex(mesh_.nCells());

    DynamicList<label> pts(16);
    DynamicList<bool> isMid(16);

    // Analyse a loop of an empty face of a cell of the given level:
    // positions of the 4 corners (level <= level) and of the mid-point of
    // each side (the lowest-level point between consecutive corners)
    auto analyseLoop = [&]
    (
        const label celli,
        const label facei,
        const labelList& loop,
        const label level,
        labelList& anchorPos,
        labelList& midPos
    )
    {
        const label n = loop.size();

        anchorPos.clear();
        forAll(loop, j)
        {
            if (newPointLevel[loop[j]] <= level)
            {
                anchorPos.append(j);
            }
        }

        if (anchorPos.size() != 4)
        {
            FatalErrorInFunction
                << "Face " << facei << " of cell " << celli
                << " (level " << level << ") at "
                << mesh_.faceCentres()[facei] << " has "
                << anchorPos.size() << " corner points instead of 4."
                << nl << "    Points after edge splitting: " << loop
                << nl << "    Inconsistent refinement levels."
                << abort(FatalError);
        }

        midPos.setSize(4);

        forAll(anchorPos, t)
        {
            const label j0 = anchorPos[t];
            const label j1 = anchorPos[(t + 1) % 4];

            label best = -1;
            label nBest = 0;

            for (label j = (j0 + 1) % n; j != j1; j = (j + 1) % n)
            {
                if
                (
                    best == -1
                 || newPointLevel[loop[j]] < newPointLevel[loop[best]]
                )
                {
                    best = j;
                    nBest = 1;
                }
                else if (newPointLevel[loop[j]] == newPointLevel[loop[best]])
                {
                    nBest++;
                }
            }

            if (best == -1 || nBest != 1)
            {
                FatalErrorInFunction
                    << "Face " << facei << " of cell " << celli
                    << " (level " << level << ") at "
                    << mesh_.faceCentres()[facei]
                    << ": no unique mid-point between corners "
                    << loop[j0] << " and " << loop[j1] << nl
                    << "    Points after edge splitting: " << loop
                    << abort(FatalError);
            }

            midPos[t] = best;
        }
    };

    // Sub-face around the corner t: from the mid-point before it to the
    // mid-point after it (loop order), closed through the face mid-point
    auto subFace = [&]
    (
        const labelList& loop,
        const labelList& midPos,
        const label t,
        const label fm
    )
    {
        const label n = loop.size();

        DynamicList<label> f;
        for (label j = midPos[(t + 3) % 4]; ; j = (j + 1) % n)
        {
            f.append(loop[j]);
            if (j == midPos[t]) break;
        }
        f.append(fm);

        return face(f);
    };

    // Record the points of loop between the mid-points around corner t as
    // belonging to child k
    auto storeChildPoints = [&]
    (
        const label celli,
        const labelList& loop,
        const labelList& midPos,
        const label t,
        const label k
    )
    {
        const label n = loop.size();

        for
        (
            label j = (midPos[(t + 3) % 4] + 1) % n;
            j != midPos[t];
            j = (j + 1) % n
        )
        {
            cellChildPoints[celli].append(loop[j]);
            cellChildIndex[celli].append(k);
        }
    };

    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const label level = cellLevel_[celli];
        const label faceA = cellEmptyFaces[celli][0];
        const label faceB = cellEmptyFaces[celli][1];

        // Face A defines the numbering of the children: child k at corner k
        extendedFace(faceA, edgeMidPoint, isInPlaneEdge, pts, isMid);
        const labelList loopA(pts);

        labelList anchorPosA;
        labelList midPosA;
        analyseLoop(celli, faceA, loopA, level, anchorPosA, midPosA);

        FixedList<label, 4> a;
        FixedList<label, 4> mA;
        FixedList<label, 4> b;
        forAll(a, k)
        {
            a[k] = loopA[anchorPosA[k]];
            mA[k] = loopA[midPosA[k]];
            b[k] = normalPartner(a[k], isInPlaneEdge);
        }

        // Face B in its own orientation
        extendedFace(faceB, edgeMidPoint, isInPlaneEdge, pts, isMid);
        const labelList loopB(pts);

        labelList anchorPosB;
        labelList midPosB;
        analyseLoop(celli, faceB, loopB, level, anchorPosB, midPosB);

        // Child of each corner of B, and the B mid-point matching mA[k]
        FixedList<label, 4> childOfB;
        forAll(anchorPosB, t)
        {
            childOfB[t] = findIndex(b, loopB[anchorPosB[t]]);

            if (childOfB[t] < 0)
            {
                FatalErrorInFunction
                    << "Corner " << loopB[anchorPosB[t]]
                    << " of face " << faceB << " of cell " << celli
                    << " is not the partner of a corner of face " << faceA
                    << ": " << a << " " << b << abort(FatalError);
            }
        }

        FixedList<label, 4> mB;
        forAll(anchorPosB, t)
        {
            // The side between B corners t and t+1 matches side k of A,
            // where {k, k+1} = {child(t), child(t+1)}
            const label c0 = childOfB[t];
            const label c1 = childOfB[(t + 1) % 4];

            label k = -1;
            if ((c0 + 1) % 4 == c1)
            {
                k = c0;
            }
            else if ((c1 + 1) % 4 == c0)
            {
                k = c1;
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot match the faces on the empty patches of cell "
                    << celli << ": " << loopA << " " << loopB
                    << abort(FatalError);
            }

            mB[k] = loopB[midPosB[t]];
        }


        // Children
        labelList& cAdded = cellAddedCells[celli];
        cAdded.setSize(4);
        cAdded[0] = celli;
        newCellLevel[celli] = level + 1;

        for (label k = 1; k < 4; k++)
        {
            cAdded[k] = meshMod.addCell(celli);
            newCellLevel(cAdded[k]) = level + 1;
        }

        const label fmA = faceMidPoint[faceA];
        const label fmB = faceMidPoint[faceB];

        // Sub-faces of A and B (orientation of the original faces kept)
        forAll(a, k)
        {
            const face newFace(subFace(loopA, midPosA, k, fmA));

            storeChildPoints(celli, loopA, midPosA, k, k);

            if (k == 0)
            {
                modifyFace(meshMod, faceA, newFace, cAdded[k], -1);
            }
            else
            {
                addFace(meshMod, faceA, newFace, cAdded[k], -1);
            }
        }

        forAll(anchorPosB, t)
        {
            const label k = childOfB[t];
            const face newFace(subFace(loopB, midPosB, t, fmB));

            storeChildPoints(celli, loopB, midPosB, t, k);

            if (t == 0)
            {
                modifyFace(meshMod, faceB, newFace, cAdded[k], -1);
            }
            else
            {
                addFace(meshMod, faceB, newFace, cAdded[k], -1);
            }
        }


        // Internal faces between child k and child k+1 (through m_k)
        forAll(a, k)
        {
            const label c1 = cAdded[k];
            const label c2 = cAdded[(k + 1) % 4];

            const label own = min(c1, c2);
            const label nei = max(c1, c2);
            const label ownAnchor = (own == c1 ? a[k] : a[(k + 1) % 4]);

            face newFace(4);
            newFace[0] = fmA;
            newFace[1] = mA[k];
            newFace[2] = mB[k];
            newFace[3] = fmB;

            // Area vector of the (planar) quad from its diagonals
            const vector area
            (
                0.5
               *(
                    (pointPos(newFace[2]) - pointPos(newFace[0]))
                  ^ (pointPos(newFace[3]) - pointPos(newFace[1]))
                )
            );
            const vector toOwner(pointPos(ownAnchor) - pointPos(mA[k]));

            if ((area & toOwner) > 0)
            {
                newFace = newFace.reverseFace();
            }

            meshMod.addFace
            (
                newFace,        // face
                own,            // owner
                nei,            // neighbour
                -1,             // master face for addition
                false,          // flux flip
                -1              // patch for face
            );
        }
    }


    // 4. All other faces: split the faces normal to the plane that span two
    //    children, re-assign owner/neighbour, insert hanging points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Refinement status of the cell on the other side of coupled faces
    labelList neiRefine;
    {
        labelList cellRef(mesh_.nCells(), 0);
        forAll(cellLabels, i)
        {
            cellRef[cellLabels[i]] = 1;
        }
        syncTools::swapBoundaryCellList(mesh_, cellRef, neiRefine);
    }

    auto childOf = [&](const label celli, const labelUList& fPts)
    {
        const labelList& cPts = cellChildPoints[celli];

        forAll(fPts, fp)
        {
            const label index = findIndex(cPts, fPts[fp]);

            if (index != -1)
            {
                return cellAddedCells[celli][cellChildIndex[celli][index]];
            }
        }

        FatalErrorInFunction
            << "Cannot find the child of refined cell " << celli
            << " for face points " << fPts
            << " (child points " << cPts << ")"
            << abort(FatalError);

        return label(-1);
    };

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(faces, facei)
    {
        const label own = faceOwner[facei];
        const bool internal = mesh_.isInternalFace(facei);
        const label nei = internal ? faceNeighbour[facei] : -1;

        const bool ownRef = refineCell[own];
        const bool neiRef = internal && refineCell[nei];

        if (isEmptyFace[facei])
        {
            // Faces of refined cells are done; others may get hanging points
            if (!ownRef)
            {
                extendedFace(facei, edgeMidPoint, isInPlaneEdge, pts, isMid);

                if (pts.size() != faces[facei].size())
                {
                    modifyFace(meshMod, facei, face(pts), own, -1);
                }
            }
            continue;
        }

        bool remoteRef = false;
        if (!internal)
        {
            const label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                remoteRef = neiRefine[facei - mesh_.nInternalFaces()] != 0;
            }
        }

        extendedFace(facei, edgeMidPoint, isInPlaneEdge, pts, isMid);

        if (!ownRef && !neiRef && !remoteRef)
        {
            if (pts.size() != faces[facei].size())
            {
                modifyFace(meshMod, facei, face(pts), own, nei);
            }
            continue;
        }

        labelList midPos;
        forAll(isMid, j)
        {
            if (isMid[j])
            {
                midPos.append(j);
            }
        }

        if (midPos.size() == 2)
        {
            // Split into two halves at the mid-points; both halves keep the
            // orientation of the original face
            const label n = pts.size();

            DynamicList<label> h1;
            for (label j = midPos[0]; ; j = (j + 1) % n)
            {
                h1.append(pts[j]);
                if (j == midPos[1]) break;
            }

            DynamicList<label> h2;
            for (label j = midPos[1]; ; j = (j + 1) % n)
            {
                h2.append(pts[j]);
                if (j == midPos[0]) break;
            }

            const face f1(h1);
            const face f2(h2);

            const label own1 = ownRef ? childOf(own, f1) : own;
            const label own2 = ownRef ? childOf(own, f2) : own;
            const label nei1 = internal ? (neiRef ? childOf(nei, f1) : nei) : -1;
            const label nei2 = internal ? (neiRef ? childOf(nei, f2) : nei) : -1;

            modifyFace(meshMod, facei, f1, own1, nei1);
            addFace(meshMod, facei, f2, own2, nei2);
        }
        else if (midPos.empty())
        {
            const face f(pts);

            const label newOwn = ownRef ? childOf(own, f) : own;
            const label newNei =
                internal ? (neiRef ? childOf(nei, f) : nei) : -1;

            modifyFace(meshMod, facei, f, newOwn, newNei);
        }
        else
        {
            FatalErrorInFunction
                << "Face " << facei << " at " << mesh_.faceCentres()[facei]
                << " normal to the plane has " << midPos.size()
                << " mid-points (expected 0 or 2): " << pts
                << abort(FatalError);
        }
    }


    // 5. Levels and history
    // ~~~~~~~~~~~~~~~~~~~~~

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    if (history_.active())
    {
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, celli)
        {
            const labelList& addedCells = cellAddedCells[celli];

            if (addedCells.size())
            {
                history_.storeSplit(celli, addedCells);
            }
        }
    }

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        refinedCells[i].transfer(cellAddedCells[cellLabels[i]]);
    }

    return refinedCells;
}


// * * * * * * * * * * * * 2D unrefinement (CDMGpuFoam)  * * * * * * * * //

Foam::labelList Foam::hexRef2D::getSplitPoints() const
{
    check2D();

    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    boolList isEmptyFace;
    boolList isInPlaneEdge;
    calcEmptyFacesAndInPlaneEdges(isEmptyFace, isInPlaneEdge);

    // In 2D the split point of a refined cell is the mid-point of its faces
    // on the empty patches: used by exactly 4 cells with the same parent
    boolList onEmptyFace(mesh_.nPoints(), false);
    forAll(isEmptyFace, facei)
    {
        if (isEmptyFace[facei])
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, fp)
            {
                onEmptyFace[f[fp]] = true;
            }
        }
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nPoints(), -1);
    labelList splitMasterLevel(mesh_.nPoints(), 0);

    for (label pointi = 0; pointi < mesh_.nPoints(); pointi++)
    {
        if (mesh_.pointCells(pointi).size() != 4 || !onEmptyFace[pointi])
        {
            splitMaster[pointi] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, celli)
    {
        const labelList& cPoints = mesh_.cellPoints(celli);

        if (visibleCells[celli] != -1 && history_.parentIndex(celli) >= 0)
        {
            const label parentIndex = history_.parentIndex(celli);

            forAll(cPoints, i)
            {
                const label pointi = cPoints[i];
                const label masterCelli = splitMaster[pointi];

                if (masterCelli == -1)
                {
                    splitMaster[pointi] = parentIndex;
                    splitMasterLevel[pointi] = cellLevel_[celli] - 1;
                }
                else if (masterCelli == -2)
                {}
                else if
                (
                    (masterCelli != parentIndex)
                 || (splitMasterLevel[pointi] != cellLevel_[celli] - 1)
                )
                {
                    splitMaster[pointi] = -2;
                }
            }
        }
        else
        {
            forAll(cPoints, i)
            {
                splitMaster[cPoints[i]] = -2;
            }
        }
    }

    // Unmark points on boundary faces other than the empty ones
    // (walls, processor boundaries: cannot merge across those)
    for
    (
        label facei = mesh_.nInternalFaces();
        facei < mesh_.nFaces();
        facei++
    )
    {
        if (!isEmptyFace[facei])
        {
            const face& f = mesh_.faces()[facei];
            forAll(f, fp)
            {
                splitMaster[f[fp]] = -2;
            }
        }
    }

    // One split point per pair of face mid-points
    DynamicList<label> splitPoints;

    forAll(splitMaster, pointi)
    {
        if (splitMaster[pointi] >= 0)
        {
            const label partner = normalPartner(pointi, isInPlaneEdge);

            if (splitMaster[partner] >= 0 && pointi < partner)
            {
                splitPoints.append(pointi);
            }
        }
    }

    return labelList(splitPoints);
}


void Foam::hexRef2D::setUnrefinement
(
    const labelList& splitPointLabels,
    polyTopoChange& meshMod
)
{
    check2D();

    if (!history_.active())
    {
        FatalErrorInFunction
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "hexRef2D::setUnrefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
    }

    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        // The 4 internal faces around each split point.  The faces on the
        // empty patches and the halves of the other faces are merged by
        // removeFaces.
        labelHashSet splitFaces(4*splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pFaces =
                mesh_.pointFaces()[splitPointLabels[i]];

            forAll(pFaces, j)
            {
                if (mesh_.isInternalFace(pFaces[j]))
                {
                    splitFaces.insert(pFaces[j]);
                }
            }
        }

        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorInFunction
                << "Initial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << " (" << facesToRemove.size() << " faces to remove, "
                << splitFaces.size() << " expected)"
                << abort(FatalError);
        }
    }

    // Check the region masters are consistent with the refinement structure
    forAll(splitPointLabels, i)
    {
        const label pointi = splitPointLabels[i];
        const labelList& pCells = mesh_.pointCells(pointi);

        if (pCells.size() != 4)
        {
            FatalErrorInFunction
                << "splitPoint " << pointi
                << " should have 4 cells using it. It has " << pCells
                << abort(FatalError);
        }

        const label masterCelli = min(pCells);

        forAll(pCells, j)
        {
            const label celli = pCells[j];
            const label region = cellRegion[celli];

            if (region == -1)
            {
                FatalErrorInFunction
                    << "Initial set of split points to unrefine does not"
                    << " seem to be consistent or not mid points"
                    << " of refined cells" << nl
                    << "cell:" << celli << " on splitPoint " << pointi
                    << " has no region to be merged into"
                    << abort(FatalError);
            }

            if (masterCelli != cellRegionMaster[region])
            {
                FatalErrorInFunction
                    << "cell:" << celli << " on splitPoint:" << pointi
                    << " in region " << region
                    << " has master:" << cellRegionMaster[region]
                    << " which is not the lowest numbered cell"
                    << " among the pointCells:" << pCells
                    << abort(FatalError);
            }
        }
    }

    // Insert all commands to combine cells
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Adapt cell levels and the history
    forAll(splitPointLabels, i)
    {
        const label pointi = splitPointLabels[i];
        const labelList& pCells = mesh_.pointCells(pointi);
        const label masterCelli = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }

        history_.combineCells(masterCelli, pCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.topoChange will take care of truncating.
}


// ************************************************************************* //
