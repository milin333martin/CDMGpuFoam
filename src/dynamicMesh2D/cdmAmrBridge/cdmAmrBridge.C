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

#include "cdmAmrBridge.H"

#include <algorithm>
#include <utility>

namespace cdmBridge
{


// * * * * * * * * * * * * * * * *  The bridge * * * * * * * * * * * * * * * //

// sidePatch gives the patch index for each of the four in-plane sides, or
// -1 where that side cannot be a boundary; frontBack is the patch for the
// two faces in the plane; cyclicMin and cyclicMax are the cyclic pair.
// The engine treats a cyclic direction as internal, which a host mesh
// cannot: an internal face needs vertices shared by both cells, and the
// two are a domain width apart. Each seam face becomes two host faces.
// This costs no refinement freedom. Where a coarse cell meets two finer
// cells across the seam it gets two half-side faces on its cyclic patch,
// which is the hanging-node treatment applied to a patch.
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
{
    const int nCells = m.nCells();
    const int nEngInt = m.nFaces();
    const int nBndSide = m.nBoundaryFaces();
    const int nFB = 2*nCells;
    const std::uint64_t NX = static_cast<std::uint64_t>(m.nx)*cFine;

    const std::uint64_t strideX =
        static_cast<std::uint64_t>(m.nx)*cFine + 1;

    const double dxf = m.Lx/(static_cast<double>(m.nx)*cFine);
    const double dyf = m.Ly/(static_cast<double>(m.ny)*cFine);

    // ---- the cells' integer boxes, used by both halves ----
    std::vector<Box> box(nCells);
    for (int c = 0; c < nCells; ++c)
    {
        box[c] = cellBoxI(m.keys[c], m.nx);
    }

    // This must be geometric. A refined cell in base column nx-1 may sit in
    // that column's left half, where its +x side is an ordinary interior face;
    // testing the base index turned those into cyclic pairs.
    const std::uint64_t NXs = NX;
    auto atSeam = [NXs](const Box& b, int d)
    {
        return (d == 0 && b.ix0 == 0)
            || (d == 1 && std::uint64_t(b.ix0) + b.s == NXs);
    };

    // Corners are generated from the cells, so none need searching for: each
    // carries a slot id through the sort and the point index it lands on is
    // scattered back into slotPoint.
    // A searched version did four binary searches per face through a
    // three-million-entry array, about seven million cache misses a rebuild.
    double t0 = wall();

    // Seam faces first: each needs four points at the far end of the domain,
    // corners of no cell. Ordinary hanging nodes need none of this, being
    // built from the finer cell whose corners they already are.
    std::vector<int> seamOf(static_cast<std::size_t>(nEngInt), -1);
    int nSeam = 0;

    if (cyclicMin >= 0)
    {
        for (int f = 0; f < nEngInt; ++f)
        {
            if (atSeam(box[m.faceCell[f]], m.faceSide[f]))
            {
                seamOf[static_cast<std::size_t>(f)] = nSeam++;
            }
        }
    }

    // slot layout: 8 per cell as 4*z + corner, corner 0..3 running
    // (x0,y0), (x1,y0), (x1,y1), (x0,y1); then 4 per seam face as
    // 2*z + e, e = 0 at iy0 and 1 at iy0+s, all at the far end
    const std::size_t cellSlots = 8*static_cast<std::size_t>(nCells);
    const std::size_t nSlots = cellSlots + 4*static_cast<std::size_t>(nSeam);

    struct Corner { std::uint64_t code; int slot; };
    std::vector<Corner> corner(nSlots);

    for (int c = 0; c < nCells; ++c)
    {
        const Box& b = box[c];
        const std::uint64_t x0 = b.ix0, x1 = b.ix0 + b.s;
        const std::uint64_t y0 = b.iy0, y1 = b.iy0 + b.s;
        const std::size_t s0 = 8*static_cast<std::size_t>(c);

        for (std::uint64_t z = 0; z < 2; ++z)
        {
            const std::size_t at = s0 + 4*static_cast<std::size_t>(z);
            corner[at + 0] =
                Corner{pointCode(x0, y0, z, strideX), int(at + 0)};
            corner[at + 1] =
                Corner{pointCode(x1, y0, z, strideX), int(at + 1)};
            corner[at + 2] =
                Corner{pointCode(x1, y1, z, strideX), int(at + 2)};
            corner[at + 3] =
                Corner{pointCode(x0, y1, z, strideX), int(at + 3)};
        }
    }

    for (int f = 0; f < nEngInt; ++f)
    {
        const int sm = seamOf[static_cast<std::size_t>(f)];
        if (sm < 0) { continue; }

        const Box& b = box[m.faceCell[f]];
        const std::uint64_t ixOther = (m.faceSide[f] == 0) ? NX : 0;
        const std::size_t s0 = cellSlots + 4*static_cast<std::size_t>(sm);

        for (std::uint64_t z = 0; z < 2; ++z)
        {
            const std::size_t at = s0 + 2*static_cast<std::size_t>(z);
            corner[at + 0] =
                Corner{pointCode(ixOther, b.iy0, z, strideX), int(at + 0)};
            corner[at + 1] = Corner
            {
                pointCode(ixOther, b.iy0 + b.s, z, strideX), int(at + 1)
            };
        }
    }

    psort
    (
        corner.begin(), corner.end(),
        [](const Corner& a, const Corner& b) { return a.code < b.code; }
    );

    std::vector<label> slotPoint(nSlots, -1);
    std::vector<std::uint64_t> code;
    code.reserve(nSlots/2 + 1);

    for (std::size_t i = 0; i < corner.size(); ++i)
    {
        if (i == 0 || corner[i].code != corner[i - 1].code)
        {
            code.push_back(corner[i].code);
        }
        slotPoint[static_cast<std::size_t>(corner[i].slot)] =
            static_cast<label>(code.size()) - 1;
    }

    const label nPoints = static_cast<label>(code.size());
    points.setSize(nPoints);

    for (label q = 0; q < nPoints; ++q)
    {
        const std::uint64_t cd = code[static_cast<std::size_t>(q)];
        const std::uint64_t iz = cd & 1;
        const std::uint64_t r = cd >> 1;

        points[q] = point
        (
            static_cast<double>(r % strideX)*dxf,
            static_cast<double>(r / strideX)*dyf,
            static_cast<double>(iz)*m.thickness
        );
    }

    pointCodeOut = code;

    tPoints = wall() - t0;

    // Which of a cell's eight corners each side's quad uses, and in what
    // order. The winding is corrected afterwards against the engine's Sf.
    static const int sideCorner[6][4] =
    {
        {0, 3, 3, 0},   // -x
        {1, 2, 2, 1},   // +x
        {0, 0, 1, 1},   // -y
        {3, 3, 2, 2},   // +y
        {0, 1, 2, 3},   // back
        {0, 1, 2, 3}    // front
    };
    static const int sideZ[6][4] =
    {
        {0, 0, 1, 1},
        {0, 0, 1, 1},
        {0, 1, 1, 0},
        {0, 1, 1, 0},
        {0, 0, 0, 0},
        {1, 1, 1, 1}
    };

    auto sideQuad = [&slotPoint](int cell, int d, label v[4])
    {
        const std::size_t s0 = 8*static_cast<std::size_t>(cell);
        for (int i = 0; i < 4; ++i)
        {
            v[i] = slotPoint
            [
                s0 + 4*static_cast<std::size_t>(sideZ[d][i])
                   + static_cast<std::size_t>(sideCorner[d][i])
            ];
        }
    };

    // The seam partner's half-face, from the four slots reserved for it
    auto seamQuad = [&slotPoint, cellSlots](int sm, label v[4])
    {
        const std::size_t s0 = cellSlots + 4*static_cast<std::size_t>(sm);
        v[0] = slotPoint[s0 + 0];   // z0, y0
        v[1] = slotPoint[s0 + 1];   // z0, y1
        v[2] = slotPoint[s0 + 3];   // z1, y1
        v[3] = slotPoint[s0 + 2];   // z1, y0
    };


    // ---- faces ----
    t0 = wall();

    std::vector<FaceRec> rec;
    rec.reserve(static_cast<std::size_t>(nEngInt + nBndSide + nFB) + 4096);

    // The geometric identity of side d of a cell, and of a z face.  sign
    // is filled by the caller once the face's orientation is settled.
    auto sideKey = [](const Box& b, int d, int cell)
    {
        const std::uint32_t x0 = b.ix0, x1 = b.ix0 + b.s;
        const std::uint32_t y0 = b.iy0, y1 = b.iy0 + b.s;

        FaceKey k;
        k.cell = cell;
        k.side = d;
        k.sign = 1;

        if (d < 2)
        {
            k.orient = 0; k.pos = (d == 1) ? x1 : x0;
            k.a0 = y0; k.a1 = y1; k.b0 = 0;
        }
        else if (d < 4)
        {
            k.orient = 1; k.pos = (d == 3) ? y1 : y0;
            k.a0 = x0; k.a1 = x1; k.b0 = 0;
        }
        else
        {
            k.orient = (d == 5) ? 3 : 2; k.pos = 0;
            k.a0 = x0; k.a1 = x1; k.b0 = y0;
        }
        return k;
    };

    // The normal of the quad as wound, for comparison against the sense the
    // face is required to have
    auto quadNormal = [&points](const label v[4]) -> vector
    {
        return 0.5*
        (
            (points[v[2]] - points[v[0]]) ^ (points[v[3]] - points[v[1]])
        );
    };

    auto orient = [&quadNormal](FaceRec& r, const vector& want)
    {
        if ((quadNormal(r.v) & want) < 0) { std::swap(r.v[1], r.v[3]); }
    };

    // internal, and the seam pairs that come out of them
    for (int f = 0; f < nEngInt; ++f)
    {
        const int c = m.faceCell[f];
        const int d = m.faceSide[f];
        const Box& b = box[c];

        const bool seam = (cyclicMin >= 0) && atSeam(b, d);

        if (!seam)
        {
            FaceRec r;
            r.patch = -1;
            r.own = m.owner[f];
            r.nei = m.neighbour[f];
            r.skey = 0;
            r.eface = f;
            r.id = sideKey(b, d, m.owner[f]);
            // Sf is owner-outward; the face's own +axis sense is d odd
            r.id.sign = (m.Sf[2*f] + m.Sf[2*f + 1] > 0) ? 1 : -1;
            sideQuad(c, d, r.v);

            // The engine's Sf is already owner-outward and sign corrected. Trust it
            // rather than re-deriving the winding: a disagreement means the face is
            // reversed, which is one comparison instead of a case analysis.
            orient(r, vector(m.Sf[2*f], m.Sf[2*f + 1], 0));
            rec.push_back(r);
            continue;
        }

        // A seam face is made from the finer of the two cells, so its y interval
        // is that cell's side. A coarser partner receives two such faces, one
        // from each finer neighbour, and they tile its side.
        const label partner = (m.owner[f] == c) ? m.neighbour[f] : m.owner[f];
        const std::uint64_t y0 = b.iy0;

        // c sits at the end its own side faces; the partner at the other
        const bool cAtMin = (d == 0);

        FaceRec rc;
        rc.patch = cAtMin ? cyclicMin : cyclicMax;
        rc.own = c;
        rc.nei = -1;
        rc.skey = static_cast<label>(y0);
        rc.eface = -1;
        rc.id = sideKey(b, d, c);
        rc.id.pos = static_cast<std::uint32_t>(cAtMin ? 0 : NX);
        rc.id.sign = cAtMin ? -1 : 1;
        // the creator's own side, which is at the seam by definition
        sideQuad(c, d, rc.v);
        orient(rc, vector(cAtMin ? -1.0 : 1.0, 0, 0));
        rec.push_back(rc);

        FaceRec rp;
        rp.patch = cAtMin ? cyclicMax : cyclicMin;
        rp.own = partner;
        rp.nei = -1;
        rp.skey = static_cast<label>(y0);
        rp.eface = -1;
        rp.id = sideKey(b, d, partner);
        rp.id.pos = static_cast<std::uint32_t>(cAtMin ? NX : 0);
        rp.id.sign = cAtMin ? 1 : -1;
        // the partner's half-face, from the slots reserved for this seam
        seamQuad(seamOf[static_cast<std::size_t>(f)], rp.v);
        orient(rp, vector(cAtMin ? 1.0 : -1.0, 0, 0));
        rec.push_back(rp);
    }

    // in-plane boundary
    for (int bi = 0; bi < nBndSide; ++bi)
    {
        const int c = m.bCell[bi];
        const int d = m.bSide[bi];

        FaceRec r;
        r.patch = sidePatch[d];
        r.own = c;
        r.nei = -1;
        r.skey = c;
        r.eface = -1;
        r.id = sideKey(box[c], d, c);
        r.id.sign = (d == 1 || d == 3) ? 1 : -1;
        sideQuad(c, d, r.v);

        // A boundary face must point out of its cell
        orient
        (
            r,
            vector
            (
                (d == 0) ? -1.0 : (d == 1) ? 1.0 : 0.0,
                (d == 2) ? -1.0 : (d == 3) ? 1.0 : 0.0,
                0.0
            )
        );
        rec.push_back(r);
    }

    // front and back
    for (int c = 0; c < nCells; ++c)
    {
        for (int k = 0; k < 2; ++k)
        {
            FaceRec r;
            r.patch = frontBack;
            r.own = c;
            r.nei = -1;
            r.skey = c;
            r.eface = -1;
            r.id = sideKey(box[c], (k == 0) ? 4 : 5, c);
            r.id.sign = (k == 0) ? -1 : 1;
            sideQuad(c, (k == 0) ? 4 : 5, r.v);
            orient(r, vector(0, 0, (k == 0) ? -1.0 : 1.0));
            rec.push_back(r);
        }
    }

    const label nFaces = static_cast<label>(rec.size());

    tFaces = wall() - t0;

    // ---- order them the way OpenFOAM requires ----
    t0 = wall();

    // Sort a permutation, not the records: a FaceRec is some sixty bytes, and
    // moving a million through a comparison sort costs far more than moving
    // sixteen-byte keys and reading the records once in the resulting order.
    // The key puts internal faces first, then groups by patch, as OpenFOAM
    // requires. Within a cyclic patch the order is position along the seam,
    // which is what makes the two halves pair up.
    std::vector<FaceSort> ord(static_cast<std::size_t>(nFaces));

    // The key is (group, a, b) in one 64-bit integer. Six bits of group leaves
    // 29 for each cell field, or 536 million cells; a 22/21/21 split capped
    // the mesh at 2,097,151 and a production run reached it and aborted.
    // The a field is a cell index for an internal face and a seam coordinate
    // for a patch face; the latter is bounded by ny*2^maxLevel, so the cell
    // bound covers both.
    if (nCells >= (1 << 29))
    {
        FatalErrorInFunction
            << nCells << " cells does not fit the 29-bit sort key"
            << abort(FatalError);
    }

    if (nPatches + 1 >= (1 << 6))
    {
        FatalErrorInFunction
            << nPatches << " patches does not fit the 6-bit group field"
            << abort(FatalError);
    }

    for (label f = 0; f < nFaces; ++f)
    {
        const FaceRec& r = rec[static_cast<std::size_t>(f)];

        const std::uint64_t grp =
            (r.patch < 0) ? 0u : static_cast<std::uint64_t>(r.patch) + 1u;
        const std::uint64_t a =
            (r.patch < 0) ? static_cast<std::uint64_t>(r.own)
                          : static_cast<std::uint64_t>(r.skey);
        const std::uint64_t b2 =
            (r.patch < 0) ? static_cast<std::uint64_t>(r.nei + 1)
                          : static_cast<std::uint64_t>(r.own);

        ord[static_cast<std::size_t>(f)] =
            FaceSort{(grp << 58) | (a << 29) | b2, f};
    }

    psort(ord.begin(), ord.end());

    // Internal faces sort first, so counting them is counting the leading run.
    // This is not the engine's internal face count: each seam face left that
    // list and became two patch faces.
    label nInt = 0;
    while
    (
        nInt < nFaces
     && rec[static_cast<std::size_t>(ord[std::size_t(nInt)].idx)].patch < 0
    )
    {
        ++nInt;
    }

    faces.setSize(nFaces);
    owner.setSize(nFaces);
    neighbour.setSize(nInt);
    faceEngine.setSize(nInt);
    faceKey.resize(static_cast<std::size_t>(nFaces));

    // Assigned, not setSize(n, 0): setSize only initialises what it adds, so
    // previous counts would survive and be incremented on top of, surfacing
    // later as a patch slicing past the end of faceAreas.
    patchSizes.setSize(nPatches);
    patchStarts.setSize(nPatches);
    patchSizes = 0;
    patchStarts = 0;

    for (label f = 0; f < nFaces; ++f)
    {
        const FaceRec& r =
            rec[static_cast<std::size_t>(ord[static_cast<std::size_t>(f)].idx)];

        face& fc = faces[f];
        fc.setSize(4);
        fc[0] = r.v[0];
        fc[1] = r.v[1];
        fc[2] = r.v[2];
        fc[3] = r.v[3];

        owner[f] = r.own;
        faceKey[static_cast<std::size_t>(f)] = r.id;
        if (r.patch < 0) { neighbour[f] = r.nei; faceEngine[f] = r.eface; }
        else             { ++patchSizes[r.patch]; }
    }

    label start = nInt;
    forAll(patchStarts, p)
    {
        patchStarts[p] = start;
        start += patchSizes[p];
    }

    // The patches must tile [nInternalFaces, nFaces) exactly. Checked here
    // because nothing downstream does: the symptom is a segmentation fault
    // inside OpenFOAM, a long way from the cause.
    if (start != nFaces)
    {
        FatalErrorInFunction
            << "patches cover " << start - nInt << " faces but there are "
            << nFaces - nInt << " boundary faces"
            << abort(FatalError);
    }

    // A cyclic pair with unequal halves cannot be matched face for face,
    // and OpenFOAM will not say so until much later.
    if (cyclicMin >= 0 && patchSizes[cyclicMin] != patchSizes[cyclicMax])
    {
        FatalErrorInFunction
            << "cyclic halves differ: " << patchSizes[cyclicMin]
            << " against " << patchSizes[cyclicMax]
            << abort(FatalError);
    }

    tSort = wall() - t0;
}



// * * * * * * * * * * * *  The topology change map  * * * * * * * * * * * * //

// Builds what polyTopoChangeMap needs so that fvMesh::mapFields can map
// every registered field with OpenFOAM's own semantics, rather than a
// subset reimplemented here.
//   cells  by Morton key            (mapRange, exact)
//   points by integer position code (a merge of two sorted arrays, exact)
//   faces  by geometry              (a range query, exact)
// cellsFromCells carries the coarsened cells, and cellMapper turns it into
// volume weights using the old cell volumes, which is how coarsening
// stays conservative.
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
{
    const label nNewCells = m.nCells();
    const label nOldCells = static_cast<label>(oldKeys.size());

    // ---- cells ----
    {
        std::vector<int> first, count;
        cdmAmr::Mesh::mapRange(oldKeys, m.keys, first, count);

        cellMap.setSize(nNewCells);
        DynamicList<objectMap> cfc(nNewCells/8 + 16);

        for (label c = 0; c < nNewCells; ++c)
        {
            const std::size_t uc = static_cast<std::size_t>(c);
            cellMap[c] = first[uc];

            if (count[uc] > 1)
            {
                labelList mo(count[uc]);
                forAll(mo, i) { mo[i] = first[uc] + i; }
                cfc.append(objectMap(c, mo));
            }
        }
        cellsFromCells.transfer(cfc);
    }

    // ---- points: a merge of two sorted, unique code arrays ----
    {
        pointMap.setSize(static_cast<label>(newCode.size()));

        std::size_t i = 0;
        for (std::size_t p = 0; p < newCode.size(); ++p)
        {
            while (i < oldCode.size() && oldCode[i] < newCode[p]) { ++i; }

            // A point refinement created has no old counterpart. polyTopoChangeMap
            // represents that as a negative entry, and pointsFromPoints stays empty:
            // no point is ever born from several old points.
            pointMap[static_cast<label>(p)] =
                (i < oldCode.size() && oldCode[i] == newCode[p])
              ? static_cast<label>(i) : -1;
        }
    }

    // Both lists are sorted by the packed key, so this is a merge rather than
    // a binary search per face. z faces are left out of both, since they
    // follow the cell map instead.
    {
        const label nNewFaces = static_cast<label>(newFace.size());
        const label nOldFaces = static_cast<label>(oldFace.size());

        faceMap.setSize(nNewFaces);
        faceMap = -1;
        label nMerged = 0, nCreated = 0, nUnmapped = 0;

        // Every old face indexed by the side of the cell it sits on. Used for z
        // faces, and to give a master to faces a refinement created from nothing.
        std::vector<int> sideFace
        (
            6*static_cast<std::size_t>(nOldCells), -1
        );
        for (label f = 0; f < nOldFaces; ++f)
        {
            const FaceKey& k = oldFace[static_cast<std::size_t>(f)];
            sideFace
            [
                6*static_cast<std::size_t>(k.cell)
              + static_cast<std::size_t>(k.side)
            ] = f;
        }

        auto ancestorFace = [&](label newCell, int orient) -> label
        {
            const label oc = (newCell >= 0) ? cellMap[newCell] : -1;
            if (oc < 0) { return -1; }
            const std::size_t base = 6*static_cast<std::size_t>(oc);
            const int want = (orient == 0) ? 0 : (orient == 1) ? 2 : 4;
            for (int k = 0; k < 6; ++k)
            {
                const int trySide = (k == 0) ? want : k - 1;
                const label f = sideFace[base + std::size_t(trySide)];
                if (f >= 0) { return f; }
            }
            return -1;
        };

        // the in-plane faces of the new mesh, sorted
        newSorted.clear();
        newSorted.reserve(static_cast<std::size_t>(nNewFaces));
        for (label f = 0; f < nNewFaces; ++f)
        {
            const FaceKey& k = newFace[static_cast<std::size_t>(f)];
            if (k.orient >= 2) { continue; }
            if (k.pos >= (1u << 20) || k.a0 >= (1u << 20)
             || k.b0 >= (1u << 20))
            {
                FatalErrorInFunction
                    << "face coordinate does not fit 20 bits:"
                    << " pos " << k.pos << " a0 " << k.a0
                    << " b0 " << k.b0 << "; base mesh too wide"
                    << abort(FatalError);
            }
            newSorted.push_back(FaceSort{packFaceKey(k), f});
        }
        psort(newSorted.begin(), newSorted.end());

        // merge against the old mesh's sorted in-plane faces
        std::size_t i = 0, oi = 0;
        while (i < newSorted.size())
        {
            const std::uint64_t pre = newSorted[i].key >> 20;

            while (oi < oldSorted.size() && (oldSorted[oi].key >> 20) < pre)
            {
                ++oi;
            }
            std::size_t oj = oi;
            while (oj < oldSorted.size() && (oldSorted[oj].key >> 20) == pre)
            {
                ++oj;
            }
            std::size_t j = i;
            while (j < newSorted.size() && (newSorted[j].key >> 20) == pre)
            {
                ++j;
            }

            // New faces sharing a plane tile it without overlapping and both groups
            // are sorted along it, so one forward sweep of the old group serves every
            // new face in this plane.
            std::size_t op = oi;

            for (std::size_t t = i; t < j; ++t)
            {
                const label f = newSorted[t].idx;
                const FaceKey& k = newFace[static_cast<std::size_t>(f)];

                while (op < oj
                    && oldFace[static_cast<std::size_t>(oldSorted[op].idx)].a1
                       <= k.a0)
                {
                    ++op;
                }

                label direct = -1, firstOverlap = -1, nOver = 0;

                for (std::size_t q = op; q < oj; ++q)
                {
                    const label of = oldSorted[q].idx;
                    const FaceKey& x = oldFace[static_cast<std::size_t>(of)];
                    if (x.a0 >= k.a1) { break; }

                    if ((x.a0 == k.a0 && x.a1 == k.a1)
                     || (x.a0 <= k.a0 && x.a1 >= k.a1))
                    {
                        direct = of;
                        break;
                    }
                    if (firstOverlap < 0) { firstOverlap = of; }
                    ++nOver;
                }

                label master = -1;
                if (direct >= 0)
                {
                    master = direct;
                }
                else if (nOver >= 1)
                {
                    // Coarsening merged several old faces into this one. facesFromFaces would
                    // express that, but faceMapper then takes its interpolative branch, which
                    // fatal-errors on any face with no addressing at all, and refinement
                    // creates exactly those. So one master is taken, keeping every face in
                    // the direct branch.
                    // cellsFromCells is untouched, so cell fields stay volume weighted and
                    // conservative; only surface fields degrade to nearest-value on a merged
                    // face, and this solver carries no surface field as state.
                    master = firstOverlap;
                    if (nOver > 1) { ++nMerged; }
                }
                bool standIn = false;
                if (direct < 0 && nOver == 0)
                {
                    // refinement created this face inside the cell it split
                    master = ancestorFace(k.cell, k.orient);
                    standIn = true;
                    ++nCreated;
                    if (master < 0) { ++nUnmapped; }
                }

                faceMap[f] = master;

                // The sign is only meaningful against a master that is really the same
                // face. A stand-in is the ancestor's face on another side of the cell, so
                // its sign would corrupt the value rather than correct it.
                if (!standIn
                 && master >= 0
                 && oldFace[static_cast<std::size_t>(master)].sign != k.sign)
                {
                    flipFaceFlux.insert(f);
                }
            }
            i = j;
        }

        // z faces follow the cell map
        for (label f = 0; f < nNewFaces; ++f)
        {
            const FaceKey& k = newFace[static_cast<std::size_t>(f)];
            if (k.orient < 2) { continue; }

            const label oc = cellMap[k.cell];
            const label of =
                (oc < 0) ? -1
              : sideFace
                [
                    6*static_cast<std::size_t>(oc)
                  + ((k.orient == 2) ? 4u : 5u)
                ];
            faceMap[f] = of;
            if (of >= 0
             && oldFace[static_cast<std::size_t>(of)].sign != k.sign)
            {
                flipFaceFlux.insert(f);
            }
            if (of < 0) { ++nUnmapped; }
        }

        // left empty on purpose; see the note on merged faces above
        facesFromFaces.clear();

        static bool reported = false;
        if (!reported)
        {
            reported = true;
            Info<< "    first map: " << nMerged << " faces merged, "
                << nCreated << " created by refinement, "
                << nUnmapped << " left unmapped" << nl << endl;
        }
    }


    // ---- reverse maps ----
    revCell.setSize(nOldCells, -1);
    revCell = -1;
    forAll(cellMap, c)
    {
        if (cellMap[c] >= 0) { revCell[cellMap[c]] = c; }
    }

    revFace.setSize(static_cast<label>(oldFace.size()), -1);
    revFace = -1;
    forAll(faceMap, f)
    {
        if (faceMap[f] >= 0) { revFace[faceMap[f]] = f; }
    }

    revPoint.setSize(static_cast<label>(oldCode.size()), -1);
    revPoint = -1;
    forAll(pointMap, p)
    {
        if (pointMap[p] >= 0) { revPoint[pointMap[p]] = p; }
    }

}

// * * * * * * * * * * * *  Fitting an existing case  * * * * * * * * * * * //

// Works out which patch each side of the base mesh belongs to from where
// its boundary faces are, rather than from patch names, which a real case
// chooses freely.
// One patch per side: a side whose faces span more than one patch is
// refused rather than quietly mis-assigned. Handling that needs a
// per-base-face table, which does not exist yet.
void deducePatchLayout
(
    const fvMesh& mesh,
    labelList& sidePatch,       // [4], patch per in-plane side, or -1
    label& frontBack,
    label& cyclicMin,
    label& cyclicMax
)
{
    const boundBox& bb = mesh.bounds();
    const scalar tol = 1.0e-6*bb.mag();

    sidePatch.setSize(4);
    sidePatch = -1;
    frontBack = -1;
    cyclicMin = -1;
    cyclicMax = -1;

    // The polyPatches, not the fvPatches. An emptyFvPatch reports size() == 0
    // because it carries no finite-volume faces, so going through the
    // fvBoundaryMesh skips the one patch this is looking for.
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        // zero-size patches say nothing about where they are; OpenFOAM-13
        // carries an "internalFaces" patch of type internal with no faces
        if (pp.size() == 0) { continue; }

        // By value. faceCentres() and faceAreas() return a vectorField::subField
        // view by value, and binding a const vectorField& to one gives a reference
        // to a temporary that dies at the end of the statement.
        const vectorField::subField Cf = pp.faceCentres();
        const vectorField::subField Sf = pp.faceAreas();

        vector nSum = Zero;
        forAll(Sf, i) { nSum += Sf[i]; }
        const vector n = nSum/(mag(nSum) + vSmall);

        int axis = 0;
        if (mag(n.y()) > mag(n.x())) { axis = 1; }
        if (mag(n.z()) > mag(n[axis])) { axis = 2; }

        const bool atMin = (n[axis] < 0);

        if (axis == 2)
        {
            // the empty patch, both planes of a one-cell-thick mesh
            if (frontBack >= 0 && frontBack != patchi)
            {
                FatalErrorInFunction
                    << "two different patches on the z planes, "
                    << mesh.boundaryMesh()[frontBack].name() << " and "
                    << pp.name()
                    << "; the bridge puts both z faces on one patch"
                    << abort(FatalError);
            }
            frontBack = patchi;
            continue;
        }

        // the patch should sit on one plane
        scalar lo = great, hi = -great;
        forAll(Cf, i)
        {
            lo = min(lo, Cf[i][axis]);
            hi = max(hi, Cf[i][axis]);
        }

        if (hi - lo > tol)
        {
            FatalErrorInFunction
                << "patch " << pp.name() << " is not planar in direction "
                << axis << ": spans " << lo << " to " << hi
                << abort(FatalError);
        }

        if (isA<cyclicPolyPatch>(pp))
        {
            if (axis != 0)
            {
                FatalErrorInFunction
                    << "cyclic patch " << pp.name() << " is on the y"
                    << " sides; the engine is only cyclic in x"
                    << abort(FatalError);
            }
            if (atMin) { cyclicMin = patchi; } else { cyclicMax = patchi; }
            continue;
        }

        const label side = 2*axis + (atMin ? 0 : 1);

        if (sidePatch[side] >= 0 && sidePatch[side] != patchi)
        {
            FatalErrorInFunction
                << "side " << side << " spans patches "
                << mesh.boundaryMesh()[sidePatch[side]].name() << " and "
                << pp.name()
                << "; a per-base-face patch table would be needed"
                << abort(FatalError);
        }
        sidePatch[side] = patchi;
    }

    if (frontBack < 0)
    {
        FatalErrorInFunction
            << "no patch found on the z planes; is this a one-cell-thick"
            << " 2D mesh?" << abort(FatalError);
    }

    if ((cyclicMin >= 0) != (cyclicMax >= 0))
    {
        FatalErrorInFunction
            << "found one half of a cyclic pair but not the other"
            << abort(FatalError);
    }

    Info<< "Patch layout deduced from the case:" << nl;
    const word sideName[4] = {"xMin", "xMax", "yMin", "yMax"};
    for (int d = 0; d < 4; ++d)
    {
        Info<< "    " << sideName[d] << "          "
            << (
                   sidePatch[d] >= 0
                 ? mesh.boundaryMesh()[sidePatch[d]].name()
                 : word("(cyclic or none)")
               ) << nl;
    }
    Info<< "    frontAndBack  " << mesh.boundaryMesh()[frontBack].name()
        << nl;
    if (cyclicMin >= 0)
    {
        Info<< "    cyclic pair   "
            << mesh.boundaryMesh()[cyclicMin].name() << " / "
            << mesh.boundaryMesh()[cyclicMax].name() << nl;
    }
    Info<< endl;
}

} // End namespace cdmBridge

// ************************************************************************* //
