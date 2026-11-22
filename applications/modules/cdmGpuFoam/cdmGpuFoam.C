/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: a chemical-diffusive-model solver
   \\    /   O peration     | for gaseous detonations in OpenFOAM-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

\*---------------------------------------------------------------------------*/

#include "cdmGpuFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"
#include "globalMeshData.H"
#include "PstreamReduceOps.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "UIPstream.H"
#include "UOPstream.H"
#include "OSspecific.H"
#include <cstdlib>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(cdmGpuFoam, 0);
    addToRunTimeSelectionTable(solver, cdmGpuFoam, fvMesh);
}
}


// The back end's C interface takes int, so a 64-bit-label build of OpenFOAM
// would silently pass the wrong thing.  Fail at compile time instead.
static_assert
(
    sizeof(Foam::label) == sizeof(int),
    "cdmGpuFoam needs a 32-bit-label OpenFOAM (WM_LABEL_SIZE=32)"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

static unsigned long long cdmTopoHash(const Foam::fvMesh& m)
{
    // Cell and face counts can return to their old values across a refine
    // and coarsen pass while the mesh is renumbered, so counts cannot detect
    // a topology change and the device would keep stale addressing.
    unsigned long long h = 1469598103934665603ull;
    const Foam::labelList& own = m.faceOwner();
    const Foam::labelList& nei = m.faceNeighbour();
    h = (h ^ (unsigned long long)m.nCells()) * 1099511628211ull;
    forAll(own, i) h = (h ^ (unsigned long long)own[i]) * 1099511628211ull;
    forAll(nei, i) h = (h ^ (unsigned long long)nei[i]) * 1099511628211ull;
    return h;
}

void Foam::solvers::cdmGpuFoam::uploadMesh()
{
    const double tEnterUpload = cdmGpuWallTime();
    nUploadMesh_++;

    // DimensionedField publicly derives from Field<scalar>, so bind the base
    // reference: mesh.V().field() resolves to OldTimeField::field(), private.
    const scalarField& Vcells = mesh.V();

    const label nCells = mesh.nCells();
    const label nIF = mesh.nInternalFaces();
    const label nFaces = mesh.nFaces();

    // ---- flattened boundary list, empty patches excluded -----------------
    // emptyFvPatch has zero size, so a 2D mesh's front and back drop out of
    // mesh.boundary() by construction: on an 8.6M cell case they would be
    // 17M faces that are not part of the discretisation.
    labelList faceToB(nFaces - nIF, -1);

    // A cyclic face is an internal face, not a boundary. Handing the device
    // a cell value on one side and an interpolated face value on the other
    // gives a half-cell Riemann problem, first order, and the two sides
    // disagree, so mass and energy leak at the seam.
    // Cyclic faces are therefore appended to the device's internal-face list
    // once each: one face, one flux, conservative by construction, with full
    // MUSCL reconstruction and no new kernel code.
    labelList faceToC(nFaces - nIF, -1);
    // The cast is load-bearing: with a bare 0 both arguments are int and
    // List's iterator-pair constructor wins over List(label, const T&).
    List<signed char> faceToCSign(nFaces - nIF, static_cast<signed char>(0));

    labelList cycPatch;
    labelList cycFace;
    {
        DynamicList<label> cp, cf;

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (fvp.size() == 0 || !fvp.coupled()) { continue; }

            // Processor patches are collected below: promoted to internal faces too,
            // but with a ghost cell as the neighbour rather than a local one.
            if (isA<processorFvPatch>(fvp)) { continue; }

            if (!isA<cyclicFvPatch>(fvp))
            {
                FatalErrorInFunction
                    << "patch " << fvp.name() << " is coupled but is a "
                    << mesh.boundaryMesh()[patchi].type()
                    << ", which cdmGpuFoam does not handle." << nl
                    << "Only cyclic and processor coupling are supported."
                    << exit(FatalError);
            }

            // owner() is cyclicFvPatch's own answer to "which side of the
            // pair am I", so the face is taken exactly once
            if (!refCast<const cyclicFvPatch>(fvp).owner()) { continue; }

            forAll(fvp, facei)
            {
                cp.append(patchi);
                cf.append(facei);
            }
        }

        cycPatch.transfer(cp);
        cycFace.transfer(cf);
    }

    const label nCyc = cycPatch.size();

    // ---- processor patches, in patch order ------------------------------
    // Every processor face becomes a device internal face with a ghost
    // neighbour. OpenFOAM matches face i on both ranks, so the exchange is
    // a contiguous swap with no index traffic.
    {
        DynamicList<label> pp, pf, hc, ps, pn;

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (fvp.size() == 0 || !isA<processorFvPatch>(fvp)) { continue; }

            ps.append(pp.size());
            pn.append
            (
                refCast<const processorFvPatch>(fvp).neighbProcNo()
            );

            const labelUList& fc = fvp.faceCells();

            forAll(fvp, facei)
            {
                pp.append(patchi);
                pf.append(facei);
                hc.append(fc[facei]);
            }
        }

        ps.append(pp.size());

        procPatch_.transfer(pp);
        procFace_.transfer(pf);
        haloCell_.transfer(hc);
        procStart_.transfer(ps);
        procNbr_.transfer(pn);
    }

    nHalo_ = procPatch_.size();

    const label nIFdev = nIF + nCyc + nHalo_;

    label nBF = 0;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];
        if (fvp.coupled()) { continue; }
        nBF += fvp.size();
    }

    bPatch_.setSize(nBF);
    bFace_.setSize(nBF);

    label bi = 0;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];
        if (fvp.size() == 0 || fvp.coupled()) { continue; }

        const label start = mesh.boundaryMesh()[patchi].start();

        forAll(fvp, facei)
        {
            faceToB[start + facei - nIF] = bi;
            bPatch_[bi] = patchi;
            bFace_[bi] = facei;
            ++bi;
        }
    }

    // ---- geometry --------------------------------------------------------
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& dc = mesh.nonOrthDeltaCoeffs();
    const volVectorField& C = mesh.C();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    // The device's internal-face arrays: the mesh's own faces first, then
    // one entry per cyclic face.
    labelList devOwn(nIFdev), devNei(nIFdev);
    scalarField SfFlat(3*nIFdev), deltaFlat(3*nIFdev);
    scalarField magSfDev(nIFdev), wDev(nIFdev), dcDev(nIFdev);

    forAll(own, facei)
    {
        devOwn[facei] = own[facei];
        devNei[facei] = nei[facei];

        const vector& s = Sf[facei];
        const vector d(C[nei[facei]] - C[own[facei]]);
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            SfFlat[3*facei + cmpt] = s[cmpt];
            deltaFlat[3*facei + cmpt] = d[cmpt];
        }

        magSfDev[facei] = magSf[facei];
        wDev[facei] = w[facei];
        dcDev[facei] = dc[facei];
    }

    // ---- cyclic faces, appended as internal faces ------------------------
    // Patch-major: fvPatch::delta() builds the whole patch's vectors at once.
    {
        label k = 0;

        while (k < nCyc)
        {
            const label patchi = cycPatch[k];

            const cyclicFvPatch& cyc =
                refCast<const cyclicFvPatch>(mesh.boundary()[patchi]);
            const cyclicFvPatch& nbr = cyc.neighbFvPatch();
            const label nbri = nbr.index();

            // Taken from OpenFOAM's own boundary geometry, which already
            // carries the cyclic transform: nothing here rebuilds it.
            const tmp<vectorField> tDelta(cyc.delta());
            const vectorField& pDelta = tDelta();

            const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
            const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
            const fvsPatchScalarField& pW = w.boundaryField()[patchi];
            const fvsPatchScalarField& pDc = dc.boundaryField()[patchi];

            const labelUList& ownFc = cyc.faceCells();
            const labelUList& neiFc = nbr.faceCells();

            const label ownStart = mesh.boundaryMesh()[patchi].start();
            const label neiStart = mesh.boundaryMesh()[nbri].start();

            while (k < nCyc && cycPatch[k] == patchi)
            {
                const label f = cycFace[k];
                const label d = nIF + k;

                devOwn[d] = ownFc[f];
                devNei[d] = neiFc[f];

                for (direction cmpt = 0; cmpt < 3; ++cmpt)
                {
                    SfFlat[3*d + cmpt] = pSf[f][cmpt];
                    deltaFlat[3*d + cmpt] = pDelta[f][cmpt];
                }

                magSfDev[d] = pMagSf[f];
                wDev[d] = pW[f];
                dcDev[d] = pDc[f];

                // Both of the pair's global faces map onto this one device
                // face, the sign saying which side each cell is on.
                faceToC[ownStart + f - nIF] = d;
                faceToCSign[ownStart + f - nIF] = 1;
                faceToC[neiStart + f - nIF] = d;
                faceToCSign[neiStart + f - nIF] = -1;

                ++k;
            }
        }
    }

    // ---- processor faces, appended as internal faces with a ghost --------
    // As the cyclic case, but the neighbour belongs to another rank, so the
    // device reads a ghost that exchangeHalo() fills each Runge-Kutta stage.
    // Patch-major for the same reason.
    forAll(procNbr_, p)
    {
        const label patchi = procPatch_[procStart_[p]];
        const fvPatch& fvp = mesh.boundary()[patchi];

        const tmp<vectorField> tDelta(fvp.delta());
        const vectorField& pDelta = tDelta();

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& pW = w.boundaryField()[patchi];
        const fvsPatchScalarField& pDc = dc.boundaryField()[patchi];

        const labelUList& ownFc = fvp.faceCells();
        const label start = mesh.boundaryMesh()[patchi].start();

        for (label h = procStart_[p]; h < procStart_[p + 1]; ++h)
        {
            const label f = procFace_[h];
            const label d = nIF + nCyc + h;

            devOwn[d] = ownFc[f];
            devNei[d] = nCells + h;          // the ghost

            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                SfFlat[3*d + cmpt] = pSf[f][cmpt];
                deltaFlat[3*d + cmpt] = pDelta[f][cmpt];
            }

            magSfDev[d] = pMagSf[f];
            wDev[d] = pW[f];
            dcDev[d] = pDc[f];

            // Only this rank's face maps onto the device face; the other side lives
            // in another rank's mesh, which is the one way this differs from cyclic.
            faceToC[start + f - nIF] = d;
            faceToCSign[start + f - nIF] = 1;
        }
    }

    // fvPatch::delta() is documented as owner-to-neighbour for a coupled
    // patch, but this cannot be compiled against those headers, so the
    // convention is checked rather than assumed: a wrong sign would degrade
    // the reconstruction at the seam quietly instead of failing.
    const label nAppended = nCyc + nHalo_;

    if (nAppended > 0)
    {
        label nBadSign = 0, nBadLength = 0;
        scalar worstRatio = 1;

        for (label k = 0; k < nAppended; k++)
        {
            const label d = nIF + k;

            const vector s(SfFlat[3*d], SfFlat[3*d + 1], SfFlat[3*d + 2]);
            const vector dl
            (
                deltaFlat[3*d], deltaFlat[3*d + 1], deltaFlat[3*d + 2]
            );

            // delta must point out of the owner, the same way as Sf
            if ((s & dl) <= 0) { nBadSign++; }

            // and its normal projection must be 1/deltaCoeffs
            const scalar proj = mag(dl & s)/max(magSfDev[d], vSmall);
            const scalar expect = 1/max(dcDev[d], vSmall);
            const scalar ratio =
                max(proj, vSmall)/max(expect, vSmall);

            if (ratio < 0.5 || ratio > 2.0)
            {
                nBadLength++;
                if (mag(ratio - 1) > mag(worstRatio - 1)) { worstRatio = ratio; }
            }
        }

        if (nBadSign > 0 || nBadLength > 0)
        {
            FatalErrorInFunction
                << "the coupled face geometry does not match the convention"
                << " this solver assumes (" << nCyc << " cyclic, "
                << nHalo_ << " processor)." << nl
                << "    faces with delta opposed to Sf: " << nBadSign
                << " of " << nAppended << nl
                << "    faces whose |delta.n| disagrees with 1/deltaCoeffs: "
                << nBadLength << " of " << nAppended
                << " (worst ratio " << worstRatio << ")" << nl
                << "fvPatch::delta() is evidently not the owner-to-neighbour"
                << " vector on this patch type." << nl
                << "Run this case with solver cdmFoam until that is"
                << " corrected; do not trust a cyclic result from the device."
                << exit(FatalError);
        }
    }

    scalarField bSfFlat(3*nBF), bMagSfFlat(nBF), bDcFlat(nBF);
    labelList bOwn(nBF);
    forAll(bPatch_, b)
    {
        const label patchi = bPatch_[b];
        const label f = bFace_[b];
        const vector s(Sf.boundaryField()[patchi][f]);
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            bSfFlat[3*b + cmpt] = s[cmpt];
        }
        bMagSfFlat[b] = magSf.boundaryField()[patchi][f];
        bDcFlat[b] = dc.boundaryField()[patchi][f];
        bOwn[b] = mesh.boundary()[patchi].faceCells()[f];
    }

    // ---- the cells behind a boundary face, de-duplicated -----------------
    // bOwn repeats at corners; marking and then collecting gives a sorted,
    // unique list in one pass and keeps the device-side gather contiguous.
    // Cells behind a cyclic patch are kept even though the device no longer
    // needs patch values there: the host still evaluates those patch fields,
    // and the diagnostics and output read them.
    {
        labelList onBoundary(nCells, 0);
        forAll(bOwn, b)
        {
            onBoundary[bOwn[b]] = 1;
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (fvp.size() == 0 || !fvp.coupled()) { continue; }

            const labelUList& fc = fvp.faceCells();
            forAll(fc, i)
            {
                onBoundary[fc[i]] = 1;
            }
        }

        label nBC = 0;
        forAll(onBoundary, celli)
        {
            if (onBoundary[celli] == 1) { ++nBC; }
        }

        bCell_.setSize(nBC);
        label i = 0;
        forAll(onBoundary, celli)
        {
            if (onBoundary[celli] == 1) { bCell_[i++] = celli; }
        }

        bcRho_.setSize(nBC);
        bcU_.setSize(nBC);
        bcP_.setSize(nBC);
        bcT_.setSize(nBC);
        bcY_.setSize(nBC);

        cdmGpuSetBoundaryCells(gpu_, nBC, bCell_.begin());
    }

    // ---- CSR cell-to-face map --------------------------------------------
    // Gathers rather than scatters: a gradient or divergence is a loop over a
    // cell's own faces, with no atomics.
    labelList cfOff(nCells + 1, 0);
    forAll(mesh.cells(), celli)
    {
        label n = 0;
        const labelList& cFaces = mesh.cells()[celli];
        forAll(cFaces, i)
        {
            const label f = cFaces[i];
            if (f < nIF) { ++n; }
            else if (faceToC[f - nIF] >= 0) { ++n; }
            else if (faceToB[f - nIF] >= 0) { ++n; }
        }
        cfOff[celli + 1] = cfOff[celli] + n;
    }

    labelList cfFace(cfOff[nCells]);
    List<signed char> cfSign(cfOff[nCells]);

    forAll(mesh.cells(), celli)
    {
        label k = cfOff[celli];
        const labelList& cFaces = mesh.cells()[celli];
        forAll(cFaces, i)
        {
            const label f = cFaces[i];
            if (f < nIF)
            {
                cfFace[k] = f;
                cfSign[k] = (own[f] == celli) ? 1 : -1;
                ++k;
            }
            else if (faceToC[f - nIF] >= 0)
            {
                // A cyclic face: an internal face on the device, and this
                // cell is on one side of it or the other.
                cfFace[k] = faceToC[f - nIF];
                cfSign[k] = faceToCSign[f - nIF];
                ++k;
            }
            else
            {
                const label b = faceToB[f - nIF];
                if (b >= 0)
                {
                    cfFace[k] = nIFdev + b;
                    cfSign[k] = 1;         // a cell always owns its boundary face
                    ++k;
                }
            }
        }
    }

    // Invariant: every internal face appears twice, once with each sign, and
    // every boundary face once with a positive sign. One pass, and it turns
    // a mis-mapped cyclic pair into a startup failure rather than a
    // conservation error nobody notices.
    {
        labelList nSeen(nIFdev + nBF, 0);
        labelList signSum(nIFdev + nBF, 0);

        forAll(cfFace, i)
        {
            nSeen[cfFace[i]]++;
            signSum[cfFace[i]] += cfSign[i];
        }

        // A processor face is the deliberate exception: a ghost owns no faces, so
        // this rank sees it once, from the owner. The other side is counted on
        // the other rank, which makes the pair conservative across the
        // decomposition rather than within this map.
        label badInternal = 0, badProc = 0, badBoundary = 0;

        for (label f = 0; f < nIF + nCyc; f++)
        {
            if (nSeen[f] != 2 || signSum[f] != 0) { badInternal++; }
        }
        for (label f = nIF + nCyc; f < nIFdev; f++)
        {
            if (nSeen[f] != 1 || signSum[f] != 1) { badProc++; }
        }
        for (label b = nIFdev; b < nIFdev + nBF; b++)
        {
            if (nSeen[b] != 1 || signSum[b] != 1) { badBoundary++; }
        }

        if (badInternal > 0 || badProc > 0 || badBoundary > 0)
        {
            FatalErrorInFunction
                << "the cell-to-face map is inconsistent." << nl
                << "    internal faces not seen exactly once from each side: "
                << badInternal << " of " << nIF + nCyc << nl
                << "    processor faces not seen exactly once as owner:     "
                << badProc << " of " << nHalo_ << nl
                << "    boundary faces not seen exactly once as owner:      "
                << badBoundary << " of " << nBF << nl
                << "A divergence built on this map would not conserve."
                << exit(FatalError);
        }
    }

    cdmGpuSetMesh
    (
        gpu_, nCells, nHalo_, nIFdev, nBF,
        devOwn.begin(), devNei.begin(),
        SfFlat.begin(), magSfDev.begin(),
        wDev.begin(), dcDev.begin(),
        deltaFlat.begin(), Vcells.begin(),
        bOwn.begin(), bSfFlat.begin(), bMagSfFlat.begin(), bDcFlat.begin(),
        cfOff.begin(), cfFace.begin(), cfSign.begin()
    );

    nCellsUploaded_ = nCells;
    nFacesUploaded_ = nFaces;
    meshTopoHash_ = cdmTopoHash(mesh);

    cycPatch_ = cycPatch;
    cycFace_ = cycFace;
    nIFdev_ = nIFdev;
    phiScratch_.setSize(nIFdev);

    // The cells this rank sends, and room for the 25 doubles a face each way.
    // Sized here rather than per stage: with AMR the halo changes size only
    // at a topology change.
    if (nHalo_ > 0)
    {
        cdmGpuSetHaloCells(gpu_, nHalo_, haloCell_.begin());
        hSend_.setSize(25*nHalo_);
        hRecv_.setSize(25*nHalo_);
        hPack_.setSize(25*nHalo_);
        hUnpack_.setSize(25*nHalo_);
    }

    if (nUploadMesh_ <= 1)
    {
        Info<< "cdmGpuFoam: uploaded " << nCells << " cells, "
            << nIFdev << " internal faces, " << nBF << " boundary faces"
            << " (" << (nFaces - nIF - nBF - 2*nCyc) << " on empty patches,"
            << " excluded)" << nl;

        if (nCyc > 0)
        {
            Info<< "    " << nCyc << " cyclic faces promoted to internal"
                << " faces, so the seam is second order and conservative"
                << nl;
        }
    }

    if (nUploadMesh_ <= 1)
    {
        Info<< "    " << bCell_.size() << " cells behind a boundary face ("
            << 100.0*scalar(bCell_.size())/scalar(nCells)
            << " % of the mesh); only these are exchanged between"
            << " Runge-Kutta stages" << endl;
    }

    tUploadMesh_ += cdmGpuWallTime() - tEnterUpload;
}


void Foam::solvers::cdmGpuFoam::uploadModel()
{
    const dictionary& schemes = mesh.schemes().dict();

    const word fluxName
    (
        schemes.lookupOrDefault<word>("fluxScheme", "HLLC")
    );

    int flux = CDM_GPU_HLLC;
    if (fluxName == "HLLCLM")  { flux = CDM_GPU_HLLCLM; }
    else if (fluxName != "HLLC")
    {
        FatalErrorInFunction
            << "fluxScheme " << fluxName << " has no GPU implementation;"
            << " cdmGpuFoam supports HLLC and HLLCLM"
            << exit(FatalError);
    }

    word limName("vanLeer");
    if (schemes.isDict("interpolationSchemes"))
    {
        const dictionary& interp = schemes.subDict("interpolationSchemes");
        limName = interp.lookupOrDefault<word>("reconstruct(rho)", "vanLeer");
    }

    int lim = -1;
    if      (limName == "upwind")     lim = CDM_GPU_UPWIND;
    else if (limName == "vanLeer")    lim = CDM_GPU_VANLEER;
    else if (limName == "Minmod")     lim = CDM_GPU_MINMOD;
    else if (limName == "ROUNDF")     lim = CDM_GPU_ROUNDF;
    else if (limName == "ROUNDA")     lim = CDM_GPU_ROUNDA;
    else if (limName == "ROUNDAplus") lim = CDM_GPU_ROUNDAPLUS;
    else if (limName == "ROUNDL")     lim = CDM_GPU_ROUNDL;
    else
    {
        FatalErrorInFunction
            << "reconstruct(rho) " << limName << " has no GPU implementation;"
            << " use upwind, vanLeer, Minmod or one of the ROUND schemes"
            << exit(FatalError);
    }

    // One limiter serves every reconstructed variable on the device, so a
    // case asking for different schemes per variable cannot be reproduced
    // faithfully. Say so rather than quietly using the first.
    if (schemes.isDict("interpolationSchemes"))
    {
        const dictionary& interp = schemes.subDict("interpolationSchemes");
        const wordList others({"reconstruct(p)", "reconstruct(Y)"});
        forAll(others, i)
        {
            if (interp.found(others[i]))
            {
                const word other
                (
                    interp.lookupOrDefault<word>(others[i], limName)
                );
                if (other != limName && other != limName + "01")
                {
                    WarningInFunction
                        << others[i] << " is " << other << " but reconstruct(rho)"
                        << " is " << limName << "; the GPU back end applies one"
                        << " limiter to every variable and will use "
                        << limName << endl;
                }
            }
        }
    }

    // Uniform for the calorically perfect CDM gas, so the average is exact.
    // Held because the boundary exchange needs it every Runge-Kutta stage.
    Cv_ = gAverage(thermo_.Cv().primitiveField());

    cdmGpuSetModel
    (
        gpu_,
        gamma_, Cv_, eOffset_,
        cdm_.A(), cdm_.Ta(), cdm_.q(),
        cdm_.kappa0(), cdm_.n(), cdm_.Pr(), cdm_.Le(),
        cdm_.maxSubstepFraction(), cdm_.maxSubsteps(),
        flux, lim, viscous_ ? 1 : 0
    );

    Info<< "cdmGpuFoam: back end " << cdmGpuBackend()
        << ", flux " << fluxName << ", limiter " << limName
        << ", viscous " << (viscous_ ? "yes" : "no") << endl;
}


void Foam::solvers::cdmGpuFoam::uploadBoundary()
{
    const label nBF = bPatch_.size();

    scalarField bRho(nBF), bP(nBF), bY(nBF), bU(3*nBF);

    forAll(bPatch_, b)
    {
        const label patchi = bPatch_[b];
        const label f = bFace_[b];

        bRho[b] = rho_.boundaryField()[patchi][f];
        bP[b]   = p_.boundaryField()[patchi][f];
        bY[b]   = Y_.boundaryField()[patchi][f];

        const vector& u = U_.boundaryField()[patchi][f];
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            bU[3*b + cmpt] = u[cmpt];
        }
    }

    cdmGpuUploadBoundary
    (
        gpu_, bRho.begin(), bU.begin(), bP.begin(), bY.begin()
    );
}


void Foam::solvers::cdmGpuFoam::exchangeBoundary(const bool fetchFromDevice)
{
    volScalarField& e = thermo_.he();

    if (fetchFromDevice)
    {
        // The device's primitives at the cells behind the patches
        cdmGpuDownloadBoundaryCells
        (
            gpu_,
            bcRho_.begin(),
            reinterpret_cast<double*>(bcU_.begin()),
            bcP_.begin(),
            bcT_.begin(),
            bcY_.begin()
        );

        scalarField& rhoI = rho_.primitiveFieldRef();
        vectorField& UI = U_.primitiveFieldRef();
        scalarField& YI = Y_.primitiveFieldRef();
        scalarField& pI = p_.primitiveFieldRef();
        scalarField& eI = e.primitiveFieldRef();

        // e = p/(rho (gamma - 1)) + eOffset and T = p/(rho R) with
        // R = Cv (gamma - 1), so the sensible energy is Cv T + eOffset.
        forAll(bCell_, i)
        {
            const label celli = bCell_[i];

            rhoI[celli] = bcRho_[i];
            UI[celli] = bcU_[i];
            YI[celli] = bcY_[i];
            pI[celli] = bcP_[i];
            eI[celli] = Cv_*bcT_[i] + eOffset_;
        }
    }

    // Evaluate the patch fields in the same order decode() does. Every
    // boundary condition in the package runs here, unported and unchanged,
    // each reading the patchInternalField just written.
    U_.correctBoundaryConditions();
    Y_.correctBoundaryConditions();
    e.correctBoundaryConditions();
    p_.correctBoundaryConditions();

    // decode() takes the patch density from psi*p. Here psi = 1/(R T) and T
    // follows from the energy the patch just evaluated, so this is the same
    // number without a whole-field thermo_.correct().
    {
        const scalar R = Cv_*(gamma_ - 1);

        volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
        const volScalarField::Boundary& pBf = p_.boundaryField();
        const volScalarField::Boundary& eBf = e.boundaryField();

        forAll(rhoBf, patchi)
        {
            // Left alone on a coupled patch: rho_ carries psiThermo's calculated
            // patch fields, so decode() assigns it from psi*p instead. Nothing reads
            // a cyclic patch's density between decodes.
            if (mesh.boundary()[patchi].coupled()) { continue; }

            fvPatchScalarField& rhop = rhoBf[patchi];
            const fvPatchScalarField& pp = pBf[patchi];
            const fvPatchScalarField& ep = eBf[patchi];

            forAll(rhop, f)
            {
                rhop[f] = pp[f]/(R*(ep[f] - eOffset_)/Cv_);
            }
        }
    }

    uploadBoundary();
}


void Foam::solvers::cdmGpuFoam::haloRound()
{
    const label n = nHalo_;

    // Where each quantity sits in the 25-doubles-a-face buffer, and how
    // many doubles it occupies per face.
    static const label off[9] = {0, 1, 4, 5, 6, 7, 10, 19, 22};
    static const label str[9] = {1, 3, 1, 1, 1, 3,  9,  3,  3};

    // A tag of our own. PstreamBuffers collided with the processor-patch
    // exchange GeometricBoundaryField::evaluate() performs later in the
    // same step.
    static const int haloTag = 1729;

    {
        scalar* s = hSend_.begin();
        cdmGpuDownloadHaloCells
        (
            gpu_,
            s,            s +    n,  s +  4*n,  s +  5*n,  s +  6*n,
            s + 7*n,      s + 10*n,  s + 19*n,  s + 22*n
        );
    }

    // Pack patch by patch, so each payload is contiguous on the wire. Both
    // ranks pack the same blocks in the same order over faces OpenFOAM has
    // already matched, so the unpack is symmetric.
    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;
        scalar* w = hPack_.begin() + 25*a;
        label q = 0;

        for (label blk = 0; blk < 9; blk++)
        {
            const scalar* src = hSend_.begin() + off[blk]*n + str[blk]*a;
            const label m = str[blk]*np;
            for (label i = 0; i < m; i++) { w[q++] = src[i]; }
        }
    }

    const label req0 = UPstream::nRequests();

    // Receives first, then sends: a posted receive cannot be mistaken for
    // an unexpected message.
    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;

        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            procNbr_[p],
            reinterpret_cast<char*>(hUnpack_.begin() + 25*a),
            25*np*sizeof(scalar),
            haloTag
        );
    }

    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;

        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            procNbr_[p],
            reinterpret_cast<const char*>(hPack_.begin() + 25*a),
            25*np*sizeof(scalar),
            haloTag
        );
    }

    UPstream::waitRequests(req0);

    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;
        const scalar* r = hUnpack_.begin() + 25*a;
        label q = 0;

        for (label blk = 0; blk < 9; blk++)
        {
            scalar* dst = hRecv_.begin() + off[blk]*n + str[blk]*a;
            const label m = str[blk]*np;
            for (label i = 0; i < m; i++) { dst[i] = r[q++]; }
        }
    }

    {
        const scalar* r = hRecv_.begin();
        cdmGpuUploadHalo
        (
            gpu_,
            r,            r +    n,  r +  4*n,  r +  5*n,  r +  6*n,
            r + 7*n,      r + 10*n,  r + 19*n,  r + 22*n
        );
    }
}


void Foam::solvers::cdmGpuFoam::exchangeHalo()
{
    if (nHalo_ == 0) { return; }

    // With a single round the gradient this rank sends was built from last
    // stage's ghost, stale by one stage, and the neighbour uses it for the
    // far side of its Riemann problem. This rank's own gradients are fine,
    // so the error is invisible locally.
    // Measured on RDE2D_T6 at 135k cells: with the seam crossing the ignition
    // region one round disagreed with serial by 5e-9 relative after six
    // steps, while a seam in quiescent gas was bit-identical.
    // So: round one makes the ghost states current, the gradient pass
    // recomputes this rank's gradients from them, and round two sends those.
    // Nothing further is needed, since a gradient depends only on cell values.
    // Two messages a stage rather than one. Round one needs only the seven
    // state values and round two only the eighteen gradients, so 25 twice
    // could become 7 and 18; worth doing if the halo shows up in a profile.
    haloRound();
    cdmGpuComputeGradients(gpu_);
    haloRound();
}


void Foam::solvers::cdmGpuFoam::downloadMassFlux()
{
    scalarField bPhi(bPatch_.size());

    // Into a scratch buffer, not straight into phi_: the device's internal
    // face list is longer than the mesh's, by one entry per cyclic face.
    cdmGpuDownloadMassFlux(gpu_, phiScratch_.begin(), bPhi.begin());

    scalarField& phiI = phi_.primitiveFieldRef();
    const label nIF = mesh.nInternalFaces();

    for (label facei = 0; facei < nIF; facei++)
    {
        phiI[facei] = phiScratch_[facei];
    }

    forAll(bPatch_, b)
    {
        phi_.boundaryFieldRef()[bPatch_[b]][bFace_[b]] = bPhi[b];
    }

    // One flux per cyclic face, written to both sides with opposite signs,
    // so phi_ stays conservative across the seam as it is everywhere else.
    forAll(cycPatch_, k)
    {
        const label patchi = cycPatch_[k];
        const label f = cycFace_[k];
        const scalar flux = phiScratch_[nIF + k];

        const label nbri =
            refCast<const cyclicFvPatch>
            (
                mesh.boundary()[patchi]
            ).neighbFvPatch().index();

        phi_.boundaryFieldRef()[patchi][f] = flux;
        phi_.boundaryFieldRef()[nbri][f] = -flux;
    }

    // Only this rank's side is written here. The other rank writes the other
    // side from its own device, bit-identically, because both solved the same
    // Riemann problem from the same reconstructed states.
    forAll(procPatch_, h)
    {
        phi_.boundaryFieldRef()[procPatch_[h]][procFace_[h]] =
            phiScratch_[nIF + cycPatch_.size() + h];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::cdmGpuFoam::cdmGpuFoam(fvMesh& mesh)
:
    cdmFoam(mesh),
    gpu_(cdmGpuCreate()),
    nCellsUploaded_(-1),
    nFacesUploaded_(-1),
    meshTopoHash_(0),
    bPatch_(),
    bFace_(),
    cycPatch_(),
    cycFace_(),
    nHalo_(0),
    nIFdev_(0),
    phiScratch_(),
    tUploadMesh_(0),
    nUploadMesh_(0),
    nFallback_(0),
    Cv_(0),
    bCell_(),
    bcRho_(),
    bcU_(),
    bcP_(),
    bcT_(),
    bcY_(),
    tUpState_(0),
    tDecode_(0),
    tUpBnd_(0),
    tBackRHS_(0),
    tDownRHS_(0),
    tDownPhi_(0),
    tOutside_(0),
    tLastExit_(0),
    nProfCalls_(0)
{
    if (Pstream::parRun())
    {
        Info<< "cdmGpuFoam: parallel run, " << Pstream::nProcs()
            << " ranks.  One device per rank; processor faces are carried"
            << " as device internal faces with ghost neighbours." << endl;

        // Per rank, because a misbinding is otherwise invisible: every rank on
        // device 0 is correct and slow, and slow in a way that looks exactly like
        // poor scaling. std::getenv rather than getEnv, which warns when unset.
        const char* vis = std::getenv("CUDA_VISIBLE_DEVICES");

        Pout<< "cdmGpuFoam: rank " << Pstream::myProcNo()
            << " of " << Pstream::nProcs()
            << " on " << hostName()
            << ", back end " << cdmGpuBackend()
            << ", CUDA_VISIBLE_DEVICES="
            << (vis != nullptr ? vis : "<unset: all ranks on device 0>")
            << ", " << mesh.nCells() << " cells"
            << endl;
    }

    if (positivityLimiter_)
    {
        FatalErrorInFunction
            << "The positivity-preserving flux limiter needs the face fluxes"
            << " on the host, which cdmGpuFoam does not produce there." << nl
            << "Set positivityLimiter { active no; } in fvSolution." << nl
            << "The first-order reconstruction fallback is still active and is"
            << " reported as usual."
            << exit(FatalError);
    }

    Info<< nl << "cdmGpuFoam: cdmFoam with the right-hand side on a"
        << " device" << nl
        << "    reconstruction, Riemann fluxes, diffusive terms and the"
        << " divergence run on the back end;" << nl
        << "    boundary conditions, chemistry and everything else are"
        << " unchanged on the host." << nl << endl;

    uploadModel();
    uploadMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::cdmGpuFoam::~cdmGpuFoam()
{
    reportProfile();
    cdmGpuDestroy(gpu_);
}


// * * * * * * * * * * * * * * * * Profiling  * * * * * * * * * * * * * * * * //

void Foam::solvers::cdmGpuFoam::reportProfile() const
{
    // v2's advance() calls the back end directly, so nProfCalls_ stays
    // zero; the mesh-rebuild count still has something to say.
    if (nProfCalls_ == 0 && nUploadMesh_ == 0)
    {
        return;
    }

    const double tBackEnd =
        tUpState_ + tDecode_ + tUpBnd_ + tBackRHS_ + tDownRHS_ + tDownPhi_
      + tUploadMesh_;
    const double tTotal = tBackEnd + tOutside_;

    if (tTotal <= 0)
    {
        return;
    }

    const double pc = 100.0/tTotal;

    Info<< nl << "cdmGpuFoam time split over " << nProfCalls_
        << " right-hand-side evaluations" << nl
        << "    (tracked wall clock " << tTotal << " s; percentages of that)"
        << nl
        << "    uploadState       " << tUpState_ << " s  "
        << tUpState_*pc << " %" << nl
        << "    decode (device)   " << tDecode_ << " s  "
        << tDecode_*pc << " %" << nl
        << "    uploadBoundary    " << tUpBnd_ << " s  "
        << tUpBnd_*pc << " %" << nl
        << "    computeRHS        " << tBackRHS_ << " s  "
        << tBackRHS_*pc << " %" << nl
        << "    downloadRHS       " << tDownRHS_ << " s  "
        << tDownRHS_*pc << " %" << nl
        << "    uploadMesh        " << tUploadMesh_ << " s  "
        << tUploadMesh_*pc << " %  (" << nUploadMesh_ << " rebuilds)" << nl
        << "    downloadMassFlux  " << tDownPhi_ << " s  "
        << tDownPhi_*pc << " %" << nl
        << "    ---- back end     " << tBackEnd << " s  "
        << tBackEnd*pc << " %" << nl
        << "    ---- host rest    " << tOutside_ << " s  "
        << tOutside_*pc << " %" << nl
        << "         (Runge-Kutta updates, host decode and boundary"
        << " evaluation," << nl
        << "          chemistry, per-step diagnostics, soot foil, writing)"
        << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::cdmGpuFoam::computeRHS
(
    const scalar deltaT,
    scalarField& dRho,
    vectorField& dRhoU,
    scalarField& dRhoE,
    scalarField& dRhoY
)
{
    // deltaT is only used by the positivity limiter, which this solver
    // refuses to run with

    // Everything since the previous return from this function was host-side:
    // the Runge-Kutta updates, decode() with its boundary evaluation, the
    // chemistry integration and the per-step diagnostics.
    const double tEnter = cdmGpuWallTime();
    if (tLastExit_ > 0)
    {
        tOutside_ += tEnter - tLastExit_;
    }

    if (cdmTopoHash(mesh) != meshTopoHash_)
    {
        uploadMesh();
    }

    // ---- state to the device ---------------------------------------------
    double t0 = cdmGpuWallTime();

    cdmGpuUploadState
    (
        gpu_,
        rho_.primitiveField().begin(),
        reinterpret_cast<const double*>(rhoU_.primitiveField().begin()),
        rhoE_.primitiveField().begin(),
        rhoY_.primitiveField().begin()
    );

    double t1 = cdmGpuWallTime();
    tUpState_ += t1 - t0;

    const int bad = cdmGpuDecode(gpu_);

    t0 = cdmGpuWallTime();
    tDecode_ += t0 - t1;

    if (bad != 0)
    {
        FatalErrorInFunction
            << "non-physical state in cell " << (bad - 1)
            << " on the device: rho = " << rho_[bad - 1]
            << ", p = " << p_[bad - 1] << ", T = " << thermo_.T()[bad - 1]
            << exit(FatalError);
    }

    uploadBoundary();

    t1 = cdmGpuWallTime();
    tUpBnd_ += t1 - t0;

    // ---- the right-hand side ---------------------------------------------
    // The halo first: the fluxes about to be formed read a ghost state and a
    // ghost gradient at every processor face. A no-op in serial.
    exchangeHalo();

    int nf = 0;
    cdmGpuComputeRHS(gpu_, &nf);
    nFallback_ += nf;

    t0 = cdmGpuWallTime();
    tBackRHS_ += t0 - t1;

    cdmGpuDownloadRHS
    (
        gpu_,
        dRho.begin(),
        reinterpret_cast<double*>(dRhoU.begin()),
        dRhoE.begin(),
        dRhoY.begin()
    );

    t1 = cdmGpuWallTime();
    tDownRHS_ += t1 - t0;

    // ---- keep phi consistent ---------------------------------------------
    scalarField bPhi(bPatch_.size());
    cdmGpuDownloadMassFlux(gpu_, phi_.primitiveFieldRef().begin(), bPhi.begin());

    forAll(bPatch_, b)
    {
        phi_.boundaryFieldRef()[bPatch_[b]][bFace_[b]] = bPhi[b];
    }

    tDownPhi_ += cdmGpuWallTime() - t1;
    nProfCalls_++;

    if (nf > 0)
    {
        Info<< "cdmGpuFoam: " << nf
            << " face reconstructions dropped to first order"
            << " (non-positive rho or p)" << endl;
    }

    tLastExit_ = cdmGpuWallTime();
}


void Foam::solvers::cdmGpuFoam::advance()
{
    const scalar deltaT = runTime.deltaTValue();

    // SSP Runge-Kutta coefficients (Shu & Osher 1988), as cdmFoam has them.
    // Kept here rather than shared so the two implementations can be compared
    // without one standing in for the other.
    static const scalar a[3][3] =
    {
        {0, 0, 0},
        {0, 0.5, 0},
        {0, 0.75, 1.0/3.0}
    };
    static const scalar b[3][3] =
    {
        {1, 0, 0},
        {1, 0.5, 0},
        {1, 0.25, 2.0/3.0}
    };
    const label s = nRKStages_ - 1;

    nLimitedFaces_ = 0;

    if (cdmTopoHash(mesh) != meshTopoHash_)
    {
        uploadMesh();
    }

    // ---- the host hands the step over ------------------------------------
    cdmGpuUploadState
    (
        gpu_,
        rho_.primitiveField().begin(),
        reinterpret_cast<const double*>(rhoU_.primitiveField().begin()),
        rhoE_.primitiveField().begin(),
        rhoY_.primitiveField().begin()
    );

    if (cdmGpuDecode(gpu_) != 0)
    {
        FatalErrorInFunction
            << "non-physical state on the device at the start of the step"
            << exit(FatalError);
    }

    cdmGpuZeroQdot(gpu_);

    // No gather: the host's primitives come from the decode() that ended the
    // previous step and the device's from the same conserved state. The patch
    // fields are still re-evaluated, for the time-dependent conditions.
    exchangeBoundary(false);

    double tMark = CDMWallTime();

    // ---- the whole step, on the device -----------------------------------
    if (chemistry_ && strang_)
    {
        cdmGpuReact(gpu_, 0.5*deltaT);

        if (cdmGpuDecode(gpu_) != 0)
        {
            FatalErrorInFunction
                << "non-physical state on the device after the first"
                << " half-step of chemistry"
                << exit(FatalError);
        }

        exchangeBoundary();
    }

    tChem_ += CDMWallTime() - tMark;

    cdmGpuStoreStage(gpu_);

    for (label stage = 0; stage < nRKStages_; stage++)
    {
        tMark = CDMWallTime();

        // Every stage, because the ghosts must carry the state the
        // neighbour holds *now*, not at the start of the step.
        exchangeHalo();

        int nf = 0;
        cdmGpuComputeRHS(gpu_, &nf);
        nFallback_ += nf;

        const double tRHSDone = CDMWallTime();
        tRHS_ += tRHSDone - tMark;

        if (nf > 0)
        {
            Info<< "cdmGpuFoam: " << nf
                << " face reconstructions dropped to first order"
                << " (non-positive rho or p)" << endl;
        }

        // phi_ belongs to the flux-dependent boundary conditions, so it has
        // to be current before the patch fields are evaluated below
        downloadMassFlux();

        cdmGpuRKStage(gpu_, deltaT, a[s][stage], b[s][stage]);

        if (cdmGpuDecode(gpu_) != 0)
        {
            FatalErrorInFunction
                << "non-physical state on the device after Runge-Kutta stage "
                << stage
                << exit(FatalError);
        }

        exchangeBoundary();

        tStage_ += CDMWallTime() - tRHSDone;
    }

    tMark = CDMWallTime();

    if (chemistry_)
    {
        cdmGpuReact(gpu_, strang_ ? 0.5*deltaT : deltaT);

        if (cdmGpuDecode(gpu_) != 0)
        {
            FatalErrorInFunction
                << "non-physical state on the device after the final"
                << " chemistry step"
                << exit(FatalError);
        }

        // No exchange: nothing on the device reads the patch values again before
        // the next step, and the host decode() below recomputes every boundary
        // field from the downloaded state, so an exchange would be undone.
    }

    tChem_ += CDMWallTime() - tMark;
    tMark = CDMWallTime();

    // ---- and hands it back -----------------------------------------------
    cdmGpuDownloadState
    (
        gpu_,
        rho_.primitiveFieldRef().begin(),
        reinterpret_cast<double*>(rhoU_.primitiveFieldRef().begin()),
        rhoE_.primitiveFieldRef().begin(),
        rhoY_.primitiveFieldRef().begin()
    );

    cdmGpuDownloadQdot(gpu_, Qdot_.primitiveFieldRef().begin());

    // Qdot accumulated over the step; the host divides by deltaT below, so
    // the device copy is left accumulated and zeroed at the next step.
    Qdot_.primitiveFieldRef() /= deltaT;
    Qdot_.correctBoundaryConditions();

    // One host decode per step rather than four, which is what makes the step
    // boundary identical to cdmFoam's: the thermo update, the fvConstraints
    // and the boundary evaluation all run exactly as they would.
    decode();

    if (frontTracking_ && mesh.time().timeIndex() % frontInterval_ == 0)
    {
        trackFronts();
    }

    // Every per-step monitor in one pass, as in cdmFoam::advance(), which this
    // method overrides and which is therefore never executed on the GPU path.
    // Was eight sweeps over the cells plus four full-size temporary fields.
    scalar pMin = great, pMax = -great;
    scalar TMin = great, TMax = -great;
    scalar QdotMax = -great;
    scalar mass = 0, energy = 0;
    label cellTMax = -1;

    {
        const scalarField& pI = p_.primitiveField();
        const scalarField& TI = thermo_.T().primitiveField();
        const scalarField& QI = Qdot_.primitiveField();
        const scalarField& rhoI = rho_.primitiveField();
        const scalarField& rhoEI = rhoE_.primitiveField();
        const scalarField& rhoYI = rhoY_.primitiveField();
        const scalarField& V = mesh.V().primitiveField();
        const scalar qCdm = cdm_.q();

        forAll(pI, celli)
        {
            const scalar pc = pI[celli];
            if (pc < pMin) { pMin = pc; }
            if (pc > pMax) { pMax = pc; }

            const scalar Tc = TI[celli];
            if (Tc < TMin) { TMin = Tc; }
            if (Tc > TMax) { TMax = Tc; cellTMax = celli; }

            if (QI[celli] > QdotMax) { QdotMax = QI[celli]; }

            const scalar Vc = V[celli];
            mass += rhoI[celli]*Vc;
            energy += (rhoEI[celli] + qCdm*rhoYI[celli])*Vc;
        }
    }

    const scalar TMaxGlobal = returnReduce(TMax, maxOp<scalar>());

    if (TMaxGlobal > TWarn_ && cellTMax >= 0 && TMax == TMaxGlobal)
    {
        Pout<< "Hot spot: T = " << TMax
            << " at " << mesh.C()[cellTMax]
            << " rho = " << rho_[cellTMax]
            << " p = " << p_[cellTMax]
            << " U = " << U_[cellTMax]
            << " V = " << mesh.V()[cellTMax] << endl;
    }

    Info<< "min/max(p) = " << returnReduce(pMin, minOp<scalar>()) << ", "
        << returnReduce(pMax, maxOp<scalar>())
        << "  min/max(T) = " << returnReduce(TMin, minOp<scalar>()) << ", "
        << TMaxGlobal
        << "  max(Qdot) = " << returnReduce(QdotMax, maxOp<scalar>()) << endl;

    Info<< "Totals: mass = " << returnReduce(mass, sumOp<scalar>())
        << "  energy = " << returnReduce(energy, sumOp<scalar>())
        << "  cells = " << mesh.globalData().nTotalCells() << endl;

    tTail_ += CDMWallTime() - tMark;
}


// ************************************************************************* //
