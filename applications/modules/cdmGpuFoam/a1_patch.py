import io, sys
C = 'cdmGpuFoam.C'
H = 'cdmGpuFoam.H'

E = []

E.append((H,
'''        labelList cycPatch_;
        labelList cycFace_;''',
'''        labelList cycPatch_;
        labelList cycFace_;

        //- Processor faces, promoted to internal faces on the device in the
        //  same way, except that the neighbour is a ghost cell rather than
        //  a local one.  The device face index is mesh.nInternalFaces()
        //  + nCyc + the position here, and the ghost cell it reads is
        //  nCells + that position.
        //
        //  A processor face exists once physically and twice in the
        //  decomposition, so unlike a cyclic face only this rank's side is
        //  mapped back into phi_.  The two ranks still agree on the flux to
        //  the last bit: both evaluate the same Riemann problem from the
        //  same reconstructed states, and HLLC satisfies
        //  F(L,R,n) = -F(R,L,-n).
        labelList procPatch_;
        labelList procFace_;

        //- The local cell behind each processor face, in that same order:
        //  what this rank sends its neighbour every Runge-Kutta stage.
        labelList haloCell_;

        //- Where each processor patch's faces begin in the lists above,
        //  with a trailing total, and the rank on the other side of each.
        labelList procStart_;
        labelList procNbr_;

        //- Halo exchange buffers, sized once per topology change.  Twenty
        //  five doubles a face: rho, U, p, T, Y, and the gradients of rho,
        //  U, p and Y.  That is the entire stencil a face reads for its
        //  neighbour, which is why one layer of halo suffices where
        //  exchanging cell values alone would need two.
        scalarField hSend_;
        scalarField hRecv_;

        //- Number of processor faces, and so of ghost cells
        label nHalo_;'''))

E.append((H,
'''        //- Keep phi_ consistent with the device's face mass flux, for the
        //  flux-dependent boundary conditions
        void downloadMassFlux();''',
'''        //- Fill the device's ghost cells from the neighbouring ranks.
        //
        //  Runs the cell-gradient pass, brings back the state and gradient
        //  of the cells behind each processor patch, swaps them with the
        //  rank on the other side, and uploads the result into the ghosts.
        //  Must happen before the fluxes are formed, which is why the
        //  gradient pass is a separate entry point: cdmGpuComputeRHS then
        //  recomputes those gradients but writes only real cells, so the
        //  ghosts it reads are the exchanged ones.
        //
        //  A no-op in serial, where there are no processor patches.
        void exchangeHalo();

        //- Keep phi_ consistent with the device's face mass flux, for the
        //  flux-dependent boundary conditions
        void downloadMassFlux();'''))

E.append((C,
'#include "cyclicFvPatch.H"',
'''#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "PstreamBuffers.H"'''))

E.append((C,
'''            if (!isA<cyclicFvPatch>(fvp))
            {
                FatalErrorInFunction
                    << "patch " << fvp.name() << " is coupled but is a "
                    << mesh.boundaryMesh()[patchi].type()
                    << ", which cdmGpuFoam does not handle." << nl
                    << "Only cyclic coupling is supported; a processor patch"
                    << " cannot occur because this solver is serial."
                    << exit(FatalError);
            }''',
'''            // Processor patches are collected separately, below: they are
            // promoted to internal faces too, but with a ghost cell as the
            // neighbour rather than a local one.
            if (isA<processorFvPatch>(fvp)) { continue; }

            if (!isA<cyclicFvPatch>(fvp))
            {
                FatalErrorInFunction
                    << "patch " << fvp.name() << " is coupled but is a "
                    << mesh.boundaryMesh()[patchi].type()
                    << ", which cdmGpuFoam does not handle." << nl
                    << "Only cyclic and processor coupling are supported."
                    << exit(FatalError);
            }'''))

E.append((C,
'''    const label nCyc = cycPatch.size();
    const label nIFdev = nIF + nCyc;''',
'''    const label nCyc = cycPatch.size();

    // ---- processor patches, in patch order ------------------------------
    // Every face of every processor patch becomes a device internal face
    // with a ghost neighbour.  OpenFOAM orders the faces of a processor
    // patch so that face i here is face i on the other rank, which is what
    // lets the exchange be a plain contiguous swap with no index traffic.
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

    const label nIFdev = nIF + nCyc + nHalo_;'''))

E.append((C,
'    // fvPatch::delta() is documented as the owner-to-neighbour vector for a',
'''    // ---- processor faces, appended as internal faces with a ghost --------
    // Identical to the cyclic case but for the neighbour: the cell on the
    // other side belongs to another rank, so the device reads a ghost that
    // exchangeHalo() fills each Runge-Kutta stage.  Patch-major for the
    // same reason - delta() builds the whole patch's vectors at once.
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

            // Only this rank's face maps onto the device face.  The other
            // side of it is in another rank's mesh, which is the single
            // way this differs from a cyclic pair.
            faceToC[start + f - nIF] = d;
            faceToCSign[start + f - nIF] = 1;
        }
    }

    // fvPatch::delta() is documented as the owner-to-neighbour vector for a'''))

E.append((C,
'''    if (nCyc > 0)
    {
        label nBadSign = 0, nBadLength = 0;
        scalar worstRatio = 1;

        for (label k = 0; k < nCyc; k++)
        {
            const label d = nIF + k;''',
'''    // Checked over the cyclic and processor faces alike: both are coupled
    // patches whose delta() this code trusts for the reconstruction.
    const label nAppended = nCyc + nHalo_;

    if (nAppended > 0)
    {
        label nBadSign = 0, nBadLength = 0;
        scalar worstRatio = 1;

        for (label k = 0; k < nAppended; k++)
        {
            const label d = nIF + k;'''))

E.append((C,
'''            FatalErrorInFunction
                << "the cyclic face geometry does not match the convention"
                << " this solver assumes." << nl
                << "    faces with delta opposed to Sf: " << nBadSign
                << " of " << nCyc << nl
                << "    faces whose |delta.n| disagrees with 1/deltaCoeffs: "
                << nBadLength << " of " << nCyc''',
'''            FatalErrorInFunction
                << "the coupled face geometry does not match the convention"
                << " this solver assumes (" << nCyc << " cyclic, "
                << nHalo_ << " processor)." << nl
                << "    faces with delta opposed to Sf: " << nBadSign
                << " of " << nAppended << nl
                << "    faces whose |delta.n| disagrees with 1/deltaCoeffs: "
                << nBadLength << " of " << nAppended'''))

E.append((C,
'        gpu_, nCells, 0 /* nHalo: serial until MPI lands */, nIFdev, nBF,',
'        gpu_, nCells, nHalo_, nIFdev, nBF,'))

E.append((C,
'''    cycPatch_ = cycPatch;
    cycFace_ = cycFace;
    nIFdev_ = nIFdev;
    phiScratch_.setSize(nIFdev);''',
'''    cycPatch_ = cycPatch;
    cycFace_ = cycFace;
    nIFdev_ = nIFdev;
    phiScratch_.setSize(nIFdev);

    // The cells this rank sends, and room for the twenty-five doubles a
    // face each way.  Sized here rather than per stage: with AMR the halo
    // changes size at every topology change and nowhere else.
    if (nHalo_ > 0)
    {
        cdmGpuSetHaloCells(gpu_, nHalo_, haloCell_.begin());
        hSend_.setSize(25*nHalo_);
        hRecv_.setSize(25*nHalo_);
    }'''))

E.append((C,
'''        phi_.boundaryFieldRef()[patchi][f] = flux;
        phi_.boundaryFieldRef()[nbri][f] = -flux;
    }
}''',
'''        phi_.boundaryFieldRef()[patchi][f] = flux;
        phi_.boundaryFieldRef()[nbri][f] = -flux;
    }

    // A processor face has only this rank's side here.  The other side is
    // written by the other rank from its own device, and the two agree to
    // the last bit because both solved the same Riemann problem from the
    // same reconstructed states.
    forAll(procPatch_, h)
    {
        phi_.boundaryFieldRef()[procPatch_[h]][procFace_[h]] =
            phiScratch_[nIF + cycPatch_.size() + h];
    }
}'''))

E.append((C,
'''    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "cdmGpuFoam runs in serial only: one process, one device, the"
            << " whole mesh." << nl
            << "Run foamRun without -parallel, on an undecomposed case."
            << exit(FatalError);
    }''',
'''    if (Pstream::parRun())
    {
        Info<< "cdmGpuFoam: parallel run, " << Pstream::nProcs()
            << " ranks.  One device per rank; processor faces are carried"
            << " as device internal faces with ghost neighbours." << endl;
    }'''))

E.append((C, '    nIFdev_(0),', '    nHalo_(0),\n    nIFdev_(0),'))

E.append((C,
'''    // ---- the right-hand side ---------------------------------------------
    int nf = 0;
    cdmGpuComputeRHS(gpu_, &nf);''',
'''    // ---- the right-hand side ---------------------------------------------
    // The halo first: the fluxes about to be formed read a ghost state and
    // a ghost gradient at every processor face.  A no-op in serial.
    exchangeHalo();

    int nf = 0;
    cdmGpuComputeRHS(gpu_, &nf);'''))

E.append((C,
'''        tMark = CDMWallTime();

        int nf = 0;
        cdmGpuComputeRHS(gpu_, &nf);''',
'''        tMark = CDMWallTime();

        // Every stage, because the ghosts must carry the state the
        // neighbour holds *now*, not at the start of the step.
        exchangeHalo();

        int nf = 0;
        cdmGpuComputeRHS(gpu_, &nf);'''))

E.append((C,
'''void Foam::solvers::cdmGpuFoam::downloadMassFlux()
{''',
'''void Foam::solvers::cdmGpuFoam::exchangeHalo()
{
    if (nHalo_ == 0) { return; }

    const label n = nHalo_;

    // Where each quantity sits in the 25-doubles-a-face buffer, and how
    // many doubles it occupies per face.
    static const label off[9] = {0, 1, 4, 5, 6, 7, 10, 19, 22};
    static const label str[9] = {1, 3, 1, 1, 1, 3,  9,  3,  3};

    // The gradients have to exist before they can be sent, and they are
    // computed over real cells only, so this pass leaves the ghosts alone.
    cdmGpuComputeGradients(gpu_);

    {
        scalar* s = hSend_.begin();
        cdmGpuDownloadHaloCells
        (
            gpu_,
            s,            s +    n,  s +  4*n,  s +  5*n,  s +  6*n,
            s + 7*n,      s + 10*n,  s + 19*n,  s + 22*n
        );
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;

        List<scalar> buf(25*np);
        label q = 0;

        for (label blk = 0; blk < 9; blk++)
        {
            const scalar* src = hSend_.begin() + off[blk]*n + str[blk]*a;
            const label m = str[blk]*np;
            for (label i = 0; i < m; i++) { buf[q++] = src[i]; }
        }

        UOPstream toNbr(procNbr_[p], pBufs);
        toNbr << buf;
    }

    pBufs.finishedSends();

    forAll(procNbr_, p)
    {
        const label a = procStart_[p];
        const label np = procStart_[p + 1] - a;

        UIPstream fromNbr(procNbr_[p], pBufs);
        List<scalar> buf(fromNbr);

        if (buf.size() != 25*np)
        {
            FatalErrorInFunction
                << "halo exchange with rank " << procNbr_[p] << " returned "
                << buf.size() << " values for " << np << " faces; expected "
                << 25*np << "." << nl
                << "  The two ranks disagree on how many faces their shared"
                << " processor patch has, which means the meshes were built"
                << " independently rather than from one decomposition."
                << exit(FatalError);
        }

        label q = 0;

        for (label blk = 0; blk < 9; blk++)
        {
            scalar* dst = hRecv_.begin() + off[blk]*n + str[blk]*a;
            const label m = str[blk]*np;
            for (label i = 0; i < m; i++) { dst[i] = buf[q++]; }
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


void Foam::solvers::cdmGpuFoam::downloadMassFlux()
{'''))


def main():
    if 'exchangeHalo' in io.open(C).read():
        sys.exit('ALREADY APPLIED')
    text = {p: io.open(p).read() for p in (C, H)}
    bad = 0
    for i, (path, old, new) in enumerate(E):
        n = text[path].count(old)
        if n != 1:
            print('FAIL  edit %-2d  %-14s anchor x%d: %s'
                  % (i, path, n, old.strip().split('\n')[0][:48]))
            bad += 1
    if bad:
        print('\n%d of %d anchors did not match exactly once. NOTHING written.'
              % (bad, len(E)))
        return 1
    for path, old, new in E:
        text[path] = text[path].replace(old, new, 1)
    for path in (C, H):
        io.open(path, 'w').write(text[path])
    print('applied %d edits to cdmGpuFoam.C and cdmGpuFoam.H' % len(E))
    return 0


if __name__ == '__main__':
    sys.exit(main())
