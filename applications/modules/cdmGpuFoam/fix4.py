import io, sys
C, H = 'cdmGpuFoam.C', 'cdmGpuFoam.H'
s = io.open(C).read()
if 'haloRound' in s: sys.exit('ALREADY APPLIED')

BODY = r'''void Foam::solvers::cdmGpuFoam::haloRound()
{
    const label n = nHalo_;

    // Where each quantity sits in the 25-doubles-a-face buffer, and how
    // many doubles it occupies per face.
    static const label off[9] = {0, 1, 4, 5, 6, 7, 10, 19, 22};
    static const label str[9] = {1, 3, 1, 1, 1, 3,  9,  3,  3};

    // A tag of our own.  PstreamBuffers was the obvious choice and was
    // wrong: its traffic collided with the processor-patch exchange that
    // GeometricBoundaryField::evaluate() performs later in the same step.
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

    // Pack patch by patch, so each patch's payload is contiguous on the
    // wire.  Both ranks pack the same blocks in the same order over faces
    // OpenFOAM has already matched, so the unpack is symmetric.
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

    // Two rounds, with the gradient pass between them.  The reason is
    // subtle enough to be worth writing down.
    //
    // A cell next to a processor face builds its gradient from its own
    // faces, one of which reads a ghost.  With a single round, the
    // gradient this rank *sends* was computed from last stage's ghost -
    // stale by one stage - and the neighbour then uses it for the far side
    // of its Riemann problem.  This rank's own gradients are fine, because
    // cdmGpuComputeRHS recomputes them after the upload.  It is only the
    // ones we send that are wrong, which is what makes the error invisible
    // from here and visible only as a disagreement with serial.
    //
    // Measured on RDE2D_T6 at 135k cells: with a seam crossing the
    // ignition region, one round disagreed with serial by 5e-9 relative
    // after six steps, while a seam in quiescent gas was bit-identical -
    // there the stale ghost happened to equal the fresh one, which is
    // exactly the kind of accident that hides a bug in the easy test.
    //
    // So: round one makes the ghost states current, the gradient pass
    // recomputes this rank's gradients from them, and round two sends
    // those.  Nothing beyond that is needed, since a gradient depends only
    // on cell values and round one has already made those exact.
    //
    // The cost is two messages a stage rather than one.  Round one only
    // needs the five state blocks and round two only the four gradient
    // blocks, so this could send 7 and 18 doubles a face instead of 25
    // twice; worth doing if the halo ever shows up in a profile.
    haloRound();
    cdmGpuComputeGradients(gpu_);
    haloRound();
}


'''

a = s.index('void Foam::solvers::cdmGpuFoam::exchangeHalo()')
b = s.index('void Foam::solvers::cdmGpuFoam::downloadMassFlux()')
s = s[:a] + BODY + s[b:]

old = '''    if (getenv("CDM_NO_QDOT_BC") == nullptr)
    {
        Qdot_.correctBoundaryConditions();
    }'''
if s.count(old) == 1:
    s = s.replace(old, '    Qdot_.correctBoundaryConditions();', 1)
io.open(C, 'w').write(s)

h = io.open(H).read()
old2 = '        void exchangeHalo();'
new2 = '''        //- One download / swap / upload of the ghost payload
        void haloRound();

        void exchangeHalo();'''
if h.count(old2) != 1: sys.exit('header anchor x%d' % h.count(old2))
io.open(H, 'w').write(h.replace(old2, new2, 1))

print('ok: two-round exchange, diagnostics removed')
print('leftover switches:', sum(s.count(k) for k in
      ('CDM_NO_HALO','CDM_HALO_NOMPI','CDM_HALO_STAGE','CDM_FORCE_GRADS','CDM_NO_QDOT_BC')))
