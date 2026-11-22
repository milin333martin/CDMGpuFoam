import io, sys
C, H = 'cdmGpuFoam.C', 'cdmGpuFoam.H'
if 'haloTag' in io.open(C).read(): sys.exit('ALREADY APPLIED')

E = [
 (H, '''        scalarField hSend_;
        scalarField hRecv_;''',
     '''        scalarField hSend_;
        scalarField hRecv_;

        //- The same payload laid out patch by patch, which is what goes on
        //  the wire: a patch's slice is then contiguous and needs no
        //  per-exchange allocation.
        scalarField hPack_;
        scalarField hUnpack_;'''),

 (C, '''#include "processorFvPatch.H"
#include "PstreamBuffers.H"''',
     '''#include "processorFvPatch.H"
#include "UIPstream.H"
#include "UOPstream.H"'''),

 (C, '''        hSend_.setSize(25*nHalo_);
        hRecv_.setSize(25*nHalo_);''',
     '''        hSend_.setSize(25*nHalo_);
        hRecv_.setSize(25*nHalo_);
        hPack_.setSize(25*nHalo_);
        hUnpack_.setSize(25*nHalo_);'''),

 (C, '''    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

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
    }''',
     '''    // Explicit non-blocking pairs on a tag of our own.
    //
    // PstreamBuffers was the obvious choice and was wrong.  Its traffic
    // collided with the processor-patch exchange that
    // GeometricBoundaryField::evaluate() performs later in the same step;
    // evaluate() then read this payload - twenty-five doubles a face -
    // into a field sized for three components, overran it, and died in a
    // libc copy with a corrupted pointer, in a place with no visible
    // connection to the halo.  A dedicated tag cannot collide however the
    // two interleave.
    static const int haloTag = 1729;

    // Pack patch by patch, so each patch's payload is contiguous on the
    // wire.  Both ranks pack the same blocks in the same order over faces
    // that OpenFOAM has already matched, so the unpack is symmetric.
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
    }'''),
]

t = {p: io.open(p).read() for p in (C, H)}
bad = [i for i, (p, o, nw) in enumerate(E) if t[p].count(o) != 1]
if bad: sys.exit('FAIL anchors %s - nothing written' % bad)
for p, o, nw in E: t[p] = t[p].replace(o, nw, 1)
for p in (C, H): io.open(p, 'w').write(t[p])
print('ok: PstreamBuffers replaced with tagged non-blocking pairs')
