import io, sys
F = 'cdmGpuFoam.C'
s = io.open(F).read()
if 'badProc' in s:
    sys.exit('ALREADY APPLIED')

OLD = '''        label badInternal = 0, badBoundary = 0;

        for (label f = 0; f < nIFdev; f++)
        {
            if (nSeen[f] != 2 || signSum[f] != 0) { badInternal++; }
        }
        for (label b = nIFdev; b < nIFdev + nBF; b++)
        {
            if (nSeen[b] != 1 || signSum[b] != 1) { badBoundary++; }
        }

        if (badInternal > 0 || badBoundary > 0)
        {
            FatalErrorInFunction
                << "the cell-to-face map is inconsistent." << nl
                << "    internal faces not seen exactly once from each side: "
                << badInternal << " of " << nIFdev << nl
                << "    boundary faces not seen exactly once as owner:      "
                << badBoundary << " of " << nBF << nl'''

NEW = '''        // A processor face is the exception, deliberately: its neighbour
        // is a ghost, and a ghost owns no faces, so this rank sees it
        // exactly once, from the owner.  The other side is counted on the
        // other rank, which is what makes the pair conservative across the
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
                << badBoundary << " of " << nBF << nl'''

if s.count(OLD) != 1:
    sys.exit('FAIL: anchor seen %d times' % s.count(OLD))
io.open(F, 'w').write(s.replace(OLD, NEW, 1))
print('ok: processor faces exempted from the twice-seen invariant')
