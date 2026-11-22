import io, sys
F = 'cdmAmrBridge.C'
s = io.open(F).read()
if '1 << 29' in s:
    sys.exit('ALREADY APPLIED')

OLD1 = '''    // 21 bits per field.  The cell count here already reaches 1.5M
    // against a limit of 2097151, so this is checked, not trusted.
    if (nCells >= (1 << 21))
    {
        FatalErrorInFunction
            << nCells << " cells does not fit the 21-bit sort key"
            << abort(FatalError);
    }'''

NEW1 = '''    // The key is (group, a, b) packed into one 64-bit integer.  The
    // group is a patch index plus one, so six bits is ample; that leaves
    // 29 bits for each of the two cell-index fields, or 536 million
    // cells.  The original 22/21/21 split capped the mesh at 2,097,151
    // and a production run reached it - 2,097,867 cells - and aborted.
    //
    // The "a" field is a cell index for an internal face and a seam
    // coordinate for a patch face; the latter is bounded by ny*2^maxLevel
    // (46,080 on this case), so the cell bound covers both.
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
    }'''

OLD2 = '            FaceSort{(grp << 42) | (a << 21) | b2, f};'
NEW2 = '            FaceSort{(grp << 58) | (a << 29) | b2, f};'

for old, new in ((OLD1, NEW1), (OLD2, NEW2)):
    if s.count(old) != 1:
        sys.exit('FAIL: anchor x%d: %s' % (s.count(old), old.strip()[:50]))
    s = s.replace(old, new, 1)

io.open(F, 'w').write(s)
print('ok: cap raised to 536,870,911 cells')
