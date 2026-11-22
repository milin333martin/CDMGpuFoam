#!/usr/bin/env python3
# Adds halo (MPI ghost-cell) support to src/gpu.  Run from src/gpu with
# halo-payload.txt beside it.  Safe to re-run: it refuses if already applied.
import io, sys, os

PAY = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'halo-payload.txt')
if not os.path.exists(PAY):
    PAY = 'halo-payload.txt'

def load_blocks(path):
    blocks, cur = {}, None
    for ln in io.open(path, encoding='utf-8').read().split('\n'):
        if ln.startswith('### '):
            cur = ln[4:].strip()
            if cur != 'END':
                blocks[cur] = []
            continue
        if cur and cur != 'END':
            blocks[cur].append(ln)
    return {k: '\n'.join(v) for k, v in blocks.items()}

B = load_blocks(PAY)

OLD_STRUCT = """    // mesh
    int nCells = 0, nIF = 0, nBF = 0, nBC = 0;
    std::vector<int> bCell;

    //- Contiguous staging for the boundary-cell exchange, so the host never
    //  touches a strided subset of a device-resident field
    std::vector<double> bcsRho, bcsU, bcsP, bcsT, bcsY;"""

OLD_ALLOC = """        rho.assign(nCells, 0.0);  rhoU.assign(3*nCells, 0.0);
        rhoE.assign(nCells, 0.0); rhoY.assign(nCells, 0.0);
        rho0 = rho; rhoU0 = rhoU; rhoE0 = rhoE; rhoY0 = rhoY;

        U.assign(3*nCells, 0.0);  p.assign(nCells, 0.0);
        T.assign(nCells, 0.0);    Y.assign(nCells, 0.0);
        c.assign(nCells, 0.0);    Qdot.assign(nCells, 0.0);

        gRho.assign(3*nCells, 0.0); gU.assign(9*nCells, 0.0);
        gP.assign(3*nCells, 0.0);   gY.assign(3*nCells, 0.0);"""

OLD_SIG = """    cdmGpuCase* c,
    int nCells, int nIF, int nBF,
    const int* owner, const int* neighbour,"""
NEW_SIG = """    cdmGpuCase* c,
    int nCells, int nHalo, int nIF, int nBF,
    const int* owner, const int* neighbour,"""

OLD_BODY = "    c->nCells = nCells; c->nIF = nIF; c->nBF = nBF;"
NEW_BODY = "    c->nCells = nCells; c->nHalo = nHalo; c->nIF = nIF; c->nBF = nBF;"

OLD_SEAM = """    std::copy(bY,   bY   + c->nBF,   c->bY.begin());
}
"""

OLD_HSIG = """    cdmGpuCase* c,
    int nCells,
    int nInternalFaces,
    int nBoundaryFaces,"""

OLD_HRHS = "/*- Evaluate the right-hand side from the current device state.  Returns the"

OLD_ZND = "        c, N, nIF, nBF, own.data(), nei.data(),"
NEW_ZND = """        c, N, 0 /* nHalo: serial, no processor faces */, nIF, nBF,
        own.data(), nei.data(),"""

EDITS = [
    ('cdmGpuKernels.cpp', OLD_STRUCT, B['STRUCT_NEW']),
    ('cdmGpuKernels.cpp', OLD_ALLOC,  B['ALLOCATE_NEW']),
    ('cdmGpuKernels.cpp', OLD_SIG,    NEW_SIG),
    ('cdmGpuKernels.cpp', OLD_BODY,   NEW_BODY),
    ('cdmGpuKernels.cpp', OLD_SEAM,   OLD_SEAM + '\n\n' + B['KERNELS_HALO_BLOCK'] + '\n'),
    ('cdmGpu.h',          OLD_HSIG,   B['HEADER_SIG_NEW']),
    ('cdmGpu.h',          OLD_HRHS,   B['HEADER_DECLS'] + '\n' + OLD_HRHS),
    ('test/testZND.cpp',  OLD_ZND,    NEW_ZND),
]

# Two of the edits above are insertions whose anchor survives them, so a
# second run would insert the block twice.  This is the stop.
GUARD = [('cdmGpuKernels.cpp', 'cdmGpuComputeGradients'),
         ('cdmGpu.h',          'cdmGpuSetHaloCells')]


def main():
    for path, marker in GUARD:
        if marker in io.open(path, encoding='utf-8').read():
            sys.exit("ALREADY APPLIED: %s already defines %s.  Nothing done."
                     % (path, marker))

    bad = 0
    for path, old, new in EDITS:
        s = io.open(path, encoding='utf-8').read()
        n = s.count(old)
        if n != 1:
            print("FAIL  %-22s anchor seen %d times: %s..."
                  % (path, n, old.split('\n')[0][:46]))
            bad += 1
            continue
        io.open(path, 'w', encoding='utf-8').write(s.replace(old, new, 1))
        print("ok    %-22s %s..." % (path, old.split('\n')[0][:46]))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
