import io, sys
F = 'cdmGpuKernels.cpp'
s = io.open(F).read()
if 'Counts, never end()' in s: sys.exit('ALREADY APPLIED')

OLD1 = '''    std::copy(c->rho.begin(),  c->rho.end(),  rho);
    std::copy(c->rhoU.begin(), c->rhoU.end(), rhoU);
    std::copy(c->rhoE.begin(), c->rhoE.end(), rhoE);
    std::copy(c->rhoY.begin(), c->rhoY.end(), rhoY);'''
NEW1 = '''    // Counts, never end().  rho is nCells + nHalo long because a ghost
    // needs a density; the host field it is being copied into is nCells
    // long and always will be.  Writing the halo into it overruns the heap
    // by nHalo doubles - invisible in serial, where nHalo is zero, and
    // fatal in parallel somewhere else entirely.
    const int n = c->nCells;

    std::copy(c->rho.begin(),  c->rho.begin()  +   n, rho);
    std::copy(c->rhoU.begin(), c->rhoU.begin() + 3*n, rhoU);
    std::copy(c->rhoE.begin(), c->rhoE.begin() +   n, rhoE);
    std::copy(c->rhoY.begin(), c->rhoY.begin() +   n, rhoY);'''

OLD2 = '''    std::copy(c->U.begin(), c->U.end(), U);
    std::copy(c->p.begin(), c->p.end(), p);
    std::copy(c->T.begin(), c->T.end(), T);
    std::copy(c->Y.begin(), c->Y.end(), Y);
    std::copy(c->Qdot.begin(), c->Qdot.end(), Qdot);'''
NEW2 = '''    // Counts, never end(): U, p, T and Y all carry nHalo ghosts past the
    // real cells, and the host fields do not.  See cdmGpuDownloadState.
    const int n = c->nCells;

    std::copy(c->U.begin(), c->U.begin() + 3*n, U);
    std::copy(c->p.begin(), c->p.begin() +   n, p);
    std::copy(c->T.begin(), c->T.begin() +   n, T);
    std::copy(c->Y.begin(), c->Y.begin() +   n, Y);
    std::copy(c->Qdot.begin(), c->Qdot.begin() + n, Qdot);'''

for o, nw in ((OLD1, NEW1), (OLD2, NEW2)):
    if s.count(o) != 1: sys.exit('FAIL anchor x%d' % s.count(o))
    s = s.replace(o, nw, 1)
io.open(F, 'w').write(s)
print('ok: five heap overruns fixed')
