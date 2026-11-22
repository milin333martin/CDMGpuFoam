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

Description
    Refinement and coarsening in the same step.

    Both flags are indexed by the cell array and refine() shifts every
    index after the first refined cell, so a coarsening flag built before
    the refinement coarsens the wrong cells.  It reads in bounds, so
    nothing warns: the only symptom is that the mesh refines about half as
    much as it was asked to, and every validity check still passes because
    none of them depends on which cells were chosen.

    This test states the property exactly: a cell of the old mesh that is
    gone, and whose parent is present, was coarsened, and it must have
    been flagged for coarsening.  Nothing else can make that true by
    accident.

\*---------------------------------------------------------------------------*/

#include "cdmAmr.h"
#include <cstdio>
#include <cmath>
#include <set>
using namespace cdmAmr;

int main()
{
    Mesh m; m.initialise(40, 24, 0.4, 0.24, true); m.buildAddressing();

    for (int step = 0; step < 12; ++step)
    {
        const double xc = 0.4*(0.15 + 0.7*step/11.0);
        const int n = m.nCells();
        std::vector<char> ref(n,0), crs(n,0);
        for (int c = 0; c < n; ++c)
        {
            const bool in = std::fabs(m.C[2*c]-xc) < 0.03;
            const int L = (int)keyLevel(m.keys[c]);
            if (in && L < 3) ref[c] = 1;
            if (!in && L > 0) crs[c] = 1;
        }

        // the keys of the cells we INTEND to coarsen, and their parents
        std::set<Key> wantParent;
        for (int c = 0; c < n; ++c) if (crs[c] && !ref[c])
        {
            const Key k = m.keys[c];
            wantParent.insert(makeKey(keyBase(k), keyLevel(k)-1, keyMorton(k)>>2));
        }

        std::vector<Key> prev = m.keys;
        std::vector<int> ni;
        m.refine(ref, &ni);

        // every old cell must be accounted for, in order
        for (std::size_t c = 1; c < ni.size(); ++c)
            if (ni[c] <= ni[c-1]) { printf("FAIL newIndex not increasing at %zu\n", c); return 1; }

        std::vector<char> crs2(m.nCells(), 0);
        for (std::size_t c = 0; c < ni.size(); ++c)
        {
            const int w = ni[c];
            if (w >= 0 && w < m.nCells() && m.keys[w] == prev[c]) crs2[w] = crs[c];
        }
        // the carried flag must sit on the same key it was built for
        for (std::size_t c = 0; c < ni.size(); ++c)
        {
            const int w = ni[c];
            if (m.keys[w] == prev[c] && crs2[w] != crs[c])
            { printf("FAIL flag lost at %zu\n", c); return 1; }
        }
        m.unrefine(crs2);
        m.balance();
        m.buildAddressing();

        // Exact statement of correctness: a cell of the old mesh that is gone and
        // whose parent is present was coarsened, and must have been flagged for
        // coarsening. Anything else means the flag landed on the wrong cell.
        std::set<Key> after(m.keys.begin(), m.keys.end());
        int wrong = 0, coarsened = 0;
        for (std::size_t c = 0; c < prev.size(); ++c)
        {
            const Key k = prev[c];
            if (keyLevel(k) == 0 || after.count(k)) continue;
            const Key par = makeKey(keyBase(k), keyLevel(k)-1, keyMorton(k)>>2);
            if (!after.count(par)) continue;        // refined away, not coarsened
            ++coarsened;
            if (!crs[c]) ++wrong;
        }

        const std::string e1 = checkKeys(m), e2 = checkMesh(m);
        printf("step %2d  cells %7d  coarsened %5d  wrongly %d  %s%s\n",
               step, m.nCells(), coarsened, wrong,
               e1.empty()?"keys ok":e1.c_str(), e2.empty()?" mesh ok":e2.c_str());
        if (!e1.empty() || !e2.empty() || wrong) return 1;
    }
    printf("OK\n");
    return 0;
}
