/*---------------------------------------------------------------------------*\
  Validation and timing for the Morton-key adaptive mesh engine.

  The engine's addressing is checked against a brute-force reference that
  knows nothing about keys: it compares every pair of cell boxes and calls
  them neighbours when they share an edge of positive length.  That is
  O(n^2) and only usable on small meshes, which is the point - it is an
  independent statement of what the answer should be, not a second copy of
  the same algorithm.
\*---------------------------------------------------------------------------*/

#include "cdmAmr.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace cdmAmr;

namespace
{

int nFail = 0;

void check(bool ok, const char* what)
{
    std::printf("  %-54s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) { ++nFail; }
}

void checkEmpty(const std::string& err, const char* what)
{
    if (err.empty())
    {
        std::printf("  %-54s ok\n", what);
    }
    else
    {
        std::printf("  %-54s FAILED: %s\n", what, err.c_str());
        ++nFail;
    }
}

double now()
{
    using clk = std::chrono::steady_clock;
    static const clk::time_point t0 = clk::now();
    return std::chrono::duration<double>(clk::now() - t0).count();
}


//- Every neighbouring pair, found by comparing boxes and nothing else
std::set<std::pair<int, int>> bruteForceFaces(const Mesh& m)
{
    const int n = m.nCells();
    std::vector<double> x0(n), y0(n), x1(n), y1(n);
    for (int c = 0; c < n; ++c) { m.cellBox(c, x0[c], y0[c], x1[c], y1[c]); }

    const double tol = 1e-12*std::max(m.Lx, m.Ly);
    std::set<std::pair<int, int>> faces;

    for (int a = 0; a < n; ++a)
    {
        for (int b = a + 1; b < n; ++b)
        {
            const double yo =
                std::min(y1[a], y1[b]) - std::max(y0[a], y0[b]);
            const double xo =
                std::min(x1[a], x1[b]) - std::max(x0[a], x0[b]);

            bool touch = false;

            if (yo > tol)
            {
                if (std::fabs(x1[a] - x0[b]) < tol) { touch = true; }
                if (std::fabs(x1[b] - x0[a]) < tol) { touch = true; }

                if (m.cyclicX)
                {
                    if (std::fabs(x1[a] - (x0[b] + m.Lx)) < tol) { touch = true; }
                    if (std::fabs(x1[b] - (x0[a] + m.Lx)) < tol) { touch = true; }
                }
            }

            if (xo > tol)
            {
                if (std::fabs(y1[a] - y0[b]) < tol) { touch = true; }
                if (std::fabs(y1[b] - y0[a]) < tol) { touch = true; }
            }

            if (touch) { faces.insert({a, b}); }
        }
    }

    return faces;
}


void compareWithBruteForce(Mesh& m, const char* what)
{
    m.buildAddressing();

    std::set<std::pair<int, int>> mine;
    for (int f = 0; f < m.nFaces(); ++f)
    {
        mine.insert({m.owner[f], m.neighbour[f]});
    }

    const std::set<std::pair<int, int>> ref = bruteForceFaces(m);

    // A coarse cell facing two finer ones gives two faces between the same
    // pair of cells, so the engine may hold duplicates that the set of pairs
    // collapses. Compare the pair sets, then account for the duplicates.
    const bool same = (mine == ref);

    std::printf("  %-54s %s  (%d faces, %zu pairs, %zu reference)\n",
        what, same ? "ok" : "FAILED",
        m.nFaces(), mine.size(), ref.size());

    if (!same)
    {
        ++nFail;
        int shown = 0;
        for (const auto& p : ref)
        {
            if (!mine.count(p) && shown < 5)
            {
                std::printf("      missing pair (%d,%d)\n", p.first, p.second);
                ++shown;
            }
        }
        for (const auto& p : mine)
        {
            if (!ref.count(p) && shown < 10)
            {
                std::printf("      spurious pair (%d,%d)\n", p.first, p.second);
                ++shown;
            }
        }
    }
}

} // End anonymous namespace


int main(int argc, char** argv)
{
    bool timingOnly = false;
    for (int a = 1; a < argc; ++a)
    {
        if (!std::strcmp(argv[a], "--timing")) { timingOnly = true; }
    }

    std::printf("Morton-key adaptive mesh engine\n\n");

    if (!timingOnly)
    {
        // ---- a uniform base mesh ------------------------------------------
        std::printf("uniform base mesh, 6 x 4, no refinement\n");
        {
            Mesh m;
            m.initialise(6, 4, 0.6, 0.4, false);
            checkEmpty(checkKeys(m), "key set valid");
            m.buildAddressing();
            checkEmpty(checkMesh(m), "mesh valid");
            check(m.nCells() == 24, "24 cells");
            // 6x4 structured: 5*4 vertical + 6*3 horizontal = 38
            check(m.nFaces() == 38, "38 internal faces");
            check(m.nBoundaryFaces() == 20, "20 boundary faces");
            compareWithBruteForce(m, "addressing matches brute force");
        }

        // ---- one cell refined, forcing hanging faces ----------------------
        std::printf("\none cell refined\n");
        {
            Mesh m;
            m.initialise(6, 4, 0.6, 0.4, false);
            std::vector<char> f(m.nCells(), 0);
            f[9] = 1;
            check(m.refine(f) == 1, "one cell refined");
            check(m.nCells() == 27, "27 cells");
            checkEmpty(checkKeys(m), "key set valid");
            m.buildAddressing();
            checkEmpty(checkMesh(m), "mesh valid");
            compareWithBruteForce(m, "addressing matches brute force");
        }

        // ---- deep refinement needing 2:1 balance --------------------------
        std::printf("\nnested refinement with 2:1 balance\n");
        {
            Mesh m;
            m.initialise(8, 8, 0.8, 0.8, false);

            for (int pass = 0; pass < 3; ++pass)
            {
                std::vector<char> f(m.nCells(), 0);
                for (int c = 0; c < m.nCells(); ++c)
                {
                    double x0, y0, x1, y1;
                    m.cellBox(c, x0, y0, x1, y1);
                    const double cx = 0.5*(x0 + x1), cy = 0.5*(y0 + y1);
                    if (std::hypot(cx - 0.4, cy - 0.4) < 0.12) { f[c] = 1; }
                }
                m.refine(f);
            }

            const int added = m.balance();
            std::printf("    %d cells, %d added by balancing\n",
                m.nCells(), added);

            checkEmpty(checkKeys(m), "key set valid");
            m.buildAddressing();
            checkEmpty(checkMesh(m), "mesh valid (includes 2:1 check)");
            compareWithBruteForce(m, "addressing matches brute force");
        }

        // ---- cyclic in x --------------------------------------------------
        std::printf("\ncyclic in x, refined against the seam\n");
        {
            Mesh m;
            m.initialise(6, 4, 0.6, 0.4, true);

            std::vector<char> f(m.nCells(), 0);
            for (int c = 0; c < m.nCells(); ++c)
            {
                double x0, y0, x1, y1;
                m.cellBox(c, x0, y0, x1, y1);
                if (x0 < 1e-12 || x1 > 0.6 - 1e-12) { f[c] = 1; }
            }
            m.refine(f);
            m.balance();

            checkEmpty(checkKeys(m), "key set valid");
            m.buildAddressing();
            checkEmpty(checkMesh(m), "mesh valid");
            compareWithBruteForce(m, "addressing matches brute force");

            // The seam must carry faces, or it is not cyclic at all
            int seam = 0;
            for (int ff = 0; ff < m.nFaces(); ++ff)
            {
                const double dx =
                    std::fabs(m.C[2*m.owner[ff]] - m.C[2*m.neighbour[ff]]);
                if (dx > 0.5*m.Lx) { ++seam; }
            }
            check(seam > 0, "faces cross the cyclic seam");
            bool onlyY = true;
            for (std::size_t b = 0; b < m.bSide.size(); ++b)
            {
                if (m.bSide[b] < 2) { onlyY = false; break; }
            }
            check(onlyY, "no boundary face on x (it is cyclic)");
        }

        // ---- refine then unrefine round trip -------------------------------
        std::printf("\nrefine then unrefine\n");
        {
            Mesh m;
            m.initialise(6, 4, 0.6, 0.4, false);
            const std::vector<Key> base = m.keys;

            std::vector<char> f(m.nCells(), 1);
            m.refine(f);
            m.refine(std::vector<char>(m.nCells(), 1));
            check(m.nCells() == 24*16, "refined twice");

            m.unrefine(std::vector<char>(m.nCells(), 1));
            m.unrefine(std::vector<char>(m.nCells(), 1));
            check(m.keys == base, "returns to the base mesh exactly");
        }

        // ---- a partial sibling group must not coarsen ----------------------
        std::printf("\npartial sibling group\n");
        {
            Mesh m;
            m.initialise(4, 4, 0.4, 0.4, false);
            std::vector<char> f(m.nCells(), 0);
            f[5] = 1;
            m.refine(f);

            const int before = m.nCells();
            std::vector<char> u(m.nCells(), 0);
            // flag only three of the four children
            int flagged = 0;
            for (int c = 0; c < m.nCells() && flagged < 3; ++c)
            {
                if (keyLevel(m.keys[c]) == 1) { u[c] = 1; ++flagged; }
            }
            check(m.unrefine(u) == 0, "three of four siblings do not coarsen");
            check(m.nCells() == before, "cell count unchanged");
        }

        // ---- field mapping -------------------------------------------------
        std::printf("\nfield mapping through refinement\n");
        {
            Mesh m;
            m.initialise(8, 8, 0.8, 0.8, false);

            std::vector<double> phi(m.nCells());
            for (int c = 0; c < m.nCells(); ++c)
            {
                double x0, y0, x1, y1;
                m.cellBox(c, x0, y0, x1, y1);
                phi[c] = 0.5*(x0 + x1) + 10.0*0.5*(y0 + y1);
            }

            const std::vector<Key> oldKeys = m.keys;
            std::vector<char> f(m.nCells(), 0);
            for (int c = 0; c < m.nCells(); c += 3) { f[c] = 1; }
            m.refine(f);

            std::vector<int> map;
            Mesh::parentMap(oldKeys, m.keys, map);

            bool allMapped = true;
            double worst = 0;
            for (int c = 0; c < m.nCells(); ++c)
            {
                if (map[c] < 0) { allMapped = false; continue; }

                double x0, y0, x1, y1;
                m.cellBox(c, x0, y0, x1, y1);
                const double exact = 0.5*(x0 + x1) + 10.0*0.5*(y0 + y1);
                // a child inherits its parent's value, so the error is
                // bounded by half the parent's diagonal in this linear field
                worst = std::max(worst, std::fabs(phi[map[c]] - exact));
            }

            check(allMapped, "every new cell has a parent");
            check(worst < 0.6, "inherited values are the parent's");
            std::printf("    worst inherited error %.4f\n", worst);
        }

        // ---- conservation across a map -----------------------------------
        //
        // The integral of a cell field over the domain must not move when the mesh
        // changes under it. Refinement cannot move it at all, since a child
        // inherits its parent's value; coarsening leaves it alone only if the four
        // children are combined by volume.
        //
        // The field is perturbed between maps deliberately. Siblings are born
        // holding their parent's value, so while nothing disturbs them taking one
        // child is arithmetically identical to averaging all four, and a
        // nearest-value map looks conservative.
        std::printf("\nconservation across a map\n");
        {
            Mesh m;
            m.initialise(60, 24, 0.15, 0.036, true);
            m.thickness = 5.0e-4;
            m.buildAddressing();

            auto integral = [](const std::vector<double>& V,
                               const std::vector<double>& phi)
            {
                // pairwise, so the reference is not limited by the
                // round-off of a long sequential sum
                std::vector<double> t(phi.size());
                for (std::size_t i = 0; i < phi.size(); ++i)
                {
                    t[i] = V[i]*phi[i];
                }
                std::size_t n = t.size();
                while (n > 1)
                {
                    const std::size_t h = (n + 1)/2;
                    for (std::size_t i = 0; i + h < n; ++i) { t[i] += t[i+h]; }
                    n = h;
                }
                return t.empty() ? 0.0 : t[0];
            };

            std::vector<double> phi(m.nCells()), phiN(m.nCells());
            for (int c = 0; c < m.nCells(); ++c)
            {
                const double x = m.C[2*c];
                phi[c] = 1.0 + 5.0*std::exp(-std::pow((x - 0.05)/0.004, 2));
            }
            phiN = phi;

            double worstCons = 0, worstNear = 0;
            int nCoarsened = 0;

            for (int step = 0; step < 30; ++step)
            {
                // the "solver": a spatially varying increment, so cells
                // that were siblings stop agreeing
                for (int c = 0; c < m.nCells(); ++c)
                {
                    const double d = 0.02*std::sin
                    (
                        400.0*m.C[2*c] + 37.0*m.C[2*c + 1] + step
                    );
                    phi[c] += d;
                    phiN[c] += d;
                }

                const double xc = 0.15*(0.15 + 0.7*step/29.0);
                std::vector<char> ref(m.nCells(), 0), crs(m.nCells(), 0);
                for (int c = 0; c < m.nCells(); ++c)
                {
                    const bool in = std::fabs(m.C[2*c] - xc) < 0.05*0.15;
                    const int L = int(keyLevel(m.keys[c]));
                    if (in && L < 3) { ref[c] = 1; }
                    if (!in && L > 0) { crs[c] = 1; }
                }

                const std::vector<Key> oldKeys = m.keys;
                const std::vector<double> oldV = m.V;
                const double I0 = integral(oldV, phi);
                const double I0n = integral(oldV, phiN);

                nCoarsened += m.unrefine(crs);
                m.refine(ref);
                m.balance();
                m.buildAddressing();

                std::vector<int> first, count;
                Mesh::mapRange(oldKeys, m.keys, first, count);
                std::vector<double> next(m.nCells());
                Mesh::mapField
                (
                    oldKeys, m.keys, first, count, phi.data(), next.data()
                );
                phi.swap(next);

                std::vector<int> pm;
                Mesh::parentMap(oldKeys, m.keys, pm);
                std::vector<double> nextN(m.nCells());
                for (int c = 0; c < m.nCells(); ++c)
                {
                    nextN[c] = phiN[pm[c]];
                }
                phiN.swap(nextN);

                worstCons = std::max
                (
                    worstCons, std::fabs(integral(m.V, phi) - I0)/std::fabs(I0)
                );
                worstNear = std::max
                (
                    worstNear,
                    std::fabs(integral(m.V, phiN) - I0n)/std::fabs(I0n)
                );
            }

            check(nCoarsened > 1000, "coarsening actually happened");
            check(worstCons < 1.0e-13, "mapField conserves the integral");
            check
            (
                worstNear > 1.0e-7,
                "parentMap does not, as documented"
            );
            std::printf("    %d groups coarsened\n", nCoarsened);
            std::printf("    mapField  worst drift %.3e\n", worstCons);
            std::printf("    parentMap worst drift %.3e\n", worstNear);
        }
    }

    // ---- timing at RDE scale ----------------------------------------------
    std::printf("\ntiming: 750 x 180 base mesh, a refined band sweeping x\n");
    {
        Mesh m;
        m.initialise(750, 180, 0.30, 0.072, true);
        m.thickness = 1e-4;

        double tRefine = 0, tUnrefine = 0, tBalance = 0, tAddressing = 0;
        const int nSteps = 20;

        for (int step = 0; step < nSteps; ++step)
        {
            const double xf = 0.30*(0.1 + 0.8*step/double(nSteps));

            // refine to three levels inside a band, as the RDE case does
            for (int lev = 0; lev < 3; ++lev)
            {
                std::vector<char> f(m.nCells(), 0);
                for (int c = 0; c < m.nCells(); ++c)
                {
                    double x0, y0, x1, y1;
                    m.cellBox(c, x0, y0, x1, y1);
                    const double cx = 0.5*(x0 + x1);
                    if (std::fabs(cx - xf) < 0.01
                     && keyLevel(m.keys[c]) == std::uint32_t(lev))
                    {
                        f[c] = 1;
                    }
                }
                double t = now();
                m.refine(f);
                tRefine += now() - t;
            }

            double t = now();
            m.balance();
            tBalance += now() - t;

            // coarsen everything outside the band
            std::vector<char> u(m.nCells(), 0);
            for (int c = 0; c < m.nCells(); ++c)
            {
                double x0, y0, x1, y1;
                m.cellBox(c, x0, y0, x1, y1);
                if (std::fabs(0.5*(x0 + x1) - xf) > 0.02) { u[c] = 1; }
            }
            t = now();
            m.unrefine(u);
            tUnrefine += now() - t;

            t = now();
            m.buildAddressing();
            tAddressing += now() - t;
        }

        std::printf("    final mesh        %d cells, %d faces\n",
            m.nCells(), m.nFaces());
        std::printf("    refine            %7.3f s  (%6.1f ms per step)\n",
            tRefine, 1000*tRefine/nSteps);
        std::printf("    unrefine          %7.3f s  (%6.1f ms per step)\n",
            tUnrefine, 1000*tUnrefine/nSteps);
        std::printf("    balance           %7.3f s  (%6.1f ms per step)\n",
            tBalance, 1000*tBalance/nSteps);
        std::printf("    buildAddressing   %7.3f s  (%6.1f ms per step)\n",
            tAddressing, 1000*tAddressing/nSteps);
        std::printf("    total             %7.3f s  (%6.1f ms per step)\n",
            tRefine + tUnrefine + tBalance + tAddressing,
            1000*(tRefine + tUnrefine + tBalance + tAddressing)/nSteps);

        checkEmpty(checkKeys(m), "key set still valid");
        checkEmpty(checkMesh(m), "mesh still valid");
    }

    std::printf("\n%s\n", nFail == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED");
    return nFail == 0 ? 0 : 1;
}

// ************************************************************************* //
