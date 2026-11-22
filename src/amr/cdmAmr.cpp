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

\*---------------------------------------------------------------------------*/

#include "cdmAmr.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <sstream>
#include <string>

#ifdef CDM_STDPAR
#include <atomic>
#include <execution>
#endif

namespace cdmAmr
{

namespace
{

//- An index range that std::for_each can walk, so the same loop body is a
//  serial loop, a host-threaded one or a device kernel depending only on
//  how this file is compiled.  Same scaffolding as the flux back end.
struct CountingIterator
{
    using value_type = int;
    using difference_type = std::ptrdiff_t;
    using pointer = const int*;
    using reference = const int&;
    using iterator_category = std::random_access_iterator_tag;

    int i;

    CountingIterator() : i(0) {}
    explicit CountingIterator(int v) : i(v) {}

    int operator*() const { return i; }
    int operator[](difference_type n) const { return i + int(n); }

    CountingIterator& operator++() { ++i; return *this; }
    CountingIterator& operator--() { --i; return *this; }
    CountingIterator operator++(int) { CountingIterator t = *this; ++i; return t; }
    CountingIterator operator--(int) { CountingIterator t = *this; --i; return t; }
    CountingIterator& operator+=(difference_type n) { i += int(n); return *this; }
    CountingIterator& operator-=(difference_type n) { i -= int(n); return *this; }

    friend CountingIterator operator+(CountingIterator a, difference_type n) { a += n; return a; }
    friend CountingIterator operator+(difference_type n, CountingIterator a) { a += n; return a; }
    friend CountingIterator operator-(CountingIterator a, difference_type n) { a -= n; return a; }
    friend difference_type operator-(CountingIterator a, CountingIterator b)
    { return difference_type(a.i) - difference_type(b.i); }

    friend bool operator==(CountingIterator a, CountingIterator b) { return a.i == b.i; }
    friend bool operator!=(CountingIterator a, CountingIterator b) { return a.i != b.i; }
    friend bool operator< (CountingIterator a, CountingIterator b) { return a.i <  b.i; }
    friend bool operator> (CountingIterator a, CountingIterator b) { return a.i >  b.i; }
    friend bool operator<=(CountingIterator a, CountingIterator b) { return a.i <= b.i; }
    friend bool operator>=(CountingIterator a, CountingIterator b) { return a.i >= b.i; }
};

template<class F>
inline void pfor(int n, F f)
{
#ifdef CDM_STDPAR
    std::for_each(std::execution::par_unseq, CountingIterator(0), CountingIterator(n), f);
#elif defined(CDM_OPENMP)
    // The shared library OpenFOAM loads cannot use stdpar: it has to be built
    // with OpenFOAM's own compiler, and mixing nvc++ objects into it caused an
    // allocator mismatch. OpenMP gives the same loops the same parallelism.
    // Most bodies write only to indices derived from their own counter.
    // balance() is the exception, where several cells may mark the same coarse
    // neighbour; it uses an atomic write.
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < n; ++k) { f(k); }
#else
    for (int k = 0; k < n; ++k) { f(k); }
#endif
}

} // End anonymous namespace

// * * * * * * * * * * * * * * * Morton coding  * * * * * * * * * * * * * * //

std::uint32_t morton2D(std::uint32_t i, std::uint32_t j)
{
    std::uint32_t m = 0;
    for (int b = 0; b < maxLevel; ++b)
    {
        m |= ((i >> b) & 1u) << (2*b);
        m |= ((j >> b) & 1u) << (2*b + 1);
    }
    return m;
}


void unMorton2D(std::uint32_t m, std::uint32_t& i, std::uint32_t& j)
{
    i = 0;
    j = 0;
    for (int b = 0; b < maxLevel; ++b)
    {
        i |= ((m >> (2*b)) & 1u) << b;
        j |= ((m >> (2*b + 1)) & 1u) << b;
    }
}


// * * * * * * * * * * * * * * * Local helpers  * * * * * * * * * * * * * * //

namespace
{

// Face directions are 0 = -x, 1 = +x, 2 = -y, 3 = +y. The step in (i, j)
// is derived from that in neighbourPos rather than read from a table: a
// namespace-scope array is not reachable from device code under -stdpar.

// Where a cell sits inside its base cell
struct Pos
{
    std::uint64_t base;
    std::uint32_t level;
    std::uint32_t i, j;
    std::uint32_t n;        // 1 << level
};

inline Pos decode(Key k)
{
    Pos p;
    p.base = keyBase(k);
    p.level = keyLevel(k);
    unMorton2D(keyMorton(k), p.i, p.j);
    p.n = 1u << p.level;
    return p;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Mesh  * * * * * * * * * * * * * * * * * * //

Mesh::Mesh() = default;
Mesh::~Mesh() = default;


namespace
{

//- Move a buffer into a heap holder that is never deleted.  The move
//  leaves the source empty with no capacity, so the source's own
//  destructor has nothing left to free.
template<class T>
void leakBuffer(std::vector<T>& v)
{
    if (v.capacity() != 0)
    {
        (void)new std::vector<T>(std::move(v));
    }
}

} // End anonymous namespace


void Mesh::releaseWithoutFreeing()
{
    leakBuffer(keys);
    leakBuffer(owner);
    leakBuffer(neighbour);
    leakBuffer(faceCell);
    leakBuffer(faceSide);
    leakBuffer(Sf);
    leakBuffer(magSf);
    leakBuffer(weights);
    leakBuffer(deltaCoeffs);
    leakBuffer(V);
    leakBuffer(C);
    leakBuffer(bCell);
    leakBuffer(bSide);
}


void Mesh::initialise(int nx_, int ny_, double Lx_, double Ly_, bool cyclicX_)
{
    nx = nx_;
    ny = ny_;
    Lx = Lx_;
    Ly = Ly_;
    cyclicX = cyclicX_;

    keys.resize(static_cast<std::size_t>(nx)*ny);
    for (std::size_t b = 0; b < keys.size(); ++b)
    {
        keys[b] = makeKey(b, 0, 0);
    }
}


int Mesh::find(Key k) const
{
    const auto it = std::lower_bound(keys.begin(), keys.end(), k);
    if (it == keys.end() || *it != k) { return -1; }
    return static_cast<int>(it - keys.begin());
}


void Mesh::cellBox
(
    int c, double& x0, double& y0, double& x1, double& y1
) const
{
    const Pos p = decode(keys[c]);
    const double dxb = Lx/nx, dyb = Ly/ny;
    const double bx = static_cast<double>(p.base % nx);
    const double by = static_cast<double>(p.base / nx);
    const double h = 1.0/p.n;

    x0 = (bx + p.i*h)*dxb;
    x1 = (bx + (p.i + 1)*h)*dxb;
    y0 = (by + p.j*h)*dyb;
    y1 = (by + (p.j + 1)*h)*dyb;
}


int Mesh::refine(const std::vector<char>& flag, std::vector<int>* newIndex)
{
    const std::size_t nOld = keys.size();

    int nRefined = 0;
    for (std::size_t c = 0; c < nOld; ++c)
    {
        if (flag[c] && keyLevel(keys[c]) + 1 < static_cast<std::uint32_t>(maxLevel))
        {
            ++nRefined;
        }
    }
    if (nRefined == 0)
    {
        if (newIndex)
        {
            newIndex->resize(nOld);
            for (std::size_t c = 0; c < nOld; ++c)
            {
                (*newIndex)[c] = static_cast<int>(c);
            }
        }
        return 0;
    }

    std::vector<Key> out;
    out.reserve(nOld + 3*nRefined);

    if (newIndex) { newIndex->assign(nOld, -1); }

    for (std::size_t c = 0; c < nOld; ++c)
    {
        const Key k = keys[c];
        const std::uint32_t L = keyLevel(k);

        if (newIndex) { (*newIndex)[c] = static_cast<int>(out.size()); }

        if (flag[c] && L + 1 < static_cast<std::uint32_t>(maxLevel))
        {
            // Appending two bits to the Morton code is exactly descending
            // one level, so the four children are 4m, 4m+1, 4m+2, 4m+3.
            const std::uint32_t m = keyMorton(k) << 2;
            for (std::uint32_t c4 = 0; c4 < 4; ++c4)
            {
                out.push_back(makeKey(keyBase(k), L + 1, m | c4));
            }
        }
        else
        {
            out.push_back(k);
        }
    }

    // No sort: the children tile their parent's z interval, so they land
    // exactly where the parent was and the array is still ordered.
    keys.swap(out);
    return nRefined;
}


int Mesh::unrefine(const std::vector<char>& flag)
{
    const std::size_t n = keys.size();

    std::vector<Key> out;
    out.reserve(n);

    int nGroups = 0;
    std::size_t c = 0;

    while (c < n)
    {
        const Key k = keys[c];
        const std::uint32_t L = keyLevel(k);
        const std::uint32_t m = keyMorton(k);

        // A coarsenable group is four siblings: same base, same level, Morton
        // codes 4p..4p+3, all present and all flagged. Sorting by key puts them
        // next to each other, so this is a local test.
        bool group = false;

        if (L > 0 && (m & 3u) == 0 && c + 3 < n)
        {
            group = true;
            for (std::uint32_t s = 0; s < 4; ++s)
            {
                const Key ks = keys[c + s];
                if
                (
                    keyBase(ks) != keyBase(k)
                 || keyLevel(ks) != L
                 || keyMorton(ks) != m + s
                 || !flag[c + s]
                )
                {
                    group = false;
                    break;
                }
            }
        }

        if (group)
        {
            out.push_back(makeKey(keyBase(k), L - 1, m >> 2));
            ++nGroups;
            c += 4;
        }
        else
        {
            out.push_back(k);
            ++c;
        }
    }

    if (nGroups == 0) { return 0; }

    // No sort, for the same reason as refine()
    keys.swap(out);
    return nGroups;
}


namespace
{

//- The position of the neighbour of p in direction d, at p's own level.
//  Returns false if it leaves the domain.
inline bool neighbourPos
(
    int nx, int ny, bool cyclicX, const Pos& p, int d,
    std::uint64_t& nbBase, std::uint32_t& ni, std::uint32_t& nj
)
{
    const long long bx = static_cast<long long>(p.base % nx);
    const long long by = static_cast<long long>(p.base / nx);

    // d = 0, 2 step backwards, d = 1, 3 step forwards; d < 2 steps in x
    const long long step = (d & 1) ? 1 : -1;
    const bool xDir = (d < 2);

    long long i = static_cast<long long>(p.i) + (xDir ? step : 0);
    long long j = static_cast<long long>(p.j) + (xDir ? 0 : step);
    long long nbx = bx, nby = by;
    const long long nn = p.n;

    if (i < 0)        { nbx = bx - 1; i = nn - 1; }
    else if (i >= nn) { nbx = bx + 1; i = 0; }

    if (j < 0)        { nby = by - 1; j = nn - 1; }
    else if (j >= nn) { nby = by + 1; j = 0; }

    if (nbx < 0 || nbx >= nx)
    {
        if (!cyclicX) { return false; }
        nbx = (nbx + nx) % nx;
    }
    if (nby < 0 || nby >= ny) { return false; }

    nbBase = static_cast<std::uint64_t>(nby*nx + nbx);
    ni = static_cast<std::uint32_t>(i);
    nj = static_cast<std::uint32_t>(j);
    return true;
}


//- The leaf covering the point (nbBase, z), whatever its level.
//
//  Leaves tile the z axis within a base cell and are sorted by z, so the
//  one covering a point is the last whose z does not exceed it.  That is a
//  single upper_bound, where walking the levels took one search per level.
inline int findCovering
(
    const Key* keys, int nKeys, std::uint64_t nbBase, std::uint32_t z,
    std::uint32_t& foundLevel
)
{
    // Hand-rolled rather than std::upper_bound: this runs inside a device
    // kernel, where a raw pointer and an integer index are the only things
    // that can be relied on.
    const Key probe = (nbBase << 20) | (static_cast<Key>(z) << 4) | 0xFu;

    int lo = 0, hi = nKeys;
    while (lo < hi)
    {
        const int mid = lo + ((hi - lo) >> 1);
        if (keys[mid] <= probe) { lo = mid + 1; } else { hi = mid; }
    }

    if (lo == 0) { return -1; }
    const int idx = lo - 1;
    const Key k = keys[idx];

    if (keyBase(k) != nbBase) { return -1; }
    if (z < keyZ(k) || z >= keyZ(k) + keyZSpan(k)) { return -1; }

    foundLevel = keyLevel(k);
    return idx;
}


//- The same, for a target expressed as a position at a given level
inline int findCoveringCell
(
    const Key* keys, int nKeys, std::uint64_t nbBase, std::uint32_t level,
    std::uint32_t ni, std::uint32_t nj, std::uint32_t& foundLevel
)
{
    const std::uint32_t z = morton2D(ni, nj) << (2*(maxLevel - level));
    return findCovering(keys, nKeys, nbBase, z, foundLevel);
}

} // End anonymous namespace


int Mesh::balance()
{
    int added = 0;

    std::vector<std::uint8_t> baseMinLevel;

    for (int sweep = 0; sweep < maxLevel + 2; ++sweep)
    {
        std::vector<char> flag(keys.size(), 0);
#ifdef CDM_STDPAR
        // nvc++ supports atomic_ref for 32- and 64-bit values, so parallel marks
        // go in 32-bit slots and convert to the char flags refine() expects once
        // the parallel pass completes.
        std::vector<std::uint32_t> atomicFlag(keys.size(), 0);
#endif
        int nFlag = 0;

        // The coarsest cell in each base cell. If that is no more than one level
        // below the cell doing the asking, nothing in that base cell can violate
        // 2:1 and the lookup is skipped entirely.
        baseMinLevel.assign(static_cast<std::size_t>(nx)*ny, 0xFFu);
        for (const Key k : keys)
        {
            std::uint8_t& mn = baseMinLevel[keyBase(k)];
            const std::uint8_t L = static_cast<std::uint8_t>(keyLevel(k));
            if (L < mn) { mn = L; }
        }

        {
            const Key* kp = keys.data();
            const int nK = static_cast<int>(keys.size());
            const std::uint8_t* bml = baseMinLevel.data();
#ifdef CDM_STDPAR
            std::uint32_t* fl = atomicFlag.data();
#else
            char* fl = flag.data();
#endif
            const int nxL = nx, nyL = ny;
            const bool cycL = cyclicX;

            // Several cells may flag the same coarse neighbour, and concurrent
            // non-atomic writes are a data race even when every writer stores the
            // same value. Relaxed atomic stores on the parallel paths; the pass
            // completes before any flag is read back.
            pfor(nK, [=](int c)
            {
                const Pos p = decode(kp[c]);

                // Levels 0 and 1 cannot have a neighbour coarser than L-1
                if (p.level < 2) { return; }

                for (int d = 0; d < 4; ++d)
                {
                    std::uint64_t nbBase;
                    std::uint32_t ni, nj;
                    if (!neighbourPos(nxL, nyL, cycL, p, d, nbBase, ni, nj))
                    {
                        continue;
                    }

                    if (std::uint32_t(bml[nbBase]) + 1u >= p.level) { continue; }

                    std::uint32_t nL;
                    const int nb = findCoveringCell
                    (
                        kp, nK, nbBase, p.level, ni, nj, nL
                    );

                    // A neighbour two or more levels coarser violates 2:1,
                    // so it is the one that has to be refined.
                    if (nb >= 0 && nL + 1 < p.level)
                    {
#if defined(CDM_STDPAR)
                        std::atomic_ref<std::uint32_t>(fl[nb]).store
                        (
                            1, std::memory_order_relaxed
                        );
#elif defined(CDM_OPENMP)
                        #pragma omp atomic write
                        fl[nb] = 1;
#else
                        fl[nb] = 1;
#endif
                    }
                }
            });

            for (std::size_t c = 0; c < flag.size(); ++c)
            {
#ifdef CDM_STDPAR
                if (atomicFlag[c]) { flag[c] = 1; ++nFlag; }
#else
                if (flag[c]) { ++nFlag; }
#endif
            }
        }

        if (nFlag == 0) { break; }
        added += refine(flag);
    }

    return added;
}


void Mesh::buildAddressing()
{
    const int n = nCells();
    const double dxb = Lx/nx, dyb = Ly/ny;

    V.assign(n, 0.0);
    C.assign(2*n, 0.0);

    {
        const Key* kp = keys.data();
        double* Vp = V.data();
        double* Cp = C.data();
        const double dxb_ = dxb, dyb_ = dyb, th = thickness;
        const int nxL = nx;

        pfor(n, [=](int c)
        {
            const Pos p = decode(kp[c]);
            const double h = 1.0/p.n;
            const double bx = double(p.base % nxL);
            const double by = double(p.base / nxL);

            const double x0 = (bx + p.i*h)*dxb_;
            const double x1 = (bx + (p.i + 1)*h)*dxb_;
            const double y0 = (by + p.j*h)*dyb_;
            const double y1 = (by + (p.j + 1)*h)*dyb_;

            Vp[c] = (x1 - x0)*(y1 - y0)*th;
            Cp[2*c] = 0.5*(x0 + x1);
            Cp[2*c + 1] = 0.5*(y0 + y1);
        });
    }

    // Counting first and filling second is what lets this run on a device:
    // every cell writes to a slot it computed itself, so there is no append
    // and no contention. The neighbour search is kept, not repeated.
    std::vector<int> nbIdx(4*n, -1);
    std::vector<std::uint8_t> nbLev(4*n, 0);
    std::vector<int> nFace(n + 1, 0), nBnd(n + 1, 0);

    {
        const Key* kp = keys.data();
        const int nK = n;
        int* nbI = nbIdx.data();
        std::uint8_t* nbL = nbLev.data();
        int* nFp = nFace.data();
        int* nBp = nBnd.data();
        const int nxL = nx, nyL = ny;
        const bool cycL = cyclicX;

        pfor(n, [=](int c)
        {
            const Pos p = decode(kp[c]);
            int nf = 0, nb_ = 0;

            for (int d = 0; d < 4; ++d)
            {
                std::uint64_t nbBase;
                std::uint32_t ni, nj;

                if (!neighbourPos(nxL, nyL, cycL, p, d, nbBase, ni, nj))
                {
                    nbI[4*c + d] = -2;        // outer boundary
                    ++nb_;
                    continue;
                }

                std::uint32_t nL;
                const int nb =
                    findCoveringCell(kp, nK, nbBase, p.level, ni, nj, nL);

                if (nb < 0) { continue; }

                // Each face is made once, by the finer side, or at equal levels by the
                // lower key. A cell facing finer neighbours makes none and they each make
                // one to it, which is how a hanging edge becomes two faces.
                if (nL > p.level) { continue; }
                if (nL == p.level && kp[c] > kp[nb]) { continue; }

                nbI[4*c + d] = nb;
                nbL[4*c + d] = static_cast<std::uint8_t>(nL);
                ++nf;
            }

            nFp[c + 1] = nf;
            nBp[c + 1] = nb_;
        });
    }

    // ---- pass 2: where each cell's faces go ----
    for (int c = 0; c < n; ++c)
    {
        nFace[c + 1] += nFace[c];
        nBnd[c + 1] += nBnd[c];
    }

    const int nF = nFace[n];
    const int nB = nBnd[n];

    owner.assign(nF, 0);
    neighbour.assign(nF, 0);
    faceCell.assign(nF, 0);
    faceSide.assign(nF, 0);
    Sf.assign(2*nF, 0.0);
    magSf.assign(nF, 0.0);
    weights.assign(nF, 0.0);
    deltaCoeffs.assign(nF, 0.0);
    bCell.assign(nB, 0);
    bSide.assign(nB, 0);

    // ---- pass 3: fill ----
    {
        const Key* kp = keys.data();
        const int* nbI = nbIdx.data();
        const int* nFp = nFace.data();
        const int* nBp = nBnd.data();
        const double* Cp = C.data();
        int* ownP = owner.data();
        int* neiP = neighbour.data();
        int* fcP = faceCell.data();
        int* fsP = faceSide.data();
        double* SfP = Sf.data();
        double* magP = magSf.data();
        double* wP = weights.data();
        double* dcP = deltaCoeffs.data();
        int* bcP = bCell.data();
        int* bsP = bSide.data();
        const double dxb_ = dxb, dyb_ = dyb, th = thickness, LxL = Lx;
        const int nxL = nx;
        const bool cycL = cyclicX;

        pfor(n, [=](int c)
        {
        const Pos p = decode(kp[c]);
        const double h = 1.0/p.n;
        const double dx = h*dxb_, dy = h*dyb_;

        const double bxx = double(p.base % nxL);
        const double byy = double(p.base / nxL);
        const double x0 = (bxx + p.i*h)*dxb_;
        const double x1 = (bxx + (p.i + 1)*h)*dxb_;
        const double y0 = (byy + p.j*h)*dyb_;
        const double y1 = (byy + (p.j + 1)*h)*dyb_;

        int f = nFp[c];
        int b = nBp[c];

        for (int d = 0; d < 4; ++d)
        {
            const int nb = nbI[4*c + d];

            if (nb == -2)
            {
                bcP[b] = c;
                bsP[b] = d;
                ++b;
                continue;
            }
            if (nb < 0) { continue; }

            const double area = ((d < 2) ? dy : dx)*th;
            const double sgn = (d == 1 || d == 3) ? 1.0 : -1.0;

            double sfx = 0, sfy = 0, fx, fy;

            if (d < 2)
            {
                sfx = sgn*area;
                fx = (d == 1) ? x1 : x0;
                fy = 0.5*(y0 + y1);
            }
            else
            {
                sfy = sgn*area;
                fx = 0.5*(x0 + x1);
                fy = (d == 3) ? y1 : y0;
            }

            // Neighbour centre, unwrapped across a cyclic seam so the
            // owner-to-neighbour vector goes the short way round
            double cnx = Cp[2*nb], cny = Cp[2*nb + 1];
            if (cycL)
            {
                if (cnx - Cp[2*c] > 0.5*LxL) { cnx -= LxL; }
                else if (Cp[2*c] - cnx > 0.5*LxL) { cnx += LxL; }
            }

            const double dvx = cnx - Cp[2*c];
            const double dvy = cny - Cp[2*c + 1];
            const double dmag = std::sqrt(dvx*dvx + dvy*dvy);

            //   SfdOwn = |Sf . (Cf - C_own)|
            //   SfdNei = |Sf . (C_nei - Cf)|
            //   w      = SfdNei/(SfdOwn + SfdNei)
            // and fall back to the ratio of plain distances only when that denominator
            // degenerates. The two agree on an orthogonal face, which is why the
            // distance ratio passed every test here; at a hanging node they differ by
            // up to 1.7% against fvMesh, and that is exactly where the refinement is.
            const double SfdOwn =
                std::fabs(sfx*(fx - Cp[2*c]) + sfy*(fy - Cp[2*c + 1]));
            const double SfdNei =
                std::fabs(sfx*(cnx - fx) + sfy*(cny - fy));
            const double SfdSum = SfdOwn + SfdNei;

            double w;
            if (SfdSum > 1e-300)
            {
                w = SfdNei/SfdSum;
            }
            else
            {
                const double dOwn = std::sqrt
                (
                    (fx - Cp[2*c])*(fx - Cp[2*c])
                  + (fy - Cp[2*c + 1])*(fy - Cp[2*c + 1])
                );
                const double dNei = std::max(dmag - dOwn, 1e-300);
                w = dNei/(dOwn + dNei);
            }
            int own = c, nei = nb;

            // OpenFOAM requires owner < neighbour
            if (own > nei)
            {
                const int t = own; own = nei; nei = t;
                sfx = -sfx;
                sfy = -sfy;
                w = 1.0 - w;
            }

            // Recorded from the creating cell, not from own/nei: the
            // swap above may have exchanged them.
            fcP[f] = c;
            fsP[f] = d;

            ownP[f] = own;
            neiP[f] = nei;
            SfP[2*f] = sfx;
            SfP[2*f + 1] = sfy;
            magP[f] = area;
            wP[f] = w;
            dcP[f] = 1.0/((dmag > 1e-300) ? dmag : 1e-300);
            ++f;
        }
        });
    }
}


void Mesh::parentMap
(
    const std::vector<Key>& oldKeys,
    const std::vector<Key>& newKeys,
    std::vector<int>& map
)
{
    map.assign(newKeys.size(), -1);

    // One covering search answers both directions: refined or unchanged, the
    // covering old leaf is the ancestor; coarsened, it is the lowest-Morton
    // child absorbed and the other three are discarded. mapField is the
    // conservative route.
    for (std::size_t c = 0; c < newKeys.size(); ++c)
    {
        const Key k = newKeys[c];
        std::uint32_t foundLevel = 0;
        map[c] = findCovering
        (
            oldKeys.data(), static_cast<int>(oldKeys.size()),
            keyBase(k), keyZ(k), foundLevel
        );
    }
}


void Mesh::mapRange
(
    const std::vector<Key>& oldKeys,
    const std::vector<Key>& newKeys,
    std::vector<int>& first,
    std::vector<int>& count
)
{
    const int nNew = static_cast<int>(newKeys.size());
    const int nOld = static_cast<int>(oldKeys.size());

    first.assign(newKeys.size(), 0);
    count.assign(newKeys.size(), 0);

    const Key* op = oldKeys.data();
    const Key* np = newKeys.data();
    int* fp = first.data();
    int* cp = count.data();

    pfor(nNew, [=](int c)
    {
        const Key k = np[c];
        const std::uint64_t base = keyBase(k);
        const std::uint32_t z0 = keyZ(k);
        const std::uint32_t z1 = z0 + keyZSpan(k);

        // The old leaf covering the new cell's first z. Under refinement it starts
        // before z0 and extends past z1; under coarsening it is the first of the
        // children. Either way it is the start of the run.
        std::uint32_t lvl = 0;
        const int lo = findCovering(op, nOld, base, z0, lvl);

        if (lo < 0)
        {
            // Can only happen if the two key sets do not describe the same base mesh,
            // which is a programming error rather than a state the caller can be in.
            fp[c] = 0;
            cp[c] = 0;
            return;
        }

        // Walk forward while the old leaves still start inside [z0, z1). At most
        // four for a single coarsening step, and the loop is the honest bound
        // rather than a hard-coded 4.
        int hi = lo + 1;
        while
        (
            hi < nOld
         && keyBase(op[hi]) == base
         && keyZ(op[hi]) < z1
        )
        {
            ++hi;
        }

        fp[c] = lo;
        cp[c] = hi - lo;
    });
}


void Mesh::mapField
(
    const std::vector<Key>& oldKeys,
    const std::vector<Key>& newKeys,
    const std::vector<int>& first,
    const std::vector<int>& count,
    const double* phiOld,
    double* phiNew
)
{
    const int nNew = static_cast<int>(newKeys.size());

    const Key* op = oldKeys.data();
    const Key* np = newKeys.data();
    const int* fp = first.data();
    const int* cp = count.data();

    pfor(nNew, [=](int c)
    {
        const std::uint32_t lNew = keyLevel(np[c]);

        double sum = 0.0;

        for (int i = 0; i < cp[c]; ++i)
        {
            const int o = fp[c] + i;
            const std::uint32_t lOld = keyLevel(op[o]);

            // The overlap of an old and a new leaf is whichever is smaller, since one
            // always contains the other in a tree. As a fraction of the new cell's
            // volume that is 4^(lNew - max(lNew, lOld)).
            const std::uint32_t lMax = (lOld > lNew) ? lOld : lNew;
            const double wgt =
                1.0/static_cast<double>(1u << (2*(lMax - lNew)));

            sum += wgt*phiOld[o];
        }

        phiNew[c] = sum;
    });
}


// * * * * * * * * * * * * * * * * Checks  * * * * * * * * * * * * * * * * * //

std::string checkKeys(const Mesh& m)
{
    std::ostringstream err;

    if (!std::is_sorted(m.keys.begin(), m.keys.end()))
    {
        err << "keys are not sorted; ";
    }

    for (std::size_t c = 1; c < m.keys.size(); ++c)
    {
        if (m.keys[c] == m.keys[c - 1])
        {
            err << "duplicate key at " << c << "; ";
            break;
        }
    }

    // Every base cell must be covered exactly once, which for a quadtree
    // means its leaves' areas summing to one base cell.
    std::vector<double> cover(static_cast<std::size_t>(m.nx)*m.ny, 0.0);
    for (const Key k : m.keys)
    {
        const double a = 1.0/static_cast<double>(1u << (2*keyLevel(k)));
        cover[keyBase(k)] += a;
    }
    for (std::size_t b = 0; b < cover.size(); ++b)
    {
        if (std::fabs(cover[b] - 1.0) > 1e-9)
        {
            err << "base cell " << b << " covered " << cover[b]
                << " times; ";
            break;
        }
    }

    return err.str();
}


std::string checkMesh(const Mesh& m)
{
    std::ostringstream err;
    const int n = m.nCells();

    // Every cell closed: its face area vectors, plus its boundary faces,
    // must sum to zero.
    std::vector<double> sx(n, 0.0), sy(n, 0.0);

    for (int f = 0; f < m.nFaces(); ++f)
    {
        sx[m.owner[f]] += m.Sf[2*f];
        sy[m.owner[f]] += m.Sf[2*f + 1];
        sx[m.neighbour[f]] -= m.Sf[2*f];
        sy[m.neighbour[f]] -= m.Sf[2*f + 1];
    }

    const double dxb = m.Lx/m.nx, dyb = m.Ly/m.ny;
    for (std::size_t b = 0; b < m.bCell.size(); ++b)
    {
        const int c = m.bCell[b];
        const int d = m.bSide[b];
        const Pos p = decode(m.keys[c]);
        const double h = 1.0/p.n;
        const double area = ((d < 2) ? h*dyb : h*dxb)*m.thickness;
        const double sgn = (d == 1 || d == 3) ? 1.0 : -1.0;

        if (d < 2) { sx[c] += sgn*area; } else { sy[c] += sgn*area; }
    }

    double worst = 0;
    int worstCell = -1;
    for (int c = 0; c < n; ++c)
    {
        const double e = std::sqrt(sx[c]*sx[c] + sy[c]*sy[c])
            /std::max(m.V[c]/m.thickness, 1e-300);
        if (e > worst) { worst = e; worstCell = c; }
    }
    if (worst > 1e-9)
    {
        err << "cell " << worstCell << " not closed, relative "
            << worst << "; ";
    }

    // Volumes sum to the domain
    double vol = 0;
    for (int c = 0; c < n; ++c) { vol += m.V[c]; }
    const double exact = m.Lx*m.Ly*m.thickness;
    if (std::fabs(vol - exact)/exact > 1e-12)
    {
        err << "volume " << vol << " against " << exact << "; ";
    }

    // No face may join a cell to itself, and all indices must be in range
    for (int f = 0; f < m.nFaces(); ++f)
    {
        if (m.owner[f] == m.neighbour[f])
        {
            err << "face " << f << " joins a cell to itself; ";
            break;
        }
        if (m.owner[f] < 0 || m.neighbour[f] >= n)
        {
            err << "face " << f << " out of range; ";
            break;
        }
        if (m.owner[f] >= m.neighbour[f])
        {
            err << "face " << f << " has owner >= neighbour; ";
            break;
        }
    }

    // 2:1 balance
    for (int f = 0; f < m.nFaces(); ++f)
    {
        const int lo = static_cast<int>(keyLevel(m.keys[m.owner[f]]));
        const int ln = static_cast<int>(keyLevel(m.keys[m.neighbour[f]]));
        if (std::abs(lo - ln) > 1)
        {
            err << "face " << f << " spans " << std::abs(lo - ln)
                << " levels; ";
            break;
        }
    }

    return err.str();
}

} // End namespace cdmAmr

// ************************************************************************* //
