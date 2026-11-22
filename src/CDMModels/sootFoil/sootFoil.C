/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam.
    GNU General Public License v3 or later; see COPYING.

\*---------------------------------------------------------------------------*/

#include "sootFoil.H"
#include "Time.H"
#include "OSspecific.H"
#include "Pstream.H"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{
    defineTypeNameAndDebug(sootFoil, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::CDM::sootFoil::outputDir() const
{
    return
        mesh_.time().globalPath()/"postProcessing"/"sootFoil"
       /mesh_.time().name();
}


void Foam::CDM::sootFoil::readRestart()
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName file(outputDir()/"sootFoil.raw");

    if (!isFile(file))
    {
        return;
    }

    std::ifstream is(file.c_str(), std::ios::binary);

    List<double> buf(foil_.size());
    is.read
    (
        reinterpret_cast<char*>(buf.begin()),
        std::streamsize(buf.size()*sizeof(double))
    );

    if (!is || is.gcount() != std::streamsize(buf.size()*sizeof(double)))
    {
        WarningInFunction
            << "Cannot read " << file << " (size mismatch with the current"
            << " sootFoil box/delta); starting a new foil" << endl;
        return;
    }

    forAll(foil_, k)
    {
        foil_[k] = buf[k];
    }

    Info<< "sootFoil: restarted from " << file << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDM::sootFoil::sootFoil(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    frameVelocity_(dict.lookupOrDefault<vector>("frameVelocity", Zero)),
    min_(dict.lookup<vector2D>("min")),
    max_(dict.lookup<vector2D>("max")),
    delta_(dict.lookup<scalar>("delta")),
    interval_(dict.lookupOrDefault<label>("interval", 1)),
    nSinceSample_(0),
    cartesian_(-1),
    rThickness_(0),
    nx_(0),
    ny_(0),
    foil_()
{
    if (mesh_.geometricD()[2] != -1)
    {
        FatalIOErrorInFunction(dict)
            << "sootFoil requires a 2D mesh in the x-y plane"
            << " (empty patches normal to z)" << exit(FatalIOError);
    }

    if (interval_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "sootFoil interval must be >= 1, is " << interval_
            << exit(FatalIOError);
    }

    if (delta_ <= 0 || max_.x() <= min_.x() || max_.y() <= min_.y())
    {
        FatalIOErrorInFunction(dict)
            << "Invalid sootFoil box " << min_ << " " << max_
            << " or delta " << delta_ << exit(FatalIOError);
    }

    nx_ = label(std::ceil((max_.x() - min_.x())/delta_));
    ny_ = label(std::ceil((max_.y() - min_.y())/delta_));

    if (scalar(nx_)*scalar(ny_) > 2e8)
    {
        FatalIOErrorInFunction(dict)
            << "sootFoil of " << nx_ << " x " << ny_
            << " pixels is too large; increase delta" << exit(FatalIOError);
    }

    foil_.setSize(nx_*ny_, 0);

    Info<< "sootFoil: " << nx_ << " x " << ny_ << " pixels of " << delta_
        << " m, lab-frame box " << min_ << " " << max_
        << ", frame velocity " << frameVelocity_ << endl;

    readRestart();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDM::sootFoil::~sootFoil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CDM::sootFoil::update(const volScalarField& p)
{
    // Sample every interval_ steps.  The first call always samples, so a
    // run never starts with a blank foil.
    if (nSinceSample_ > 0 && nSinceSample_ < interval_)
    {
        ++nSinceSample_;
        return;
    }
    nSinceSample_ = 1;

    // Lab coordinate = frame coordinate + frameVelocity*t
    const vector shift(frameVelocity_*mesh_.time().value());

    const vectorField& C = mesh_.C().primitiveField();
    const scalarField& V = mesh_.V();
    const scalarField& pI = p.primitiveField();

    const scalar rDelta = 1/delta_;

    // One-off: does the centre-and-volume box agree with the real point bounds
    // for every cell? If it does, the cell-to-point addressing is never
    // touched again.
    if (cartesian_ < 0)
    {
        const boundBox& bb = mesh_.bounds();
        rThickness_ = 1/max(bb.span().z(), vSmall);

        const pointField& points = mesh_.points();
        const labelListList& cellPoints = mesh_.cellPoints();

        cartesian_ = 1;
        forAll(cellPoints, celli)
        {
            const labelList& cPts = cellPoints[celli];
            scalar a0 = great, a1 = -great, b0 = great, b1 = -great;
            forAll(cPts, i)
            {
                const point& pt = points[cPts[i]];
                a0 = min(a0, pt.x());  a1 = max(a1, pt.x());
                b0 = min(b0, pt.y());  b1 = max(b1, pt.y());
            }

            const scalar h = std::sqrt(V[celli]*rThickness_);
            const scalar tol = 1e-6*h;

            if
            (
                mag(a1 - a0 - h) > tol || mag(b1 - b0 - h) > tol
             || mag(0.5*(a0 + a1) - C[celli].x()) > tol
             || mag(0.5*(b0 + b1) - C[celli].y()) > tol
            )
            {
                cartesian_ = 0;
                break;
            }
        }

        Info<< "sootFoil: cell boxes from centre and volume: "
            << (cartesian_ ? "exact" : "not exact, using point bounds")
            << endl;
    }

    const pointField* pointsPtr = nullptr;
    const labelListList* cellPointsPtr = nullptr;
    if (!cartesian_)
    {
        pointsPtr = &mesh_.points();
        cellPointsPtr = &mesh_.cellPoints();
    }

    forAll(C, celli)
    {
        scalar x0, x1, y0, y1;

        if (cartesian_)
        {
            const scalar h = 0.5*std::sqrt(V[celli]*rThickness_);
            x0 = C[celli].x() - h;  x1 = C[celli].x() + h;
            y0 = C[celli].y() - h;  y1 = C[celli].y() + h;
        }
        else
        {
            const labelList& cPts = (*cellPointsPtr)[celli];
            x0 = great; x1 = -great; y0 = great; y1 = -great;
            forAll(cPts, i)
            {
                const point& pt = (*pointsPtr)[cPts[i]];
                x0 = min(x0, pt.x());  x1 = max(x1, pt.x());
                y0 = min(y0, pt.y());  y1 = max(y1, pt.y());
            }
        }

        x0 += shift.x() - min_.x();
        x1 += shift.x() - min_.x();
        y0 += shift.y() - min_.y();
        y1 += shift.y() - min_.y();

        // Pixels whose centres (i + 0.5)*delta lie in [x0, x1)
        label i0 = label(std::ceil(x0*rDelta - 0.5));
        label i1 = label(std::ceil(x1*rDelta - 0.5)) - 1;
        label j0 = label(std::ceil(y0*rDelta - 0.5));
        label j1 = label(std::ceil(y1*rDelta - 0.5)) - 1;

        // Cells smaller than a pixel: the pixel containing the centre
        if (i0 > i1)
        {
            i0 = i1 = label
            (
                std::floor((C[celli].x() + shift.x() - min_.x())*rDelta)
            );
        }
        if (j0 > j1)
        {
            j0 = j1 = label
            (
                std::floor((C[celli].y() + shift.y() - min_.y())*rDelta)
            );
        }

        i0 = max(i0, label(0));
        j0 = max(j0, label(0));
        i1 = min(i1, nx_ - 1);
        j1 = min(j1, ny_ - 1);

        const scalar pc = pI[celli];

        for (label j = j0; j <= j1; j++)
        {
            scalar* row = foil_.begin() + j*nx_;

            for (label i = i0; i <= i1; i++)
            {
                if (pc > row[i])
                {
                    row[i] = pc;
                }
            }
        }
    }
}


void Foam::CDM::sootFoil::write()
{
    // Global maximum on the master
    Pstream::listCombineGather(foil_, maxEqOp<scalar>());

    if (!Pstream::master())
    {
        return;
    }

    const fileName dir(outputDir());
    mkDir(dir);

    // Restart data (native doubles)
    {
        List<double> buf(foil_.size());
        forAll(foil_, k)
        {
            buf[k] = foil_[k];
        }

        std::ofstream os
        (
            (dir/"sootFoil.raw").c_str(),
            std::ios::binary
        );
        os.write
        (
            reinterpret_cast<const char*>(buf.begin()),
            std::streamsize(buf.size()*sizeof(double))
        );
    }

    // Legacy VTK image (big-endian float32), pixel centres as points
    {
        std::ofstream os
        (
            (dir/"sootFoil.vtk").c_str(),
            std::ios::binary
        );

        os  << "# vtk DataFile Version 3.0\n"
            << "CDMGpuFoam sootFoil (lab frame) t = "
            << mesh_.time().value() << "\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_POINTS\n"
            << "DIMENSIONS " << nx_ << " " << ny_ << " 1\n"
            << "ORIGIN " << min_.x() + 0.5*delta_ << " "
            << min_.y() + 0.5*delta_ << " 0\n"
            << "SPACING " << delta_ << " " << delta_ << " " << delta_ << "\n"
            << "POINT_DATA " << nx_*ny_ << "\n"
            << "SCALARS sootFoil float 1\n"
            << "LOOKUP_TABLE default\n";

        const uint16_t one = 1;
        const bool littleEndian =
            *reinterpret_cast<const uint8_t*>(&one) == 1;

        forAll(foil_, k)
        {
            const float v = float(foil_[k]);
            uint32_t u;
            std::memcpy(&u, &v, sizeof(u));

            if (littleEndian)
            {
                u = __builtin_bswap32(u);
            }

            os.write(reinterpret_cast<const char*>(&u), sizeof(u));
        }

        os << "\n";
    }

    Info<< "sootFoil: written " << dir/"sootFoil.vtk" << endl;
}


// ************************************************************************* //
