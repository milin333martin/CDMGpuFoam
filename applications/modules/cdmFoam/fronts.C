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

Description
    Tracking of the leading shock and the leading flame, globally and inside
    each channel of an arrestor, with the reduction over all the processors.

    The leading shock is the furthest cell (along direction) with
    p > pThreshold; the leading flame the furthest cell with Y < YFlame,
    i.e. the tip of the burnt gas.  Speeds are differences over a window of
    samples.  Results are written to

        postProcessing/fronts/<startTime>/fronts.dat

    with one row per sample: time, then for the whole domain and for each
    band, the shock position and speed, the flame position and speed, and
    the shock-flame separation.

\*---------------------------------------------------------------------------*/

#include "cdmFoam.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::cdmFoam::trackFronts()
{
    if (!frontTracking_) return;

    const label nBands = frontBands_.size();
    const label n = nBands + 1;                 // band 0 is the whole domain

    scalarField xShock(n, -great);
    scalarField xFlame(n, -great);

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& pI = p_.primitiveField();
    const scalarField& YI = Y_.primitiveField();

    forAll(pI, celli)
    {
        const scalar x = C[celli] & frontDirection_;
        const bool shock = pI[celli] > frontPThreshold_;
        const bool flame = YI[celli] < YFlame_;

        if (!shock && !flame) continue;

        if (shock) xShock[0] = max(xShock[0], x);
        if (flame) xFlame[0] = max(xFlame[0], x);

        const scalar y = C[celli].y();

        forAll(frontBands_, bandi)
        {
            if (y > frontBands_[bandi].first() && y < frontBands_[bandi].second())
            {
                if (shock)
                {
                    xShock[bandi + 1] = max(xShock[bandi + 1], x);
                }
                if (flame)
                {
                    xFlame[bandi + 1] = max(xFlame[bandi + 1], x);
                }
                break;
            }
        }
    }

    // One reduction for all the bands
    Pstream::listCombineGather(xShock, maxEqOp<scalar>());
    Pstream::listCombineScatter(xShock);
    Pstream::listCombineGather(xFlame, maxEqOp<scalar>());
    Pstream::listCombineScatter(xFlame);

    // History for the speeds
    frontTimes_.append(runTime.value());
    frontShock_.append(xShock);
    frontFlame_.append(xFlame);

    if (frontTimes_.size() > frontWindow_)
    {
        const label nDrop = frontTimes_.size() - frontWindow_;

        for (label i = 0; i < frontWindow_; i++)
        {
            frontTimes_[i] = frontTimes_[i + nDrop];
            frontShock_[i] = frontShock_[i + nDrop];
            frontFlame_[i] = frontFlame_[i + nDrop];
        }

        frontTimes_.setSize(frontWindow_);
        frontShock_.setSize(frontWindow_);
        frontFlame_.setSize(frontWindow_);
    }

    const label nS = frontTimes_.size();
    const scalar dt = frontTimes_[nS - 1] - frontTimes_[0];

    scalarField DShock(n, 0), DFlame(n, 0);

    if (nS > 1 && dt > vSmall)
    {
        forAll(DShock, i)
        {
            if (frontShock_[0][i] > -great && frontShock_[nS - 1][i] > -great)
            {
                DShock[i] = (frontShock_[nS - 1][i] - frontShock_[0][i])/dt;
            }
            if (frontFlame_[0][i] > -great && frontFlame_[nS - 1][i] > -great)
            {
                DFlame[i] = (frontFlame_[nS - 1][i] - frontFlame_[0][i])/dt;
            }
        }
    }

    // Report and write
    if (Pstream::master())
    {
        if (!frontFilePtr_.valid())
        {
            // the directory is named by the time the run started from,
            // like the other postProcessing output
            const fileName dir
            (
                runTime.globalPath()/"postProcessing"/"fronts"/frontStartName_
            );
            mkDir(dir);

            frontFilePtr_.set(new OFstream(dir/"fronts.dat"));

            OFstream& os = frontFilePtr_();
            os  << "# Leading shock (p > " << frontPThreshold_
                << " Pa) and leading flame (Y < " << YFlame_ << ")" << nl
                << "# columns: time";

            forAll(xShock, i)
            {
                const word name
                (
                    i == 0 ? word("all") : frontNames_[i - 1]
                );
                os  << "  xShock_" << name << "  DShock_" << name
                    << "  xFlame_" << name << "  DFlame_" << name
                    << "  sep_" << name;
            }
            os  << endl;
        }

        OFstream& os = frontFilePtr_();
        os  << runTime.value();

        forAll(xShock, i)
        {
            const scalar xs = xShock[i] > -great ? xShock[i] : scalar(0);
            const scalar xf = xFlame[i] > -great ? xFlame[i] : scalar(0);
            const scalar sep =
                (xShock[i] > -great && xFlame[i] > -great) ? xs - xf : 0;

            os  << ' ' << xs << ' ' << DShock[i]
                << ' ' << xf << ' ' << DFlame[i] << ' ' << sep;
        }
        os  << endl;
    }

    if (xShock[0] > -great)
    {
        Info<< "Leading shock: x = " << xShock[0]
            << "  D(" << nS << " samples) = " << DShock[0] << " m/s";

        if (xFlame[0] > -great)
        {
            Info<< "   flame: x = " << xFlame[0]
                << "  separation = " << xShock[0] - xFlame[0] << " m";
        }

        Info<< endl;
    }
}


// ************************************************************************* //
