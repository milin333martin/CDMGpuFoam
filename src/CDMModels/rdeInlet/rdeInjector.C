/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: a chemical-diffusive-model solver
   \\    /   O peration     | for gaseous detonations in OpenFOAM-13
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, distributed under the GNU General
    Public License version 3 or later.

\*---------------------------------------------------------------------------*/

#include "rdeInjector.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rdeInjector::rdeInjector(const dictionary& dict)
:
    p0_(dict.lookup<scalar>("p0")),
    T0_(dict.lookup<scalar>("T0")),
    gamma_(dict.lookup<scalar>("gamma")),
    R_(dict.lookup<scalar>("R")),
    slotWidth_(dict.lookup<scalar>("slotWidth")),
    slotSpacing_(dict.lookup<scalar>("slotSpacing")),
    slotOffset_(dict.lookupOrDefault<scalar>("slotOffset", 0)),
    tangent_(dict.lookupOrDefault<vector>("tangent", vector(1, 0, 0))),
    direction_(dict.lookupOrDefault<vector>("direction", Zero)),
    saturating_(true)
{
    const word chokedModel =
        dict.lookupOrDefault<word>("chokedModel", "saturating");

    if (chokedModel == "saturating")
    {
        saturating_ = true;
    }
    else if (chokedModel == "unbounded")
    {
        saturating_ = false;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown chokedModel " << chokedModel << nl
            << "    valid choices: saturating, unbounded"
            << exit(FatalIOError);
    }

    if (p0_ <= 0 || T0_ <= 0 || R_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "p0, T0 and R must all be positive: "
            << "p0 = " << p0_ << ", T0 = " << T0_ << ", R = " << R_
            << exit(FatalIOError);
    }

    if (gamma_ <= 1)
    {
        FatalIOErrorInFunction(dict)
            << "gamma must exceed 1: " << gamma_
            << exit(FatalIOError);
    }

    if (slotWidth_ <= 0 || slotSpacing_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "slotWidth and slotSpacing must both be positive: "
            << "slotWidth = " << slotWidth_
            << ", slotSpacing = " << slotSpacing_
            << exit(FatalIOError);
    }

    if (slotWidth_ > slotSpacing_)
    {
        FatalIOErrorInFunction(dict)
            << "slotWidth (" << slotWidth_ << ") exceeds slotSpacing ("
            << slotSpacing_ << "), so the slots overlap and the injector"
            << " is continuous.  Set slotWidth == slotSpacing if that is"
            << " what you intend." << exit(FatalIOError);
    }

    const scalar magTangent = mag(tangent_);

    if (magTangent < small)
    {
        FatalIOErrorInFunction(dict)
            << "tangent must be a non-zero vector" << exit(FatalIOError);
    }

    tangent_ /= magTangent;

    if (mag(direction_) > small)
    {
        direction_ /= mag(direction_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rdeInjector::pCrit() const
{
    return p0_*pow(2/(gamma_ + 1), gamma_/(gamma_ - 1));
}


Foam::scalar Foam::rdeInjector::vMax() const
{
    return sqrt(2*gamma_*R_*T0_/(gamma_ - 1));
}


bool Foam::rdeInjector::inSlot(const scalar x) const
{
    // fmod is negative for negative arguments, so fold into [0, slotSpacing)
    scalar s = fmod(x, slotSpacing_);

    if (s < 0)
    {
        s += slotSpacing_;
    }

    return s < slotWidth_;
}


Foam::rdeInjector::regime Foam::rdeInjector::state(const scalar p) const
{
    if (p >= p0_)
    {
        return regime::blocked;
    }

    if (!saturating_)
    {
        // the reference behaviour: the choked branch is never taken
        return regime::subsonic;
    }

    return (p <= pCrit()) ? regime::choked : regime::subsonic;
}


Foam::scalar Foam::rdeInjector::T(const scalar p) const
{
    switch (state(p))
    {
        case regime::blocked:
        {
            return T0_;
        }

        case regime::choked:
        {
            // sonic throat
            return T0_*2/(gamma_ + 1);
        }

        case regime::subsonic:
        {
            const scalar r = max(p/p0_, small);
            return T0_*pow(r, (gamma_ - 1)/gamma_);
        }
    }

    return T0_;
}


Foam::scalar Foam::rdeInjector::V(const scalar p) const
{
    switch (state(p))
    {
        case regime::blocked:
        {
            return 0;
        }

        case regime::choked:
        {
            // sonic: the exit speed is the local speed of sound
            return sqrt(gamma_*R_*T(p));
        }

        case regime::subsonic:
        {
            const scalar r = max(p/p0_, small);
            const scalar term = 1 - pow(r, (gamma_ - 1)/gamma_);
            return sqrt(2*gamma_*R_*T0_*max(term, scalar(0))/(gamma_ - 1));
        }
    }

    return 0;
}


void Foam::rdeInjector::report(Ostream& os, const word& patchName) const
{
    const scalar Tc = T0_*2/(gamma_ + 1);

    os  << "RDE injector on patch " << patchName << nl
        << "    plenum          p0 = " << p0_ << " Pa, T0 = " << T0_ << " K" << nl
        << "    critical        p  = " << pCrit() << " Pa ("
        << pCrit()/p0_ << " p0)" << nl
        << "    sonic exit      T  = " << Tc << " K, V = "
        << sqrt(gamma_*R_*Tc) << " m/s" << nl
        << "    slots           " << slotWidth_ << " m wide every "
        << slotSpacing_ << " m (open fraction "
        << slotWidth_/slotSpacing_ << ")" << nl
        << "    choked model    "
        << (saturating_ ? "saturating" : "unbounded") << endl;

    if (!saturating_)
    {
        WarningInFunction
            << "chokedModel 'unbounded' reproduces a reference"
               " implementation in which the choked regime is never"
               " entered.  The injection velocity is then unbounded below"
               " the critical pressure, approaching " << vMax()
            << " m/s against the physical sonic value "
            << sqrt(gamma_*R_*Tc) << " m/s.  Use 'saturating' unless you"
               " are reproducing results from that implementation." << endl;
    }
}


void Foam::rdeInjector::write(Ostream& os) const
{
    os.writeKeyword("p0")          << p0_          << token::END_STATEMENT << nl;
    os.writeKeyword("T0")          << T0_          << token::END_STATEMENT << nl;
    os.writeKeyword("gamma")       << gamma_       << token::END_STATEMENT << nl;
    os.writeKeyword("R")           << R_           << token::END_STATEMENT << nl;
    os.writeKeyword("slotWidth")   << slotWidth_   << token::END_STATEMENT << nl;
    os.writeKeyword("slotSpacing") << slotSpacing_ << token::END_STATEMENT << nl;
    os.writeKeyword("slotOffset")  << slotOffset_  << token::END_STATEMENT << nl;
    os.writeKeyword("tangent")     << tangent_     << token::END_STATEMENT << nl;

    if (mag(direction_) > small)
    {
        os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("chokedModel")
        << word(saturating_ ? "saturating" : "unbounded")
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
