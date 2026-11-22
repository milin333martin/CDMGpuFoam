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

Class
    Foam::rdeInjector

Description
    The plenum-fed slot injector of a rotating detonation engine, shared by
    the four \c rdeInlet* patch fields so that the regime logic exists once
    rather than four times.

    The injector face is an array of slots of width \c slotWidth repeating
    every \c slotSpacing along the patch's tangential coordinate; the gaps
    between them are solid wall. Each slot is fed from a plenum held at total
    pressure \c p0 and total temperature \c T0, and discharges isentropically
    into whatever chamber pressure \f$p\f$ the solution presents to it. Three
    regimes follow, from the critical pressure

    \f[
        p_{cr} = p_0 \left(\frac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}
    \f]

    \table
        Regime      | Condition          | Behaviour \\
        blocked     | \f$p \ge p_0\f$    | no injection: the slot is a wall \\
        subsonic    | \f$p_{cr} < p < p_0\f$ | isentropic discharge to the
                                           chamber pressure \\
        choked      | \f$p \le p_{cr}\f$ | the throat is sonic: the exit state
                                           saturates and stops responding to
                                           the chamber \\
    \endtable

    In the subsonic regime the exit state follows from isentropic expansion to
    the chamber pressure,

    \f[
        T = T_0 \left(\frac{p}{p_0}\right)^{(\gamma-1)/\gamma}, \qquad
        V = \sqrt{\frac{2\gamma}{\gamma-1} R T_0
                  \left[1 - \left(\frac{p}{p_0}\right)^{(\gamma-1)/\gamma}
                  \right]} .
    \f]

    In the choked regime the throat is sonic and the exit state is fixed at
    the critical state, independent of the chamber:

    \f[
        p = p_{cr}, \qquad
        T = T_0 \frac{2}{\gamma+1}, \qquad
        V = \sqrt{\gamma R T} .
    \f]

    \par The \c chokedModel switch

    Setting \c chokedModel to \c unbounded reverts to the behaviour of the
    reference implementation this was ported from, in which the subsonic
    relations are applied at all pressures below \f$p_0\f$ and the choked
    branch is never taken. That lets the injection velocity grow towards the
    vacuum-discharge limit \f$\sqrt{2\gamma R T_0/(\gamma-1)}\f$, which for a
    typical plenum is several times the sonic value, and it is reached exactly
    in the low-pressure fill region behind the detonation. It is provided only
    for reproducing results computed with that implementation and should not
    be used otherwise.

Usage
    \table
        Property     | Description                        | Required | Default
        p0           | plenum total pressure [Pa]         | yes      |
        T0           | plenum total temperature [K]        | yes      |
        gamma        | ratio of specific heats             | yes      |
        R            | specific gas constant [J/(kg K)]    | yes      |
        slotWidth    | width of one injector slot [m]      | yes      |
        slotSpacing  | slot-to-slot pitch [m]              | yes      |
        slotOffset   | position of the first slot [m]      | no       | 0
        tangent      | direction the slot pattern runs in  | no       | (1 0 0)
        direction    | injection direction                 | no       | patch normal
        chokedModel  | saturating \| unbounded             | no       | saturating
    \endtable

SourceFiles
    rdeInjector.C

\*---------------------------------------------------------------------------*/

#ifndef rdeInjector_H
#define rdeInjector_H

#include "dictionary.H"
#include "scalarField.H"
#include "vectorField.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rdeInjector Declaration
\*---------------------------------------------------------------------------*/

class rdeInjector
{
public:

    //- The discharge regime at a face
    enum class regime
    {
        blocked,        //- p >= p0: no injection
        subsonic,       //- pCrit < p < p0
        choked          //- p <= pCrit
    };


private:

    // Private Data

        //- Plenum total pressure [Pa]
        scalar p0_;

        //- Plenum total temperature [K]
        scalar T0_;

        //- Ratio of specific heats
        scalar gamma_;

        //- Specific gas constant [J/(kg K)]
        scalar R_;

        //- Width of one slot [m]
        scalar slotWidth_;

        //- Slot-to-slot pitch [m]
        scalar slotSpacing_;

        //- Tangential position of the first slot [m]
        scalar slotOffset_;

        //- Direction the slot pattern repeats along
        vector tangent_;

        //- Injection direction; zero means "use the inward patch normal"
        vector direction_;

        //- Saturate the exit state at the sonic condition once choked.
        //  False reproduces the reference implementation, in which the
        //  subsonic relations are applied at every pressure below p0.
        bool saturating_;


public:

    // Constructors

        //- Construct from a dictionary
        explicit rdeInjector(const dictionary& dict);


    // Member Functions

        // Access

            inline scalar p0() const { return p0_; }
            inline scalar T0() const { return T0_; }
            inline scalar gamma() const { return gamma_; }
            inline scalar R() const { return R_; }
            inline const vector& direction() const { return direction_; }

        //- The critical (choking) pressure
        scalar pCrit() const;

        //- The vacuum-discharge velocity: the bound the unbounded model
        //  approaches, reported so that it can be compared with the sonic value
        scalar vMax() const;

        //- The tangential coordinate of a face centre: its distance along
        //  \c tangent, measured from \c slotOffset
        inline scalar coordinate(const vector& Cf) const
        {
            return (Cf & tangent_) - slotOffset_;
        }

        //- Is the tangential coordinate x inside a slot?
        bool inSlot(const scalar x) const;

        //- The discharge regime at chamber pressure p
        regime state(const scalar p) const;

        //- Static temperature of the injected gas at chamber pressure p
        scalar T(const scalar p) const;

        //- Speed of the injected gas at chamber pressure p
        scalar V(const scalar p) const;

        //- Static pressure the injector imposes when choked
        //  (meaningless in the other regimes)
        inline scalar pChoked() const { return pCrit(); }

        //- Report the model to a stream, once, at construction
        void report(Ostream& os, const word& patchName) const;

        //- Write the dictionary entries
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
