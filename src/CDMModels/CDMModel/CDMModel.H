/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM.

    CDMGpuFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CDMGpuFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CDMGpuFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::CDM::CDMModel

Description
    Chemical-diffusive model (CDM) of Kaplan, Oezgen & Oran (2019) as used by
    Lu, Kaplan & Oran, Combust. Flame 230 (2021) 111417.

    Single irreversible reaction R -> P of a calorically perfect gas with
    identical gamma and molecular weight for reactant and product:

    \verbatim
        dY/dt   = -A rho Y exp(-Ta/T)          (Y: reactant mass fraction)
        Qdot    = -q rho dY/dt                 (heat release, W/m^3)

        K       = kappa0 Cp T^n                (thermal conductivity)
        mu      = Pr K/Cp   = Pr kappa0 T^n
        rho D   = K/(Cp Le) = kappa0 T^n / Le
    \endverbatim

    The reaction is integrated per cell at constant volume and constant
    velocity with an exponential (positivity preserving, exactly energy
    conserving) sub-stepped update in which the temperature is updated with
    the heat released in each sub-step.

    Dictionary constant/CDMProperties.  Parameters may be given either in SI:

    \verbatim
        A        1.064e10;     // m^3/(kg s)
        Ta       16375.1;      // K           (Ea/R)
        q        4.7154e6;     // J/kg
        kappa0   2.367e-6;     // kg/(m s K^n)
    \endverbatim

    or exactly as printed in the tables of Lu et al. (2021):

    \verbatim
        A_cgs       1.064e13;  // cm^3/(g s)
        EaByRT0     54.95;
        qMwByRT0    46.36;
        kappa0_cgs  2.367e-5;  // g/(s cm K^n)
        T0          298;
    \endverbatim

    in which case the specific gas constant (Cp - Cv) of the thermo is used to
    convert q.  Common optional entries:

    \verbatim
        n                   0.7;    // transport temperature exponent
        Pr                  1;
        Le                  1;
        maxSubstepFraction  0.2;    // max k*h per chemistry sub-step
        maxSubsteps         1000;
        cpuLoad             no;     // yes: per-cell chemistry cost for the
                                    // loadBalancer distributor
    \endverbatim

SourceFiles
    CDMModel.C

\*---------------------------------------------------------------------------*/

#ifndef CDMModel_H
#define CDMModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CDM
{

/*---------------------------------------------------------------------------*\
                          Class CDMModel Declaration
\*---------------------------------------------------------------------------*/

class CDMModel
:
    public IOdictionary
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Specific gas constant Cp - Cv [J/(kg K)]
        const scalar Rspecific_;

        //- Record the chemistry CPU time per cell for load balancing
        Switch cpuLoad_;

        //- Pre-exponential factor [m^3/(kg s)]
        scalar A_;

        //- Activation temperature Ea/R [K]
        scalar Ta_;

        //- Heat release per unit mass of reactant [J/kg]
        scalar q_;

        //- Transport constant [kg/(m s K^n)]
        scalar kappa0_;

        //- Transport temperature exponent
        scalar n_;

        //- Prandtl number
        scalar Pr_;

        //- Lewis number
        scalar Le_;

        //- Maximum k*h per chemistry sub-step
        scalar maxSubstepFraction_;

        //- Maximum number of chemistry sub-steps per call
        label maxSubsteps_;


public:

    //- Runtime type information
    TypeName("CDMModel");


    // Constructors

        //- Construct from mesh and the specific gas constant Cp - Cv
        //  (needed only to convert q from the non-dimensional form)
        CDMModel(const fvMesh& mesh, const scalar Rspecific);

        //- Disallow default bitwise copy construction
        CDMModel(const CDMModel&) = delete;


    //- Destructor
    virtual ~CDMModel();


    // Member Functions

        // Access

            scalar A() const
            {
                return A_;
            }

            scalar Ta() const
            {
                return Ta_;
            }

            scalar q() const
            {
                return q_;
            }

            scalar kappa0() const
            {
                return kappa0_;
            }

            scalar n() const
            {
                return n_;
            }

            scalar Pr() const
            {
                return Pr_;
            }

            scalar Le() const
            {
                return Le_;
            }


        // Transport

            //- Return kappa0*T^n = K/Cp  [kg/(m s)]
            tmp<volScalarField> kappaByCp(const volScalarField& T) const;

            //- Return the largest diffusivity-like coefficient factor
            //  max(Pr, 1, 1/Le) used for the explicit diffusion limit
            scalar maxDiffusivityFactor() const
            {
                return max(max(Pr_, scalar(1)), 1/Le_);
            }


        // Chemistry

            //- Largest k*h permitted in a chemistry sub-step
            scalar maxSubstepFraction() const
            {
                return maxSubstepFraction_;
            }

            //- Sub-step ceiling, after which the rate is frozen
            label maxSubsteps() const
            {
                return maxSubsteps_;
            }

            //- Reaction rate constant k = A rho exp(-Ta/T) [1/s]
            inline scalar k(const scalar rho, const scalar T) const
            {
                return A_*rho*exp(-Ta_/T);
            }

            //- Integrate the reaction over deltaT in every cell at constant
            //  density and velocity.  Updates rhoY and rhoE and adds the
            //  energy released per unit volume to heatReleased.
            void react
            (
                const scalar deltaT,
                const scalarField& rho,
                const scalarField& T,
                const scalarField& Cv,
                scalarField& rhoY,
                scalarField& rhoE,
                scalarField& heatReleased
            ) const;


        //- Print a summary of the model parameters
        void printSummary() const;

        //- Read
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CDMModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CDM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
