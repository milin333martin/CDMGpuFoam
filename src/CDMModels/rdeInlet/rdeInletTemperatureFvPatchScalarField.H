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
    Foam::rdeInletTemperatureFvPatchScalarField

Description
    Temperature at the slot injector of a rotating detonation engine.

    The regime is decided from the chamber pressure in the cell adjoining each
    face (see Foam::rdeInjector). Between the slots, and in a slot the chamber
    has pressurised shut, the patch is an adiabatic wall and the temperature is
    zero-gradient. In a flowing slot the temperature is the plenum total
    temperature expanded isentropically: to the chamber pressure while
    subsonic, and to the fixed sonic value \f$2T_0/(\gamma+1)\f$ once choked.

    Pair with the other three \c rdeInlet* conditions on the same patch;
    all four take the same injector entries.

Usage
    \verbatim
    inlet
    {
        type            rdeInletTemperature;
        p0              5e5;        // plenum total pressure [Pa]
        T0              300;        // plenum total temperature [K]
        gamma           1.125;
        R               365.1;      // Ru/Mw [J/(kg K)]
        slotWidth       2e-3;
        slotSpacing     1e-2;
        value           uniform 300;
    }
    \endverbatim

See also
    Foam::rdeInjector

SourceFiles
    rdeInletTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rdeInletTemperatureFvPatchScalarField_H
#define rdeInletTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "rdeInjector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rdeInletTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rdeInletTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- The injector model
        rdeInjector injector_;


public:

    //- Runtime type information
    TypeName("rdeInletTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        rdeInletTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rdeInletTemperatureFvPatchScalarField
        (
            const rdeInletTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rdeInletTemperatureFvPatchScalarField
        (
            const rdeInletTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rdeInletTemperatureFvPatchScalarField
        (
            const rdeInletTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new rdeInletTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
