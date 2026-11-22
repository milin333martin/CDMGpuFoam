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
    Foam::rdeInletPressureFvPatchScalarField

Description
    Pressure at the slot injector of a rotating detonation engine.

    The regime is decided from the chamber pressure in the cell adjoining each
    face (see Foam::rdeInjector). Between the slots, and in a slot the chamber
    has pressurised shut, the patch is a wall and the pressure is
    zero-gradient. In a slot flowing subsonically the chamber still sets the
    pressure, so it is again zero-gradient. Only when the slot chokes does the
    injector impose a pressure of its own, the critical pressure, because the
    sonic throat stops the chamber communicating upstream.

    Pair with \c rdeInletTemperature, \c rdeInletVelocity and \c rdeInletY on
    the same patch; all four take the same injector entries.

Usage
    \verbatim
    inlet
    {
        type            rdeInletPressure;
        p0              5e5;        // plenum total pressure [Pa]
        T0              300;        // plenum total temperature [K]
        gamma           1.125;
        R               365.1;      // Ru/Mw [J/(kg K)]
        slotWidth       2e-3;
        slotSpacing     1e-2;
        value           uniform 101325;
    }
    \endverbatim

See also
    Foam::rdeInjector

SourceFiles
    rdeInletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rdeInletPressureFvPatchScalarField_H
#define rdeInletPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "rdeInjector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rdeInletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rdeInletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- The injector model
        rdeInjector injector_;


public:

    //- Runtime type information
    TypeName("rdeInletPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        rdeInletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rdeInletPressureFvPatchScalarField
        (
            const rdeInletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rdeInletPressureFvPatchScalarField
        (
            const rdeInletPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rdeInletPressureFvPatchScalarField
        (
            const rdeInletPressureFvPatchScalarField&,
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
                new rdeInletPressureFvPatchScalarField(*this, iF)
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
