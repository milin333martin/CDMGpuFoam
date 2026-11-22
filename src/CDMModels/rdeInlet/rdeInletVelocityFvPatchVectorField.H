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
    Foam::rdeInletVelocityFvPatchVectorField

Description
    Velocity at the slot injector of a rotating detonation engine.

    The regime is decided from the chamber pressure in the cell adjoining each
    face (see Foam::rdeInjector). Between the slots, and in a slot the chamber
    has pressurised shut, the patch is a solid wall and the velocity is fixed
    at zero. In a flowing slot the speed follows from isentropic discharge:
    growing as the chamber pressure falls while subsonic, then fixed at the
    sonic speed once choked.

    Injection is along the inward patch normal unless \c direction is given.

    Pair with the other three \c rdeInlet* conditions on the same patch;
    all four take the same injector entries.

Usage
    \verbatim
    inlet
    {
        type            rdeInletVelocity;
        p0              5e5;        // plenum total pressure [Pa]
        T0              300;        // plenum total temperature [K]
        gamma           1.125;
        R               365.1;      // Ru/Mw [J/(kg K)]
        slotWidth       2e-3;
        slotSpacing     1e-2;
        direction       (0 1 0);   // optional: default is the inward normal
        value           uniform (0 0 0);
    }
    \endverbatim

See also
    Foam::rdeInjector

SourceFiles
    rdeInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rdeInletVelocityFvPatchVectorField_H
#define rdeInletVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"
#include "rdeInjector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rdeInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rdeInletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private Data

        //- The injector model
        rdeInjector injector_;


public:

    //- Runtime type information
    TypeName("rdeInletVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        rdeInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rdeInletVelocityFvPatchVectorField
        (
            const rdeInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rdeInletVelocityFvPatchVectorField
        (
            const rdeInletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rdeInletVelocityFvPatchVectorField
        (
            const rdeInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new rdeInletVelocityFvPatchVectorField(*this, iF)
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
