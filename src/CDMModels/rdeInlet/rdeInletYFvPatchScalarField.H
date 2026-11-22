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
    Foam::rdeInletYFvPatchScalarField

Description
    Reactant mass fraction at the slot injector of a rotating
    detonation engine.

    The regime is decided from the chamber pressure in the cell adjoining each
    face (see Foam::rdeInjector). A flowing slot, choked or subsonic, carries
    the injected mixture at \c Yinj. Between the slots, and in a slot the
    chamber has pressurised shut, the patch is a wall and the mass fraction is
    zero-gradient.

    Pair with the other three \c rdeInlet* conditions on the same patch;
    all four take the same injector entries.

Usage
    \verbatim
    inlet
    {
        type            rdeInletY;
        p0              5e5;        // plenum total pressure [Pa]
        T0              300;        // plenum total temperature [K]
        gamma           1.125;
        R               365.1;      // Ru/Mw [J/(kg K)]
        slotWidth       2e-3;
        slotSpacing     1e-2;
        Yinj            1;          // reactant mass fraction injected
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::rdeInjector

SourceFiles
    rdeInletYFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rdeInletYFvPatchScalarField_H
#define rdeInletYFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "rdeInjector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rdeInletYFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rdeInletYFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- The injector model
        rdeInjector injector_;

        //- Reactant mass fraction of the injected mixture
        scalar Yinj_;


public:

    //- Runtime type information
    TypeName("rdeInletY");


    // Constructors

        //- Construct from patch, internal field and dictionary
        rdeInletYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rdeInletYFvPatchScalarField
        (
            const rdeInletYFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rdeInletYFvPatchScalarField
        (
            const rdeInletYFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rdeInletYFvPatchScalarField
        (
            const rdeInletYFvPatchScalarField&,
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
                new rdeInletYFvPatchScalarField(*this, iF)
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
