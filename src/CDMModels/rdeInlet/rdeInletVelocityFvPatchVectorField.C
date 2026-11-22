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

#include "rdeInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rdeInletVelocityFvPatchVectorField::rdeInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false),
    injector_(dict)
{
    // mixedFvPatchField is constructed with valuesRequired = false, so
    // refValue, refGradient and valueFraction are sized but not read from the
    // dictionary: the injector computes all three in updateCoeffs before they
    // are used.
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(refValue());
    }
}


Foam::rdeInletVelocityFvPatchVectorField::rdeInletVelocityFvPatchVectorField
(
    const rdeInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    injector_(ptf.injector_)
{}


Foam::rdeInletVelocityFvPatchVectorField::rdeInletVelocityFvPatchVectorField
(
    const rdeInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    injector_(ptf.injector_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rdeInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& Cf = patch().Cf();

    // The inward normal: the default injection direction
    const vectorField nIn(-patch().nf());

    const fvPatchScalarField& pp =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const scalarField pChamber(pp.patchInternalField());

    const bool fixedDirection = (mag(injector_.direction()) > small);

    forAll(Cf, facei)
    {
        const scalar pc = pChamber[facei];

        if
        (
            !injector_.inSlot(injector_.coordinate(Cf[facei]))
         || injector_.state(pc) == rdeInjector::regime::blocked
        )
        {
            // wall: no penetration
            valueFraction()[facei] = 1;
            refValue()[facei] = Zero;
            refGrad()[facei] = Zero;
        }
        else
        {
            const vector dir =
                fixedDirection ? injector_.direction() : nIn[facei];

            valueFraction()[facei] = 1;
            refValue()[facei] = injector_.V(pc)*dir;
            refGrad()[facei] = Zero;
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::rdeInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    injector_.write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rdeInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
