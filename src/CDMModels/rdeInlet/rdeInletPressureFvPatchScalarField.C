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

#include "rdeInletPressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rdeInletPressureFvPatchScalarField::rdeInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    injector_(dict)
{
    // mixedFvPatchField is constructed with valuesRequired = false, so
    // refValue, refGradient and valueFraction are sized but not read from the
    // dictionary: the injector computes all three in updateCoeffs before they
    // are used.
    refValue() = injector_.pCrit();
    refGrad() = Zero;
    valueFraction() = 0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }

    injector_.report(Info, p.name());
}


Foam::rdeInletPressureFvPatchScalarField::rdeInletPressureFvPatchScalarField
(
    const rdeInletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    injector_(ptf.injector_)
{}


Foam::rdeInletPressureFvPatchScalarField::rdeInletPressureFvPatchScalarField
(
    const rdeInletPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    injector_(ptf.injector_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rdeInletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& Cf = patch().Cf();

    // The chamber pressure is the value in the adjoining cell, not the value
    // currently on the patch
    const scalarField pChamber(patchInternalField());

    const scalar pc = injector_.pCrit();

    forAll(Cf, facei)
    {
        if (!injector_.inSlot(injector_.coordinate(Cf[facei])))
        {
            // between the slots: wall
            valueFraction()[facei] = 0;
            refGrad()[facei] = 0;
            continue;
        }

        if (injector_.state(pChamber[facei]) == rdeInjector::regime::choked)
        {
            // sonic throat: the injector sets the pressure
            valueFraction()[facei] = 1;
            refValue()[facei] = pc;
            refGrad()[facei] = 0;
        }
        else
        {
            // blocked or subsonic: the chamber sets the pressure
            valueFraction()[facei] = 0;
            refGrad()[facei] = 0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::rdeInletPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    injector_.write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rdeInletPressureFvPatchScalarField
    );
}


// ************************************************************************* //
