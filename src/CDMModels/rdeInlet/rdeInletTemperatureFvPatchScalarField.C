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

#include "rdeInletTemperatureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rdeInletTemperatureFvPatchScalarField::
rdeInletTemperatureFvPatchScalarField
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
    refValue() = injector_.T0();
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
}


Foam::rdeInletTemperatureFvPatchScalarField::
rdeInletTemperatureFvPatchScalarField
(
    const rdeInletTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    injector_(ptf.injector_)
{}


Foam::rdeInletTemperatureFvPatchScalarField::
rdeInletTemperatureFvPatchScalarField
(
    const rdeInletTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    injector_(ptf.injector_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rdeInletTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& Cf = patch().Cf();

    // The chamber pressure is the value in the adjoining cell
    const fvPatchScalarField& pp =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const scalarField pChamber(pp.patchInternalField());

    forAll(Cf, facei)
    {
        const scalar pc = pChamber[facei];

        if
        (
            !injector_.inSlot(injector_.coordinate(Cf[facei]))
         || injector_.state(pc) == rdeInjector::regime::blocked
        )
        {
            // wall: adiabatic
            valueFraction()[facei] = 0;
            refGrad()[facei] = 0;
        }
        else
        {
            // flowing: isentropic expansion from the plenum
            valueFraction()[facei] = 1;
            refValue()[facei] = injector_.T(pc);
            refGrad()[facei] = 0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::rdeInletTemperatureFvPatchScalarField::write(Ostream& os) const
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
        rdeInletTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
