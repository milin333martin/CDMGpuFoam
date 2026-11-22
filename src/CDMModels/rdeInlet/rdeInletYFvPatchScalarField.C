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

#include "rdeInletYFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rdeInletYFvPatchScalarField::rdeInletYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    injector_(dict),
    Yinj_(dict.lookupOrDefault<scalar>("Yinj", 1))
{
    if (Yinj_ < 0 || Yinj_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "Yinj must lie in [0, 1]: " << Yinj_
            << exit(FatalIOError);
    }

    // mixedFvPatchField is constructed with valuesRequired = false, so
    // refValue, refGradient and valueFraction are sized but not read from the
    // dictionary: the injector computes all three in updateCoeffs before they
    // are used.
    refValue() = Yinj_;
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


Foam::rdeInletYFvPatchScalarField::rdeInletYFvPatchScalarField
(
    const rdeInletYFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    injector_(ptf.injector_),
    Yinj_(ptf.Yinj_)
{}


Foam::rdeInletYFvPatchScalarField::rdeInletYFvPatchScalarField
(
    const rdeInletYFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    injector_(ptf.injector_),
    Yinj_(ptf.Yinj_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rdeInletYFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& Cf = patch().Cf();

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
            // wall: nothing is being injected, so let the chamber decide
            valueFraction()[facei] = 0;
            refGrad()[facei] = 0;
        }
        else
        {
            // flowing: the injected mixture
            valueFraction()[facei] = 1;
            refValue()[facei] = Yinj_;
            refGrad()[facei] = 0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::rdeInletYFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    injector_.write(os);
    writeEntry(os, "Yinj", Yinj_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rdeInletYFvPatchScalarField
    );
}


// ************************************************************************* //
