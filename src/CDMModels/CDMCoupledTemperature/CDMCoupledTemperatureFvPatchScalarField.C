/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam.
    GNU General Public License v3 or later; see COPYING.

\*---------------------------------------------------------------------------*/

#include "CDMCoupledTemperatureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::CDMCoupledTemperatureFvPatchScalarField::kappaCDM() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>("K"))
    {
        FatalErrorInFunction
            << "Field K (the CDM thermal conductivity) is not registered in "
            << "region " << mesh.name() << nl
            << "    Set" << nl
            << "        PIMPLE { CHT yes; }" << nl
            << "    in the fvSolution of this region so that CDMGpuFoam "
            << "publishes it." << exit(FatalError);
    }

    const volScalarField& K = mesh.lookupObject<volScalarField>("K");

    return tmp<scalarField>(new scalarField(K.boundaryField()[patch().index()]));
}


void Foam::CDMCoupledTemperatureFvPatchScalarField::getThis
(
    tmp<scalarField>& kappa,
    tmp<scalarField>& sumKappaTByDelta,
    tmp<scalarField>& sumKappaByDelta,
    scalarField& sumq,
    tmp<scalarField>& qByKappa
) const
{
    kappa = kappaCDM();

    qByKappa = sumq/kappa();

    sumq = 0;
}


void Foam::CDMCoupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& sumKappaTByDeltaNbr,
    tmp<scalarField>& sumKappaByDeltaNbr,
    tmp<scalarField>& qNbr
) const
{
    sumKappaByDeltaNbr = kappaCDM()*patch().deltaCoeffs();

    sumKappaTByDeltaNbr = sumKappaByDeltaNbr()*patchInternalField();
}


void Foam::CDMCoupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& TrefNbr,
    tmp<scalarField>& qNbr
) const
{
    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );

    TrefNbr = tmp<scalarField>(new scalarField(Tp));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDMCoupledTemperatureFvPatchScalarField::
CDMCoupledTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledTemperatureFvPatchScalarField(p, iF, dict)
{}


Foam::CDMCoupledTemperatureFvPatchScalarField::
CDMCoupledTemperatureFvPatchScalarField
(
    const CDMCoupledTemperatureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    coupledTemperatureFvPatchScalarField(psf, p, iF, mapper)
{}


Foam::CDMCoupledTemperatureFvPatchScalarField::
CDMCoupledTemperatureFvPatchScalarField
(
    const CDMCoupledTemperatureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledTemperatureFvPatchScalarField(psf, iF)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        CDMCoupledTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
