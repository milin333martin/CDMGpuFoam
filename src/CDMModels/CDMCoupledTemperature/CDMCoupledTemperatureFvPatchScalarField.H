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

Class
    Foam::CDMCoupledTemperatureFvPatchScalarField

Description
    coupledTemperature for the fluid side of a conjugate heat transfer case
    run with CDMGpuFoam.

    Identical to coupledTemperature (the same mixed condition, the same
    mapping to the neighbouring region), except that the patch thermal
    conductivity is the chemical-diffusive one,

        K = kappa0 Cp T^n,

    taken from the field the solver registers (K in the mesh database)
    rather than from a thermophysicalTransportModel, which CDMGpuFoam does
    not construct.  The neighbouring (solid) side uses the standard
    coupledTemperature.

Usage
    \verbatim
    <patchName>
    {
        type            CDMCoupledTemperature;
        value           uniform 298;
    }
    \endverbatim
    on a mapped(Wall) patch, with

    \verbatim
    PIMPLE { CHT yes; }
    \endverbatim

    in the fluid region's fvSolution so that the solver registers K.

SourceFiles
    CDMCoupledTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CDMCoupledTemperatureFvPatchScalarField_H
#define CDMCoupledTemperatureFvPatchScalarField_H

#include "coupledTemperatureFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class CDMCoupledTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CDMCoupledTemperatureFvPatchScalarField
:
    public coupledTemperatureFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- The CDM thermal conductivity on this patch
        tmp<scalarField> kappaCDM() const;

        //- As coupledTemperature, with the CDM conductivity
        virtual void getThis
        (
            tmp<scalarField>& kappa,
            tmp<scalarField>& sumKappaTByDelta,
            tmp<scalarField>& sumKappaByDelta,
            scalarField& qTot,
            tmp<scalarField>& qByKappa
        ) const;

        //- As coupledTemperature, with the CDM conductivity
        virtual void getNbr
        (
            tmp<scalarField>& sumKappaTByDeltaNbr,
            tmp<scalarField>& sumKappaByDeltaNbr,
            tmp<scalarField>& qNbr
        ) const;

        //- As coupledTemperature (wall resistance form)
        virtual void getNbr
        (
            tmp<scalarField>& TwNbr,
            tmp<scalarField>& qNbr
        ) const;


public:

    //- Runtime type information
    TypeName("CDMCoupledTemperature");


    // Constructors

        //- Construct from patch, internal field and dictionary
        CDMCoupledTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        CDMCoupledTemperatureFvPatchScalarField
        (
            const CDMCoupledTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        CDMCoupledTemperatureFvPatchScalarField
        (
            const CDMCoupledTemperatureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        CDMCoupledTemperatureFvPatchScalarField
        (
            const CDMCoupledTemperatureFvPatchScalarField&,
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
                new CDMCoupledTemperatureFvPatchScalarField(*this, iF)
            );
        }
};


} // End namespace Foam

#endif

// ************************************************************************* //
