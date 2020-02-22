/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "omegaViscosityRatioFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchField<scalar>(p, iF),
    ratio_(p.size(), 1.0)
{}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    ratio_(ptf.ratio_)
{}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    ratio_("ratio", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
            (
                scalarField("value", dict, p.size())
            );
    }
    else
    {
        this->operator==(patchInternalField());
    }
}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    ratio_(owfpsf.ratio_)
{}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    ratio_(owfpsf.ratio_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void omegaViscosityRatioFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const label patchi = this->patch().index();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kp = k.boundaryField()[patchi];
    
    const tmp<scalarField> tnu = turbModel.nu(patchi);
    const scalarField& nu = tnu();

    operator==(kp/(ratio_*nu));
}




void omegaViscosityRatioFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "ratio", ratio_);
    writeEntry(os, "value", *this);
}

void omegaViscosityRatioFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<scalar>::autoMap(m);
    m(ratio_, ratio_);
}


void omegaViscosityRatioFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<scalar>::rmap(ptf, addr);

    const omegaViscosityRatioFvPatchScalarField& tiptf =
        refCast<const omegaViscosityRatioFvPatchScalarField>(ptf);

    ratio_.rmap(tiptf.ratio_, addr);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaViscosityRatioFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
