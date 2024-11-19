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
    inletOutletFvPatchScalarField(p, iF),
    ratio_(1.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    ratio_(ptf.ratio_)
{}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    ratio_(readScalar(dict.lookup("ratio")))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& owfpsf
)
:
    inletOutletFvPatchScalarField(owfpsf),
    ratio_(owfpsf.ratio_)
{}


omegaViscosityRatioFvPatchScalarField::omegaViscosityRatioFvPatchScalarField
(
    const omegaViscosityRatioFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(owfpsf, iF),
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
    
    //operator==(kp/(ratio_*nu));
    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    this->refValue() = max(kp/(ratio_*nu), SMALL);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();

}




void omegaViscosityRatioFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#if (OPENFOAM >= 1812)
    os.writeEntry("ratio", ratio_);
    //ratio_.writeEntry("ratio", os);
    this->writeEntry("value", os);
#else
    writeEntry(os, "ratio", ratio_);
    writeEntry(os, "value", *this);
#endif
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
