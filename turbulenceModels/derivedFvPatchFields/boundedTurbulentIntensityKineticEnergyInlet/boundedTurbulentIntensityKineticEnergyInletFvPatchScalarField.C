/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    intensity_(0.0),
    k_min_(0.0),
    k_max_(GREAT),
    UName_("U")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    intensity_(ptf.intensity_),
    k_min_(ptf.k_min_),
    k_max_(ptf.k_max_),
    UName_(ptf.UName_)
{}

Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    intensity_(dict.get<scalar>("intensity")),
    k_min_(dict.getOrDefault<scalar>("k_min", 0.0)),
    k_max_(dict.getOrDefault<scalar>("k_max", GREAT)),
    UName_(dict.getOrDefault<word>("U", "U"))
{
    fvPatchFieldBase::readDict(dict);
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    if (intensity_ < 0 || intensity_ > 1)
    {
        FatalErrorInFunction
            << "Turbulence intensity should be specified as a fraction 0-1 "
               "of the mean velocity\n"
               "    value given is " << intensity_ << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    intensity_(ptf.intensity_),
    k_min_(ptf.k_min_),
    k_max_(ptf.k_max_),
    UName_(ptf.UName_)
{}


Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    intensity_(ptf.intensity_),
    k_min_(ptf.k_min_),
    k_max_(ptf.k_max_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& Up =
        patch().lookupPatchField<volVectorField>(UName_);

    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    this->refValue() = max( min(1.5*sqr(intensity_)*magSqr(Up), k_max_), k_min_);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("intensity", intensity_);
    os.writeEntryIfDifferent<scalar>("k_min", k_min_, 0.0);
    os.writeEntryIfDifferent<scalar>("k_max", k_max_, GREAT);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
