/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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

#include "roughnessAmplificationLangelFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::roughnessAmplificationLangelFvPatchScalarField::
roughnessAmplificationLangelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    CAr1_(2000.0),
    CAr2_(1.0),
    CAr3_(13.5),
    ks_(p.size(), Zero)
{
}


Foam::roughnessAmplificationLangelFvPatchScalarField::
roughnessAmplificationLangelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    CAr1_(dict.getOrDefault<scalar>("CAr1", 2000.0)),
    CAr2_(dict.getOrDefault<scalar>("CAr2", 1.0)),
    CAr3_(dict.getOrDefault<scalar>("CAr3", 13.5)),
    ks_("ks", dict, p.size())
{
}


Foam::roughnessAmplificationLangelFvPatchScalarField::
roughnessAmplificationLangelFvPatchScalarField
(
    const roughnessAmplificationLangelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    CAr1_(ptf.CAr1_),
    CAr2_(ptf.CAr2_),
    CAr3_(ptf.CAr3_),
    ks_(ptf.ks_, mapper)
{}


Foam::roughnessAmplificationLangelFvPatchScalarField::
roughnessAmplificationLangelFvPatchScalarField
(
    const roughnessAmplificationLangelFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    CAr1_(ptf.CAr1_),
    CAr2_(ptf.CAr2_),
    CAr3_(ptf.CAr3_),
    ks_(ptf.ks_)
{}


Foam::roughnessAmplificationLangelFvPatchScalarField::
roughnessAmplificationLangelFvPatchScalarField
(
    const roughnessAmplificationLangelFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    CAr1_(ptf.CAr1_),
    CAr2_(ptf.CAr2_),
    CAr3_(ptf.CAr3_),
    ks_(ptf.ks_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::roughnessAmplificationLangelFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    ks_.autoMap(m);
}


void Foam::roughnessAmplificationLangelFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const roughnessAmplificationLangelFvPatchScalarField& tiptf =
        refCast<const roughnessAmplificationLangelFvPatchScalarField>(ptf);

    ks_.rmap(tiptf.ks_, addr);
}


void Foam::roughnessAmplificationLangelFvPatchScalarField::updateCoeffs()
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

    const label patchi = patch().index();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField uTau = sqrt(nuw*magUp/y);   // Valid for yPlus < 11 !

    fixedValueFvPatchScalarField::operator==
    (
        CAr1_ / (1 + exp(CAr3_ - CAr2_*uTau*ks_/nuw) )
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::roughnessAmplificationLangelFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("CAr1", CAr1_);
    os.writeEntry("CAr2", CAr2_);
    os.writeEntry("CAr3", CAr3_);
    ks_.writeEntry("ks", os);
    fvPatchScalarField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        roughnessAmplificationLangelFvPatchScalarField
    );
}

// ************************************************************************* //
