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

#include "roughnessAmplificationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::roughnessAmplificationFvPatchScalarField::
roughnessAmplificationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    cr1_(8.0),
    ks_(p.size(), Zero)
{
}


Foam::roughnessAmplificationFvPatchScalarField::
roughnessAmplificationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    cr1_(dict.getOrDefault<scalar>("cr1", 8.0)),
    ks_("ks", dict, p.size())
{
}


Foam::roughnessAmplificationFvPatchScalarField::
roughnessAmplificationFvPatchScalarField
(
    const roughnessAmplificationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    cr1_(ptf.cr1_),
    ks_(ptf.ks_, mapper)
{}


Foam::roughnessAmplificationFvPatchScalarField::
roughnessAmplificationFvPatchScalarField
(
    const roughnessAmplificationFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    cr1_(ptf.cr1_),
    ks_(ptf.ks_)
{}


Foam::roughnessAmplificationFvPatchScalarField::
roughnessAmplificationFvPatchScalarField
(
    const roughnessAmplificationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    cr1_(ptf.cr1_),
    ks_(ptf.ks_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::roughnessAmplificationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    ks_.autoMap(m);
}


void Foam::roughnessAmplificationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const roughnessAmplificationFvPatchScalarField& tiptf =
        refCast<const roughnessAmplificationFvPatchScalarField>(ptf);

    ks_.rmap(tiptf.ks_, addr);
}


void Foam::roughnessAmplificationFvPatchScalarField::updateCoeffs()
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
        cr1_*uTau*ks_/nuw 
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::roughnessAmplificationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("cr1", cr1_);
    ks_.writeEntry("ks", os);
    fvPatchScalarField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        roughnessAmplificationFvPatchScalarField
    );
}

// ************************************************************************* //
