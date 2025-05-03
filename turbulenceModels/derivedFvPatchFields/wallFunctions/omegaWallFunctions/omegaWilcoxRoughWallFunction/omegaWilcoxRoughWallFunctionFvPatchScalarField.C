/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

#include "omegaWilcoxRoughWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::
omegaWilcoxRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    beta1_(0.075),
    ks_(p.size(), Zero)
{
}


Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::
omegaWilcoxRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    ks_("ks", dict, p.size())
{
}


Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::
omegaWilcoxRoughWallFunctionFvPatchScalarField
(
    const omegaWilcoxRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    beta1_(ptf.beta1_),
    ks_(ptf.ks_, mapper)
{}


Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::
omegaWilcoxRoughWallFunctionFvPatchScalarField
(
    const omegaWilcoxRoughWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    beta1_(ptf.beta1_),
    ks_(ptf.ks_)
{}


Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::
omegaWilcoxRoughWallFunctionFvPatchScalarField
(
    const omegaWilcoxRoughWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    beta1_(ptf.beta1_),
    ks_(ptf.ks_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    ks_.autoMap(m);
}


void Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const omegaWilcoxRoughWallFunctionFvPatchScalarField& tiptf =
        refCast<const omegaWilcoxRoughWallFunctionFvPatchScalarField>(ptf);

    ks_.rmap(tiptf.ks_, addr);
}


void Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::updateCoeffs()
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

    const scalarField ksPlus = ks_*uTau/nuw;

    scalarField& omega = *this;

    scalar smoothWall = sqrt(2500*beta1_/60);
    
    forAll(omega, facei)
    {
        if (ks_[facei] < smoothWall*y[facei])
        {
            omega[facei] = 10*6*nuw[facei]/(beta1_*sqr(y[facei]));
        }
        else if (ksPlus[facei] <= 25)
        {
            //scalar Sr = sqr(50/ksPlus[facei]);
            //omega[facei] = sqr(uTau[facei])*Sr/nuw[facei];
            omega[facei] = 2500*nuw[facei]/sqr(ks_[facei]);
        }
        else
        {
            //scalar Sr = 100/ksPlus[facei];
            //omega[facei] = sqr(uTau[facei])*Sr/nuw[facei];
            omega[facei] = 100*uTau[facei]/ks_[facei];
        }
    }
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("beta1", beta1_);
    ks_.writeEntry("ks", os);
    fvPatchScalarField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaWilcoxRoughWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
