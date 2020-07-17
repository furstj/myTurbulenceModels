/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Felix Langfeldt
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "ReThetatInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ReThetatInletFvPatchScalarField::ReThetatInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

ReThetatInletFvPatchScalarField::ReThetatInletFvPatchScalarField
(
    const ReThetatInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_)
{}

ReThetatInletFvPatchScalarField::ReThetatInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

ReThetatInletFvPatchScalarField::ReThetatInletFvPatchScalarField
(
    const ReThetatInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    kName_(ptf.kName_)
{}


ReThetatInletFvPatchScalarField::ReThetatInletFvPatchScalarField
(
    const ReThetatInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    kName_(ptf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar ReThetatInletFvPatchScalarField::ReThetatInlet(scalar Tu) const
{
    if (Tu<1.3)
    {
        return 1173.51 - 589.428*Tu + 0.2196/sqr(Tu);
    }
    else
    {
        return 331.50*pow(Tu - 0.5658, -0.671);
    }
}

void ReThetatInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    forAll(kp, faceI)
    {
        scalar Tu = scalar(100.0)*sqrt(kp[faceI]/scalar(1.5))/mag(Up[faceI]);
        this->refValue()[faceI] = ReThetatInlet(Tu);
    }

    this->valueFraction() = 1.0 - pos(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void ReThetatInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
#if (OPENFOAM_PLUS > 1712  || OPENFOAM >= 1912)
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("k", "k", kName_);
    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
#else
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
#endif
#if (OPENFOAM >= 1912)
    this->writeEntry("value", os);
#else
    writeEntry(os, "value", *this);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    ReThetatInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
