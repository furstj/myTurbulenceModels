/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "nutUKnoppWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

#if (OPENFOAM >= 1812)
tmp<scalarField> nutUKnoppWallFunctionFvPatchScalarField::calcNut() const
#else
tmp<scalarField> nutUKnoppWallFunctionFvPatchScalarField::nut() const
#endif
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


tmp<scalarField> nutUKnoppWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();


    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    forAll(uTau, faceI)
    {
        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);
        scalar ut0 = ut;
        if (ut > ROOTVSMALL)
        {
            // First calculate log law solution
            int iter = 0;
            scalar err = GREAT;
            scalar Flog = 0.0;

            do
            {
                scalar f =
                    - ut*y[faceI]/nuw[faceI]        // yPlus (LHS)
                    + (1/kappa)*log(E*(ut*y[faceI]/nuw[faceI])); // RHS

                scalar df =                         // df/du_ut
                    y[faceI]/nuw[faceI]
                    + 1/(kappa*ut);

                scalar uTauNew = ut + f/df;         // Newton iteration
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            Flog = max(0.0, ut);
            //Info << "Flog " << Flog << endl;

            // Then calculate Reichardt law solution
            iter = 0;
            err = GREAT;
            ut = ut0;
            scalar FRei = 0.0;

            do
            {
                scalar f =
                    - magUp[faceI]/ut               // uPlus (LHS)
                    + log(1 + 0.4*(ut*y[faceI]/nuw[faceI]))/kappa
                    + 7.8*(1 - exp(-(ut*y[faceI])/(nuw[faceI]*11.0)) -
                      ((ut*y[faceI])/(nuw[faceI]*11.0))*
                      exp(-ut*y[faceI]/(nuw[faceI]*3.0)));  // RHS

                scalar df =
                    magUp[faceI]/sqr(ut)
                    + (1/kappa)*((0.4*y[faceI])/(nuw[faceI] + 0.4*y[faceI]*ut))
                    + 7.8*((y[faceI]/(nuw[faceI]*11.0))*exp(-y[faceI]*ut/
                    (nuw[faceI]*11.0)) + (y[faceI]/(nuw[faceI]*11.0))*
                    exp(-y[faceI]*ut/(nuw[faceI]*3.0))*
                    (y[faceI]*ut/(nuw[faceI]*3.0) - 1.0));  // df/d_ut

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            FRei = max(0.0, ut);

            //Info << "FRei " << FRei << endl;

            // Reichardt blending
            scalar phib1 = tanh(pow4(y[faceI]*FRei/(nuw[faceI]*27.0)));
            scalar FReim = (1 - phib1)*FRei + phib1*Flog;

            // Finally solution to Spalding's law
            iter = 0;
            err = GREAT;
            ut = ut0;
            scalar FSp = 0.0;

            do
            {
                scalar kUu = min(kappa*magUp[faceI]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[faceI]/nuw[faceI]
                    + magUp[faceI]/ut
                    + 1/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[faceI]/nuw[faceI]
                    + magUp[faceI]/sqr(ut)
                    + 1/E*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            FSp = max(0.0, ut);

            //Info << "FSp "<< FSp << endl;

            // Final blending
            scalar phiko = tanh(sqr(y[faceI]*FRei/(nuw[faceI]*50.0)));
            scalar Fko = (1 - phiko)*FSp + phiko*FReim;

            uTau[faceI] = max(0.0, Fko);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUKnoppWallFunctionFvPatchScalarField::
nutUKnoppWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutUKnoppWallFunctionFvPatchScalarField::
nutUKnoppWallFunctionFvPatchScalarField
(
    const nutUKnoppWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutUKnoppWallFunctionFvPatchScalarField::
nutUKnoppWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutUKnoppWallFunctionFvPatchScalarField::
nutUKnoppWallFunctionFvPatchScalarField
(
    const nutUKnoppWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutUKnoppWallFunctionFvPatchScalarField::
nutUKnoppWallFunctionFvPatchScalarField
(
    const nutUKnoppWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUKnoppWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutUKnoppWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
#if (OPENFOAM >= 1812)
    this->writeEntry("value", os);
#else
    writeEntry(os,"value", *this);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUKnoppWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
