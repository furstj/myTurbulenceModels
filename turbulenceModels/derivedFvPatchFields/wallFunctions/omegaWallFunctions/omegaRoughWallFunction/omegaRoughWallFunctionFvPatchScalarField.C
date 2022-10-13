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

#include "omegaRoughWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaRoughWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaRoughWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void omegaRoughWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<scalar>("beta1", 0.075, beta1_);
    wallCoeffs_.writeEntries(os);
    roughnessHeight_.writeEntry("roughnessHeight", os);
}


void omegaRoughWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaRoughWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label celli = faceCells[i];
                weights[celli]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


omegaRoughWallFunctionFvPatchScalarField&
omegaRoughWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaRoughWallFunctionFvPatchScalarField& opf =
        refCast<const omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaRoughWallFunctionFvPatchScalarField&>(opf);
}


void omegaRoughWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void omegaRoughWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();
    
    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

	const scalar uTau = Cmu25 * sqrt(k[celli]);
	
	const scalar ksPlus = roughnessHeight_[facei] * uTau / nuw[facei];
        
        scalar omegaVis;
        
	if (ksPlus > 25.0)
        {
            const scalar SR = 100.0/ksPlus;
            omegaVis = sqr(uTau)*SR/nuw[facei];
        }
        else if (ksPlus > 5.0)
        {
            const scalar ksPlusMin = min(4.3*pow(yPlus,0.85), 8.0);
            const scalar SR = sqr(50/max(ksPlus,ksPlusMin));
            omegaVis = sqr(uTau)*SR/nuw[facei];
        }
        else
        {
            omegaVis = 6.0*nuw[facei]/(beta1_*sqr(y[facei]));
        }
        
        const scalar omegaLog = sqrt(k[celli])/(Cmu25*kappa*y[facei]);

        switch (blender_)
        {
            case blenderType::STEPWISE:
            {
                if (yPlus > yPlusLam)
                {
                    omega0[celli] += w*omegaLog;
                }
                else
                {
                    omega0[celli] += w*omegaVis;
                }
                break;
            }

            case blenderType::BINOMIAL:
            {
                omega0[celli] +=
                    w*pow
                    (
                        pow(omegaVis, n_) + pow(omegaLog, n_),
                        scalar(1)/n_
                    );
                break;
            }

            case blenderType::MAX:
            {
                // (PH:Eq. 27)
                omega0[celli] += max(omegaVis, omegaLog);
                break;
            }

            case blenderType::EXPONENTIAL:
            {
                // (PH:Eq. 31)
                const scalar Gamma = 0.01*pow4(yPlus)/(1 + 5*yPlus);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                omega0[celli] +=
                    w*(omegaVis*exp(-Gamma) + omegaLog*exp(-invGamma));
                break;
            }
            case blenderType::TANH:
            {
                // (KAS:Eqs. 33-34)
                const scalar phiTanh = tanh(pow4(0.1*yPlus));
                const scalar b1 = omegaVis + omegaLog;
                const scalar b2 =
                    pow(pow(omegaVis, 1.2) + pow(omegaLog, 1.2), 1.0/1.2);

                omega0[celli] += phiTanh*b1 + (1 - phiTanh)*b2;
                break;
            }
        }

        if (!(blender_ == blenderType::STEPWISE) || yPlus > yPlusLam)
        {
            G0[celli] +=
                w
                *(nutw[facei] + nuw[facei])
                *magGradUw[facei]
                *Cmu25*sqrt(k[celli])
                /(kappa*y[facei]);
        }
        
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    wallFunctionBlenders(),
    wallCoeffs_(),
    beta1_(0.075),
    roughnessHeight_(p.size(), Zero),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    wallCoeffs_(ptf.wallCoeffs_),
    beta1_(ptf.beta1_),
    roughnessHeight_(ptf.roughnessHeight_, mapper),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::BINOMIAL, scalar(2)),
    wallCoeffs_(dict),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    roughnessHeight_("roughnessHeight", dict, p.size()),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    wallFunctionBlenders(owfpsf),
    wallCoeffs_(owfpsf.wallCoeffs_),
    beta1_(owfpsf.beta1_),
    roughnessHeight_(owfpsf.roughnessHeight_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    wallFunctionBlenders(owfpsf),
    wallCoeffs_(owfpsf.wallCoeffs_),
    beta1_(owfpsf.beta1_),
    roughnessHeight_(owfpsf.roughnessHeight_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& omegaRoughWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& omegaRoughWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaRoughWallFunctionFvPatchScalarField::updateCoeffs()
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

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
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

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintomega.append(omega[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        //scalarField(constraintomega.xfer())
	constraintomega
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
