/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "transitionPrihoda24.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

namespace Foam
{

defineTypeNameAndDebug(transitionPrihoda24, 0);
addToRunTimeSelectionTable(transitionModel, transitionPrihoda24, params);

transitionPrihoda24::transitionPrihoda24(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
):

    transitionModel(dict, U, k, omega),

    CSS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CSS",
            dict,
            2.75
        )
    ),

    A_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "A",
            dict,
            620.0
        )
    ),

    B_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "B",
            dict,
            453.0
        )
    ),

    Csep_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Csep",
            dict,
            2.0
        )
    ),

    Av_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Av",
            dict,
            550.0
        )
    ),


    Clim_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Clim",
            dict,
            0.875
        )
    ),

    a1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "a1",
            dict,
            0.3
        )
    ),

    a2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "a2",
            dict,
            0.45
        )
    ),

    fSS_
    (
        IOobject
        (
            IOobject::groupName("fSS", U.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),

    nus_
    (
        IOobject
        (
            IOobject::groupName("nus", U.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimViscosity, 0.0)
    ),

    nul_
    (
        IOobject
        (
            IOobject::groupName("nul", U.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimViscosity, 0.0)
    ),

    gammaInt_
    (
        IOobject
        (
            IOobject::groupName("gammaInt", U.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),

    Psep_
    (
        IOobject
        (
            IOobject::groupName("Psep", U.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimVelocity)/dimTime, 0.0)
    ),

    y_(wallDist::New(U.mesh()).y())

{
}


void transitionPrihoda24::correct(
    const volScalarField& nu, 
    const volScalarField& S, 
    const volScalarField& W
    )
{
    const dimensionedScalar kMin("kMin", sqr(dimVelocity), 1.e-10);
    const dimensionedScalar omegaMin("omegaMin", inv(dimTime), 1.e-10);

    const volScalarField ReW = max(k_, kMin) / (nu * max(W, omegaMin));
    const volScalarField ReV = sqr(y_) * W / nu;

    const volScalarField Rey = sqrt(max(this->k_, kMin)) * y_ / nu;
    
    fSS_ = exp(-sqr(CSS_ / ReW));
    
    nus_ = fSS_ * this->k_ / max(this->omega_, Clim_ * S / a1_);
    nul_ = (1 - fSS_) * this->k_ / max(this->omega_, Clim_ * S / a2_);

    const volScalarField zetaT = max(ReV - B_, 0.0);

    gammaInt_ = min(zetaT / A_, 1.0); 

    //ReV = sqr(y_) * S / nu;
    tmp<volScalarField> Fsep = min(max(ReV / (2.2 * Av_) - 1.0, 0.0), 1.0);
    Psep_ = Csep_ * Fsep * nu * sqr(S);
}


tmp<volScalarField> transitionPrihoda24::nut() const
{
    return nus_ + nul_;
}


tmp<volScalarField> transitionPrihoda24::Pk(
    const volScalarField& S, 
    const volScalarField& W
    ) const
{
    return nus_ * sqr(S);
}

tmp<fvScalarMatrix> transitionPrihoda24::kSource() const
{
    return transitionModel::kSource() + (1.0 - gammaInt_) * Psep_;
}

} // End namespace Foam

// ************************************************************************* //
