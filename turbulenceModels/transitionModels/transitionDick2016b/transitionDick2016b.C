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

#include "transitionDick2016b.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

namespace Foam
{

defineTypeNameAndDebug(transitionDick2016b, 0);
addToRunTimeSelectionTable(transitionModel, transitionDick2016b, params);

transitionDick2016b::transitionDick2016b(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
):

    transitionModel(dict, U, k, omega),

    Agamma_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Agamma",
            dict,
            12.0
        )
    ),

    CS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CS",
            dict,
            21.0
        )
    ),

    CA_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CA",
            dict,
            1.0
        )
    ),

    CKH_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CKH",
            dict,
            10.0
        )
    ),

    Ck_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ck",
            dict,
            6.0
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


void transitionDick2016b::correct(
    const volScalarField& nu, 
    const volScalarField& S, 
    const volScalarField& W
    )
{
    const dimensionedScalar kMin("kMin", sqr(dimVelocity), 1.e-10);
    const dimensionedScalar omegaMin("omegaMin", inv(dimTime), 1.e-10);

    const volScalarField Rey = sqrt(max(this->k_, kMin)) * y_ / nu;
    
    const scalar betaStar = 0.09;

    volScalarField X = tanh( - W * (S - W) / (CKH_ * sqr(betaStar*omega_)));
    volScalarField FA = (1.0 - tanh(k_/(Ck_*nu*max(omega_, omegaMin))));
    
    volScalarField CSS_ = CS_ * (1.0 + CA_*FA);
    
    fSS_ = exp(-sqr(CSS_ / Rey));

    nus_ = fSS_ * this->k_ / max(this->omega_, Clim_ * S / a1_);
    nul_ = (1 - fSS_) * this->k_ / max(this->omega_, Clim_ * S / a2_);

    gammaInt_ = min(max(Rey / Agamma_ - 1.0, 0.0), 1.0); 

    tmp<volScalarField> ReV = sqr(y_) * S / nu;
    tmp<volScalarField> Fsep = min(max(ReV / (2.2 * Av_) - 1.0, 0.0), 1.0);
    Psep_ = Csep_ * Fsep * nu * sqr(S);
}


tmp<volScalarField> transitionDick2016b::nut() const
{
    return nus_ + nul_;
}


tmp<volScalarField> transitionDick2016b::Pk(
    const volScalarField& S, 
    const volScalarField& W
    ) const
{
    return nus_ * sqr(S);
}

tmp<fvScalarMatrix> transitionDick2016b::kSource() const
{
    return transitionModel::kSource() + (1.0 - gammaInt_) * Psep_;
}

} // End namespace Foam

// ************************************************************************* //
