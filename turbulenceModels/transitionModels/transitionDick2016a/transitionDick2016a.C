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

#include "transitionDick2016a.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

namespace Foam
{

defineTypeNameAndDebug(transitionDick2016a, 0);
addToRunTimeSelectionTable(transitionModel, transitionDick2016a, params);

transitionDick2016a::transitionDick2016a(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
):

    transitionModel(dict, U, k, omega),

    AT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "AT",
            dict,
            10.0
        )
    ),

    CT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CT",
            dict,
            15.5
        )
    ),

    CSS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CSS",
            dict,
            2.0
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

    y_(wallDist::New(U.mesh()).y())

{
}


void transitionDick2016a::correct(
    const volScalarField& nu, 
    const volScalarField& S, 
    const volScalarField& W
    )
{
    const dimensionedScalar kMin("kMin", sqr(dimVelocity), 1.e-10);

    fSS_ = exp(-sqr(CSS_ * nu * W / max(this->k_, kMin)));

    nus_ = fSS_ * this->k_ / max(this->omega_, Clim_ * S / a1_);
    nul_ = (1 - fSS_) * this->k_ / max(this->omega_, Clim_ * S / a2_);

    volScalarField Rey = sqrt(this->k_) * y_ / nu;
    volScalarField zetaT = max(Rey - CT_, 0.0);
 
    gammaInt_ = min(zetaT / AT_, 1.0); 
}


tmp<volScalarField> transitionDick2016a::nut() const
{
    return nus_ + nul_;
}


tmp<volScalarField> transitionDick2016a::Pk(
    const volScalarField& S, 
    const volScalarField& W
    ) const
{
    return nus_ * sqr(S);
}


} // End namespace Foam

// ************************************************************************* //
