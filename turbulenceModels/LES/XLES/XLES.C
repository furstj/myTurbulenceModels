
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "XLES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> XLES<BasicTurbulenceModel>::Lt() const
{
    return sqrt(this->k_())/(this->betaStar_*this->omega_());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> XLES<BasicTurbulenceModel>::FDES() const
{
  return max(Lt()/(CDES_*this->delta()()), scalar(1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> XLES<BasicTurbulenceModel>::epsilonByk() const
{
    return FDES() * this->betaStar_ * this->omega_;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> XLES<BasicTurbulenceModel>::omegaSource() const
{
    return fvm::Su(
        max(
            this->alphaD_/this->omega_ * (fvc::grad(this->k_) & fvc::grad(this->omega_)),
            dimensionedScalar("0",inv(sqr(dimTime)), 0.0)
        ),
        this->omega_
    );
}

template<class BasicTurbulenceModel>
void XLES<BasicTurbulenceModel>::correctNut()
{
  this->nut_ = min
    (
     this->k_/this->omega_,
     this->betaStar_*CDES_*this->delta()*sqrt(this->k_)
     );

  this->nut_.correctBoundaryConditions();
  fv::options::New(this->mesh_).correct(this->nut_);
  BasicTurbulenceModel::correctNut();

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
XLES<BasicTurbulenceModel>::XLES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmega<
        LESeddyViscosity<BasicTurbulenceModel>,
        BasicTurbulenceModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    
    alphaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD",
            this->coeffDict_,
            0.5
        )
    ),

    productionLimiter_
    (
        Switch::lookupOrAddToDict
        (
            "productionLimiter",
            this->coeffDict_,
            true
        )
    ),

    shockLimiter_
    (
        Switch::lookupOrAddToDict
        (
            "shockLimiter",
            this->coeffDict_,
            true
        )
    ),

    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.61
        )
    )
{
    this->coeffDict_.set("alphaK", 2.0/3.0);
    this->coeffDict_.set("beta", 0.75);
    this->coeffDict_.set("gamma", 0.55);

    this->read();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool XLES<BasicTurbulenceModel>::read()
{
    if
    (
        kOmega<LESeddyViscosity<BasicTurbulenceModel>, BasicTurbulenceModel>
        ::read()
    )
    {
        alphaD_.readIfPresent(this->coeffDict());
        productionLimiter_.readIfPresent("productionLimiter", this->coeffDict());
        shockLimiter_.readIfPresent("shockLimiter", this->coeffDict());
        CDES_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
