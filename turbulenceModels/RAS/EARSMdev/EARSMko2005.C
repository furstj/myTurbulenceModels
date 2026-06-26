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

#include "EARSMko2005.H"
#include "fvOptions.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EARSMko2005<BasicTurbulenceModel>::EARSMko2005
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
    EARSMBase<kOmegaSST<BasicTurbulenceModel>>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    kInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            sqr(dimVelocity),
            1.e-10
        )
    )
{
    if (type == typeName)
    {
        auto originalCoeffs = this->getOriginalCoeffs();

        if (!originalCoeffs.found("alphaK1")) this->coeffDict_.set("alphaK1", 1.1);
        if (!originalCoeffs.found("alphaK2")) this->coeffDict_.set("alphaK2", 1.1);
        if (!originalCoeffs.found("alphaOmega1")) this->coeffDict_.set("alphaOmega1", 0.53);
        if (!originalCoeffs.found("alphaOmega2")) this->coeffDict_.set("alphaOmega2", 1.0);
        if (!originalCoeffs.found("beta1")) this->coeffDict_.set("beta1", 0.0747);
        if (!originalCoeffs.found("beta2")) this->coeffDict_.set("beta2", 0.0828);
        if (!originalCoeffs.found("gamma1")) this->coeffDict_.set("gamma1", 0.518);
        if (!originalCoeffs.found("gamma2")) this->coeffDict_.set("gamma2", 0.44);

        this->read();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> EARSMko2005<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    
    // Compute gradKgradOmegaByOmega from CDkOmega
    // SST defines CDkOmega = 2*alphaOmega2_ * (gradK & gradOmega) / omega
    // We need to remove the 2*alphaOmega2_ factor
    volScalarField gradKgradOmegaByOmega = CDkOmega / (2.0 * this->alphaOmega2_);

    // EARSM fMix blending function (Hellsten 2005)
    tmp<volScalarField> Gamma = min
    (
        min
        (
            max
            (
                sqrt(this->k_) / (this->betaStar_ * this->omega_ * this->y_),
                scalar(500) * this->nu() / (this->omega_ * sqr(this->y_))
            ),
            scalar(20) * this->k_ / max(sqr(this->y_) * gradKgradOmegaByOmega, scalar(200) * kInf_)
        ),
        scalar(10)
    );

    return tanh(scalar(1.5) * pow4(Gamma));
}


template<class BasicTurbulenceModel>
bool EARSMko2005<BasicTurbulenceModel>::read()
{
    return EARSMBase<kOmegaSST<BasicTurbulenceModel>>::read(typeName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
