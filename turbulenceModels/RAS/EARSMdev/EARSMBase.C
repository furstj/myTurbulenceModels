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

#include "EARSMBase.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ScaleModel>
EARSMBase<ScaleModel>::EARSMBase
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
    ScaleModel(
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    CTau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTau_",
            this->coeffDict_,
            6.0
        )
    ),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            -0.72
        )
    ),

    nearWallDamping_
    (
        Switch::lookupOrAddToDict
        (
            "nearWallDamping",
            this->coeffDict_,
            false
        )
    ),

    B2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B2",
            this->coeffDict_,
            1.8
        )
    ),

    CyA_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CyA",
            this->coeffDict_,
            0.092
        )
    ),

    CyB_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CyB",
            this->coeffDict_,
            1.2e-4
        )
    )
{
}
 

template<class ScaleModel>
void EARSMBase<ScaleModel>::correctNut(const volScalarField& S2)
{
    Info << "********** EARSMBase: Correcting nut" << endl;
  ScaleModel::correctNut(S2);
}

template<class ScaleModel>
bool EARSMBase<ScaleModel>::read(const word& type)
{
    return ScaleModel::read();
}

} // End namespace RASModels
} // End namespace Foam
