/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "kOmegaSSTCCM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTCCM<BasicTurbulenceModel>::kOmegaSSTCCM
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
    kOmegaSST
    <
        BasicTurbulenceModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTCCM<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S("S", sqrt(2.0)*mag(symm(tgradU())));
    const volScalarField W("W", sqrt(2.0)*mag(skew(tgradU())));
    dimensionedScalar eps("eps", dimless/dimTime, 1.e-10);
    
    const volScalarField rs = S/max(W,eps);
    const volScalarField rt = W/max(S,eps)*(W/max(S,eps) - 1.0);

    const scalar cr1 = 1.0;
    const scalar cr2 = 2.0;
    const scalar cr3 = 1.0;
    tmp<volScalarField> frot = (1 + cr1)*(2*rs)/(rs + 1)*(1 - cr3*atan(cr2*rt)) - cr1;
    const scalar cScale = 1.0;
    tmp<volScalarField> frt = max(0.0, min(frot, 1.25));
    volScalarField fr  = max(0.0, 1.0 + cScale*(frt - 1.0));
    
    return kOmegaSST<BasicTurbulenceModel>::Pk(G)*fr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
