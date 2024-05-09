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

#include "transitionNone.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

namespace Foam
{

defineTypeNameAndDebug(transitionNone, 0);
addToRunTimeSelectionTable(transitionModel, transitionNone, params);

transitionNone::transitionNone(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
):
    transitionModel(dict, U, k, omega)
{
}


void transitionNone::correct(
    const volScalarField& nu, 
    const volScalarField& S, 
    const volScalarField& W
    )
{
}


tmp<volScalarField> transitionNone::nut() const
{
    const dimensionedScalar omegaMin("omegaMin", inv(dimTime), 1.e-10);
    return k_ / max(omega_, omegaMin);
}


tmp<volScalarField> transitionNone::Pk(
    const volScalarField& S, 
    const volScalarField& W
    ) const
{
    const dimensionedScalar omegaMin("omegaMin", inv(dimTime), 1.e-10);
    return k_/max(omega_, omegaMin) * sqr(S);
}


tmp<volScalarField> transitionNone::gammaInt() const
{
    return tmp<volScalarField>::New(
        IOobject
        (
            IOobject::groupName("gammaInt", U_.group()),
            this->time_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    );
    
} // End namespace Foam

} // End namespace Foam
// ************************************************************************* //
