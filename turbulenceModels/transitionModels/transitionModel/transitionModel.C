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

#include "transitionModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

namespace Foam
{

defineTypeNameAndDebug(transitionModel, 0);
defineRunTimeSelectionTable(transitionModel, params);


autoPtr<transitionModel>
transitionModel::New
(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
)
{
    word transitionModelName
    (
        dict.lookupOrAddDefault(word("transitionModel"), word("none"))
    );

    Info<< "Selecting transition model "
        << transitionModelName << endl;

    auto cstrIter = paramsConstructorTablePtr_->cfind(transitionModelName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown transitionModel "
            << transitionModelName << endl << endl
            << "Valid transitionModels are : " << endl
            << paramsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<transitionModel>(cstrIter()(dict, U, k, omega));
}


transitionModel::transitionModel(
    dictionary& dict,
    const volVectorField& U,
    const volScalarField& k,
    const volScalarField& omega
):
    dict_(dict),
    mesh_(U.mesh()),
    time_(mesh_.time()),
    U_(U),
    k_(k),
    omega_(omega)
{}



}
// ************************************************************************* //
