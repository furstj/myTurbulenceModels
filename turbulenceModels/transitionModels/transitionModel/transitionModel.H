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

Class
    Foam::transitionModel

Description

SourceFiles
    transitionModel.C

\*---------------------------------------------------------------------------*/

#ifndef transitionModel_H
#define transitionModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class transitionModel Declaration
\*---------------------------------------------------------------------------*/

class transitionModel
{
protected:
    // Protected data
    dictionary& dict_;
    
    const fvMesh& mesh_;

    const Time& time_;

    const volVectorField& U_;
    
    const volScalarField& k_;
    
    const volScalarField& omega_;

    // Private Member Functions

private:
    //- Disallow default bitwise copy construct
    transitionModel(const transitionModel&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const transitionModel&) = delete;


public:
    TypeName("transitionModel");

    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        transitionModel,
        params,
        (
            dictionary& dict,
            const volVectorField& U,
            const volScalarField& k,
            const volScalarField& omega
        ),
        (dict, U, k, omega)
    );

    
    // Constructors

    //- Construct from components
    transitionModel(
            dictionary& dict,
            const volVectorField& U,
            const volScalarField& k,
            const volScalarField& omega
    );
       

    // Selectors

    
    //- Return a reference to the selected transition model
    static autoPtr<transitionModel> New(
        dictionary& dict,
        const volVectorField& U,
        const volScalarField& k,
        const volScalarField& omega
    );
    

    //- Destructor
    virtual ~transitionModel() = default;
    

    // Member Functions
   
    virtual void correct(const volScalarField& nu, const volScalarField& S, const volScalarField& W) = 0;

    virtual tmp<volScalarField> nut() const = 0;

    virtual tmp<volScalarField> Pk(
        const volScalarField& S, 
        const volScalarField& W
    ) const = 0;
   
   virtual tmp<volScalarField> gammaInt() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
