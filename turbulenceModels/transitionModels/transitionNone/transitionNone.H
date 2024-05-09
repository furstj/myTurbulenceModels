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
    Foam::transitionNone

Description

    Model for fully turbulent case

SourceFiles
    transitionModel.C

\*---------------------------------------------------------------------------*/

#ifndef transitionNone_H
#define transitionNone_H

#include "transitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                         Class transitionModel Declaration
\*---------------------------------------------------------------------------*/

class transitionNone
:
    public transitionModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    transitionNone(const transitionNone&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const transitionNone&) = delete;

public:

    TypeName("none");
        
    // Constructors

    //- Construct from components
    transitionNone(
            dictionary& dict,
            const volVectorField& U,
            const volScalarField& k,
            const volScalarField& omega
    );
    
    //- Destructor
    virtual ~transitionNone() {};
    

    // Member Functions

    //- Main update/correction routine
    virtual void correct(const volScalarField& nu, const volScalarField& S, const volScalarField& W);
    
    virtual tmp<volScalarField> nut() const;

    virtual tmp<volScalarField> Pk(
        const volScalarField& S, 
        const volScalarField& W
    ) const;

    virtual tmp<volScalarField> gammaInt() const;
     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
