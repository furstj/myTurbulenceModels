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
    Foam::transitionDick2016a

Description

   Kubacki, Dick: An algebraic model for bypass transition in turbomachinery boundary layer flows

SourceFiles
    transitionModel.C

\*---------------------------------------------------------------------------*/

#ifndef transitionDick2016a_H
#define transitionDick2016a_H

#include "transitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                         Class transitionModel Declaration
\*---------------------------------------------------------------------------*/

class transitionDick2016a
:
    public transitionModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    transitionDick2016a(const transitionDick2016a&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const transitionDick2016a&) = delete;

protected:

    dimensionedScalar AT_;

    dimensionedScalar CT_;

    dimensionedScalar CSS_;

    dimensionedScalar Clim_;

    dimensionedScalar a1_;

    dimensionedScalar a2_;

    volScalarField fSS_;

    volScalarField nus_;

    volScalarField nul_;

    volScalarField gammaInt_;

    const volScalarField& y_;

public:

    TypeName("Dick2016a");
        
    // Constructors

    //- Construct from components
    transitionDick2016a(
            dictionary& dict,
            const volVectorField& U,
            const volScalarField& k,
            const volScalarField& omega
    );
    
    //- Destructor
    virtual ~transitionDick2016a() {};
    

    // Member Functions

    //- Main update/correction routine
    virtual void correct(const volScalarField& nu, const volScalarField& S, const volScalarField& W);
    
    virtual bool hasShearStressLimiter() const
    {
        return true;
    }   

    virtual tmp<volScalarField> nut() const;

    virtual tmp<volScalarField> Pk(
        const volScalarField& S, 
        const volScalarField& W
    ) const;

    virtual tmp<volScalarField> gammaInt() const
    {
        return gammaInt_;
    }
     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
