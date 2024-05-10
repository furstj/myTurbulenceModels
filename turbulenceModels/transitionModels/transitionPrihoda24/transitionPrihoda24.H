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
    Foam::transitionPrihoda24

Description

    Modification of the DIck & Kubacki model proposed by J. Prihoda

   Kubacki, S., & Dick, E. (2016). An algebraic intermittency model for bypass, 
     separation-induced and wake-induced transition. 
     International Journal of Heat and Fluid Flow, 62, 344–361. 
     https://doi.org/10.1016/j.ijheatfluidflow.2016.09.013
     
SourceFiles
    transitionModel.C

\*---------------------------------------------------------------------------*/

#ifndef transitionPrihoda24_H
#define transitionPrihoda24_H

#include "transitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                         Class transitionModel Declaration
\*---------------------------------------------------------------------------*/

class transitionPrihoda24
:
    public transitionModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    transitionPrihoda24(const transitionPrihoda24&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const transitionPrihoda24&) = delete;

protected:

    dimensionedScalar CSS_;

    dimensionedScalar A_;

    dimensionedScalar B_;

    dimensionedScalar Csep_;

    dimensionedScalar Av_;

    dimensionedScalar Clim_;

    dimensionedScalar a1_;

    dimensionedScalar a2_;

    volScalarField fSS_;

    volScalarField nus_;

    volScalarField nul_;

    volScalarField gammaInt_;

    volScalarField Psep_;

    const volScalarField& y_;

public:

    TypeName("Prihoda24");
        
    // Constructors

    //- Construct from components
    transitionPrihoda24(
            dictionary& dict,
            const volVectorField& U,
            const volScalarField& k,
            const volScalarField& omega
    );
    
    //- Destructor
    virtual ~transitionPrihoda24() {};
    

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
     
    virtual tmp<fvScalarMatrix> kSource() const;
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
