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

Class
    Foam::RASModels::mykkLOmegaDev

Group
    grpRASTurbulence

Description
    Low Reynolds-number k-kl-omega turbulence model for
    incompressible and compressible flows.

    This turbulence model is described in:
    \verbatim
        Walters, D. K., & Cokljat, D. (2008).
        A three-equation eddy-viscosity model for Reynolds-averaged
        Navier–Stokes simulations of transitional flow.
        Journal of Fluids Engineering, 130(12), 121401.
    \endverbatim

    however the paper contains several errors which must be corrected for the
    model to operation correctly as explained in

    \verbatim
        Furst, J. (2013).
        Numerical simulation of transitional flows with laminar kinetic energy.
        Engineering MECHANICS, 20(5), 379-388.
    \endverbatim

    All these corrections and updates are included in this implementation.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! WARNING: this is just a sandbox model for testing different modifications !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SourceFiles
    mykkLOmegaDev.C

\*---------------------------------------------------------------------------*/

#ifndef mykkLOmegaDev_H
#define mykkLOmegaDev_H

#include "mykkLOmega.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class mykkLOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class mykkLOmegaDev
:
    public mykkLOmega<BasicTurbulenceModel>
{
    // Private memmber functions

    // Disallow default bitwise copy construct and assignment
    mykkLOmegaDev(const mykkLOmegaDev&);
    mykkLOmegaDev& operator=(const mykkLOmegaDev&);

protected:
    
    // Protected data
    
    // Model coefficients
    
    // Fields
    
    virtual tmp<volScalarField> lambdaT() const;

    virtual tmp<volScalarField> fTaul
    (
        const volScalarField& lambdaEff,
        const volScalarField& ktL,
        const volScalarField& omega
    ) const;

    
public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mykkLOmegaDev");
    
    
    // Constructors
    
    //- Construct from components
    mykkLOmegaDev
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );
    
    
    //- Destructor
    virtual ~mykkLOmegaDev()
        {}


    // Member Functions

    //- Re-read model coefficients if they have changed
    virtual bool read();

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "mykkLOmegaDev.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
