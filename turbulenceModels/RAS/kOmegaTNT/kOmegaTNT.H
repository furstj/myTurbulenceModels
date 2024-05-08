/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::compressible::RASModels::mykOmegaTNT

Description
    Implementation of the k-omega-TNT turbulence model.

    References:
    \verbatim
        J.C. Kok
        Resolving the dependence on free-stream values for the k-omega turbulence model
        NLR-TP-99295 report, 1999
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            beta        0.075;
            gamma       0.55;  // alpha in the article
            alphaK      0.66;  // alphaK in the article
            alphaOmega  0.5;   // alphaOmega in the article
            alphaD      0.5;   // sigmaD in the article
            productionLimiter  true;
            shockLimiter       true;
        }
    \endverbatim


\*---------------------------------------------------------------------------*/

#ifndef kOmegaTNT_H
#define kOmegaTNT_H

#include "kOmega.H"
#include "eddyViscosity.H"
#include "transitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaTNT Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>  
class kOmegaTNT
:
    public kOmega<BasicTurbulenceModel>
{
    // Private Member Functions
    
    // Disallow default bitwise copy construct and assignment
    
    kOmegaTNT(const kOmegaTNT&);
    kOmegaTNT& operator=(const kOmegaTNT&);

protected:
    // Protected data
    
    // Model coefficients
    dimensionedScalar alphaD_;

    Switch productionLimiter_;
    
    Switch shockLimiter_;

    autoPtr<transitionModel> transitionModel_;

    virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaTNT");


    // Constructors
    
    //- Construct from components
    kOmegaTNT
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
    virtual ~kOmegaTNT()
    {}

    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
    
    //- Calculate the turbulence viscosity for first step
    
    //- Read RASProperties dictionary
    virtual bool read();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaTNT.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
