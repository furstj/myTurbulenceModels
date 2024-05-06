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
    Foam::compressible::RASModels::mykOmegaTNTTrans

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
        kOmegaTNTTransCoeffs
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

#ifndef kOmegaTNTTrans_H
#define kOmegaTNTTrans_H

#include "kOmega.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaTNTTrans Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>  
class kOmegaTNTTrans
:
    public kOmega<BasicTurbulenceModel>
{
    // Private Member Functions
    
    // Disallow default bitwise copy construct and assignment
    
    kOmegaTNTTrans(const kOmegaTNTTrans&);
    kOmegaTNTTrans& operator=(const kOmegaTNTTrans&);

    
protected:
    // Protected data
    
    // Model coefficients
    dimensionedScalar alphaD_;

    Switch productionLimiter_;
    
    Switch shockLimiter_;
    
    //- Transition model parameters
    dimensionedScalar CSS_;
    dimensionedScalar A_;
    dimensionedScalar B_;
    dimensionedScalar a1_;
    dimensionedScalar a2_;
    dimensionedScalar Clim_;

    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaTNTTrans");


    // Constructors
    
    //- Construct from components
    kOmegaTNTTrans
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
    virtual ~kOmegaTNTTrans()
    {}

    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
    
    //- Read RASProperties dictionary
    virtual bool read();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaTNTTrans.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
