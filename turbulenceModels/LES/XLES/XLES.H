/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::LESModels::XLES

Description
    Implementation of the XLES turbulence model for
    incompressible and compressible flows.

    DES model described in:
    \verbatim
        KOK, Johan C., H. S. DOL, B. OSKAM a Harmen VAN DER VEN. 
	EXTRA-LARGE EDDY SIMULATION OF MASSIVELY SEPARATED FLOWS. 
	AIAA Paper. 2004, (2004–264), 1–12. doi:10.2514/6.2004-264
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

Group
    grpLESTurbulence

See also

SourceFiles
    XLES.C

\*---------------------------------------------------------------------------*/

#ifndef XLES_H
#define XLES_H

#include "kOmegaBase.H"
#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class XLES
:
    public Foam::kOmega
    <
        LESeddyViscosity<BasicTurbulenceModel>,
        BasicTurbulenceModel
    >
{
    // Private Member Functions
    
    // Disallow default bitwise copy construct and assignment
    
    XLES(const XLES&);
    XLES& operator=(const XLES&);

protected:

    // Protected data
    
    // Model constants
    dimensionedScalar alphaD_;

    Switch productionLimiter_;
    
    Switch shockLimiter_;
 
    //- DES coefficient
    dimensionedScalar CDES_;


    // Protected Member Functions

    //- Return the turbulent length-scale
    tmp<volScalarField::Internal> Lt() const;

    //- Return the DES dissipation rate multiplier
    tmp<volScalarField::Internal> FDES() const;

    tmp<volScalarField::Internal> epsilonByk() const;

    virtual tmp<fvScalarMatrix> omegaSource() const;

    virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("XLES");


    // Constructors

    //- Construct from components
    XLES
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
    virtual ~XLES()
    {}

    // Member Functions
 
    //- Read model coefficients if they have changed
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "XLES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
