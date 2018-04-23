/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::incompressible::RASModels::gammaSST

Description
    Implementation of the gamma-SST turbulence model for predicting
    laminar-turbulent transition on arbitrary geometries and grid topo-
    logies.


    The default model coefficients correspond to the following:
    @verbatim
        gammaSSTCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

    
                Flength       100;
	        ca2             0.06;
		ce2             50;
		sigmaGamma      1;
                CPG1            14.68;
                CPG1lim          1.5;
                CPG2            -7.34;
                CPG3             0;
                CPG2lim         3.0;
                CTU1          100.0;
                CTU2         1000.0;
                CTU3            1.0;
        }
    @endverbatim

SourceFiles
    gammaSST.C

\*---------------------------------------------------------------------------*/

#ifndef gammaSST_H
#define gammaSST_H

#include "fvOptions.H"
#include "kOmegaSST.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class gammaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class gammaSST
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions
    
    // Disallow default bitwise copy construct and assignment
    gammaSST(const gammaSST&);
    gammaSST& operator=(const gammaSST&);
    
protected:
    // Protected data
        
    // Model coefficients
    
    dimensionedScalar Flength_;
    
    dimensionedScalar ca2_;
    dimensionedScalar ce2_;
    dimensionedScalar sigmaGamma_;
    dimensionedScalar CPG1_;
    dimensionedScalar CPG1lim_;
    dimensionedScalar CPG2_;
    dimensionedScalar CPG3_;
    dimensionedScalar CPG2lim_;
    dimensionedScalar CTU1_;
    dimensionedScalar CTU2_;
    dimensionedScalar CTU3_;
    dimensionedScalar ReThetacLim_;
    dimensionedScalar Ck_;
    dimensionedScalar CSEP_;

    // Fields
    volScalarField gammaInt_;


    // Protected member functions

    // Model functions of gammaSST model
    tmp<volScalarField> F1(const volScalarField& CDkOmega) const;	
    tmp<volScalarField> ReThetac() const;
    tmp<volScalarField> Fonset(const volScalarField& S) const;
    tmp<volScalarField> Fonset1(const volScalarField& S) const;
    tmp<volScalarField> Fturb() const;
    tmp<volScalarField> TuL() const;
    tmp<volScalarField> FPG() const;

    
public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("gammaSST");


    // Constructors

        //- Construct from components
        gammaSST
        (            
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName  = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~gammaSST()
        {}
    
    
    // Member Functions

    
    //- Return the effective diffusivity for gamma
    tmp<volScalarField> DgammaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DgammaEff", (this->nut_/sigmaGamma_) + this->nu())
                );
        }
    
    //- Return the intermittency
    virtual tmp<volScalarField> gammaInt() const
        {
            return gammaInt_;
        }
    
    virtual tmp<volScalarField> Rt() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("Rt", this->k_/(this->nu()*this->omega_)) 
                    // LANGTRY's definition
                );
        }
    
    
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
#   include "gammaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
