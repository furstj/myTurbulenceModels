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

\*---------------------------------------------------------------------------*/

#include "mykkLOmegaDev.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaDev<BasicTurbulenceModel>::lambdaT() const
{
  
  volScalarField fINT2("fINT2", min(scalar(1), 0.9*sqr(this->kl_/this->kt_)) );

  volScalarField lambdaF("lambdaF", sqrt(this->kt_ + this->kl_ * fINT2) / this->omega_);
  
  return tmp<volScalarField>(new volScalarField(
	"lambdaT",
	sqrt(this->kt_ + this->kl_ * fINT2) / this->omega_
  ));
}
  
template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaDev<BasicTurbulenceModel>::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return tmp<volScalarField>(new volScalarField(
        "fTaul",
        scalar(1)
        - exp
        (
            -this->CtauL_ * ktL
            /
            (
                sqr
                (
                    lambdaEff * this->omega_
                    + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    ));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykkLOmegaDev<BasicTurbulenceModel>::mykkLOmegaDev
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    mykkLOmega<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    )
{
    
    {
      Info << "================================================" << endl;
      Info << "| WARNING: This is unstable development model! |" << endl;
      Info << "================================================" << endl;
      
      this->CnatCrit_ = 1800;
      this->CrNat_ = 0.04;

        // Evaluating nut_ is complex so start from the field read from file
        this->nut_.correctBoundaryConditions();
        
        this->printCoeffs(type);

	if (debug) 
            Info << "Debug switch is on!" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mykkLOmegaDev<BasicTurbulenceModel>::read()
{
    if (  mykkLOmega<BasicTurbulenceModel>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void mykkLOmegaDev<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    mykkLOmega<BasicTurbulenceModel>::correct();

    if (debug && this->runTime_.outputTime()) {
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
