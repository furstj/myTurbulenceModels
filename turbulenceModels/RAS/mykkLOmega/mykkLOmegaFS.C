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

#include "mykkLOmegaFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::L(const volScalarField& ReOmega) const
{
    const volScalarField& p = this->mesh_.objectRegistry::lookupObject<volScalarField>("p");
    const volVectorField& U_ = this->U_;

    dimensionedScalar pTot("pTot", p.dimensions(), 
    gMax( volScalarField( p + 0.5*magSqr(U_) ) ) );
    
    dimensionedScalar uMin("uMin", dimVelocity, VSMALL);
    
    volScalarField Ue = sqrt( 2.0 * (pTot - p) ); 
    
    volScalarField dpdx = (fvc::grad(p) & U_) / max( mag(U_), uMin); 
  
    volScalarField K = - this->nu() / pow3(max(Ue,uMin)) * dpdx;

    return ( sqr(ReOmega) * K );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::RevByReTheta(const volScalarField& L) const
{
    volScalarField Lapg = max( min(L, 0.0), -1.5);
    return  2.1884 * ( 1.0 - 0.95419*Lapg - 0.13183*sqr(Lapg) );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::lambdaTheta(const volScalarField& L) const
{
    return tmp<volScalarField>(new volScalarField(
        "lambdaTheta",
        L / sqr(RevByReTheta(L))
    ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::ReThetac(const volScalarField& lambda) const
{
    return tmp<volScalarField>(new volScalarField(
        "ReThetac",
        200.69 / (1.0 - 55.287 * lambda + 3.4992e+05 * pow4(lambda) ) 
    ));   
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::CtsCrit(const volScalarField& L) const
{
    return tmp<volScalarField>(new volScalarField(
        "CtsCrit",
        ReThetac(lambdaTheta(L)) * RevByReTheta(L)
    ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::CnatCrit(const volScalarField& L) const
{
    return tmp<volScalarField>(new volScalarField(
            "CnatCrit",
            this->CnatCrit_ * CtsCrit(L) / 439.19
        ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    volScalarField L_ = L(ReOmega);

    return tmp<volScalarField>(new volScalarField(
            "BetaTS",
            scalar(1) - exp(-sqr(max(ReOmega - CtsCrit(L_), scalar(0)))/this->Ats_)
        ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaFS<BasicTurbulenceModel>::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    volScalarField L_ = L(ReOmega);

    return tmp<volScalarField>(new volScalarField(
        "phiNAT",
        max
        (
            ReOmega
            - CnatCrit(L_)
            / (
                fNatCrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    ));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykkLOmegaFS<BasicTurbulenceModel>::mykkLOmegaFS
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
        // Evaluating nut_ is complex so start from the field read from file
        this->nut_.correctBoundaryConditions();
        
        this->printCoeffs(type);

	if (debug) 
            Info << "Debug switch is on!" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mykkLOmegaFS<BasicTurbulenceModel>::read()
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
void mykkLOmegaFS<BasicTurbulenceModel>::correct()
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
