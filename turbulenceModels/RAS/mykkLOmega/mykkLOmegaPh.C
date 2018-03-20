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

#include "mykkLOmegaPh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::Ue(const volScalarField& p, const volVectorField& U) const
{

    if ( p.dimensions() == dimensionSet(0, 2, -2, 0, 0) ) 
    { 
        dimensionedScalar pTot("pTot", p.dimensions(), 
        gMax( volScalarField( p + 0.5*magSqr(U) ) ) );
        
        return tmp<volScalarField>(new volScalarField(
            "Ue",
            sqrt( 2.0 * (pTot - p) )
        ));
    } 
    else 
    {
        const basicThermo& thermo =
            this->mesh_.objectRegistry::lookupObject<basicThermo>("thermophysicalProperties");

	const volScalarField gamma("gamma", thermo.Cp() / thermo.Cv());
	const volScalarField a("a", sqrt( gamma * p / thermo.rho() ) );

	dimensionedScalar pTot( "pTot", p.dimensions(), gMax( volScalarField(
            p * pow( 1.0 + (gamma-1.0)/2 * magSqr(U)/sqr(a), gamma/(gamma-1) ) 
        )));
        
        return tmp<volScalarField>(new volScalarField(
            "Ue",
            sqrt(2/(gamma-1) * (pow( p/pTot, (1-gamma)/gamma ) - 1.0)) * a
        ));
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::L(const volScalarField& ReOmega) const
{
    const volScalarField& p = this->mesh_.objectRegistry::lookupObject<volScalarField>("p");
    const volVectorField& U_ = this->U_;

    dimensionedScalar uMin("uMin", dimVelocity, VSMALL);
    
    volScalarField dpdx = (fvc::grad(p) & U_) / max( mag(U_), uMin); 
  
    tmp<volScalarField> K;
    if ( p.dimensions() == dimensionSet(0, 2, -2, 0, 0) ) 
        K = - this->nu() / pow3(max(Ue(p,U_),uMin)) * dpdx;
    else
    {
        const basicThermo& thermo =
            this->mesh_.objectRegistry::lookupObject<basicThermo>("thermophysicalProperties");
        K = - this->nu() / pow3(max(Ue(p,U_),uMin)) * dpdx / thermo.rho();
    }

    return ( sqr(ReOmega) * K );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::CtsCrit(const volScalarField& L) const
{
    return tmp<volScalarField>(new volScalarField(
        "CtsCrit",
        this->CtsCrit0_ / ( scalar(1) - this->CnatApg_ * min(L, scalar(0)) )
    ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::CnatCrit(const volScalarField& L) const
{
    return tmp<volScalarField>(new volScalarField(
            "CnatCrit",
            this->CnatCrit_ / ( scalar(1) - this->CnatApg_ * min(L, scalar(0)) )
        ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    volScalarField L_ = L(ReOmega);

    return tmp<volScalarField>(new volScalarField(
            "BetaTS",
            scalar(1) - exp(-sqr(max(ReOmega - CtsCrit(L_), scalar(0)))/this->Ats_)
        ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmegaPh<BasicTurbulenceModel>::phiNAT
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
mykkLOmegaPh<BasicTurbulenceModel>::mykkLOmegaPh
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
     ),
    CtsCrit0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit0",
            this->coeffDict_,
            536.40
        )
    ),
    CnatApg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatApg",
            this->coeffDict_,
            8.963
        )
    )

{
    
    {
        // Evaluating nut_ is complex so start from the field read from file
        this->nut_.correctBoundaryConditions();
        if (type == typeName)
        {
            this->printCoeffs(type);
            if (debug) 
            {
                Info << "Debug switch is on!" << endl;
            }
        }

    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mykkLOmegaPh<BasicTurbulenceModel>::read()
{
    if (  mykkLOmega<BasicTurbulenceModel>::read())
    {
        CtsCrit0_.readIfPresent(this->coeffDict());
        CnatApg_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void mykkLOmegaPh<BasicTurbulenceModel>::correct()
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
