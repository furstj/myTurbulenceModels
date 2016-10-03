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

#include "kOmegaTNT.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaTNT<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaTNT<BasicTurbulenceModel>::kOmegaTNT
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
    eddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

   alphaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaStar",
            this->coeffDict_,
            1.0
        )
    ),

    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            5.0/9.0
        )
    ),

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            3.0/40.0
        )
    ),

    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            2.0/3.0
        )
    ),

    sigmaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega",
            this->coeffDict_,
            0.5
        )
    ),

    sigmaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaD",
            this->coeffDict_,
            0.5
        )
    ),

    vortexModification_
    (
     Switch::lookupOrAddToDict
     (
      "vortexModification",
      this->coeffDict_,
      false
      )
    ),

    productionLimiter_
    (
     Switch::lookupOrAddToDict
     (
      "productionLimiter",
      this->coeffDict_,
      true
      )
    ),

    shockLimiter_
    (
     Switch::lookupOrAddToDict
     (
      "shockLimiter",
      this->coeffDict_,
      true
      )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaTNT<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        alphaStar_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());
        vortexModification_.readIfPresent("vortexModification", this->coeffDict());
        productionLimiter_.readIfPresent("productionLimiter", this->coeffDict());
        shockLimiter_.readIfPresent("shockLimiter", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaTNT<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField S2("S2", 2.0*magSqr(dev(symm(tgradU()))) );
    volScalarField O2("O2", 2.0*magSqr(skew(tgradU())) );

    volScalarField GbyK = alphaStar_/omega_*S2 - 2./3.*divU;

    // Limiter based on e 
    if (productionLimiter_)
      GbyK = min(GbyK, 20*betaStar_*omega_);

    // Limiter according to Wallin (PhD. thesis, paper 6)
    if (shockLimiter_)
      GbyK = min(GbyK,  sqrt(S2/2.0));

    volScalarField G
    (
        this->GName(),
        k_ * GbyK
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega = max(
        sigmaD_*rho/omega_*(fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar("0", rho.dimensions()/sqr(dimTime), 0.0)
    );

    // Turbulent frequency equation 
    // source term modified according to NLR-TP-2001-238
    volScalarField GOmega("RASModel::GOmega",
			  vortexModification_ ?
			  alphaStar_ * alphaOmega_ * rho * max(S2, O2) :
			  alphaOmega_*rho*(alphaStar_*S2 - 2./3.*divU*omega_)
			  );

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        GOmega
      - fvm::Sp(beta_*alpha*rho*omega_, omega_)
      + CDkOmega
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
