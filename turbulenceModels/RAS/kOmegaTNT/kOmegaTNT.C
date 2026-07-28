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

\*---------------------------------------------------------------------------*/

#include "kOmegaTNT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaTNT<BasicTurbulenceModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& /* S2 not used */
) const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        gradU() && devTwoSymm(gradU())
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaTNT<BasicTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return GbyNu0;
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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

    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.66
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),
    gamma_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.55
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
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
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    sigmaD_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaD",
            this->coeffDict_,
            0.5
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
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
bool kOmegaTNT<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());
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
void kOmegaTNT<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaTNT<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return G;
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
    volScalarField& nut = this->nut_;
    volScalarField& k = k_;
    volScalarField& omega = omega_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    
    const volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    tgradU.clear();


    // Limiter based on epsilon
    if (productionLimiter_)
        G = min(G,  20 * betaStar_ * k() * omega());

    // Limiter according to Wallin (PhD. thesis, paper 6)
    if (shockLimiter_)
    {
        volScalarField::Internal S2 = 2.0 * magSqr(dev(symm(tgradU())));
        G = min(G,  k() * sqrt(S2/2.0));
    }
    // Update omega and G at the wall
    omega.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    volScalarField::Internal CDkOmega = max(
        sigmaD_ / max(omega(), this->omegaMin_)() * (fvc::grad(k) & fvc::grad(omega))(),
        dimensionedScalar("0", inv(sqr(dimTime)), 0.0)
    );

    
    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega)
      + fvm::div(alphaRhoPhi, omega)
      - fvm::laplacian(alpha * rho * this->DomegaEff(), omega)
     ==
        //gamma_ * alpha * rho() * GbyNu0
        gamma_ * alpha * rho() * G * omega_/k_
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha() * rho() * divU, omega)
      - fvm::Sp(beta_ * alpha() * rho() * omega(), omega)
      + alpha * rho * CDkOmega
      + fvOptions(alpha, rho, omega)
    );


    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega);
    bound(omega, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k)
      + fvm::div(alphaRhoPhi, k)
      - fvm::laplacian(alpha * rho * this->DkEff(), k)
     ==
        alpha() * rho() * Pk(G)
      - fvm::SuSp((2.0/3.0) * alpha() * rho() * divU, k)
      - fvm::Sp(betaStar_ * alpha() * rho() * omega(), k)
      + fvOptions(alpha, rho, k)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k);
    bound(k, this->kMin_);

    this->correctNut();
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
