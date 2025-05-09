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
):
            
    kOmega<BasicTurbulenceModel>
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
    
    alphaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD",
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
    this->alphaK_ = 2.0/3.0;

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
bool kOmegaTNT<BasicTurbulenceModel>::read()
{
    if (kOmega<BasicTurbulenceModel>::read())
    {
        alphaD_.readIfPresent(this->coeffDict());
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
    volScalarField& nut = this->nut_;
    volScalarField& k_ = this->k_;
    volScalarField& omega_ = this->omega_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    
    volScalarField::Internal S2 = 2.0 * magSqr(dev(symm(tgradU())));
    volScalarField::Internal O2 = 2.0 * magSqr(skew(tgradU()));
    
    const volScalarField::Internal GbyNu
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        tgradU()  && devTwoSymm(tgradU())
    );
    tgradU.clear();

    volScalarField::Internal G(this->GName(), nut()*GbyNu);

    // Limiter based on epsilon
    if (productionLimiter_)
        G = min(G,  20 * this->Cmu_ * k_() * omega_());
    
    // Limiter according to Wallin (PhD. thesis, paper 6)
    if (shockLimiter_)
        G = min(G,  k_() * sqrt(S2/2.0));
    
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    volScalarField::Internal CDkOmega = max(
        this->alphaD_ / max(omega_, this->omegaMin()) * (fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar("0",inv(sqr(dimTime)), 0.0)
    );


    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha * rho * this->DomegaEff(), omega_)
     ==
        this->gamma_ * alpha * rho() * GbyNu
      - fvm::SuSp(((2.0/3.0)*this->gamma_)*alpha() * rho() * divU, omega_)
      - fvm::Sp(this->beta_ * alpha() * rho() * omega_(), omega_)
      + alpha * rho * CDkOmega
      + fvOptions(alpha, rho, omega_)
    );


    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin());

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha * rho * this->DkEff(), k_)
     ==
        alpha() * rho() * G
      - fvm::SuSp((2.0/3.0) * alpha() * rho() * divU, k_)
      - fvm::Sp(this->Cmu_ * alpha() * rho() * omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    this->correctNut();
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
