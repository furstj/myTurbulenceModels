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

#include "kOmegaTrans.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaTrans<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    this->nut_ = k_ / this->omegaTilde(S2);
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void kOmegaTrans<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(dev(symm(fvc::grad(this->U_)))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> 
kOmegaTrans<BasicTurbulenceModel>::omegaTilde(const volScalarField& S2) const
{
    return tmp<volScalarField>
        (
            new volScalarField
            (
                "omegaBar",
                max(omega_, Clim_*sqrt(S2/Cmu_))
        )
        );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> 
kOmegaTrans<BasicTurbulenceModel>::beta(const volTensorField& gradU) const
{
    volTensorField Omega(skew(gradU));
    volSymmTensorField Shat(symm(gradU) - 0.5*tr(gradU)*I);
    volScalarField Xomega( mag( (Omega & Omega) && Shat ) / pow(Cmu_*omega_,3) );

    return tmp<volScalarField>
        (
            new volScalarField
            (
                "beta",
                beta0_ * (1 + 85*Xomega) / (1 + 100*Xomega)
            )
        );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> 
kOmegaTrans<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                k_,
                dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
            )
        );
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> 
kOmegaTrans<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                omega_,
                dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
            )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaTrans<BasicTurbulenceModel>::kOmegaTrans
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
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    beta0_ 
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta0",
            this->coeffDict_,
            0.0708
        )
    ),
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.875
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.6
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
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
            0.125
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
bool kOmegaTrans<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        Clim_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaTrans<BasicTurbulenceModel>::correct()
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

    volScalarField G(
        this->GName(),
        this->nut_ * ( dev(twoSymm(tgradU())) && tgradU() ) 
    );

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    tmp<volScalarField> beta_ = this->beta(tgradU);
    tgradU.clear();

    volScalarField CDkOmega = max(
        sigmaD_/omega_*(fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar("0", sqr(inv(dimTime)), 0.0)
    );

    // Turbulent frequency equation 
    // source term modified according to NLR-TP-2001-238
    dimensionedScalar kMin("kMin", sqr(dimVelocity), VSMALL);
    tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
            + fvm::div(alphaRhoPhi, omega_)
            - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
            ==
            gamma_ * alpha * rho * G * omega_ / max(k_, kMin)
            - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rho*divU, omega_)
            - fvm::Sp(beta_*alpha*rho*omega_, omega_)
            + CDkOmega
            + omegaSource()
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
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Cmu_*alpha*rho*omega_, k_)
      + kSource()
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
