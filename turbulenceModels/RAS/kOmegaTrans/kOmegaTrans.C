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
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kOmegaTrans<BasicMomentumTransportModel>::correctNut(const volScalarField& nus, const volScalarField& nul)
{
    this->nut_ = nus + nul;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}

template<class BasicMomentumTransportModel>
void kOmegaTrans<BasicMomentumTransportModel>::correctNut()
{
    volTensorField gradU(fvc::grad(this->U_));
    volScalarField S( sqrt(2*magSqr(dev(symm(gradU)))) );
    volScalarField W( sqrt(2*magSqr(skew(gradU))) );
    volScalarField fSS_(this->fSS(S,W));

    correctNut(this->nus(S, fSS_), this->nul(S,fSS_));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> 
kOmegaTrans<BasicMomentumTransportModel>::intermittency() const
{
    return tmp<volScalarField>
        (
            new volScalarField
            (
                "intermittency",
                min( max( sqrt(k_)*y_/(Agamma_*this->nu()) - 1.0, 0.0), 1.0)
            )
        );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> 
kOmegaTrans<BasicMomentumTransportModel>::nus(const volScalarField& S, const volScalarField& fSS) const
{
    return tmp<volScalarField>
        (
            new volScalarField
            (
                "nus",
                fSS * k_ / max(omega_, Clim_*S / sqrt(Cmu_))
            )
        );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> 
kOmegaTrans<BasicMomentumTransportModel>::nul(const volScalarField& S, const volScalarField& fSS) const
{
    return tmp<volScalarField>
        (
            new volScalarField
            (
                "nul",
                (1-fSS) * k_ / max(omega_, Clim_*S / al_)
            )
        );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> 
kOmegaTrans<BasicMomentumTransportModel>::fSS(const volScalarField& S, const volScalarField& W) const
{
    volScalarField fW( 1 - tanh( k_ / (CW_ * this->nu() * omega_) ) );
    volScalarField psi( tanh( - W * (S - W) / ( Cpsi_ * sqr(Cmu_*omega_) ) ) );
    volScalarField CSS( CS_ * (1 + CA_*fW*psi) );
    dimensionedScalar kMin("kMin", sqr(dimVelocity), ROOTVSMALL);
    
    return tmp<volScalarField>
        (
            new volScalarField
            (
                "fSS",
                exp( -sqr(CSS * this->nu() / y_) / max(k_,kMin) )
            )
        );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> 
kOmegaTrans<BasicMomentumTransportModel>::beta(const volTensorField& gradU) const
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


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> 
kOmegaTrans<BasicMomentumTransportModel>::kSource() const
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

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> 
kOmegaTrans<BasicMomentumTransportModel>::omegaSource() const
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

template<class BasicMomentumTransportModel>
kOmegaTrans<BasicMomentumTransportModel>::kOmegaTrans
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
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
    al_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "al",
            this->coeffDict_,
            0.45
        )
    ),
    Agamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Agamma",
            this->coeffDict_,
            12.0
        )
    ),
    CS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CS",
            this->coeffDict_,
            21.0
        )
    ),
    CA_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CA",
            this->coeffDict_,
            1.0
        )
    ),
    Cpsi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cpsi",
            this->coeffDict_,
            10.0
        )
    ),
    CW_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CW",
            this->coeffDict_,
            6.0
        )
    ),
    Csep_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Csep",
            this->coeffDict_,
            2.0
        )
    ),
    AV_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "AV",
            this->coeffDict_,
            550.0
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

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

template<class BasicMomentumTransportModel>
bool kOmegaTrans<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel> >::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        Clim_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());
        al_.readIfPresent(this->coeffDict());
        CS_.readIfPresent(this->coeffDict());
        CA_.readIfPresent(this->coeffDict());
        Cpsi_.readIfPresent(this->coeffDict());
        CW_.readIfPresent(this->coeffDict());
        Csep_.readIfPresent(this->coeffDict());
        AV_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kOmegaTrans<BasicMomentumTransportModel>::correct()
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

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicMomentumTransportModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S( sqrt(2*magSqr(dev(symm(tgradU())))) );
    volScalarField W( sqrt(2*magSqr(skew(tgradU()))) );
    volScalarField fSS_(this->fSS(S,W));

    volScalarField nus_( this->nus(S, fSS_) );
    
    volScalarField G(
        this->GName(),
        nus_ * ( dev(twoSymm(tgradU())) && tgradU() ) 
        // ( nus_ * dev(twoSymm(tgradU())) - 2./3.*k_*I ) && tgradU()
        // nus_ * sqr(S)
    );

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

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
            + rho * CDkOmega
            + omegaSource()
	    + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    volScalarField gammaInt = this->intermittency();
    tmp<volScalarField> Rv = sqr(y_) * S / this->nu();
    tmp<volScalarField> Fsep = min( max( Rv / (2.2*AV_) - 1.0, 0.0), 1.0);

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
        + fvm::div(alphaRhoPhi, k_)
        - fvm::laplacian(alpha*rho*DkEff(), k_)
        ==
        gammaInt * alpha*rho*G
        + (1.0 - gammaInt) * Csep_ * Fsep * this->nu() * sqr(S)
        - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
        - fvm::Sp(Cmu_*alpha*rho*omega_, k_)
        + kSource()
	+ fvOptions(alpha, rho, omega_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(nus_, this->nul(S, fSS_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
