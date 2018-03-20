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

#include "EARSM.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> EARSM<BasicTurbulenceModel>::fMix
(
    const volScalarField& gradKgradOmegaByOmega
) const
{
    tmp<volScalarField> Gamma = min
        (
            max
            (
                sqrt(k_) / (betaStar_ * omega_ * y_),
                500.0 * this->nu() / (omega_ * sqr(y_))
            ),
            20.0 * k_ /
            max( sqr(y_) * gradKgradOmegaByOmega, 200.0 * kInf_)
        );
    
    
    return tanh(1.5 * pow4(Gamma));
}


template<class BasicTurbulenceModel>
void EARSM<BasicTurbulenceModel>::correctNut()
{
  correctNonlinearStress(fvc::grad(this->U_));
}


template<class BasicTurbulenceModel>
volScalarField EARSM<BasicTurbulenceModel>::N
(
    const volScalarField& A3p,
    const volScalarField& P1,
    const volScalarField& P2
) const
{
    volScalarField N = A3p / 3.0;

    forAll(N, i)
    {
        if (P2[i] < 0)
        {
            N[i] += 2*pow(sqr(P1[i]) - P2[i], 1./6.)
                * cos( 1./3.*acos( P1[i]/sqrt(sqr(P1[i]) - P2[i]))) ;
        }
        else
        {
            scalar a = max(P1[i] + sqrt(P2[i]), 0.0);
            scalar b = P1[i] - sqrt(P2[i]);
            N[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
        }
    };

    forAll(N.boundaryField(), patchi)
    {
        fvPatchScalarField& pN = N.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pP1 = P1.boundaryField()[patchi];
        const fvPatchScalarField& pP2 = P2.boundaryField()[patchi];

        forAll(pN, i)
        {
            if (pP2[i] < 0)
            {
                pN[i] += 2*pow(sqr(pP1[i]) - pP2[i], 1./6.)
                    * cos( 1./3.*acos( pP1[i]/sqrt(sqr(pP1[i]) - pP2[i]))) ;
            }
            else
            {
                scalar a = max(pP1[i] + sqrt(pP2[i]), 0.0);
                scalar b = pP1[i] - sqrt(pP2[i]);
                pN[i] += pow(a, 1./3.) + sign(b) * pow(fabs(b), 1./3.);
            }
        };
        
    };

    return N;
}


template<class BasicTurbulenceModel>
void EARSM<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
    scalar Ctau = 6.0;
    volScalarField tau(
        max
        (
            1.0 / (this->betaStar_ * this->omega_),
            Ctau * sqrt(this->nu() / (this->betaStar_ * max(this->k_, this->kMin_) * this->omega_))
        ));
    
    volSymmTensorField S(tau * dev(symm(gradU)));
    volTensorField     W(tau * skew(gradU));

    volScalarField IIS  = tr(S & S);
    volScalarField IIW  = tr(W & W);
    // volScalarField IIIS = tr(S & S & S);
    volScalarField IV   = tr(S & W & W);
    
    scalar Neq = 81.0 / 20.0;
    scalar CDiff = 2.2;
    volScalarField beta1eq = - 6.0/5.0 * Neq / (sqr(Neq) - 2*IIW);
    volScalarField A3p = 9.0/5.0 + 9.0/4.0 * CDiff * max(1 + beta1eq*IIS, 0.0);
    volScalarField P1 = (sqr(A3p)/27 + (9.0/20.0)*IIS - (2.0/3.0)*IIW) * A3p;
    volScalarField P2 = sqr(P1) - pow3(sqr(A3p)/9 + 0.9*IIS + (2.0/3.0)*IIW);
    
    volScalarField N = this->N(A3p, P1, P2);

    volScalarField Q = 5.0/6.0*(sqr(N) - 2*IIW)*(2*sqr(N)-IIW);

    volScalarField beta1 = -N*(2.0*sqr(N) - 7.0*IIW) / Q;
    volScalarField beta3 = -12.0 * IV / (N * Q);
    volScalarField beta4 = -2.0 * (sqr(N)  - 2.0*IIW) / Q;
    volScalarField beta6 = -6.0 * N / Q;
    volScalarField beta9 =  6.0 / Q;

    
    volScalarField Cmu = - 0.5 * (beta1 + IIW * beta6);

    this->nut_ = Cmu * this->k_ * tau;
    this->nut_.correctBoundaryConditions();

    
    this->nonlinearStress_ = this->k_ * symm(
        beta3 * ( (W & W) - (1.0/3.0) * IIW * I )
        + beta4 * ( (S & W) - (W & S) )
        + beta6 * ( (S & W & W) + (W & W & S) - IIW * S - (2.0/3.0) * IV * I)
        + beta9 * ( (W & S & W & W) - (W & W & S & W) )
    );

    this->nonlinearStress_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EARSM<BasicTurbulenceModel>::EARSM
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
    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >
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
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.518
        )
    ),

    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),

    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.0747
        )
    ),

    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            1.1
        )
    ),

    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.1
        )
    ),

    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.53
        )
    ),

    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            1.0
        )
    ),

    alphaD1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD1",
            this->coeffDict_,
            1.0
        )
    ),

    alphaD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD2",
            this->coeffDict_,
            0.4
        )
    ),

    kInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kInf_",
            this->coeffDict_,
            sqr(dimVelocity),
            1.e-10
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
    ),
    y_(wallDist::New(this->mesh_).y())
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
bool EARSM<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {    
        betaStar_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        alphaD1_.readIfPresent(this->coeffDict());
        alphaD2_.readIfPresent(this->coeffDict());
        kInf_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void EARSM<BasicTurbulenceModel>::validate()
{
    this->correctNut();
}


template<class BasicTurbulenceModel>
void EARSM<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField>   tgradU = fvc::grad(U);

    volScalarField G
    (
        this->GName(),
        (nut * dev(twoSymm(tgradU())) - this->nonlinearStress_) && tgradU()
    );
    

    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField gradKgradOmegaByOmega
    (
        (fvc::grad(k_) & fvc::grad(omega_)) / omega_
    );

    volScalarField fMix( this->fMix(gradKgradOmegaByOmega) );

    {
        volScalarField gamma( this->gamma(fMix) );
        volScalarField beta( this->beta(fMix) );
        volScalarField alphaD( this->alphaD(fMix) );

        tmp<volScalarField> CDOmega = alphaD * alpha * rho *
            max( gradKgradOmegaByOmega, dimensionedScalar("zero", inv(sqr(dimTime)), 0.0));
        
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(fMix), omega_)
         ==
            alpha*rho*gamma * omega_ / max(k_, this->kMin_) * G
            - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
            - fvm::Sp(alpha*rho*beta*omega_, omega_)
            + CDOmega
            + fvOptions(alpha, rho, omega_)            
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }
    
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(fMix), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
        - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNonlinearStress(tgradU());
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
