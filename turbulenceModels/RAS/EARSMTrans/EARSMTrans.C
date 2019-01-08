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

#include "EARSMTrans.H"
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
void EARSMTrans<BasicTurbulenceModel>::correctNut()
{
  correctNonlinearStress(fvc::grad(this->U_));
}


template<class BasicTurbulenceModel>
volScalarField EARSMTrans<BasicTurbulenceModel>::N
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
void EARSMTrans<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
    volScalarField tau(
        max
        (
            1.0 / (this->betaStar_ * this->omega_),
            Ctau_ * sqrt(this->nu() / (this->betaStar_ * max(this->k_, this->kMin_) * this->omega_))
        ));
    
    volSymmTensorField S(tau * dev(symm(gradU)));
    volTensorField     W(-tau * skew(gradU));

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

    BasicTurbulenceModel::correctNut();

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EARSMTrans<BasicTurbulenceModel>::EARSMTrans
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

    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            2.0/3.0
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

    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            3.0/40.0
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

    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            5.0/9.0
        )
    ),

    Ctau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau",
            this->coeffDict_,
            6.0
        )
    ),

    CSS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSS",
            this->coeffDict_,
            2.75 //3.25
        )
    ),

    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            14.5/8.
        )
    ),

    AT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "AT",
            this->coeffDict_,
            1.0
        )
    ),

    productionLimiter_
    (
        Switch::lookupOrAddToDict
        (
            "productionLimiter",
            this->coeffDict_,
            false
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
bool EARSMTrans<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {    
        betaStar_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        Ctau_.readIfPresent(this->coeffDict());
        CSS_.readIfPresent(this->coeffDict());
        CT_.readIfPresent(this->coeffDict());
        AT_.readIfPresent(this->coeffDict());
        productionLimiter_.readIfPresent("productionLimiter", this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void EARSMTrans<BasicTurbulenceModel>::validate()
{
    this->correctNut();
}


template<class BasicTurbulenceModel>
void EARSMTrans<BasicTurbulenceModel>::correct()
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

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField>   tgradU = fvc::grad(U);
    volScalarField  W( sqrt(2*magSqr(skew(tgradU()))) );
    volScalarField fSS = exp( -sqr(CSS_*this->nu()*W/max(k_,this->kMin_)) );
    volScalarField zetaT = max(k_/(this->nu()*W) - CT_, 0.0);
    volScalarField gammaInt = min(zetaT / AT_, 1.0);

    volScalarField G
    (
        this->GName(),
        (fSS * nut * dev(twoSymm(tgradU())) - this->nonlinearStress_) && tgradU()
    );

    if (productionLimiter_)
    {
        G = min(G, 10*betaStar_*k_*omega_);
    }
    
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega = max(
        this->sigmaD_ / max(omega_, this->omegaMin()) * (fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar("0",inv(sqr(dimTime)), 0.0)
    );

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha * rho * this->DomegaEff(), omega_)
     ==
        this->gamma_ * alpha * rho * G * omega_/max(k_,this->kMin())
      - fvm::SuSp(((2.0/3.0)*this->alphaOmega_)*alpha * rho * divU, omega_)
      - fvm::Sp(this->beta_ * alpha * rho * omega_, omega_)
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
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        gammaInt * alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
        - fvm::Sp(this->betaStar_*alpha*rho*omega_, k_)
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
