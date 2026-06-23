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

#include "EARSMBase.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ScaleModel>
EARSMBase<ScaleModel>::EARSMBase
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
    ScaleModel(
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    CTau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTau_",
            this->coeffDict_,
            6.0
        )
    ),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            -0.72
        )
    ),

    nearWallDamping_
    (
        Switch::lookupOrAddToDict
        (
            "nearWallDamping",
            this->coeffDict_,
            false
        )
    ),

    B2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B2",
            this->coeffDict_,
            1.8
        )
    ),

    CyA_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CyA",
            this->coeffDict_,
            0.092
        )
    ),

    CyB_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CyB",
            this->coeffDict_,
            1.2e-4
        )
    ),

    aEx_
    (
        IOobject
        (
            "aEx",
            this->U_.time().timeName(),
            this->U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->U_.mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor(Foam::Zero))
    )
{
}


template<class ScaleModel>
volScalarField EARSMBase<ScaleModel>::N(
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

template<class ScaleModel>
void EARSMBase<ScaleModel>::correctAnisotropy(const volTensorField& gradU)
{
    volScalarField tau(
        max
        (
            1.0 / (this->betaStar_ * this->omega_),
            this->CTau_ * sqrt(this->nu() / (this->betaStar_ * max(this->k_, this->kMin_) * this->omega_))
        ));

    volSymmTensorField S("S", tau * dev(symm(gradU)));
    volTensorField     W("W", -tau * skew(gradU));
    // NOTE: Wij = 1/2(dui/dxj - duj/dxi) = - skew(grad(U))

     volScalarField IIS  = tr(S & S);
     
     /*
     if (this->curvatureCorrection_)
    {
        const volVectorField& U = this->U_;
        const surfaceScalarField& phi = this->phi_;
        const rhoField& rho = this->rho_;

        volSymmTensorField DSDt = tau*dev(symm(fvc::grad(fvc::ddt(U))))
                + fvc::div(phi/fvc::interpolate(rho),S,"div(phiv,S)");
        volVectorField SDeps = *(skew(S & DSDt))*2;
        volScalarField IIIS = tr(S & S & S);

        volTensorField B = (pow(IIS,2)*I + 12*IIIS*S + 6*IIS*(S&S))/
            max(2*pow(IIS,3) - 12*pow(IIIS,2), 1.e-10);

        volVectorField BSDeps = B & SDeps;
        volTensorField Wr = *(BSDeps);
        W -= (tau/this->A0_)*Wr;
    }
        */

    volScalarField IIW  = tr(W & W);
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

    
    this->aEx_ =  symm(
          beta3 * ( (W & W) - (1.0/3.0) * IIW * I )
        + beta4 * ( (S & W) - (W & S) )
        + beta6 * ( (S & W & W) + (W & W & S) - IIW * S - (2.0/3.0) * IV * I)
        + beta9 * ( (W & S & W & W) - (W & W & S & W) )
    );

    // this->nonlinearStress_.correctBoundaryConditions();
    // BasicTurbulenceModel::correctNut();

}


template<class ScaleModel>
void EARSMBase<ScaleModel>::correctNut()
{
    Info << "********** EARSMBase: Correcting nut()" << endl;
    const volVectorField& U = this->U_;
    tmp<volTensorField>   tgradU = fvc::grad(U);
    this->correctAnisotropy(tgradU);
}


template<class ScaleModel>
void EARSMBase<ScaleModel>::correctNut(const volScalarField& S2)
{
    Info << "********** EARSMBase: Correcting nut(S2)" << endl;
    const volVectorField& U = this->U_;
    tmp<volTensorField>   tgradU = fvc::grad(U);
    this->correctAnisotropy(tgradU);
}


template<class ScaleModel>
tmp<volSymmTensorField> EARSMBase<ScaleModel>::R() const
{
    Info << "********** EARSMBase: Computing R" << endl;
    return ScaleModel::R();
}


template<class ScaleModel>
tmp<volSymmTensorField> EARSMBase<ScaleModel>::devRhoReff() const
{    
    return devRhoReff(this->U_);
}


template<class ScaleModel>
tmp<volSymmTensorField> EARSMBase<ScaleModel>::devRhoReff
(
    const volVectorField& U
) const
{
    Info << "********** EARSMBase: Computing devRhoReff with U" << endl;
 
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ScaleModel::devRhoReff(U) + this->alpha_*this->rho_*this->k()*this->aEx_
        )
    );

    return ScaleModel::devRhoReff(U);
}


template<class ScaleModel>
tmp<fvVectorMatrix> EARSMBase<ScaleModel>::divDevRhoReff(volVectorField& U) const
{
    Info << "********** EARSMBase: Computing divDevRhoReff(U)" << endl;
    return ScaleModel::divDevRhoReff(U) + this->alpha_*this->rho_*this->k()*fvc::div(this->aEx_);
}


template<class ScaleModel>
tmp<fvVectorMatrix> EARSMBase<ScaleModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    Info << "********** EARSMBase: Computing divDevRhoReff(rho, U)" << endl;
    return ScaleModel::divDevRhoReff(rho, U) + this->alpha_*rho*this->k()*fvc::div(this->aEx_);
}


template<class ScaleModel>
void EARSMBase<ScaleModel>::validate()
{
    Info << "********** EARSMBase: Validating" << endl; 
    this->correctNut();
    ScaleModel::validate();
}


template<class ScaleModel>
void EARSMBase<ScaleModel>::correct()
{
    Info << "********** EARSMBase: Correcting" << endl;
    ScaleModel::correct();
}


template<class ScaleModel>
bool EARSMBase<ScaleModel>::read(const word& type)
{
    Info << "********** EARSMBase: Reading" << endl;
    return ScaleModel::read();
}


} // End namespace RASModels
} // End namespace Foam
