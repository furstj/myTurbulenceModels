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

#include "kv2Omega.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::fv(const volScalarField& Ret) const
{
    return tmp<volScalarField>(new volScalarField(
        "fv",
        1.0 - exp(-sqrt(Ret)/Anu_)
    ));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::fINT() const
{
    return tmp<volScalarField>(new volScalarField(
        "fINT",
        (
            min
            (
                v2_ / (CINT_* ( k_ + this->kMin_)),
                dimensionedScalar("1.0", dimless, 1.0)
            )
        )
    ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::fSS(const volScalarField& W) const
{
    return tmp<volScalarField>(new volScalarField(
        "fSS",
        exp( -sqr(CSS_ * this->nu() * W / (v2_ + this->kMin_) ) )
    ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::Cmu(const volScalarField& S) const
{
    return tmp<volScalarField>(new volScalarField(
        "Cmu",
        1.0/(A0_ + AS_*(S/max(omega_,this->omegaMin_)))
    ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::betaTS(const volScalarField& ReW) const
{
    return tmp<volScalarField>(new volScalarField(
            "betaTS",
            scalar(1) - exp( -sqr( max(ReW - CTScrit_, scalar(0)) ) /ATS_ )
        ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::lambdaT() const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaT",
            sqrt(this->v2_) / max(this->omega_, this->omegaMin_)
    ));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::lambdaEff(const volScalarField& lambdaT) const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaEff",
            min( this->Clambda_ * y_, lambdaT)
        ));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& v2l,
    const volScalarField& W
) const
{
    const dimensionedScalar vMin("ROOTVSMALL",  dimLength*inv(dimTime), ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "fTaul",
        scalar(1)
        - exp
        (
            -Ctau1_ * v2l
            /
            sqr( max( lambdaEff * W, vMin ) )
        )
    ));
};
        
        
template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& v2s
) const
{
    return tmp<volScalarField>(new volScalarField(
        "alphaT",
        fv * betaStar_ * sqrt(v2s) * lambdaEff
    ));
}


// TODO
/*
template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return tmp<volScalarField>(new volScalarField(
        "fOmega",        
        scalar(1)
        - exp
        (
            -0.41
            *pow4
            (
                lambdaEff
                / (
                    lambdaT
                    + dimensionedScalar
                    (
                        "ROTVSMALL",
                        lambdaT.dimensions(),
                        ROOTVSMALL
                    )
                )
            )
        )
    ));
}

    */

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::phiBP(const volScalarField& W) const
{
    const dimensionedScalar wMin("ROOTVSMALL", inv(dimTime), ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "phiBP",
        min
        (
            max
            (
                v2_/ ( this->nu() * max(W, wMin) ) - CBPcrit_, 
                scalar(0)
            ),
            scalar(50.0)
        )
    ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::phiNAT
(
    const volScalarField& ReW,
    const volScalarField& fNATcrit
) const
{
    const dimensionedScalar small("ROTVSMALL", dimless, ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "phiNAT",
        max
        (
            ReW - CNATcrit_ / max(fNATcrit, small),
            scalar(0)
        )
    ));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::D(const volScalarField& k) const
{
    return 2.0*this->nu()*magSqr(fvc::grad(sqrt(k)));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kv2Omega<BasicMomentumTransportModel>::F1() const
{
    const volScalarField CDkOmega(
        "CDkOmega",
        max(
            (2 * this->sigmaW2_)*
            (fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        )
    );

    const volScalarField arg1(
        "arg1",
        min(
            max(
                sqrt(this->v2_)/(this->omega_*this->y_),
                500.0 * this->nu() * this->betaStar_ / (sqr(this->y_) * this->omega_)
            ),
            4.0 * this->sigmaW2_ * this->k_ / (CDkOmega * sqr(this->y_))
            )
    );

    return tmp<volScalarField>( new volScalarField(
        "F1",
        tanh(pow4(arg1))
    ));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kv2Omega<BasicMomentumTransportModel>::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kv2Omega<BasicMomentumTransportModel>::kv2Omega
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

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    AS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "AS",
            this->coeffDict_,
            2.12
        )
    ),
    Anu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anu",
            this->coeffDict_,
            3.8
        )
    ),
    ABP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ABP",
            this->coeffDict_,
            0.2
        )
    ),
    ANAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ANAT",
            this->coeffDict_,
            200
        )
    ),
    ATS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ATS",
            this->coeffDict_,
            200
        )
    ),
    CBPcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CBPcrit",
            this->coeffDict_,
            1.5
        )
    ),
    CNC_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNC",
            this->coeffDict_,
            0.1
        )
    ),
    CNATcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNATcrit",
            this->coeffDict_,
            1450
        )
    ),
    CINT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CINT",
            this->coeffDict_,
            0.95
        )
    ),
    CTScrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTScrit",
            this->coeffDict_,
            1000
        )
    ),
    CRNAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CRNAT",
            this->coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            this->coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            this->coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            this->coeffDict_,
            0.32
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            this->coeffDict_,
            0.035
        )
    ),
    CSS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSS",
            this->coeffDict_,
            3.0
        )
    ),
    Ctau1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau1",
            this->coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            this->coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.92
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            this->coeffDict_,
            1.15
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            2.495
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
    PrTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "PrTheta",
            this->coeffDict_,
            0.85
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1
        )
    ),
    sigmaW_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaW",
            this->coeffDict_,
            1.17
        )
    ),
    sigmaW2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaW2",
            this->coeffDict_,
            1.856
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
    v2_
    (
        IOobject
        (
            IOobject::groupName("v2", U.group()),
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
    bound(v2_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        // Evaluating nut_ is complex so start from the field read from file
        this->nut_.correctBoundaryConditions();

        this->printCoeffs(type);

	if (debug) 
	  Info << "Debug switch is on!" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kv2Omega<BasicMomentumTransportModel>::read()
{
    if (  eddyViscosity<RASModel<BasicMomentumTransportModel> >::read())
    {    
        A0_.readIfPresent(this->coeffDict());
        AS_.readIfPresent(this->coeffDict());
        Anu_.readIfPresent(this->coeffDict());
        ABP_.readIfPresent(this->coeffDict());
        ANAT_.readIfPresent(this->coeffDict());
        ATS_.readIfPresent(this->coeffDict());
        CBPcrit_.readIfPresent(this->coeffDict());
        CNC_.readIfPresent(this->coeffDict());
        CNATcrit_.readIfPresent(this->coeffDict());
        CINT_.readIfPresent(this->coeffDict());
        CTScrit_.readIfPresent(this->coeffDict());
        CRNAT_.readIfPresent(this->coeffDict());
        C11_.readIfPresent(this->coeffDict());
        C12_.readIfPresent(this->coeffDict());
        CR_.readIfPresent(this->coeffDict());
        CalphaTheta_.readIfPresent(this->coeffDict());
        CSS_.readIfPresent(this->coeffDict());
        Ctau1_.readIfPresent(this->coeffDict());
        Cw1_.readIfPresent(this->coeffDict());
        Cw2_.readIfPresent(this->coeffDict());
        CwR_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        PrTheta_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaW_.readIfPresent(this->coeffDict());
        sigmaW2_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kv2Omega<BasicMomentumTransportModel>::validate()
{}


// TODO
template<class BasicMomentumTransportModel>
void kv2Omega<BasicMomentumTransportModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    eddyViscosity<RASModel<BasicMomentumTransportModel> >::correct();

    // Local references
    const alphaField& alpha_ = this->alpha_;
    const rhoField& rho_ = this->rho_;
    const surfaceScalarField& alphaRhoPhi_ = this->alphaRhoPhi_;
    const volVectorField& U_ = this->U_;
    const dimensionedScalar& kMin_ = this->kMin_;
    const dimensionedScalar& omegaMin_ = this->omegaMin_;
    volScalarField& nut_ = this->nut_;
    volScalarField& omega_ = this->omega_;
    volScalarField& k_ = this->k_;
    volScalarField& v2_ = this->v2_;

    const dimensionedScalar small("ROTVSMALL", dimless, ROOTVSMALL);

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();
    const volScalarField W(sqrt(2.0)*mag(skew(gradU)));
    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

    const volScalarField lambdaT_ = lambdaT();
  
    const volScalarField lambdaEff_ = lambdaEff(lambdaT_);

    const volScalarField fW
      (  "fW",
        pow
        (
            lambdaEff_
            /max(lambdaT_,dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
      );
    


    const volScalarField v2s("v2s", fSS(W) * fW * v2_);

    tmp<volScalarField> ReT( sqr(fW) * v2_ / (this->nu() * max(omega_, omegaMin_)) );
    
    const volScalarField fv_ = fv(ReT);

    const volScalarField nuTs(  
        "nuTs",
        fv_ * fINT() * Cmu(sqrt(S2)) * sqrt(v2s) * lambdaEff_
    );


    const volScalarField v2l("v2l", v2_ - v2s);
    const volScalarField ReW("ReW", sqr(y_) * W / this->nu() );
    const volScalarField nuTl
      (  "nuTl",
      min
      (
          C11_* fTaul(lambdaEff_, v2l, W) * W * sqr(lambdaEff_) 
          * sqrt(v2l) * lambdaEff_ /this->nu()
          + C12_ * betaTS(ReW) * pow4(lambdaEff_/Clambda_) * sqr(W) / this->nu()
          ,
          0.5*(k_ - v2s)/max(sqrt(S2), omegaMin_)
      )
      );

    const volScalarField F1star("F1star", scalar(1) - (scalar(1)-F1())*fSS(W) );

    const volScalarField alphaT_ = alphaT(lambdaEff_, fv_, v2s);

    // Bypass source term divided by (k-v2)
    const volScalarField RBP
      (   "RBP",
        CR_ * (1.0 - exp(-phiBP(W)/ABP_))*omega_
      /max(fW, small)
    );

    const volScalarField fNATcrit("fNATcrit", 1.0 - exp(-CNC_*sqrt(k_)*y_/this->nu()));

    // Natural source term divided by (k-v2)
    const volScalarField RNAT
      (   "RNAT",
        CRNAT_ * (1.0 - exp(-phiNAT(ReW, fNATcrit)/ANAT_)) * W
    );


    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha_, rho_, omega_)
            + fvm::div(alphaRhoPhi_, omega_)
            - fvm::laplacian(alpha_*rho_*DomegaEff(alphaT_), omega_)
            ==
            alpha_ * rho_ * (
                Cw1_ * omega_ / max(v2_, kMin_) * nuTs * S2
                - fvm::SuSp(
                    (1.0 - CwR_/max(fW,small)) * (k_ - v2_) * (RBP + RNAT)
                    /max(v2_, kMin_)
                    , omega_)
                - fvm::Sp(Cw2_*sqr(fW)*omega_, omega_)
                + betaStar_ * 2 * (1.0 - F1star) * sigmaW2_ / max(omega_, omegaMin_) *
                (fvc::grad(k_) & fvc::grad(omega_) )
            )
        );
    
    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> v2Eqn
        (
            fvm::ddt(alpha_, rho_, v2_)
            + fvm::div(alphaRhoPhi_, v2_)
            - fvm::laplacian(alpha_*rho_*DkEff(alphaT_), v2_)
            ==
            alpha_ * rho_ * (
                nuTs * S2 
                + (RBP + RNAT) * k_ 
                - fvm::Sp(RBP + RNAT, v2_)
                - fvm::Sp(omega_ + D(v2_)/max(v2_,kMin_), v2_)
            )
    );

    v2Eqn.ref().relax();
    v2Eqn.ref().boundaryManipulate(v2_.boundaryFieldRef());

    solve(v2Eqn);
    bound(v2_, kMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha_, rho_, k_)
            + fvm::div(alphaRhoPhi_, k_)
            - fvm::laplacian(alpha_*rho_*DkEff(alphaT_), k_)
            ==
            alpha_*rho_*(
                nut_ * S2
                - fvm::Sp((omega_*min(k_,v2_) + D(k_))/max(k_,kMin_), k_)
            )
        );
    
    kEqn.ref().relax();
    kEqn.ref().boundaryManipulate(k_.boundaryFieldRef());

    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate turbulent viscosity
    nut_ = nuTs + nuTl;
    nut_.correctBoundaryConditions();

#ifdef HAVE_ALPHAT
    // Re-calculate thermal diffusivity
    this->alphat_ = alpha_ * rho_ * 
      ( 
       fW * v2_ / max(k_,kMin_) * nuTs / this->PrTheta_ 
       + (scalar(1.0) - fW) * CalphaTheta_ * sqrt(v2_) * lambdaEff_
	);
#endif    

    /*
    if (debug && this->runTime_.outputTime()) {
      lambdaEff_.write();
      fw.write();
      ktS.write();
      nuts.write();
      Pkt.write();
      ktL.write();
      ReOmega.write();
      nutl.write();
      Pkl.write();
      alphaTEff.write();
      Rbp.write();
      fNatCrit.write();
      Rnat.write();
      BetaTS(ReOmega)().write(); 
      fINT()().write();
      fSS(Omega)().write();
      Cmu(sqrt(S2))().write();
      fTaul(lambdaEff_,ktL,Omega)().write(); 
      fOmega(lambdaEff_,lambdaT_)().write();
      phiBP(Omega)().write();
      phiNAT(ReOmega,fNatCrit)().write();
      y_.write();
    }
        */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
