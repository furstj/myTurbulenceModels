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

#include "mykkLOmega.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::fv(const volScalarField& Ret) const
{
    return tmp<volScalarField>(new volScalarField(
        "fv",
        1.0 - exp(-sqrt(Ret)/Av_)
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::fINT() const
{
    return tmp<volScalarField>(new volScalarField(
        "fINT",
        (
            min
            (
                kt_/(Cint_*(kl_ + kt_ + this->kMin_)),
                dimensionedScalar("1.0", dimless, 1.0)
            )
        )
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
{
    return tmp<volScalarField>(new volScalarField(
        "fSS",
        exp(-sqr(Css_*this->nu()*Omega/(kt_ + this->kMin_)))
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
{
    return tmp<volScalarField>(new volScalarField(
        "Cmu",
        1.0/(A0_ + As_*(S/(omega_ + this->omegaMin_)))
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    return tmp<volScalarField>(new volScalarField(
            "BetaTS",
            scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_)
        ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::lambdaT() const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaT",
            sqrt(this->kt_) / (this->omega_ + this->omegaMin_)
        ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::lambdaEff(const volScalarField& lambdaT) const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaEff",
            min( this->Clambda_ * y_, lambdaT)
        ));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return tmp<volScalarField>(new volScalarField(
        "fTaul",
        scalar(1)
        - exp
        (
            -CtauL_*ktL
            /
            (
                sqr
                (
                    lambdaEff*Omega
                    + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return tmp<volScalarField>(new volScalarField(
        "alphaT",
        fv*CmuStd_*sqrt(ktS)*lambdaEff
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::fOmega
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


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
{
    return tmp<volScalarField>(new volScalarField(
        "phiBP",
        min
        (
            max
            (
                kt_/this->nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        Omega.dimensions(),
                        ROOTVSMALL
                    )
                )
              - CbpCrit_,
                scalar(0)
            ),
            scalar(50.0)
        )
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return tmp<volScalarField>(new volScalarField(
        "phiNAT",
        max
        (
            ReOmega
          - CnatCrit_
            / (
                fNatCrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    ));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mykkLOmega<BasicTurbulenceModel>::D(const volScalarField& k) const
{
    return this->nu()*magSqr(fvc::grad(sqrt(k)));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void mykkLOmega<BasicTurbulenceModel>::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykkLOmega<BasicTurbulenceModel>::mykkLOmega
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

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            this->coeffDict_,
            2.12
        )
    ),
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            this->coeffDict_,
            6.75
        )
    ),
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            this->coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            this->coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            this->coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            this->coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            this->coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            this->coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            this->coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            this->coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
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
            0.12
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
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            this->coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
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
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            this->coeffDict_,
            1.5
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
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
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
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            this->coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            this->coeffDict_,
            1.17
        )
    ),
    kt_
    (
        IOobject
        (
            IOobject::groupName("kt", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kl_
    (
        IOobject
        (
            IOobject::groupName("kl", U.group()),
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_
        ),
        kt_*omega_ + D(kl_) + D(kt_)
    ),
    y_(wallDist::New(this->mesh_).y())
{
    bound(kt_, this->kMin_);
    bound(kl_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(epsilon_, this->epsilonMin_);

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

template<class BasicTurbulenceModel>
bool mykkLOmega<BasicTurbulenceModel>::read()
{
    if (  eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        A0_.readIfPresent(this->coeffDict());
        As_.readIfPresent(this->coeffDict());
        Av_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Anat_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Ats_.readIfPresent(this->coeffDict());
        CbpCrit_.readIfPresent(this->coeffDict());
        Cnc_.readIfPresent(this->coeffDict());
        CnatCrit_.readIfPresent(this->coeffDict());
        Cint_.readIfPresent(this->coeffDict());
        CtsCrit_.readIfPresent(this->coeffDict());
        CrNat_.readIfPresent(this->coeffDict());
        C11_.readIfPresent(this->coeffDict());
        C12_.readIfPresent(this->coeffDict());
        CR_.readIfPresent(this->coeffDict());
        CalphaTheta_.readIfPresent(this->coeffDict());
        Css_.readIfPresent(this->coeffDict());
        CtauL_.readIfPresent(this->coeffDict());
        Cw1_.readIfPresent(this->coeffDict());
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        CwR_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        CmuStd_.readIfPresent(this->coeffDict());
        Sigmak_.readIfPresent(this->coeffDict());
        Sigmaw_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void mykkLOmega<BasicTurbulenceModel>::validate()
{}


template<class BasicTurbulenceModel>
void mykkLOmega<BasicTurbulenceModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U_ = this->U_;
    const dimensionedScalar& kMin_ = this->kMin_;
    const dimensionedScalar& omegaMin_ = this->omegaMin_;
    const dimensionedScalar& epsilonMin_ = this->epsilonMin_;
    volScalarField& nut_ = this->nut_;
    volScalarField& omega_ = this->omega_;
    volScalarField& kt_ = this->kt_;
    volScalarField& kl_ = this->kl_;
    
    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    const volScalarField lambdaT_ = lambdaT();
  
    const volScalarField lambdaEff_ = lambdaEff(lambdaT_);

    const volScalarField fw
      (  "fw",
        pow
        (
            lambdaEff_
           /(lambdaT_ + dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
    );

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

    const volScalarField ktS("ktS", fSS(Omega)*fw*kt_);

    const volScalarField nuts
      (  "nuts",
        fv(sqr(fw)*kt_/this->nu()/(omega_ + omegaMin_))
       *fINT()
       *Cmu(sqrt(S2))*sqrt(ktS)*lambdaEff_
    );
    const volScalarField Pkt("Pkt", nuts*S2);

    const volScalarField ktL("ktL", kt_ - ktS);
    const volScalarField ReOmega("ReOmega", sqr(y_)*Omega/this->nu());
    const volScalarField nutl
      (  "nutl",
        min
        (
            C11_*fTaul(lambdaEff_, ktL, Omega)*Omega*sqr(lambdaEff_)
           *sqrt(ktL)*lambdaEff_/this->nu()
          + C12_*BetaTS(ReOmega)*ReOmega*sqr(y_)*Omega
        ,
            0.5*(kl_ + ktL)/(sqrt(S2) + omegaMin_)
        )
    );

    const volScalarField Pkl("Pkl", nutl*S2);

    const volScalarField alphaTEff
      ( "alphaTEff",
        alphaT(lambdaEff_, fv(sqr(fw)*kt_/this->nu()/(omega_ + omegaMin_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);

    const volScalarField Rbp
      (   "Rbp",
        CR_*(1.0 - exp(-phiBP(Omega)()/Abp_))*omega_
       /(fw + fwMin)
    );

    const volScalarField fNatCrit("fNatCrit", 1.0 - exp(-Cnc_*sqrt(kl_)*y_/this->nu()));

    // Natural source term divided by kl_
    const volScalarField Rnat
      (   "Rnat",
        CrNat_*(1.0 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );


    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
     fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(alphaTEff), omega_)
     ==
        alpha*rho*Cw1_*Pkt*omega_/(kt_ + kMin_)
      - fvm::SuSp
        (
            alpha*rho*(1.0 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)/(kt_ + kMin_)
          , omega_
        )
      - fvm::Sp(alpha*rho*Cw2_*sqr(fw)*omega_, omega_)
      + alpha*rho*(
            Cw3_*fOmega(lambdaEff_, lambdaT_)*alphaTEff*sqr(fw)*sqrt(kt_)
        )().dimensionedInternalField()/pow3(y_.dimensionedInternalField())
    );

    omegaEqn().relax();
    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    const volScalarField Dl(D(kl_));

    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
     fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha*rho*this->nu(), kl_)
     ==
        alpha*rho*Pkl
      - fvm::Sp(alpha*rho*(Rbp + Rnat + Dl/(kl_ + kMin_)), kl_)
    );

    klEqn().relax();
    klEqn().boundaryManipulate(kl_.boundaryField());

    solve(klEqn);
    bound(kl_, kMin_);


    const volScalarField Dt(D(kt_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
     fvm::ddt(alpha, rho, kt_)
      + fvm::div(alphaRhoPhi, kt_)
      - fvm::laplacian(alpha*rho*DkEff(alphaTEff), kt_)
     ==
        alpha*rho*Pkt
      + alpha*rho*(Rbp + Rnat)*kl_
      - fvm::Sp(alpha*rho*(omega_ + Dt/(kt_+ kMin_)), kt_)
    );

    ktEqn().relax();
    ktEqn().boundaryManipulate(kt_.boundaryField());

    solve(ktEqn);
    bound(kt_, kMin_);


    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = kt_*omega_ + Dl + Dt;
    bound(epsilon_, epsilonMin_);


    // Re-calculate turbulent viscosity
    nut_ = nuts + nutl;
    nut_.correctBoundaryConditions();

#ifdef HAVE_ALPHAT
    // Re-calculate thermal diffusivity
    this->alphat_ = alpha * rho * 
      ( 
       fw * kt_ / max(kt_ + kl_,kMin_) * nuts / this->PrTheta_ 
       + (scalar(1.0) - fw) * CalphaTheta_ * sqrt(kt_) * lambdaEff_
	);
#endif    

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
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
