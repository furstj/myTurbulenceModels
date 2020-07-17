/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "gammaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::F1(const volScalarField& CDkOmega) const
{
  return max(
	     kOmegaSST<BasicTurbulenceModel>::F1(CDkOmega),
	     exp(-sqr(pow4(this->y_*sqrt(this->k_)/(scalar(120)*this->nu()))))
	     );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::ReThetac() const
{
  return CTU1_ + CTU2_*exp(-CTU3_*TuL()*FPG() );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::Fonset(const volScalarField& S) const
{
    return tmp<volScalarField>
	(
	    new volScalarField
     	    (
                IOobject
       	        (
	            "Fonset",
    	            this->runTime_.timeName(),
	            this->mesh_,
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
		max
		(
                    min(Fonset1(S), 2.0) - max(1.0-pow3(Rt()/3.5),0.0),
                    0.0
		)
   	    )
	);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::Fonset1(const volScalarField& S) const
{
    return sqr(this->y_)*S/this->nu() / (2.2*ReThetac());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::Fturb() const
{
    return exp(-pow4(Rt()/2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::TuL() const
{
    return min(100 * sqrt(2.0/3.0*this->k_) / (this->omega_ * this->y_), 100.0);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSST<BasicTurbulenceModel>::FPG() const
{
    volVectorField n = fvc::grad(this->y_);
    volScalarField lambdaThetaL = 
        min( 1.0, max( -1.0, 
        -7.57e-3 * ( fvc::grad(this->U_ & n) & n) * sqr(this->y_) / this->nu() + 0.0128));

    tmp<volScalarField> tFPG(new volScalarField("FPG", lambdaThetaL));
 
    volScalarField& FPG_ = tFPG.ref();
    forAll(FPG_, i) {
        if (lambdaThetaL[i]>=0) 
            FPG_[i] = min(1 + CPG1_.value()*lambdaThetaL[i], CPG1lim_.value());
        else
            FPG_[i] = min(1 + CPG2_.value()*lambdaThetaL[i] + 
            CPG3_.value()*min(lambdaThetaL[i]+0.0681,0), 
            CPG2lim_.value());
        FPG_[i] = max(FPG_[i], 0.0);
    }

    return tFPG;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
gammaSST<BasicTurbulenceModel>::gammaSST
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
            
            kOmegaSST<BasicTurbulenceModel>
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

            Flength_
            (
               dimensioned<scalar>::lookupOrAddToDict
              (
            "Flength",
            this->coeffDict_,
            100.0
        )
    ),
    ca2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50.0
        )
    ),
    sigmaGamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaGamma",
            this->coeffDict_,
            1.0
        )
    ),
    CPG1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1",
            this->coeffDict_,
            14.68
        )
    ),
    CPG1lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1lim",
            this->coeffDict_,
            1.5
        )
    ),
    CPG2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2",
            this->coeffDict_,
            -7.34
        )
    ),
    CPG3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG3",
            this->coeffDict_,
            0.0
        )
    ),
    CPG2lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2lim",
            this->coeffDict_,
            3.0
        )
    ),
    CTU1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU1",
            this->coeffDict_,
            100.
        )
    ),
    CTU2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU2",
            this->coeffDict_,
            1000.
        )
    ),
    CTU3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU3",
            this->coeffDict_,
            1.
        )
    ),
    ReThetacLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ReThetacLim",
            this->coeffDict_,
            1100.
        )
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            1.
        )
    ),
    CSEP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSEP",
            this->coeffDict_,
            1.
        )
    ),

    gammaInt_
    (
        IOobject
        (
            "gamma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	this->mesh_
    )
{    
    if (type == typeName)
    {
        this->correctNut();
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool gammaSST<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
	Flength_.readIfPresent(this->coeffDict());
	ca2_.readIfPresent(this->coeffDict());
	ce2_.readIfPresent(this->coeffDict());
	sigmaGamma_.readIfPresent(this->coeffDict());
        CPG1_.readIfPresent(this->coeffDict());
        CPG1lim_.readIfPresent(this->coeffDict());
        CPG2_.readIfPresent(this->coeffDict());
        CPG3_.readIfPresent(this->coeffDict());
        CPG2lim_.readIfPresent(this->coeffDict());
        CTU1_.readIfPresent(this->coeffDict());
        CTU2_.readIfPresent(this->coeffDict());
        CTU3_.readIfPresent(this->coeffDict());
        ReThetacLim_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        CSEP_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void gammaSST<BasicTurbulenceModel>::correct()
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
    volScalarField& omega_ = this->omega_;
    volScalarField& k_ = this->k_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(2*magSqr(symm(tgradU())));
    const volScalarField S("S", sqrt(S2));
    const volScalarField W("Omega", sqrt(2*magSqr(skew(tgradU()))));

    volScalarField G(this->GName(), nut*S*W);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    const volScalarField CDkOmega
        ( "CD",
        (2*this->alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
        );
    
    const volScalarField F1("F1", this->F1(CDkOmega));

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega_)
         ==
            alpha*rho*gamma*S*W
          - fvm::Sp(alpha*rho*beta*omega_, omega_)
          - fvm::SuSp
            (
                alpha*rho*(F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
        );

        omegaEqn.ref().relax();

        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

        solve(omegaEqn);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    const volScalarField FonLim(
        "FonLim",
        min( max(sqr(this->y_)*S/this->nu() / (
            2.2*ReThetacLim_) - 1., 0.), 3.)
    );
    const volScalarField PkLim(
        "PkLim",
        5*Ck_ * max(gammaInt()-0.2,0.) * (1-gammaInt()) * FonLim * 
        max(3*CSEP_*this->nu() - this->nut_, 0.*this->nut_) * S * W
    );

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), k_)
     ==
        alpha*rho*(G*gammaInt() + PkLim)
        - fvm::Sp(max(gammaInt(),scalar(0.1)) * alpha*rho*this->betaStar_*omega_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

#if (OPENFOAM_PLUS >= 1712 || OPENFOAM >= 1912)
    this->correctNut(S2);
#else
    this->correctNut(S2, this->F23());
#endif

   // Intermittency equation (2)
    volScalarField Pgamma1 = Flength_ * S * gammaInt_ * Fonset(S);
    volScalarField Pgamma2 = ca2_ * W * gammaInt_ * Fturb();
    tmp<fvScalarMatrix> gammaEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
            + fvm::div(alphaRhoPhi, gammaInt_)
            - fvm::laplacian(alpha*rho*this->DgammaEff(), gammaInt_)
            ==
            alpha*rho*Pgamma1 - fvm::Sp(alpha*rho*Pgamma1, gammaInt_) +
            alpha*rho*Pgamma2 - fvm::Sp(alpha*rho*ce2_*Pgamma2, gammaInt_)
        ); 
    
    gammaEqn.ref().relax();
    solve(gammaEqn);

    bound(gammaInt_,scalar(0));

    if (debug && this->runTime_.outputTime()) {
        S.write();
        W.write();
        F1.write();
        CDkOmega.write();
        const volScalarField Pgamma("Pgamma", Pgamma1*(scalar(1)-gammaInt_));
        Pgamma.write();
        const volScalarField Egamma("Egamma", Pgamma2*(scalar(1)-ce2_*gammaInt_));
        Egamma.write();
        FonLim.write();
        PkLim.write();
        Fonset(S)().write();
        FPG()().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
