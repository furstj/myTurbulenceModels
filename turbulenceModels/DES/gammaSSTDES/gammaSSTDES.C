/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "gammaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void gammaSSTDES<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::correctNut(S2);

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void gammaSSTDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::r
(
    const volScalarField& nur,
    const volScalarField& magGradU
) const
{
    const dimensionedScalar eps(magGradU.dimensions(), SMALL);

    tmp<volScalarField> tr =
        min(nur/(max(magGradU, eps)*sqr(this->kappa_*this->y_)), scalar(10));

    tr.ref().boundaryFieldRef() == 0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::S2
(
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tS2 =
        kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::S2(gradU);

    if (this->useSigma_)
    {
        volScalarField& S2 = tS2.ref();

        const volScalarField& k = this->k_;
        const volScalarField& omega = this->omega_;

        const volScalarField CDkOmega
        (
            (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
        );

        const volScalarField F1(this->F1(CDkOmega));

        const volScalarField CDES(this->CDES(F1));
        const volScalarField dTilda(this->dTilda(mag(gradU), CDES));
        const volScalarField lengthScaleRAS(this->lengthScaleRAS());
        const volScalarField Ssigma(this->Ssigma(gradU));

        S2 =
            pos(dTilda - lengthScaleRAS)*S2
          + (scalar(1) - pos(dTilda - lengthScaleRAS))*sqr(Ssigma);
    }

    return tS2;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    return min(lengthScaleLES(CDES), lengthScaleRAS());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gammaSSTDES<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    volScalarField CDES(this->CDES(F1));
    return sqrt(this->k_())/dTilda(mag(gradU), CDES)()();
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gammaSSTDES<BasicTurbulenceModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& S2
) const
{
    if (this->useSigma_)
    {
        return S2();
    }

    return
        kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::GbyNu0(gradU, S2);
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gammaSSTDES<BasicTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return GbyNu0; // Unlimited
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::F1(const volScalarField& CDkOmega) const
{
  return max(
             kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::F1(CDkOmega),
             exp(-sqr(pow4(this->y_*sqrt(this->k_)/(scalar(120)*this->nu()))))
             );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::ReThetac() const
{
  return CTU1_ + CTU2_*exp(-CTU3_*TuL()*FPG() );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::Fonset(const volScalarField& S) const
{
    tmp<volScalarField> Fons(
        max
        (
            min(Fonset1(S), 2.0) - max(1.0-pow3(Rt()/3.5),0.0),
            0.0
        )
    );

    if (this->crossFlow_)
    {
        Fons = max(Fons, this->FonsetCF());
    }

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
                Fons
            )
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::Fonset1(const volScalarField& S) const
{
    return sqr(this->y_)*S/this->nu() / (2.2*ReThetac());
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::FonsetCF() const
{
    tmp<volVectorField> w(fvc::curl(this->U_));
    const dimensionedScalar wMin("VSMALL", inv(dimTime), VSMALL);
    tmp<volVectorField> ew( w() / max(mag(w()), wMin));

    tmp<volVectorField> n(-fvc::grad(this->y_));
    tmp<volScalarField> Psi(mag(n() & fvc::grad(ew))*this->y_);

    tmp<volScalarField> lambda (
        min( 0.0477, max( 0.0, 
        -7.57e-3 * ( fvc::grad(this->U_ & n()) & n()) * sqr(this->y_) / this->nu() + 0.0174))
    );
    
    tmp<volScalarField> gLambda (
        min(2.3, max( 1.0, ((27864.0*lambda()-1962.0)*lambda()+54.3)*lambda() + 1.0))
    );
    lambda.clear();
 
    tmp<volScalarField> Rev(sqr(this->y_) * mag(w) / this->nu());

    tmp<volScalarField> TC1(this->CRSF_/150.8*0.684/gLambda*Psi*Rev);

    return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "FonsetCF",
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                min
                 (
                    max(100.0*(TC1 - 1.0), 0.0),
                    1.0
                )
            )
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::Fturb() const
{
    return exp(-pow4(Rt()/2));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::TuL() const
{
    return min(100 * sqrt(2.0/3.0*this->k_) / (this->omega_ * this->y_), 100.0);
}

   
template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::FPG() const
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
gammaSSTDES<BasicTurbulenceModel>::gammaSSTDES
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
    kOmegaSSTBase<DESModel<BasicTurbulenceModel>>
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

    useSigma_
    (
        Switch::getOrAddToDict
        (
            "useSigma",
            this->coeffDict_,
            false
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CDESkom_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDESkom",
            this->coeffDict_,
            0.82
        )
    ),
    CDESkeps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDESkeps",
            this->coeffDict_,
            0.60
        )
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
    crossFlow_(
        Switch::lookupOrAddToDict
        (
            "crossFlow",
            this->coeffDict_,
            false
        )
    ),
    CRSF_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CRSF",
            this->coeffDict_,
            1.
        )
    ),
    gammaInt_
    (
        IOobject
        (
            "gammaInt",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )


{
    // Note: Ctrans coeff is model specific; for k-w = 60
    this->Ctrans_ =
        dimensioned<scalar>::getOrAddToDict("Ctrans", this->coeffDict_, 60.0);

    if (type == typeName)
    {
        WarningInFunction
            << "This model is not recommended and will be deprecated in future "
            << "releases. Please consider using the DDES/IDDES versions instead"
            << endl;

        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool gammaSSTDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTBase<DESModel<BasicTurbulenceModel>>::read())
    {
        useSigma_.readIfPresent("useSigma", this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        CDESkom_.readIfPresent(this->coeffDict());
        CDESkeps_.readIfPresent(this->coeffDict());

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
        crossFlow_.readIfPresent("crossFlow", this->coeffDict());
        CRSF_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
gammaSSTDES<BasicTurbulenceModel>::lengthScaleRAS() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    return sqrt(k)/(this->betaStar_*omega);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
gammaSSTDES<BasicTurbulenceModel>::lengthScaleLES
(
    const volScalarField& CDES
) const
{
    return CDES*this->delta();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gammaSSTDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volVectorField& U = this->U_;

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    const volScalarField F1(this->F1(CDkOmega));

    return tmp<volScalarField>::New
    (
        IOobject::scopedName("DES", "LESRegion"),
        neg(dTilda(mag(fvc::grad(U)), CDES(F1)) - lengthScaleRAS())
    );
}

template<class BasicTurbulenceModel>
void gammaSSTDES<BasicTurbulenceModel>::correct()
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
    fv::options& fvOptions(fv::options::New(this->mesh_));

    DESModel<BasicTurbulenceModel>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );
    
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    const volScalarField S("S", sqrt(S2));
    const volScalarField W("Omega", sqrt(2*magSqr(skew(tgradU()))));

    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S*W));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);


    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

        //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    this->omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //omega_.correctBoundaryConditions();

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    //const volScalarField F1_(F1(CDkOmega));
    const volScalarField F23(this->F23());

        {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S*W);

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
          + alpha()*rho()*beta*sqr(this->omegaInf_)
          + this->omegaSource()
          + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    {
        const volScalarField FonLim(
            "FonLim",
            min( max(sqr(this->y_)*S/this->nu() / (
                2.2*this->ReThetacLim_) - 1., 0.), 3.)
        );

        const volScalarField PkLim(
            "PkLim",
            5*Ck_ * max(gammaInt()-0.2,0.) * (1-gammaInt()) * FonLim * 
            max(3*this->CSEP_*this->nu() - this->nut_, 0.*this->nut_) * S * W
        );

        const volScalarField::Internal& gammaIn(this->gammaInt_.internalField());

        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, this->k_)
          + fvm::div(alphaRhoPhi, this->k_)
          - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
         ==
            alpha()*rho()*(this->Pk(G)*gammaIn + PkLim)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU*gammaIn, this->k_)
          - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU())*max(gammaIn,scalar(0.1)), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(this->k_);
        bound(this->k_, this->kMin_);
    }


    {
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
    }

    correctNut(S*W);

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
