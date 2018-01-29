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

#include "kOmegaWilcox06.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
kOmegaWilcox06<BasicTurbulenceModel>::kOmegaWilcox06
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
  
  alphaOmega_
  (
   dimensioned<scalar>::lookupOrAddToDict
   (
    "alphaOmega",
    this->coeffDict_,
    13./25.
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
  beta_
  (
   dimensioned<scalar>::lookupOrAddToDict
   (
    "beta",
    this->coeffDict_,
    0.0708
    )
   ),
  sigmaK_
  (
   dimensioned<scalar>::lookupOrAddToDict
   (
    "sigmaK",
    this->coeffDict_,
    3.0/5.0
    )
   ),
  sigmaOmega_
  (
   dimensioned<scalar>::lookupOrAddToDict
   (
    "sigmaOmega",
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
      1./8.
      )
     ),
  Clim_
  (
   dimensioned<scalar>::lookupOrAddToDict
   (
    "Clim",
    this->coeffDict_,
    7.0/8.0
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

    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();

    this->printCoeffs(type);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicTurbulenceModel>
bool kOmegaWilcox06<BasicTurbulenceModel>::read()
{
  if (  eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
      alphaOmega_.readIfPresent(this->coeffDict());
      betaStar_.readIfPresent(this->coeffDict());
      beta_.readIfPresent(this->coeffDict());
      sigmaK_.readIfPresent(this->coeffDict());
      sigmaOmega_.readIfPresent(this->coeffDict());
      sigmaD_.readIfPresent(this->coeffDict());
      Clim_.readIfPresent(this->coeffDict());
      
      return true;
    }
  else
    {
      return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaWilcox06<BasicTurbulenceModel>::omegaBar(const volSymmTensorField& Sbar) const
{
    return  max(
        omega_,
        Clim_*sqrt(2*magSqr(Sbar)/betaStar_)
    );
}

template<class BasicTurbulenceModel>
void kOmegaWilcox06<BasicTurbulenceModel>::correct()
{
  
  if (!this->turbulence_)
    {
      return;
    }
  
  volTensorField gradU(fvc::grad(this->U_));
  volSymmTensorField Sbar(dev(symm(gradU)));
    
  eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();
  
  volScalarField GbyNu = 2*(Sbar && gradU);
  volScalarField G(this->GName(), this->nut_*GbyNu);
  
  tmp<volTensorField> Omega(skew(gradU));
  tmp<volSymmTensorField> Shat(symm(gradU) - 0.5*tr(gradU)*I);
  
  tmp<volScalarField> Xomega( 
			     mag( (Omega() & Omega()) && Shat() ) /
			     pow(betaStar_*omega_,3) 
			      );
  volScalarField fBeta( (1+85*Xomega())/(1+100*Xomega()) );
  
  Xomega.clear();
  Shat.clear();
  Omega.clear();
  
  volScalarField CDkOmega = max(
				sigmaD_/omega_*(fvc::grad(k_) & fvc::grad(omega_)),
				dimensionedScalar("0", inv(sqr(dimTime)), 0.0)
				);
  
  // Turbulent frequency equation 
  // source term modified according to NLR-TP-2001-238
  tmp<fvScalarMatrix> omegaEqn
    (
     fvm::ddt(omega_)
     + fvm::div(this->phi_, omega_)
     - fvm::laplacian(DomegaEff(), omega_)
     ==
     alphaOmega_*GbyNu
     - fvm::Sp(beta_*fBeta*omega_, omega_)
     + CDkOmega
     );
  
    omegaEqn.ref().relax();
    
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);
       
    // Turbulent kinetic energy equation
    
    tmp<fvScalarMatrix> kEqn
      (
       fvm::ddt(k_)
       + fvm::div(this->phi_, k_)
       - fvm::laplacian(DkEff(), k_)
       ==
       G 
       - fvm::Sp(betaStar_*omega_, k_)
       );
    
    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);
    
    this->correctNut(Sbar);
}

template<class BasicTurbulenceModel>
void kOmegaWilcox06<BasicTurbulenceModel>::correctNut(const volSymmTensorField& Sbar)
{
  // Re-calculate viscosity
  this->nut_ = k_/omegaBar(Sbar);
  this->nut_.correctBoundaryConditions();

  BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void kOmegaWilcox06<BasicTurbulenceModel>::correctNut()
{
  volSymmTensorField Sbar(dev(symm(fvc::grad(this->U_))));
  this->correctNut(Sbar);
}

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
