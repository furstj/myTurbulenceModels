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

Class
    Foam::RASModels::EARSM

Group
    grpRASTurbulence

Description
   Implementation of the EARSM turbulence model.

    This turbulence model is described in:
    \verbatim
        Hellsten, Antti K. (2005). 
        "New Advanced K-W Turbulence Model for High-Lift Aerodynamics." 
        AIAA Journal 43 (9): 1857–69. doi:10.2514/1.13754.
    \endverbatim

    The default model coefficients are
    \verbatim
        EARSMCoeffs
        {
            betaStar       0.09
	    gamma1         0.518
            gamma2         0.44
	    beta1          0.0747
	    beta2          0.0828
	    alphaK1        1.1
	    alphaK2        1.1
	    alphaOmega1    0.53
	    alphaOmega2    1.0
	    alphaD1        1.0
	    alphaD2        0.4
            kInf           1.e-10
            curvatureCorrection false
            A0             -0.72
        }
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    EARSM.C

\*---------------------------------------------------------------------------*/

#ifndef EARSM_H
#define EARSM_H

#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class EARSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class EARSM
:
    public nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel> >
{
    // Private memmber functions
    
    // Disallow default bitwise copy construct and assignment
    EARSM(const EARSM&);
    EARSM& operator=(const EARSM&);
    
    
protected:
    
    // Protected data
    
    // Model coefficients

    dimensionedScalar betaStar_;

    dimensionedScalar gamma1_;
    dimensionedScalar gamma2_;

    dimensionedScalar beta1_;
    dimensionedScalar beta2_;

    dimensionedScalar alphaK1_;
    dimensionedScalar alphaK2_;
    
    dimensionedScalar alphaOmega1_;
    dimensionedScalar alphaOmega2_;
    
    dimensionedScalar alphaD1_;
    dimensionedScalar alphaD2_;

    dimensionedScalar kInf_;

    Switch curvatureCorrection_;
    
    dimensionedScalar A0_;
    
    // Fields
    
    volScalarField k_;
    volScalarField omega_;

    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;


    // Protected Member Functions

    virtual tmp<volScalarField> fMix(const volScalarField& gradKgradOmegaByOmega) const;
    
    tmp<volScalarField> blend
    (
        const volScalarField& fMix,
        const dimensionedScalar& psi1,
        const dimensionedScalar& psi2
    ) const
    {
        return fMix*(psi1 - psi2) + psi2;
    }

    tmp<volScalarField> alphaK(const volScalarField& fMix) const
    {
        return blend(fMix, alphaK1_, alphaK2_);
    }

    tmp<volScalarField> alphaOmega(const volScalarField& fMix) const
    {
        return blend(fMix, alphaOmega1_, alphaOmega2_);
    }

    tmp<volScalarField> alphaD(const volScalarField& fMix) const
    {
        return blend(fMix, alphaD1_, alphaD2_);
    }

    tmp<volScalarField> beta(const volScalarField& fMix) const
    {
        return blend(fMix, beta1_, beta2_);
    }

    tmp<volScalarField> gamma(const volScalarField& fMix) const
    {
        return blend(fMix, gamma1_, gamma2_);
    }

    virtual void correctNut();

    virtual volScalarField N(
        const volScalarField& A3p,
        const volScalarField& P1,
        const volScalarField& P2
    ) const;
    
    virtual void correctNonlinearStress(const volTensorField& gradU);

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("EARSM");


    // Constructors

        //- Construct from components
        EARSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~EARSM()
    {}


    // Member Functions
    
    //- Re-read model coefficients if they have changed
    virtual bool read();
    
    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff(const volScalarField& fMix) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", alphaK(fMix)*this->nut_ + this->nu())
                );
        }
    
    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff(const volScalarField& fMix) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DomegaEff", alphaOmega(fMix)*this->nut_ + this->nu())
                );
        }
    
    //- Return the total fluctuation kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
        }
    
    
    //- Return the turbulence specific dissipation rate
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    //- Return the total fluctuation kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                    )
                );
            
        }
    
    //- Validate the turbulence fields after construction
    //  Update turbulence viscosity and other derived fields as requires
    virtual void validate();
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "EARSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
