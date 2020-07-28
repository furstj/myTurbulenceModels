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

Class
    Foam::compressible::RASModels::mykOmegaWilcox06

Description
    Implementation of the updated k-omega model by Wilcox.


\*---------------------------------------------------------------------------*/

#ifndef kOmegaWilcox06_H
#define kOmegaWilcox06_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaWilcox06 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaWilcox06
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel> >
{
    // Private data
    
    // Model coefficients
    
    dimensionedScalar alphaOmega_;
    dimensionedScalar betaStar_;
    dimensionedScalar beta_;
    dimensionedScalar sigmaK_;
    dimensionedScalar sigmaOmega_;
    dimensionedScalar sigmaD_;
    
    dimensionedScalar Clim_; 
    
    // Fields
    
    volScalarField k_;
    volScalarField omega_;

    const volScalarField& y_;
    
    // Private member functions
    
    tmp<volScalarField> omegaBar(const volSymmTensorField& Sbar) const;

    // Disallow copy constructor and assignment
    kOmegaWilcox06(const kOmegaWilcox06&);
    kOmegaWilcox06& operator=(const kOmegaWilcox06&);
    
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;

    //- Runtime type information
    TypeName("kOmegaWilcox06");
    
    
    // Constructors

    //- Construct from components
    kOmegaWilcox06
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
    virtual ~kOmegaWilcox06()
        {}
    

    // Member Functions
    
    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", sigmaK_*this->nut() + this->nu())
                );
        }

    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DomegaEff", sigmaOmega_*this->nut() + this->nu())
                );
        }


    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
        }
    
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    //- Return the turbulence kinetic energy dissipation rate
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
    //virtual void validate();
    
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

    //- Update the turbulent viscosity
    virtual void correctNut();
    virtual void correctNut(const volSymmTensorField& Sbar);
    
    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "kOmegaWilcox06.C"
#endif


#endif

// ************************************************************************* //
