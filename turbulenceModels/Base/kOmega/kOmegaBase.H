/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::kOmegaBase

Group
    grpTurbulence

Description
    Implementation of the standard k-omega turbulence model for
    incompressible and compressible flows.

    \verbatim
        kOmegaCoeffs
        {
            alphaK      0.5;
            alphaOmega  0.5;
            beta        0.072;
            betaStar    0.09;
            gamma       0.52;
        }
    \endverbatim

SourceFiles
    kOmegaBase.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaBase_H
#define kOmegaBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kOmega Declaration
\*---------------------------------------------------------------------------*/

template<class MomentumTransportModel, class BasicMomentumTransportModel>
class kOmega
:
    public MomentumTransportModel
{
    // Private Member Functions

    // Disallow default bitwise copy construct and assignment
    //kOmega(const kOmega&);
    //void operator=(const kOmega&);
    

protected:

    // Protected data
    
    // Model coefficients
    
    dimensionedScalar alphaK_;
    dimensionedScalar alphaOmega_;
    dimensionedScalar gamma_;
    dimensionedScalar beta_;
    dimensionedScalar betaStar_;

    // Fields
    
    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;
    
    volScalarField k_;
    volScalarField omega_;
    
    
    // Protected Member Functions
    virtual void correctNut();

    //- Return k production rate
    virtual tmp<volScalarField::Internal> Pk
    (
        const volScalarField::Internal& G
    ) const;
    
    //- Return epsilon/k which for standard RAS is betaStar*omega
    virtual tmp<volScalarField::Internal> epsilonByk() const;

    virtual tmp<fvScalarMatrix> kSource() const;
    
    virtual tmp<fvScalarMatrix> omegaSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        kOmega
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    //- Destructor
    virtual ~kOmega()
    {}

    
    // Member Functions
    
    //- Re-read model coefficients if they have changed
    virtual bool read();
    
    //- Return the effective diffusivity for k
    virtual tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", this->alphaK_*this->nut_ + this->nu())
                );
        }
    
    //- Return the effective diffusivity for omega
    virtual tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DomegaEff",
                        this->alphaOmega_*this->nut_ + this->nu()
                    )
                );
        }
    
    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
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
    
    //- Return the turbulence kinetic energy dissipation rate
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
