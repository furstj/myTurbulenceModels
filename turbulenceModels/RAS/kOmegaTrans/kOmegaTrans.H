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
    Foam::RASModels::kOmega

Group
    grpRASTurbulence

Description
    An algebraic transition model for incompressible and compressible flows.

    References:    
    \verbatim
        Kubacki, S. & Dick, E. (2016)
	An algebraic model for prediction of bypass and separation-induced transition 
        in turbomachinery boundary layer flows
	ETMM11 conference proceedings, Palermo
    \endverbatim

    \verbatim
        Wilcox, D. C. (2008).
        Formulation of the k-omega Turbulence Model Revisited
        AIAA Journal, Vol. 46, No. 11, pp. 2823-2838. 
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaTransCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            gamma       0.52;  // Denoted by alpha in the paper
            beta0       0.0708; 
            alphak      0.6;   // Denoted by sigmaK in the paper
            alphaOmega  0.5;   // Denoted by sigmaOmega in the paper
            sigmaD      0.125; // Denoted by sigmaD0 in the paper
	    Clim        0.875;
            al          0.45;
            Agamma     12.0;
            CS         21.0;
            CA          1.0;
            Cpsi       10.0;
            CW          6.0;
            Csep        2.0;
            Av        550.0;
        }
    \endverbatim

SourceFiles
    kOmegaTrans.C

\*---------------------------------------------------------------------------*/

    #ifndef kOmegaTrans_H
    #define kOmegaTrans_H

    #include "RASModel.H"
    #include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
    Class kOmegaTrans Declaration
    \*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaTrans
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel> >
{
    
protected:
    
    // Protected data
    
    // Model coefficients
    
    dimensionedScalar Cmu_;
    dimensionedScalar beta0_;
    dimensionedScalar Clim_;
    dimensionedScalar gamma_;
    dimensionedScalar alphaK_;
    dimensionedScalar alphaOmega_;
    dimensionedScalar sigmaD_;
    dimensionedScalar al_;
    dimensionedScalar Agamma_;
    dimensionedScalar CS_;
    dimensionedScalar CA_;
    dimensionedScalar Cpsi_;
    dimensionedScalar CW_;
    dimensionedScalar Csep_;
    dimensionedScalar AV_;

    // Fields
    const volScalarField& y_;
    volScalarField k_;
    volScalarField omega_;
    
    
    // Protected Member Functions
    
    virtual void correctNut();
    virtual void correctNut(const volScalarField& nus, const volScalarField& nul);

    virtual tmp<volScalarField> intermittency() const;
    virtual tmp<volScalarField> nus(const volScalarField& S, const volScalarField& fSS) const;
    virtual tmp<volScalarField> nul(const volScalarField& S, const volScalarField& fSS) const;
    virtual tmp<volScalarField> beta(const volTensorField& gradU) const;
    virtual tmp<volScalarField> fSS(const volScalarField& S, const volScalarField& W) const;
    virtual tmp<fvScalarMatrix> kSource() const;
    virtual tmp<fvScalarMatrix> omegaSource() const;

public:
    
    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaTrans");


    // Constructors
    
    //- Construct from components
    kOmegaTrans
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
    virtual ~kOmegaTrans()
        {}
    
    
    // Member Functions
    
    //- Read RASProperties dictionary
    virtual bool read();
    
    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DkEff",
                        alphaK_*this->k_/this->omega_ + this->nu()
                    )
                );
        }

    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DomegaEff",
                        alphaOmega_*this->k_/this->omega_ + this->nu()
                    )
                );
        }

    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> k() const
        {
            return k_;
        }

    //- Return the turbulence specific dissipation rate
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
                        Cmu_*k_*omega_,
                        omega_.boundaryField().types()
                    )
                );
        }

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #ifdef NoRepository
    #   include "kOmegaTrans.C"
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
