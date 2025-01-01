/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::LESModels::gammaSSTDES

Group
    grpDESTurbulence

Description
    k-omega-SST DES transition and turbulence model for incompressible and compressible flows.

    Reference:
    \verbatim
        Strelets, M. (2001).
        Detached Eddy Simulation of Massively Separated Flows.
        39th AIAA Aerospace Sciences Meeting and Exhibit, Reno, NV.

        Menter, Florian R., Pavel E. Smirnov, Tao Liu a Ravikanth Avancha. 
        A one-equation local correlation-based transition model. 
        Flow, Turbulence and Combustion. 2015, 95(4), 583–619. ISSN 1386-6184. 
        doi:10.1007/s10494-015-9622-4

    \endverbatim

Note
    The default values of the DES constants implemented are code-specific
    values calibrated for OpenFOAM using decaying isotropic turbulence, and
    hence deviate slightly from the values suggested in the reference
    publication.

SourceFiles
    gammaSSTDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gammaSSTDES_H
#define Foam_gammaSSTDES_H

#include "DESModel.H"
#include "kOmegaSSTBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class gammaSSTDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class gammaSSTDES
:
    public kOmegaSSTBase<DESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        gammaSSTDES(const gammaSSTDES&) = delete;

        //- No copy assignment
        void operator=(const gammaSSTDES&) = delete;


protected:

    // Protected Data

        //- Switch to activate grey-area enhanced sigma-(D)DES
        Switch useSigma_;

        // Model coefficients

        dimensionedScalar kappa_;

        //- DES coefficients
        dimensionedScalar CDESkom_;
        dimensionedScalar CDESkeps_;

        // Transiton model coefficients
        dimensionedScalar Flength_;
    
        dimensionedScalar ca2_;
        dimensionedScalar ce2_;
        dimensionedScalar sigmaGamma_;
        dimensionedScalar CPG1_;
        dimensionedScalar CPG1lim_;
        dimensionedScalar CPG2_;
        dimensionedScalar CPG3_;
        dimensionedScalar CPG2lim_;
        dimensionedScalar CTU1_;
        dimensionedScalar CTU2_;
        dimensionedScalar CTU3_;
        dimensionedScalar ReThetacLim_;
        dimensionedScalar Ck_;
        dimensionedScalar CSEP_;

        Switch crossFlow_;
        dimensionedScalar CRSF_;

        // Intermittency field
        volScalarField gammaInt_;


    // Protected Member Functions

        //- Blending for CDES parameter
        virtual tmp<volScalarField> CDES(const volScalarField& F1) const
        {
            return this->blend(F1, CDESkom_, CDESkeps_);
        }

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();

        tmp<volScalarField> r
        (
            const volScalarField& nur,
            const volScalarField& magGradU
        ) const;

        //- Return square of strain rate
        virtual tmp<volScalarField> S2
        (
            const volTensorField& gradU
        ) const;

        //- Return length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;

        //- Return epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField::Internal> GbyNu0
        (
            const volTensorField& gradU,
            const volScalarField& S2
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;

        // Model functions of gammaSST model
        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;       
        tmp<volScalarField> ReThetac() const;
        tmp<volScalarField> Fonset(const volScalarField& S) const;
        tmp<volScalarField> Fonset1(const volScalarField& S) const;
        tmp<volScalarField> FonsetCF() const;
        tmp<volScalarField> Fturb() const;
        tmp<volScalarField> TuL() const;
        tmp<volScalarField> FPG() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("gammaSSTDES");


    // Constructors

        //- Construct from components
        gammaSSTDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~gammaSSTDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- RAS length scale
        virtual tmp<volScalarField> lengthScaleRAS() const;

        //- LES length scale
        virtual tmp<volScalarField> lengthScaleLES
        (
            const volScalarField& CDES
        ) const;

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;

        //- Return the effective diffusivity for gamma
        tmp<volScalarField> DgammaEff() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DgammaEff", (this->nut_/sigmaGamma_) + this->nu())
                );
        }

            //- Return the intermittency
        virtual tmp<volScalarField> gammaInt() const
        {
            return this->gammaInt_;
        }
    
        virtual tmp<volScalarField> Rt() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("Rt", this->k_/(this->nu()*this->omega_)) 
                    // LANGTRY's definition
                );
        }
    
    
        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "gammaSSTDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
