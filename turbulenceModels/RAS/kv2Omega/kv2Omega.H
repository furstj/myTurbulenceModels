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
    Foam::RASModels:kv2Omega

Group
    grpRASTurbulence

Description
   Implementation of the kv2Omega turbulence model for predicting
    laminar-turbulent transition on arbitrary geometries and grid topo-
    logies.

    This turbulence model is described in:
    \verbatim
        Lopez, M.,  Walters, D. K. (2016). 
	Prediction of transitional and fully turbulent flow using an alternative 
	to the laminar kinetic energy approach. 
	Journal of Turbulence, 17(3), 253–273. 
	doi:10.1080/14685248.2015.1062509
    \endverbatim

    The default model coefficients are
    \verbatim
        kv2OmegaCoeffs
        {
            A0             4.04
            AS             2.12
            Anu            3.8
            ABP            0.6
            ANAT           200
            ATS            200
            CBPcrit        1.5
            CNC            0.1
            CNATcrit       1450
            CINT           0.95
            CTScrit        1000
            CRNAT          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.32
            CalphaTheta    0.035
            CSS            3.0
            Ctau1          4360
            Cw1            0.44
            Cw2            0.92
            CwR            1.15
            Clambda        2.495
            betaStar       0.09
            Prt            0.85
            sigmaK         1
            sigmaW         1.17
            sigmaW2        1.856
        }
    \endverbatim

SourceFiles
    kv2Omega.C

\*---------------------------------------------------------------------------*/

#ifndef kv2Omega_H
#define kv2Omega_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kv2Omega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kv2Omega
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel> >
{
    // Private memmber functions
    
    // Disallow default bitwise copy construct and assignment
    kv2Omega(const kv2Omega&);
    kv2Omega& operator=(const kv2Omega&);
    
    
protected:
    
    // Protected data
    
    // Model coefficients
    
    dimensionedScalar A0_;
    dimensionedScalar AS_;
    dimensionedScalar Anu_;
    dimensionedScalar ABP_;
    dimensionedScalar ANAT_;
    dimensionedScalar ATS_;
    dimensionedScalar CBPcrit_;
    dimensionedScalar CNC_;
    dimensionedScalar CNATcrit_;
    dimensionedScalar CINT_;
    dimensionedScalar CTScrit_;
    dimensionedScalar CRNAT_;
    dimensionedScalar C11_;
    dimensionedScalar C12_;
    dimensionedScalar CR_;
    dimensionedScalar CalphaTheta_;
    dimensionedScalar CSS_;
    dimensionedScalar Ctau1_;
    dimensionedScalar Cw1_;
    dimensionedScalar Cw2_;
    dimensionedScalar CwR_;
    dimensionedScalar Clambda_;
    dimensionedScalar betaStar_;
    dimensionedScalar PrTheta_;
    dimensionedScalar sigmaK_;
    dimensionedScalar sigmaW_;
    dimensionedScalar sigmaW2_;
    
    // Fields
    
    volScalarField k_;
    volScalarField v2_;
    volScalarField omega_;

    //- Wall distance
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;


    // Protected Member Functions

    virtual void correctNut();

    virtual tmp<volScalarField> fv(const volScalarField& Ret) const;

    virtual tmp<volScalarField> fINT() const;

    virtual tmp<volScalarField> fSS(const volScalarField& W) const;
    
    virtual tmp<volScalarField> Cmu(const volScalarField& S) const;
    
    virtual tmp<volScalarField> betaTS(const volScalarField& ReW) const;
    
    virtual tmp<volScalarField> lambdaT() const;

    virtual tmp<volScalarField> lambdaEff(const volScalarField& lambdaT) const;

    virtual tmp<volScalarField> fTaul
    (
        const volScalarField& lambdaEff,
        const volScalarField& v2l,
        const volScalarField& W
    ) const;
    
    virtual tmp<volScalarField> alphaT
    (
        const volScalarField& lambdaEff,
        const volScalarField& fv,
        const volScalarField& v2s
    ) const;
    
  /*
    virtual tmp<volScalarField> fOmega
    (
        const volScalarField& lambdaEff,
        const volScalarField& lambdaT
    ) const;
  */
    
    virtual tmp<volScalarField> phiBP(const volScalarField& W) const;
    
    virtual tmp<volScalarField> phiNAT
    (
        const volScalarField& ReOmega,
        const volScalarField& fNatCrit
    ) const;
    
    virtual tmp<volScalarField> D(const volScalarField& k) const;

    virtual tmp<volScalarField> F1() const;
    
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kv2Omega");


    // Constructors

        //- Construct from components
        kv2Omega
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
    virtual ~kv2Omega()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaT/sigmaK_ + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaT/sigmaW_ + this->nu())
            );
        }

        //- Return the total fluctuation kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the wall normal fluctuation kinetic energy
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
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
                    new volScalarField(
                        "epsilon", 
                        omega_ * min(k_, v2_) + D(k_)
                    ));
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
#   include "kv2Omega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
