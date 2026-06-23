/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::EARSMBase

Description
    Base class implementation of the EARSM stress model

SourceFiles
    EARSMBase.C

\*---------------------------------------------------------------------------*/

#ifndef EARSMBase_H
#define EARSMBase_H

namespace Foam
{
namespace RASModels
{

template <class ScaleModel>
class EARSMBase 
: 
public ScaleModel
{

    using alphaField = typename ScaleModel::alphaField;
    using rhoField = typename ScaleModel::rhoField;
    using transportModel = typename ScaleModel::transportModel;

protected:

    virtual void correctNut(const volScalarField& S2);

    virtual void correctNut();

    virtual volScalarField N(
        const volScalarField& A3p,
        const volScalarField& P1,
        const volScalarField& P2
    ) const;

    virtual void correctAnisotropy(const volTensorField& gradU);

    virtual bool read(const word& typeName);

    dimensionedScalar CTau_;

    dimensionedScalar A0_;
    
    Switch nearWallDamping_;

    dimensionedScalar B2_;

    dimensionedScalar CyA_;

    dimensionedScalar CyB_;

    volSymmTensorField aEx_;
    
public:

    EARSMBase<ScaleModel>
    (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = word::null
    );
    

    virtual ~EARSMBase<ScaleModel>()
    {}

    //- Return the Reynolds stress tensor
    virtual tmp<volSymmTensorField> R() const;

    //- Return the effective stress tensor
    virtual tmp<volSymmTensorField> devRhoReff() const;

    //- Return the effective stress tensor based on a given velocity field
    virtual tmp<volSymmTensorField> devRhoReff
    (
        const volVectorField& U
    ) const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevRhoReff
    (
        const volScalarField& rho,
        volVectorField& U
    ) const;

    //- Validate the turbulence fields after construction
    //  Update turbulence viscosity and other derived fields as requires
    virtual void validate();

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

};

} // End namespace RASModels
} // End namespace Foam

#include "EARSMBase.C"

#endif

