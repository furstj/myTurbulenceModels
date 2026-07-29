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
    Foam::RASModels::EARSMgammaSST

Group
    grpRASTurbulence

Description
   Implementation of the EARSM turbulence model.

    This turbulence model combines the EARSM nonlinear stress formulation
    with the k-omega gammaSST scale-determining equation.

    The scale determining model is described in:
    \verbatim
        MENTER, Florian R., Pavel E. SMIRNOV, Tao LIU a Ravikanth AVANCHA. 
        A one-equation local correlation-based transition model. 
        Flow, Turbulence and Combustion. 2015, 95(4), 583–619. ISSN 1386-6184. 
        doi:10.1007/s10494-015-9622-4
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    EARSMgammaSST.C

\*---------------------------------------------------------------------------*/

#ifndef EARSMgammaSST_H
#define EARSMgammaSST_H

#include "EARSMBase.H"
#include "gammaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class EARSMgammaSST Declaration
\*---------------------------------------------------------------------------*/


template<class BasicTurbulenceModel>
class EARSMgammaSST
:
    public EARSMBase<gammaSST<BasicTurbulenceModel>> 
{
    // Disallow default bitwise copy construct and assignment
    EARSMgammaSST(const EARSMgammaSST&);
    EARSMgammaSST& operator=(const EARSMgammaSST&);

protected:

    // Protected data


    // Member Functions

    //- Return the first blending function (EARSM fMix formula)
    virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

public:

    using scaleModel = kOmegaSSTLM<BasicTurbulenceModel>;
    using alphaField = typename scaleModel::alphaField;
    using rhoField = typename scaleModel::rhoField;
    using transportModel = typename scaleModel::transportModel;


    //- Runtime type information
    TypeName("EARSMgammaSST");


    // Constructors

        //- Construct from components
        EARSMgammaSST
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
    virtual ~EARSMgammaSST() = default;
 
    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EARSMgammaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
