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
    Foam::RASModels::EARSMlm2009

Group
    grpRASTurbulence

Description
   Implementation of the EARSM turbulence model (Hellsten 2005).

    This turbulence model combines the EARSM nonlinear stress formulation
    with the k-omega SST scale-determining equation.

    The model is described in:
    \verbatim
        Hellsten, Antti K. (2005).
        "New Advanced K-W Turbulence Model for High-Lift Aerodynamics."
        AIAA Journal 43 (9): 1857–69. doi:10.2514/1.13754.
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    EARSMlm2009.C

\*---------------------------------------------------------------------------*/

#ifndef EARSMlm2009_H
#define EARSMlm2009_H

#include "EARSMBase.H"
#include "kOmegaSSTLM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class EARSMlm2009 Declaration
\*---------------------------------------------------------------------------*/


template<class BasicTurbulenceModel>
class EARSMlm2009
:
    public EARSMBase<kOmegaSSTLM<BasicTurbulenceModel>> 
{
    // Disallow default bitwise copy construct and assignment
    EARSMlm2009(const EARSMlm2009&);
    EARSMlm2009& operator=(const EARSMlm2009&);

protected:

   // Protected data

    //- Minimum value of k to avoid division by zero
    dimensionedScalar kInf_;


    // Member Functions

    //- Return the first blending function (EARSM fMix formula)
    virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

public:

    using scaleModel = kOmegaSSTLM<BasicTurbulenceModel>;
    using alphaField = typename scaleModel::alphaField;
    using rhoField = typename scaleModel::rhoField;
    using transportModel = typename scaleModel::transportModel;


    //- Runtime type information
    TypeName("EARSMlm2009");


    // Constructors

        //- Construct from components
        EARSMlm2009
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
    virtual ~EARSMlm2009() = default;
 
    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EARSMlm2009.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
