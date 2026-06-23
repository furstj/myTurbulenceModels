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
    Foam::RASModels::EARSMko2005

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
    EARSMko2005.C

\*---------------------------------------------------------------------------*/

#ifndef EARSMko2005_H
#define EARSMko2005_H

#include "EARSMBase.H"
#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class EARSMko2005 Declaration
\*---------------------------------------------------------------------------*/


template<class BasicTurbulenceModel>
class EARSMko2005
:
    public EARSMBase<kOmegaSST<BasicTurbulenceModel>> 
{
    // Disallow default bitwise copy construct and assignment
    EARSMko2005(const EARSMko2005&);
    EARSMko2005& operator=(const EARSMko2005&);

public:

    using scaleModel = kOmegaSST<BasicTurbulenceModel>;
    using alphaField = typename scaleModel::alphaField;
    using rhoField = typename scaleModel::rhoField;
    using transportModel = typename scaleModel::transportModel;


    //- Runtime type information
    TypeName("EARSMko2005");


    // Constructors

        //- Construct from components
        EARSMko2005
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
    virtual ~EARSMko2005()
    {}


    // Member Functions

    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EARSMko2005.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
