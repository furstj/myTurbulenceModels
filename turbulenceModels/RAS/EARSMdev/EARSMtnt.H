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
    Foam::RASModels::EARSMtnt

Group
    grpRASTurbulence

Description
   Implementation of the EARSM turbulence model with k-omega-TNT scale determination.

    This turbulence model combines the EARSM nonlinear stress formulation
    with the k-omega-TNT scale-determining equation (Kok 1999).

    The model is described in:
    \verbatim
        Kok, J.C. (1999).
        "Resolving the dependence on free-stream values for the k-omega turbulence model."
        NLR-TP-99295 report.
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

SourceFiles
    EARSMtnt.C

\*---------------------------------------------------------------------------*/

#ifndef EARSMtnt_H
#define EARSMtnt_H

#include "EARSMBase.H"
#include "kOmegaTNT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


/*---------------------------------------------------------------------------*\
                       Class EARSMtnt Declaration
\*---------------------------------------------------------------------------*/


template<class BasicTurbulenceModel>
class EARSMtnt
:
    public EARSMBase<kOmegaTNT<BasicTurbulenceModel>>
{
    // Disallow default bitwise copy construct and assignment
    EARSMtnt(const EARSMtnt&);
    EARSMtnt& operator=(const EARSMtnt&);

public:

    using scaleModel = kOmegaTNT<BasicTurbulenceModel>;
    using alphaField = typename scaleModel::alphaField;
    using rhoField = typename scaleModel::rhoField;
    using transportModel = typename scaleModel::transportModel;


    //- Runtime type information
    TypeName("EARSMtnt");


    // Constructors

        //- Construct from components
        EARSMtnt
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
    virtual ~EARSMtnt()
    {}


    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EARSMtnt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


