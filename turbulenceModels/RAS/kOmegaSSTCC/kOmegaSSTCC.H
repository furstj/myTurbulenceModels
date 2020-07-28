/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::RASModels::kOmegaSST

Description
    Specialisation for RAS of the generic kOmegaSSTBase base class.
    For more information, see Description of kOmegaSSTBase.H

    The model includes modified curvature correction proposed by Menter and Smirnov
    \verbatim
       P. E. Smirnov and F. R. Menter, 
       Sensitization of the SST turbulence model to rotation
       and curvature by applying the Spalart-Shur correction term,
       Journal of Turbomachinery, vol. 131, pp. 041010–1, 2009.
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

See also
    Foam::kOmegaSSTCC

SourceFiles
    kOmegaSSTCC.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTCC_H
#define kOmegaSSTCC_H

#include "kOmegaSST.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSSTCC
:
    public kOmegaSST<BasicMomentumTransportModel>
{

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTCC");


    // Constructors

        //- Construct from components
        kOmegaSSTCC
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
    virtual ~kOmegaSSTCC()
    {}

    //- Return k production rate
    virtual tmp<volScalarField::Internal> Pk
    (
        const volScalarField::Internal& G
    ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTCC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
