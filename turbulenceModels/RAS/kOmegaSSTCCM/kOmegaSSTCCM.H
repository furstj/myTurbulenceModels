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

    The model includes modified curvature correction proposed in
    \verbatim
        Quiang, Z., Yong, Y. (2013)
        A new simpler rotation/curvature correction method for
        Spalart–Allmaras turbulence model
        Chinese Journal of Aeronautics, 26(2)
    \endverbatim

Author
    Jiri Furst <Jiri.Furst@gmail.com>

See also
    Foam::kOmegaSSTCCM

SourceFiles
    kOmegaSSTCCM.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTCCM_H
#define kOmegaSSTCCM_H

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
class kOmegaSSTCCM
:
    public kOmegaSST<BasicMomentumTransportModel>
{

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTCCM");


    // Constructors

        //- Construct from components
        kOmegaSSTCCM
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
    virtual ~kOmegaSSTCCM()
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
    #include "kOmegaSSTCCM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
