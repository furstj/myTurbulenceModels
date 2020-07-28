/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IncompressibleMomentumTransportModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "kinematicMomentumTransportModels.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
#define makeRASModel(Type)                                              \
    makeTemplatedMomentumTransportModel                                               \
    (transportModelIncompressibleMomentumTransportModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                               \
    (transportModelIncompressibleMomentumTransportModel, LES, Type)
    */

// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "gammaSST.H"
makeRASModel(gammaSST);

#include "kOmegaSSTCC.H"
makeRASModel(kOmegaSSTCC);

#include "kOmegaSSTCCM.H"
makeRASModel(kOmegaSSTCCM);

#include "kOmegaTNT.H"
makeRASModel(kOmegaTNT);

#include "kOmegaWilcox06.H"
makeRASModel(kOmegaWilcox06);

#include "kv2Omega.H"
makeRASModel(kv2Omega);

#include "mykkLOmega.H"
makeRASModel(mykkLOmega);

#include "mykkLOmegaFS.H"
makeRASModel(mykkLOmegaFS);

#include "mykkLOmegaPh.H"
makeRASModel(mykkLOmegaPh);

#include "kOmegaTrans.H"
makeRASModel(kOmegaTrans);

#include "EARSM.H"
makeRASModel(EARSM);

#include "EARSMWallin.H"
makeRASModel(EARSMWallin);

#include "EARSMTrans.H"
makeRASModel(EARSMTrans);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "XLES.H"
makeLESModel(XLES)

// ************************************************************************* //
