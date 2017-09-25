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

#include "CompressibleTurbulenceModel.H"
#include "compressibleTransportModel.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
  
  typedef TurbulenceModel<
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    fluidThermo
    > fluidThermocompressibleTurbulenceModel;
  
  typedef ThermalDiffusivity< CompressibleTurbulenceModel<fluidThermo> >
  fluidThermoCompressibleTurbulenceModel;

  typedef RASModel< EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
  RASfluidThermoCompressibleTurbulenceModel;

}

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoCompressibleTurbulenceModel, LES, Type)


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "gammaSST.H"
makeRASModel(gammaSST);

#include "kOmegaTNT.H"
makeRASModel(kOmegaTNT);

#include "kOmegaWilcox06.H"
makeRASModel(kOmegaWilcox06);

#include "EARSM.H"
makeRASModel(EARSM);

#include "kOmegaTrans.H"
makeRASModel(kOmegaTrans);

#define HAVE_ALPHAT

#include "kv2Omega.H"
makeRASModel(kv2Omega);

#include "mykkLOmega.H"
makeRASModel(mykkLOmega);

#include "mykkLOmegaPh.H"
makeRASModel(mykkLOmegaPh);

#include "mykkLOmegaFS.H"
makeRASModel(mykkLOmegaFS);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //


// ************************************************************************* //
