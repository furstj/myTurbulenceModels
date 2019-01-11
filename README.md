# myTurbulenceModels

The repository contains several additional models for OpenFOAM.

## Quick Start
### Compilation

There is an **Allwmake** script which should compile a dynamic library. 

Please be sure to use correct branch:

* **OF30** compatible with OpenFOAM 3.0
* **master** development branch oriented to latest version of OpenFOAM

### Running
Put

>> libs ("libmyCompressibleTurbulenceModels.so");

or 

>> libs ("libmyIncompressibleTurbulenceModels.so");

to your controlDict.

Then you can select additional turbulence models via turbulenceProperties, namely:
* **mykkLOmega** - the k-kl-omega model of Walters & Cokljat 
* **mykkLOmegaPh** - modification of the k-kl-omega for APG using Pohlhausen profiles
* **mykkLOmegaFS** - modification of the k-kl-omega for APG using Falkner-Skan profiles
* **mykkLOmegaDev** - experimental development version
* **kv2Omega** - k-v2-omega model of Walters & Lopez
* **gammaSST** - gamma-SST model of Menter, Smirnov, Liu, Avancha
* **EARSM** - Explicit algebraic Reynolds stress model with optional curvature correction

There are some test cases in **testCases** directory.

## Documentation

For current doccumentation see please 
[wiki pages](https://github.com/furstj/myTurbulenceModels/wiki/User-guide).

