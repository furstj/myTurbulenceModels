# myTurbulenceModels

The repository contains several additional models for OpenFOAM.

## Quick Start
### Compilation

There is an **Allwmake** script which should compile a dynamic library. 

Please be sure to use correct branch:

* **OF3** compatible with OpenFOAM 3.0
* **default** development branch oriented to latest version of OpenFOAM

### Running
Put

>> libs ("libmyCompressibleTurbulenceModels.so");

or 

>> libs ("libmyIncompressibleTurbulenceModels.so");

to your controlDict.

Then you can select additional turbulence models via turbulenceProperties.

## Documentation

For current doccumentation see please 
[wiki pages](https://github.com/furstj/myTurbulenceModels/wiki/User-guide).

