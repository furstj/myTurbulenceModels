/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a turbulent kinetic energy condition,
    based on user-supplied turbulence intensity, defined as a fraction of the
    mean velocity. Moreover, the kinetic energy is bounde from up and bellow:

        \f[
            k_p = max( min( 1.5 (I |U|)^2, k_max), k_min)
        \f]

    where

    \vartable
        k_p     | kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
        k_max   | upper bound for k_p
        k_min   | lower bound for k_p
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

Usage
    \table
        Property     | Description             | Required    | Default value
        intensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        k_min        | lower bound for k_p     | no          | 0
        k_max        | upper bound for k_p     | no          | GREAT
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        boundedTurbulentIntensityKineticEnergyInlet;
        intensity   0.05;           // 5% turbulence
        k_min       1.e-5;          // bound k for below (allows phi=0)
        value       uniform 1;      // placeholder
    }
    \endverbatim

See also
    Foam::inletOutletFvPatchField

SourceFiles
    boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField_H
#define boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        scalar intensity_;

        //- Lower bound for turbulent kinetic energy
        scalar k_min_;

        //- Upper bound for turbulent kinetic energy
        scalar k_max_;
    
        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("boundedTurbulentIntensityKineticEnergyInlet");


    // Constructors

        //- Construct from patch and internal field
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        //  onto a new patch
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new boundedTurbulentIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
