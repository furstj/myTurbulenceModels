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
    Foam::nutUKnoppWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on velocity,  using blend of Spalding's law
    and Reichardt's law to give a continuous nut profile to the wall (y+ = 0)

        \f[
          F_{Sp,3} = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
                - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
          F_{Rei} = \frac{ln(1 + 0.4 y^+)}{\kappa} + 7.8 \left( 1 -\,
                - exp(-y^+/11.0) - \frac{y^+}{11.0} exp(-y^+/3.0)\right)
          F_{log} = \frac{1}{\kappa} ln(E y^+)
          F_{Rei,m} = (1 - \phi_{b1}) F_{Rei} + \phi_{b1} F_{log}
          \phi_{b1} = tanh((y^+/27)^4)
          F_{k \omega} = (1 - \phi_{k \omega}) F_{Sp,3} +
               \phi_{k \omega} F_{Rei,m}
          \phi_{k \omega} = tanh((y^+/50)^2)
        \f]

    where
    \vartable
        y^+     | non-dimensional position
        u^+     | non-dimensional velocity
        \kappa  | Von Karman constant
    \endvartable

    Model described in:
    \verbatim
        Knopp, T., Alrutz, T. and Schwamborn, D. (2006)
        "A Grid and flow adaptive wall-function method for RANS turbulence
        modelling"
        Journal of Computational Physics,
        Volume 220, 19-40
    \endverbatim

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutUKnoppWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUKnoppWallFunctionFvPatchScalarField.C

Author
    Pekka Pasanen (https://github.com/zordiack/foam-dev)

\*---------------------------------------------------------------------------*/

#ifndef nutUKnoppWallFunctionFvPatchScalarField_H
#define nutUKnoppWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUKnoppWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUKnoppWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;


public:

    //- Runtime type information
    TypeName("nutUKnoppWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUKnoppWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUKnoppWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUKnoppWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUKnoppWallFunctionFvPatchScalarField
        (
            const nutUKnoppWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUKnoppWallFunctionFvPatchScalarField
        (
            const nutUKnoppWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUKnoppWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUKnoppWallFunctionFvPatchScalarField
        (
            const nutUKnoppWallFunctionFvPatchScalarField&,
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
                new nutUKnoppWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
