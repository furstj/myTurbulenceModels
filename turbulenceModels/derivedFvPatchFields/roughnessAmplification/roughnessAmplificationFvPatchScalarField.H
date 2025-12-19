/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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
    Foam::roughnessAmplificationFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a roughnessAmplification condition,
    calculated as:

        \f[
            Ar = cr1*ksPlus
        \f]
    The boundary condition is compatible with Kuzulovic2022 and Langel2017a model.

Usage
    \table
        Property   | Description                     | Req'd | Default
        ks         | equivalent sand grain roughness | yes   | 0.0
        cr1        | multiplicative factor           | no    | 8.0          
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       roughnessAmplification;
        cr1        8.0;
        ks         uniform 10e-6;
        value      uniform 0; // optional initial value
    }
    \endverbatim

SourceFiles
    roughnessAmplificationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef roughnessAmplificationFvPatchScalarField_H
#define roughnessAmplificationFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class roughnessAmplificationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class roughnessAmplificationFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data

        //- Multiplicative factor for BC
        scalar cr1_;

        //- Field of equivalent sand grain roughness
        scalarField ks_;


public:

    //- Runtime type information
    TypeName("roughnessAmplification");


    // Constructors

        //- Construct from patch and internal field
        roughnessAmplificationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        roughnessAmplificationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        roughnessAmplificationFvPatchScalarField
        (
            const roughnessAmplificationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        roughnessAmplificationFvPatchScalarField
        (
            const roughnessAmplificationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new roughnessAmplificationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        roughnessAmplificationFvPatchScalarField
        (
            const roughnessAmplificationFvPatchScalarField&,
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
                new roughnessAmplificationFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
