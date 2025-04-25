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
    Foam::roughnessAmplificationLangelFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a roughnessAmplificationLangel condition,
    calculated as:

        \f[
            Ar = CAr1 / (1 + exp(-(CAr2*ksPlus - CAr3)) )
        \f]
    The boundary condition is compatible with Langel2017 model.

    Note that the formula (6.67) in Langel's paper probably contains a typo. 
    There should be minus sign in front of CAr3, see fig. 6.1 in the paper.

Usage
    \table
        Property   | Description                     | Req'd | Default
        ks         | equivalent sand grain roughness | yes   | 0.0
        CAr1       | correlation coeff               | no    | 2000.0
        CAr2       | correlation coeff               | no    | 1.0
        CAr3       | correlation coeff               | no    | 13.5
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       roughnessAmplificationLangel;
        CAr        8.0;
        ks         uniform 10e-6;
        value      uniform 0; // optional initial value
    }
    \endverbatim

SourceFiles
    roughnessAmplificationLangelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef roughnessAmplificationLangelFvPatchScalarField_H
#define roughnessAmplificationLangelFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class roughnessAmplificationLangelFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class roughnessAmplificationLangelFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

    // Private Data

        //- Correlation coeffs
        scalar CAr1_;
        scalar CAr2_;
        scalar CAr3_;

        //- Field of equivalent sand grain roughness
        scalarField ks_;


public:

    //- Runtime type information
    TypeName("roughnessAmplificationLangel");


    // Constructors

        //- Construct from patch and internal field
        roughnessAmplificationLangelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        roughnessAmplificationLangelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        roughnessAmplificationLangelFvPatchScalarField
        (
            const roughnessAmplificationLangelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        roughnessAmplificationLangelFvPatchScalarField
        (
            const roughnessAmplificationLangelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new roughnessAmplificationLangelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        roughnessAmplificationLangelFvPatchScalarField
        (
            const roughnessAmplificationLangelFvPatchScalarField&,
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
                new roughnessAmplificationLangelFvPatchScalarField
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
