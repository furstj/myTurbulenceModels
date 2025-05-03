/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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
    Foam::omegaWilcoxRoughWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition implements Wilcox's rough wall condition for
    omega.

    If the surface roughness $k_s^+ = 0$ then 
    \f[
         \omega_w = 10\frac{6\nu}{\beta y_1^2} 
    \f]
    otherwise
    \f[
         \omega_w = \frac{u_{\tau}^2 S_r}{\nu}
    \f]

    where $S_r$ is ether $S_r = (50/k_s^+)^2$ if $k_s^+ < 25$ or $S_r = 100/k+s^+$

    This boundary condition does not modify source term in the k-omega model and can
    be therefore used only in low Reynolds mode!

Usage

    \table
        Property     | Description             | Required    | Default value
        beta1        | Model coefficient       | no          | 0.075
        ks           | equivalent roughness    | no          | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            omegaRoughWallFunction;
    }
    \endverbatim


See also
    Foam::fixedInternalValueFvPatchField
    Foam::epsilonWallFunctionFvPatchScalarField


SourceFiles
    omegaWilcoxRoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaWilcoxRoughWallFunctionFvPatchScalarField_H
#define omegaWilcoxRoughWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class omegaWilcoxRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaWilcoxRoughWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- beta1 coefficient
        scalar beta1_;

        //- Distributed equivalent sand grain roughness height
        scalarField ks_;

public:

    //- Runtime type information
    TypeName("omegaWilcoxRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaWilcoxRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaWilcoxRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        omegaWilcoxRoughWallFunctionFvPatchScalarField
        (
            const omegaWilcoxRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        omegaWilcoxRoughWallFunctionFvPatchScalarField
        (
            const omegaWilcoxRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaWilcoxRoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaWilcoxRoughWallFunctionFvPatchScalarField
        (
            const omegaWilcoxRoughWallFunctionFvPatchScalarField&,
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
                new omegaWilcoxRoughWallFunctionFvPatchScalarField
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
