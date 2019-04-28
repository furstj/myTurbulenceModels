/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Felix Langfeldt
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::incompressible::ReThetatInletFvPatchScalarField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a inlet boundary condition for the
    ReThetat variable, which is used in the gammaReThetatSST transitional
    turbulence model. The inlet condition is based upon the turbulence
    intensity, defined as a fraction of the mean velocity:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | turbulence kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    The inlet value of ReThetat is a function of the turbulence intensity:

        \f[
            ReThetat_p = ReThetat(I)
        \f]

    where

    \vartable
        ReThetat_p     | ReThetat at the patch
        I                   | turbulence intensity
        ReThetat()          | empirical correlation
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type        ReThetatInlet;
        value       uniform 1;      // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    ReThetatInletFvPatchScalarField.C

Author
    Felix Langfeldt - original implementation
    Jiri Furst - made compatible with the kOmegaSSTLM model

\*---------------------------------------------------------------------------*/

#ifndef ReThetatInletFvPatchScalarField_H
#define ReThetatInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
   Class ReThetatInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ReThetatInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        
        //- Name of the velocity field
        word UName_;

        //- Name of the turbulent kinetic energy field
        word kName_;


    scalar ReThetatInlet(scalar Tu) const;
    
public:

    //- Runtime type information
    TypeName("ReThetatInlet");


    // Constructors

        //- Construct from patch and internal field
        ReThetatInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ReThetatInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ReThetatInletFvPatchScalarField
        //  onto a new patch
        ReThetatInletFvPatchScalarField
        (
            const ReThetatInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ReThetatInletFvPatchScalarField
        (
            const ReThetatInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ReThetatInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        ReThetatInletFvPatchScalarField
        (
            const ReThetatInletFvPatchScalarField&,
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
                new ReThetatInletFvPatchScalarField
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
