/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::omegaViscosityRatioFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides an inlet boundary condition for omega based on
    given turbulent-laminar viscosity ratio.

    \f[
        \omega_{bc} = k/(ratio*\nu)
    \f]

    where

    \vartable
        k       | TKE
        \nu     | molecular viscosity
        ratio   | \nu_t/\nu ratio
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        ratio        | viscosity ratio         | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type  omegaViscosityRatio;
        ratio uniform 100;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    omegaViscosityRatioFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaViscosityRatioFvPatchScalarField_H
#define omegaViscosityRatioFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class omegaViscosityRatioFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaViscosityRatioFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
protected:

    // Private data

    //- Viscosity ratio
    scalar ratio_;

public:

    //- Runtime type information
    TypeName("omegaViscosityRatio");

    
    // Constructors

    //- Construct from patch and internal field
    omegaViscosityRatioFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );
    
    //- Construct from patch, internal field and dictionary
    omegaViscosityRatioFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );
    
    //- Construct by mapping given
    // omegaViscosityRatioFvPatchScalarField
    //  onto a new patch
    omegaViscosityRatioFvPatchScalarField
    (
        const omegaViscosityRatioFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );
    
    //- Construct as copy
    omegaViscosityRatioFvPatchScalarField
    (
        const omegaViscosityRatioFvPatchScalarField&
    );
    
    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                    new omegaViscosityRatioFvPatchScalarField(*this)
                );
        }
    
    //- Construct as copy setting internal field reference
    omegaViscosityRatioFvPatchScalarField
    (
        const omegaViscosityRatioFvPatchScalarField&,
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
                    new omegaViscosityRatioFvPatchScalarField(*this, iF)
                );
        }
    
    
    // Member functions
    
    // Access
    
    //- Return reference to the viscosity ratio field
    const scalar& ratio() const
        {
            return ratio_;
        }
    
        
    // Evaluation functions
    
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();


    // I-O
    
    //- Write
    virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
