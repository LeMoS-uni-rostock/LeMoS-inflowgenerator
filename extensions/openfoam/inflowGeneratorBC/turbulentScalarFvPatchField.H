/*
 * This file is part of Insight CAE, a workbench for Computer-Aided Engineering 
 * Copyright (C) 2014  Hannes Kroeger <hannes@kroegeronline.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#ifndef turbulentScalarFvPatchField_H
#define turbulentScalarFvPatchField_H


#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class extendedFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class turbulentScalarFvPatchField
:
    public fixedValueFvPatchScalarField
{

protected:
    word UName_;
    
    autoPtr<FieldDataProvider<scalar> >  mean_;    
    autoPtr<FieldDataProvider<scalar> >  variance_;    
    autoPtr<FieldDataProvider<vector> >  turbulentFlux_;
    
    autoPtr<volVectorField> w123_;
    autoPtr<volScalarField> w4radic_, w4_, flucraw_, flucunlim_, flucmean_, flucc1_, flucc2_;

public:

    //- Runtime type information
    TypeName("turbulentScalar");


    // Constructors

        //- Construct from patch and internal field
        turbulentScalarFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentScalarFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentScalarFvPatchField
        //  onto a new patch
        turbulentScalarFvPatchField
        (
            const turbulentScalarFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentScalarFvPatchField
        (
            const turbulentScalarFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new turbulentScalarFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentScalarFvPatchField
        (
            const turbulentScalarFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new turbulentScalarFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

        
            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

