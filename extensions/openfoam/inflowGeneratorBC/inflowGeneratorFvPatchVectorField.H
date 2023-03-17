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

#ifndef inflowGeneratorFvPatchVectorField_H
#define inflowGeneratorFvPatchVectorField_H

#include "fixedValueFvPatchField.H"

#include "boost/shared_ptr.hpp"
#include "boost/ptr_container/ptr_set.hpp"
#include <set>
#include "Random.H"

#include "SLList.H"
#include "wallDist.H"
#include "cloud.H"

#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "treeDataFace.H"

#include "structurebasedinflowgenerator.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{




/*---------------------------------------------------------------------------*\
                     Class inflowGeneratorFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/




template<class TurbulentStructure>
class inflowGeneratorFvPatchVectorField
: public fixedValueFvPatchField<vector>,
  //public StructureBasedInflowGenerator<TurbulentStructure,face,SubList,const pointField&>,
  public cloud
{

public:
    typedef
        StructureBasedInflowGenerator<TurbulentStructure,face,SubList,const pointField&>
        structureBasedInflowGenerator;

    
protected:

    autoPtr<structureBasedInflowGenerator> inflowGenerator_;
            
    /**
     * output of lagrangian cloud leads to errors during reconstructPar
     * => make it optional
     */
    Switch lagrangianOutput_;

    label curTimeIndex_;

    Switch scaleToMassflow_;

//    VortonList filterVortons
//    (
//        const inflowGeneratorFvPatchVectorField<TurbulentStructure>& ptf,
//        const fvPatchFieldMapper& mapper,
//        const VortonList& vlist
//    ) const;
    

public:

    //- Runtime type information
    TypeName("inflowGenerator");


    // Constructors

        //- Construct from patch and internal field
        inflowGeneratorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inflowGeneratorFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inflowGeneratorFvPatchVectorField
        //  onto a new patch
        inflowGeneratorFvPatchVectorField
        (
            const inflowGeneratorFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        inflowGeneratorFvPatchVectorField
        (
            const inflowGeneratorFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new inflowGeneratorFvPatchVectorField(*this)
            );
        }


        //- Construct as copy setting internal field reference
        inflowGeneratorFvPatchVectorField
        (
            const inflowGeneratorFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new inflowGeneratorFvPatchVectorField(*this, iF)
            );
        }


    // Member functions


    // Mapping functions

	//- Map (and resize as needed) from self given a mapping object
	virtual void autoMap
	(
	    const fvPatchFieldMapper&
	);

	virtual void autoMap(const mapPolyMesh&) {} // dummy for cloud
	
	//- Reverse map the given fvPatchField onto this fvPatchField
	virtual void rmap
	(
	    const fvPatchField<vector>&,
	    const labelList&
	);

    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

    inline const objectRegistry & db() const { return fixedValueFvPatchVectorField::db(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "inflowGeneratorFvPatchVectorField.C"

#endif

// ************************************************************************* //
