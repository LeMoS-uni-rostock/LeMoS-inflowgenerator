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

#ifndef inflowGeneratorBaseFvPatchVectorField_H
#define inflowGeneratorBaseFvPatchVectorField_H

#include "boostRandomGen.H"
#include "fixedValueFvPatchFields.H"

//#include "Time.H"
#include "SLList.H"
#include "Switch.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"

#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "treeDataFace.H"

#include "PatchTools.H"

#include "fielddataprovider.h"

#include "globalPatch.H"
#include "recursiveapply.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{




/**
 * represents the spatial input data for turbulence generation
 * Is a pointField and holds spatial base for the input fields
 */
class inflowInputDataField
{
    pointField x_;
    const Time& time_;
    
    /**
     * Prescribed mean velocity field on which the turbulent fluctuations will be superimposed.
     */
    autoPtr<FieldDataProvider<vector> >  Umean_;
    mutable label curUmeanAccessTimeIndex_ = -1;
    mutable tmp<Field<vector> > Umeancur_;
    
    /**
     * Prescribed reynolds stress tensor
     */
    autoPtr<FieldDataProvider<symmTensor> > R_;
    mutable label curRAccessTimeIndex_ = -1;
    mutable tmp<Field<symmTensor> > Rcur_;

    /**
     * Prescribed length scale field, 
     * one length scale per principal direction of the reynolds stress tensor.
     * The first component belongs to the direction of the largest reynolds stress component.
     */
    autoPtr<FieldDataProvider<vector> > L_;
    mutable label curLAccessTimeIndex_ = -1;
    mutable tmp<Field<vector> > Lcur_;

    /**
     * Prescribed excess parameter
     */
    autoPtr<FieldDataProvider<scalar> >  c_;
    mutable label curcAccessTimeIndex_ = -1;
    mutable tmp<Field<scalar> > ccur_;
    
    
public:
    inflowInputDataField(const pointField& x, const Time& time);
    inflowInputDataField(const pointField& x, const Time& time, const dictionary& dict);
    inflowInputDataField(const inflowInputDataField& o);
    
    inline label num() const { return x_.size(); }
    
    const Field<point>& x() const;
    const Field<vector>& Umean() const;
    const Field<symmTensor>& R() const;
    const Field<vector>& L() const;
    const scalarField& c() const;
    
    tmp<Field<vector> > flowDirection() const;
    vector flowDirection(label I) const;

    void write(Ostream& os) const;

};



class inflowGeneratorBaseFvPatchVectorField
:
    public fixedValueFvPatchField<vector>,
    public inflowInputDataField
{
  
public:
  struct ProcessStepInfo
  {
    /**
     * # structures currently in simulation
     */
    int n_total;

    /**
     * # structures generated in current step
     */
    int n_generated;

    /**
     * # structures deleted in current step
     */
    int n_removed;

//     /**
//      * # structures that induced fluctuations in patch
//      */
//     int n_induced;
  };
  
private:

    mutable label curflucAccessTimeIndex_ = -1;
    mutable tmp<Field<vector> > curfluctuations_;

    mutable label cursflucAccessTimeIndex_ = -1;
    mutable tmp<Field<scalar> > cursfluctuations_;

protected:
    BoostRandomGen ranGen_;
    
    
    label curTimeIndex_;
    Switch scaleToMassflow_;

    virtual void writeStateVisualization
    (  
      int i,
      const vectorField& u,
      const vectorField* uMean=NULL,
      const symmTensorField* uPrime2Mean=NULL,
      const scalarField* mult = NULL
    ) const =0;
        
    vector randomTangentialDeflection(label fi);

    typedef autoPtr< indexedOctree<treeDataFace> > FaceSearcher;
    mutable FaceSearcher boundaryTree_;
    
    virtual tmp<scalarField> calculateScalarFluctuations() =0;
    
public:
    virtual void computeConditioningFactor(int writeInterval=1000, int nsteps=10000);

    //- Runtime type information
    TypeName("inflowGeneratorBase");


    // Constructors

        //- Construct from patch and internal field
        inflowGeneratorBaseFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inflowGeneratorBaseFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inflowGeneratorBaseFvPatchVectorField
        //  onto a new patch
        inflowGeneratorBaseFvPatchVectorField
        (
            const inflowGeneratorBaseFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        inflowGeneratorBaseFvPatchVectorField
        (
            const inflowGeneratorBaseFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        inflowGeneratorBaseFvPatchVectorField
        (
            const inflowGeneratorBaseFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );


	//- Map (and resize as needed) from self given a mapping object
	virtual void autoMap
	(
	    const fvPatchFieldMapper&
	);

	//- Reverse map the given fvPatchField onto this fvPatchField
	virtual void rmap
	(
	    const fvPatchField<vector>&,
	    const labelList&
	);
	
	// Member functions
	const Field<vector>& fluctuations() const;
	const scalarField& normalizedScalarFluctuations() const;

	
	/**
	 * return the max edge lengths for all faces
	 */
	tmp<scalarField> edgeLengths(bool maxL=false) const;

    
	inline vector averageMeanVelocity() const
	{ 
        return 
            gSum(Umean()*patch().magSf())
            /
            gSum(patch().magSf())
            ; 
    }

    // Evaluation functions

    virtual tmp<vectorField> continueFluctuationProcess(scalar t, ProcessStepInfo *info=NULL, double multidensity=1.0) =0;
    
	virtual void restartFluctuationProcess() =0;
	
	virtual void updateCoeffs();
	
    //- Write
    virtual void write(Ostream&) const;
};

}

#endif
