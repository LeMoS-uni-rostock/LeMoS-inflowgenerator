
/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
     Foam::DSRFG1FvPatchVectorField
\*--------------------------------------------------------------------------*/
#ifndef DSRFG1FvPatchVectorField_H
#define DSRFG1FvPatchVectorField_H

#include "Random.H"
#include "IFstream.H"
#include "labelList.H"
#include "edgeList.H"
#include "point.H"
#include "DynamicField.H"
#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fileName.H"
#define PI 3.14159

namespace Foam{
/*---------------------------------------------------------------------------*\
                   Class DSRFG1FvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class DSRFG1FvPatchVectorField
:
   public fixedValueFvPatchField<vector>
{
   //-private Members
	//- Object to Random class, rand_ is used to generate
	//- Randon variables
	Random rand_;

	//- Number of Spectrum Parts
	 label SpectrumParts_;

	//- Number of Modes
	 label NModes_;

	//-file path to Umean
	fileName Ufile;

	//- Mean Velocity
	vectorField U_;
 
	scalar UMean_;

	//- Current time Index (used for Updating)
	label curTimeIndex_;

	//- Wave Number
	scalar k0_; //Either user or dictionary

	//- Turbulence Lenght Scale
	scalar Ls_; 	
	fileName Lfile;
	scalarField Le_;
	
	//-file path to R stresses
	fileName Rfile;
	//Reynolds Stresses
	symmTensorField R_;

	scalar Theta1;
	scalar Theta2;	

	label DataCalIndex_;

	//-Deltax
	scalarField deltax;
  	//-Deltay
  	scalarField deltay;
	//-Deltaz
  	scalarField deltaz;

	scalarField Distancey;
	
	

	//Energy Spectrum

	Field<scalarField> E_;
	

	//Wave Number for each Energy Spectrum Part

	scalarField k_m;

	//random values for each spectrum part and nodes are
	// stored in zeta,zi,P,Q,K,Omega

	struct P{
        	Field<scalarField> mn;
    	} p[3];

    	struct Q{
        	Field<scalarField> mn;
    	} q[3];

    	struct K{
        	Field<scalarField> mn;
    	} k[3];
	
	
	Field<scalarField> Omega;
  

  protected:

  //- Protected Data Menbers

	//- Adamjan Experimental Coefficients
	    static const scalar C_L;
    	static const scalar C_tau;
   	  static const scalar Beta;
    	static const scalar Alpha;
    	static const scalar L_eta;

  public:

	//- Runtime type information
	 TypeName("DSRFG");

  //-Constructors

	//- Construct from patch and internal field
	DSRFG1FvPatchVectorField
	(
	  const fvPatch&,
          const DimensionedField<vector, volMesh>&
	);

	//- Construct from patch, internal field and dictionary
         DSRFG1FvPatchVectorField
         (
             const fvPatch&,
             const DimensionedField<vector, volMesh>&,
             const dictionary&
         );

	//- Construct by mapping given DSRFG1FvPatchVectorField
         //  onto a new patch
         DSRFG1FvPatchVectorField
         (
             const DSRFG1FvPatchVectorField&,
             const fvPatch&,
             const DimensionedField<vector, volMesh>&,
             const fvPatchFieldMapper&
         );

         //- Construct as copy
         DSRFG1FvPatchVectorField
         (
             const DSRFG1FvPatchVectorField&
         );

	//- Construct and return a clone
         virtual tmp<fvPatchField<vector> > clone() const
        {
             return tmp<fvPatchField<vector> >
            (
                new DSRFG1FvPatchVectorField(*this)
             );
        }

         //- Construct as copy setting internal field reference
         DSRFG1FvPatchVectorField
         (
            const DSRFG1FvPatchVectorField&,
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
                 new DSRFG1FvPatchVectorField(*this, iF)
             );
         }

//Destructor
//  ~DSRFG1FvPatchVectorField();

// Member functions

	//Calculating private Members

		//-Calculate Deltay's and Distances from wall

		void CalculateDeltas();

		//-Calculate k0 and Ls

		void CalcK0AndLs();

		//-Calculates K parts of Spectrum
		void EnergySpectrum();

		//-Calculate P, Q, K, Omega Fields
		void RandomValuesToPQK();

		//- Input from file to R

		void InputR();

		//- Input from file to U
		void InputU();

        // Access
		//- Return the Wave Number
		inline scalar k0() const
		{
			return k0_;
		}

		//- Return the Turbulence Number
		inline scalar Ls() const
		{
			return Ls_;
		}

	// Mapping functions

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



           // Evaluation functions



               //- Update the coefficients associated with the patch field
               virtual void updateCoeffs();


           	//- Write
           	virtual void write(Ostream&) const;


  };

} // End namespace Foam

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
 //************************************************************************* //
