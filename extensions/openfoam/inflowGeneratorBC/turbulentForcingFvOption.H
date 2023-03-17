
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

#ifndef turbulentForcingFvOption_H
#define turbulentForcingFvOption_H



#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Random.H"
#if OF_VERSION>=040000 //(defined(OFplus)||defined(OFdev)||defined(OF301))||defined(OFesi1806)
#include "cellSetOption.H"
#else
#include "fvOption.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class turbulentForcing Declaration
\*---------------------------------------------------------------------------*/

class turbulentForcing
:
#if OF_VERSION>=040000 //(defined(OFplus)||defined(OFdev)||defined(OF301))||defined(OFesi1806)
    public cellSetOption
#else
    public option
#endif
{
    // Private data

        label curTimeIndex_;
        
        //- drag coefficients along the coordinate directions
        scalar intensity_, relax_;

        vectorField fluctuation_;
        
        Random rand_;
	
        //- Disallow default bitwise copy construct
        turbulentForcing(const turbulentForcing&);

        //- Disallow default bitwise assignment
        void operator=(const turbulentForcing&);


public:

    //- Runtime type information
    TypeName("turbulentForcing");


    // Constructors

        //- Construct from explicit source name and mesh
        turbulentForcing
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);

        // I-O

//             //- Write the source properties
//             virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

