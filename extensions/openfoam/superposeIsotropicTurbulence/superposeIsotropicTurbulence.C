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

#include "pipe.h"
#include "channel.h"

#include "fvCFD.H"
#include "OFstream.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fixedGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"
#include "wallDist.H"
#include "interpolationTable.H"

#include <boost/concept_check.hpp>
#include <boost/assign.hpp>

#include "base/boost_include.h"
#include "boostRandomGen.H"
#include "TurbulentSpotVolumeInserter.H"
#include "meshSearch.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validArgs.append("structure type");
    argList::validArgs.append("k");
    argList::validArgs.append("l");
    
    using Foam::vector; // fix issue in OF5
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    
    wallDist y(mesh);
    volVectorField gradyw=fvc::grad(y.y());
    gradyw/=mag(gradyw)+dimensionedScalar("", gradyw.dimensions(), Foam::SMALL);
    gradyw*=y.y();
  
    Foam::word structureType(UNIOF_ADDARG(args,0));
    Foam::scalar k=readScalar(IStringStream(UNIOF_ADDARG(args,1))());
    Foam::scalar l=readScalar(IStringStream(UNIOF_ADDARG(args,2))());
    
    boundBox bb(mesh.points());
    
    double c=2.;
    Foam::vector Lbb=bb.max()-bb.min();
    double Vbb=Lbb.x()*Lbb.y()*Lbb.z();
    int ns = ceil( c*Vbb / pow(l,3) );
    
    Info<<"Inserting "<<ns<<" turbulent structures."<<endl;

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
	IOobject
	(
	    "U",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );
    
    dictionary id;
    id.add("Umean", "uniform ( 1e-6 1e-6 1e-6 )" );
    double r=2.*k/3.;
    id.add("R", boost::str(boost::format("uniform (%g 0 0 %g 0 %g)")%r%r%r).c_str() );
    id.add("L", boost::str(boost::format("uniform (%g %g %g)")%l%l%l).c_str() );
    id.add("c", "uniform 1");
    
    Info<<id;
    
    inflowInputDataField input(mesh.points(), runTime, id);
    autoPtr<turbulentSpotVolumeInserter> ins = turbulentSpotVolumeInserter::New(structureType, mesh, "", input, id);
    meshSearch ms(mesh);
    
    Random ran(0);
    for (int i=0; i<ns; i++)
    {
        point p = bb.min() + cmptMultiply( bb.max()-bb.min(), ran.
#if OF_VERSION>=060000 //defined(OFesi1806)
                                           sample01<vector>()
#else
                                           vector01()
#endif
                                           );
        label j=ms.findCell(p);
        if (j>=0)
        {
            ins->insertSpot(p, j); 
        }
    }
    
    Info<<"Spots generated, computing fluctuations."<<endl;
    
    U+=ins->uPrime();
        
    U.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
