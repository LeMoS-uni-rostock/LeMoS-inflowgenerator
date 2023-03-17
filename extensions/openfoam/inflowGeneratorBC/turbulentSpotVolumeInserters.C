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

#include "TurbulentSpotVolumeInserter.H"
#include "turbulentStructure.H"
//#include "hatSpot.H"
// #include "gaussianSpot.H"
// #include "decayingTurbulenceSpot.H"
#include "decayingTurbulenceVorton.H"
#include "anisotropicVorton.H"
// #include "anisotropicVorton2.H"

#include "className.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(turbulentSpotVolumeInserter, 0);
defineRunTimeSelectionTable(turbulentSpotVolumeInserter, components);

autoPtr<turbulentSpotVolumeInserter> turbulentSpotVolumeInserter::New
(
  const word& typeName,
  const fvMesh& mesh, 
  const word& cellZoneName,
  const inflowInputDataField& input,
  const dictionary& dict
)
{

  componentsConstructorTable::iterator cstrIter =
      componentsConstructorTablePtr_->find(typeName);

  if (cstrIter == componentsConstructorTablePtr_->end())
  {
      FatalErrorIn
      (
	  "turbulentSpotVolumeInserter::New(const word& typeName, const fvMesh& mesh, const word& cellZoneName)"
      )   << "Unknown type "
	  << typeName << nl << nl
	  << "Valid coordinateSystem types are :" << nl
	  << componentsConstructorTablePtr_->sortedToc()
	  << exit(FatalIOError);
  }

  return autoPtr<turbulentSpotVolumeInserter>(cstrIter()(mesh, cellZoneName, input, dict)); 
}

turbulentSpotVolumeInserter::~turbulentSpotVolumeInserter()
{
}

#define makeSpotInserter(TS) \
typedef TurbulentSpotVolumeInserterImpl<TS> TS##VolumeInserter;\
defineTemplateTypeNameAndDebugWithName \
( \
    TS##VolumeInserter,\
    #TS,\
    0\
);\
addToRunTimeSelectionTable(turbulentSpotVolumeInserter, TS##VolumeInserter, components);

//makeSpotInserter(hatSpot);
// makeSpotInserter(gaussianSpot);
// makeSpotInserter(decayingTurbulenceSpot);
makeSpotInserter(decayingTurbulenceVorton);
makeSpotInserter(anisotropicVorton_PseudoInv);
// makeSpotInserter(anisotropicVorton2);

}
