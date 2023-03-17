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

#ifndef TurbulentSpotVolumeInserter_H
#define TurbulentSpotVolumeInserter_H

#include "boostRandomGen.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"

#include "boost/shared_ptr.hpp"
#include <set>

namespace Foam
{
  
class turbulentSpotVolumeInserter
{
  
public:
  //- Runtime type information
  TypeName("turbulentSpotVolumeInserter");
  
  // Declare run-time constructor selection table
  declareRunTimeSelectionTable
  (
      autoPtr,
      turbulentSpotVolumeInserter,
      components,
      (
        const fvMesh& mesh, 
        const word& cellZoneName,
        const inflowInputDataField& input,
        const dictionary& dict
      ),
      (mesh, cellZoneName, input, dict)
  );
  
  //- Select constructed from dictionary and objectRegistry
  static autoPtr<turbulentSpotVolumeInserter> New
  (
    const word& typeName,
    const fvMesh& mesh, 
    const word& cellZoneName,
    const inflowInputDataField& input,
    const dictionary& dict
  );

  virtual ~turbulentSpotVolumeInserter();
  virtual void insertSpot
  (
    const point& pt,
    label I
  ) =0;
  virtual tmp<volVectorField> uPrime() const =0;
  
};




template<class TS>
class TurbulentSpotVolumeInserterImpl
: public turbulentSpotVolumeInserter
{
  BoostRandomGen rg_;
  
  const fvMesh& mesh_;
//   const cellZone& cz_;

  const inflowInputDataField& input_;
  typename TS::Parameters sp_;
  
  typedef std::shared_ptr<TS> TurbulentStructurePtr;
  typedef std::set<TurbulentStructurePtr> TurbulentStructureSet;
  TurbulentStructureSet vortons_;
  
  void createSpots();
  
public:
  //- Runtime type information
  TypeName("TurbulentSpotVolumeInserterImpl");
  
  TurbulentSpotVolumeInserterImpl
  (
      const fvMesh& mesh, 
      const word& cellZoneName,
      const inflowInputDataField& input,
      const dictionary& dict
  );
  
  virtual void insertSpot
  (
    const point& pt,
    label I
  );
  
  virtual tmp<volVectorField> uPrime() const;
  
};
  
}

#include "TurbulentSpotVolumeInserter.C"

#endif

