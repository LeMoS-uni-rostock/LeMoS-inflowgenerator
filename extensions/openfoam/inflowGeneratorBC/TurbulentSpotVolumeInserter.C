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
#include "fvCFD.H"

#include "uniof.h"

namespace Foam
{
  


const cellZone& findZoneOrFail(const fvMesh& mesh, const word& cellZoneName)
{
  label i = mesh.cellZones().findZoneID(cellZoneName);
  
  if (i<0)
  {
    FatalErrorIn("TurbulentSpotVolumeInserterImpl")
    <<"cell zone "<<cellZoneName<<" not found in mesh!"
    <<endl<<abort(FatalError);
  }
  
  return mesh.cellZones()[i];
}

template<class TS>
TurbulentSpotVolumeInserterImpl<TS>::TurbulentSpotVolumeInserterImpl
(
    const fvMesh& mesh, 
    const word& cellZoneName,
    const inflowInputDataField& input,
    const dictionary& dict
)
: mesh_(mesh),
  //cz_(findZoneOrFail(mesh_, cellZoneName))
  input_(input),
  sp_(input_, dict)
{
}

template<class TS>
void TurbulentSpotVolumeInserterImpl<TS>::insertSpot
(
  const point& pt,
  label I
)
{
  TurbulentStructurePtr snew;
  snew.reset
  (
    new TS
    (
      sp_, rg_, pt, 
      Foam::vector::zero,
      I
    )
  );
  vortons_.insert(snew);
}






template<class TurbulentStructure, class Type=Foam::vector>
class RecursiveVolApply
{
public:
    typedef 
        Type (TurbulentStructure::*FluctuationFunctionPtr)(const vector& x) const 
        ;
  
protected:
  const fvMesh& mesh_;
  Field<Type>& fluctuations_;
  FluctuationFunctionPtr func_;
  
  /**
   * if applied value falls below this tolerance, recursive search will be stopped
   */
  static constexpr const scalar tol_ = 1e-3;
  
  void recursiveapply
  (
    const TurbulentStructure& v, 
    label I, 
    labelHashSet& visited,
    label depth,
    scalar sf
  )
  {
    depth++;
    if (I<0) return;
    
//     vector u = v.fluctuation(patch_.faceCentres()[I]);
    Type u = (v.*func_)(mesh_.C()[I]);
    
    if (mag(u)>0.0/*tol_*/)
    {
      
      fluctuations_[I] += u/* *sf*/;
      visited.insert(I);
      
      labelList lneigh=mesh_.cellCells()[I];
      
      // visit local neighbours
      forAll(lneigh, j)
      {
        label nfi=lneigh[j];
        if (!visited.found(nfi))
        {
            recursiveapply(v, nfi, visited, depth, sf);
        }
      }
    }
  }

public:
  RecursiveVolApply
  (
    const fvMesh& mesh,
    Field<Type>& fluctuations,
    FluctuationFunctionPtr func = &TurbulentStructure::fluctuation
  )
  : mesh_(mesh),
    fluctuations_(fluctuations),
    func_(func)
  {
  }
  
  
  label apply
  (
    const TurbulentStructure& v,
    label startfaceI,
    scalar sf
  )
  {
    labelHashSet visited;
    recursiveapply(v, startfaceI, visited, 0, sf);
    return visited.size();
  }
  
};


template<class TS>
tmp<volVectorField> TurbulentSpotVolumeInserterImpl<TS>::uPrime() const
{
  tmp<volVectorField> uPrime
  (
    new volVectorField
    (
      IOobject
      (
        "uPrime",
        mesh_.time().timeName(),
        mesh_.time(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh_,
      dimensionedVector("uPrimeDefault", dimVelocity, vector::zero)
    )
  );

  
  RecursiveVolApply<TS> apl(uPrime().mesh(), UNIOF_TMP_NONCONST(uPrime));
  labelList n_affected(vortons_.size(), 0);
  size_t j=0;
  BOOST_FOREACH(const TurbulentStructurePtr& vorton, vortons_)
  {
    n_affected[j++]=apl.apply(*vorton, vorton->creaFace(), 1.);
  }
  
  return uPrime;
}

}
