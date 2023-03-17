#ifndef RECURSIVEAPPLY_H
#define RECURSIVEAPPLY_H

#include "point.H"
#include "Field.H"
#include "HashSet.H"

namespace Foam
{


template<class TurbulentStructure, class PatchType, class Type=Foam::vector>
class RecursiveApply
{
public:
    typedef
        Type (TurbulentStructure::*FluctuationFunctionPtr)(const vector& x) const
        ;

protected:
  const PatchType& patch_;
  Field<Type>& fluctuations_;
  FluctuationFunctionPtr func_;

  /**
   * if applied value falls below this tolerance, recursive search will be stopped
   */
  static constexpr const scalar tol_ = 1e-3;

  void recursiveApply
  (
    const TurbulentStructure& v,
    label faceI,
    labelHashSet& visited,
    label depth,
    scalar scaleFactor,
    scalar Lmax
  )
  {
    depth++;
    if (faceI<0) return;

//     vector u = v.fluctuation(patch_.faceCentres()[faceI]);
    point x = patch_.faceCentres()[faceI];
    Type u = (v.*func_)(x);

    if (mag(v - x)<Lmax)
    {

      fluctuations_[faceI] += u * scaleFactor;
      visited.insert(faceI);

      labelList lneigh=patch_.faceFaces()[faceI];

      // visit local neighbours
      forAll(lneigh, j)
      {
        label nfi=lneigh[j];
        if (!visited.found(nfi))
        {
            recursiveApply(v, nfi, visited, depth, scaleFactor, Lmax);
        }
      }
    }
  }

public:
  RecursiveApply
  (
    const PatchType& patch,
    Field<Type>& fluctuations,
    FluctuationFunctionPtr func = &TurbulentStructure::fluctuation
  )
  : patch_(patch),
    fluctuations_(fluctuations),
    func_(func)
  {
  }


  label apply
  (
    const TurbulentStructure& v,
    label startfaceI,
    scalar scaleFactor,
    scalar Lmax
  )
  {
    labelHashSet visited;
    recursiveApply(v, startfaceI, visited, 0, scaleFactor, Lmax);
    return visited.size();
  }

};

}

#endif // RECURSIVEAPPLY_H
