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

#ifndef combinedVorton_H
#define combinedVorton_H

#include "inflowGeneratorBaseFvPatchVectorField.H"
#include "turbulentStructure.H"
#include "fielddataprovider.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class combinedVorton;

Ostream& operator<<(Ostream&, const combinedVorton&);
Istream& operator>>(Istream&, combinedVorton&);

/*---------------------------------------------------------------------------*\
                           Class combinedVorton Declaration
\*---------------------------------------------------------------------------*/

class combinedVorton
: public turbulentStructure
{

public:

  class Parameters
  : public turbulentStructure::Parameters
  {
  public:
      fileName sigmaFileName_, gammaFileName_;
      VectorSpaceBase coordinateSystem_;
      scalarField rx1_, ry1_, rz1_, rx2_, ry2_, rz2_;
      scalarField sx1_, sy1_, sz1_, sx2_, sy2_, sz2_;
      scalarField residual_sigma_, error_gamma_1_, error_gamma_2_;

      Parameters(const inflowInputDataField& ifp, const dictionary&);

      void write(Ostream&) const;

      virtual int nParamFields() const;
      virtual string paramFieldName(int i) const;
      virtual ParamField paramField(int i) const;
  };

protected:
    const Parameters& p_;
    scalar epsilon1_, epsilon2_;

public:

    // Constructors
    combinedVorton
    (
        const Parameters& p,
        Istream&
    );
    combinedVorton
    (
        const Parameters& p,
        BoostRandomGen&,
        const point& footPoint,
        const vector& initialDelta,
        label creaface
    );
    
    combinedVorton(const combinedVorton& p);
    
    ~combinedVorton();
    
    tensor parameters() const
    {
      return tensor
      (
        epsilon1_, epsilon2_, 0,
        0, 0, 0,
        0, 0, 0
      );
    }

    vector fluctuation(const vector&) const;
    
    // Check

    // Edit
    void randomize(BoostRandomGen&);

    void operator=(const combinedVorton&);
//     bool operator!=(const combinedVorton&) const;

    autoPtr<combinedVorton> clone() const;
    
    virtual vector randomParameter() const
    {
      return vector(epsilon1_, epsilon2_, 0);
    }

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&, const combinedVorton&);
    friend Istream& operator>>(Istream&, combinedVorton&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
