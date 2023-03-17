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

#include "gaussianSpot.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"

namespace Foam
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gaussianSpot::gaussianSpot
(
    const Parameters& p,
    Istream& s
)
: turbulentStructure(p, s),
  epsilon_(s)
{}


gaussianSpot::gaussianSpot
(
    const Parameters& p,
    BoostRandomGen& r, 
    const point& footPoint, 
    const vector& initialDelta,
    label creaface
)
: turbulentStructure(p, r, footPoint, initialDelta, creaface),
  epsilon_(pTraits<vector>::one)
{}


gaussianSpot::gaussianSpot(const gaussianSpot& o)
: turbulentStructure(o),
  epsilon_(o.epsilon_)
{}


vector gaussianSpot::fluctuation(const vector& x) const
{
    vector delta_x = x - location();

    vector L=p_.inputData().L()[creaFace_];
    vector e1=p_.e1()[creaFace_], e2=p_.e2()[creaFace_], e3=p_.e3()[creaFace_];
    
    if 
        (
            (mag(delta_x&e1)  < L[0]) &&
            (mag(delta_x&e2)  < L[1]) &&
            (mag(delta_x&e3)  < L[2])
        )
    {
      vector f=
           exp(- 4.0*magSqr(delta_x&e1)  / sqr(L[0]) )
          *exp(- 4.0*magSqr(delta_x&e2)  / sqr(L[1]) )
          *exp(- 4.0*magSqr(delta_x&e3)  / sqr(L[2]) )
          * pTraits<vector>::one;

      return p_.Lund(creaFace_) & 
      (
        cmptMultiply(epsilon_, f) / sqrt( 0.0820292 )
      );
    }
  else
    return pTraits<vector>::zero;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //



void gaussianSpot::randomize(BoostRandomGen& rand)
{
  turbulentStructure::randomize(rand);
  epsilon_.x() = 2.0*(rand() - 0.5);
  epsilon_.y() = 2.0*(rand() - 0.5);
  epsilon_.z() = 2.0*(rand() - 0.5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gaussianSpot::~gaussianSpot()
{}


autoPtr<gaussianSpot> gaussianSpot::clone() const
{
    return autoPtr<gaussianSpot>
        (
            new gaussianSpot(*this)
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void gaussianSpot::operator=(const gaussianSpot& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("gaussianSpot::operator=(const gaussianSpot&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    turbulentStructure::operator=(rhs);
    epsilon_=rhs.epsilon_;
}

bool gaussianSpot::operator!=(const gaussianSpot& o) const
{
    return 
        (location()!=o.location())
        ||
        (epsilon_!=o.epsilon_);
}

Ostream& operator<<(Ostream& s, const gaussianSpot& ht)
{
    s << *static_cast<const turbulentStructure*>(&ht);
    writeComponents(s, ht.epsilon_);
    return s;
}

Istream& operator>>(Istream& s, gaussianSpot& ht)
{
    s >> *static_cast<turbulentStructure*>(&ht);
    readComponents(s, ht.epsilon_);
    return s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
