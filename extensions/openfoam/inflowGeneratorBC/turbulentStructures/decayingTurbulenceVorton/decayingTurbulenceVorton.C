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

#include "decayingTurbulenceVorton.H"
#include "inflowGeneratorBaseFvPatchVectorField.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


decayingTurbulenceVorton::decayingTurbulenceVorton
(
    const Parameters& p,
    Istream& s
)
: turbulentStructure(p, s),
  epsilon_(s)
{}

decayingTurbulenceVorton::decayingTurbulenceVorton
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


decayingTurbulenceVorton::decayingTurbulenceVorton(const decayingTurbulenceVorton& o)
: turbulentStructure(o),
  epsilon_(o.epsilon_)
{}

vector decayingTurbulenceVorton::fluctuation(const vector& x) const
{
    vector delta_x = x - location();

    scalar L3=
        ( p_.inputData().L()[creaFace_].x()
         * p_.inputData().L()[creaFace_].y()
         * p_.inputData().L()[creaFace_].z()
        );
    scalar L = pow(L3, 1./3.);

    scalar c=p_.inputData().c()[creaFace_];

    vector f=
       (M_PI/L) * exp( -0.5*M_PI*magSqr(delta_x)/pow(L,2) ) * (epsilon_^delta_x);

    return
        p_.Lund(creaFace_) &
        (
          f
          /
          sqrt( c*M_PI/3. )
        );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


void decayingTurbulenceVorton::randomize(BoostRandomGen& rand)
{
  turbulentStructure::randomize(rand);
  epsilon_.x() = 2.0*(rand() - 0.5);
  epsilon_.y() = 2.0*(rand() - 0.5);
  epsilon_.z() = 2.0*(rand() - 0.5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

decayingTurbulenceVorton::~decayingTurbulenceVorton()
{}


autoPtr<decayingTurbulenceVorton> decayingTurbulenceVorton::clone() const
{
    return autoPtr<decayingTurbulenceVorton>
        (
            new decayingTurbulenceVorton(*this)
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void decayingTurbulenceVorton::operator=(const decayingTurbulenceVorton& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("decayingTurbulenceVorton::operator=(const decayingTurbulenceVorton&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    turbulentStructure::operator=(rhs);
    epsilon_=rhs.epsilon_;
}

bool decayingTurbulenceVorton::operator!=(const decayingTurbulenceVorton& o) const
{
    return 
        (location()!=o.location())
        ||
        (epsilon_!=o.epsilon_);
}

Ostream& operator<<(Ostream& s, const decayingTurbulenceVorton& ht)
{
    s << *static_cast<const turbulentStructure*>(&ht);
    writeComponents(s, ht.epsilon_);
    return s;
}

Istream& operator>>(Istream& s, decayingTurbulenceVorton& ht)
{
    s >> *static_cast<turbulentStructure*>(&ht);
    readComponents(s, ht.epsilon_);
    return s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
