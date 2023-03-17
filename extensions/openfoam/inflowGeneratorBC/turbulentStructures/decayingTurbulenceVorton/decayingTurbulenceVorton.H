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

#ifndef decayingTurbulenceVorton_H
#define decayingTurbulenceVorton_H

#include "turbulentStructure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class decayingTurbulenceVorton;

Ostream& operator<<(Ostream&, const decayingTurbulenceVorton&);
Istream& operator>>(Istream&, decayingTurbulenceVorton&);

/*---------------------------------------------------------------------------*\
                           Class decayingTurbulenceVorton Declaration
\*---------------------------------------------------------------------------*/

class decayingTurbulenceVorton
: public turbulentStructure
{

protected:
    vector epsilon_;

public:

    // Constructors
    decayingTurbulenceVorton(const Parameters& p, Istream&);
    decayingTurbulenceVorton
    (
        const Parameters& p,
        BoostRandomGen&, 
        const point& footPoint, 
        const vector& initialDelta,
        label creaface        
    );
    decayingTurbulenceVorton(const decayingTurbulenceVorton& p);
    
    ~decayingTurbulenceVorton();
    
    tensor parameters() const
    {
      return tensor
      (
        epsilon_[0], epsilon_[1], epsilon_[2],
        0., 0., 0.,
        0., 0., 0.
      );
    }

    vector fluctuation(const vector&) const;
    
    // Check

    // Edit
    void randomize(BoostRandomGen&);

    void operator=(const decayingTurbulenceVorton&);
    bool operator!=(const decayingTurbulenceVorton&) const;

    // Required for storage in SLList
    autoPtr<decayingTurbulenceVorton> clone() const;

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&, const decayingTurbulenceVorton&);
    friend Istream& operator>>(Istream&, decayingTurbulenceVorton&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
