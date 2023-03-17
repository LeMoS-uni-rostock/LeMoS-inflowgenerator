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

#ifndef gaussianSpot_H
#define gaussianSpot_H

#include "turbulentStructure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gaussianSpot;

Ostream& operator<<(Ostream&, const gaussianSpot&);
Istream& operator>>(Istream&, gaussianSpot&);

/*---------------------------------------------------------------------------*\
                           Class gaussianSpot Declaration
\*---------------------------------------------------------------------------*/

class gaussianSpot
: public turbulentStructure
{

protected:
    vector epsilon_;

public:

    // Constructors
    gaussianSpot(const Parameters& p, Istream&);
    gaussianSpot
    (
        const Parameters& p,
        BoostRandomGen&, 
        const point& footPoint, 
        const vector& initialDelta,
        label creaface
    );
    gaussianSpot(const gaussianSpot& p);
    
    ~gaussianSpot();
    
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

    void operator=(const gaussianSpot&);
    bool operator!=(const gaussianSpot&) const;


    autoPtr<gaussianSpot> clone() const;

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&, const gaussianSpot&);
    friend Istream& operator>>(Istream&, gaussianSpot&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
