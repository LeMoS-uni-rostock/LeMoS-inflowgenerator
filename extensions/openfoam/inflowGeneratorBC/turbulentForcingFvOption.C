
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

#if OF_VERSION>=010700 //ndef OF16ext

#include "turbulentForcingFvOption.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulentForcing, 0);

    addToRunTimeSelectionTable
    (
        option,
        turbulentForcing,
        dictionary
    );
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulentForcing::turbulentForcing
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
#if OF_VERSION>=040000 //(defined(OFplus)||defined(OFdev)||defined(OF301))||defined(OFesi1806)
    cellSetOption(sourceName, modelType, dict, mesh),
#else
    option(sourceName, modelType, dict, mesh),
#endif
    curTimeIndex_(-1),
    intensity_(readScalar(coeffs_.lookup("intensity"))),
    relax_(readScalar(coeffs_.lookup("relax"))),
    fluctuation_(cells_.size(), vector::zero),
    rand_(-1)
{
    fieldNames_.setSize(1);
    fieldNames_[0]="U";
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    Setting turbulent forcing with intensity = " << intensity_ << " and relaxation = " << relax_ << nl << endl;
}

bool Foam::fv::turbulentForcing::read(const Foam::dictionary& dict)
{
  if (option::read(dict))
  {
    coeffs_.lookup("intensity") >> intensity_;
    coeffs_.lookup("relax") >> relax_;
    
    return true;
  }
  else
  {
    return false;
  }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbulentForcing::correct(volVectorField& U)
{
    if ( curTimeIndex_ != this->mesh().time().timeIndex()  ) // update only once per time step
    {
        scalar magUmean=mag(gAverage(UIndirectList<vector>(U, cells_)()));
        forAll(fluctuation_, i)
        {
            fluctuation_[i] = (1.-relax_)*fluctuation_[i] + relax_* (2.*
                                                         #if OF_VERSION>=060000 //defined(OFesi1806)
                                                                     rand_.sample01<vector>()
                                                         #else
                                                                     rand_.vector01()
                                                         #endif
                                                                     -vector::one)*intensity_*magUmean;
        }
        
        Info<<"fluctuation: min/max/avg = "<<gMin(fluctuation_)<<" / "<<gMax(fluctuation_)<<" / "<<gAverage(fluctuation_)<<endl;
        curTimeIndex_ = this->mesh().time().timeIndex();
    }
}

void Foam::fv::turbulentForcing::constrain(fvMatrix<vector>& eqn, const label fieldi)
{
    eqn.setValues(cells_, UIndirectList<vector>(eqn.psi(), cells_)() + fluctuation_);
}


// ************************************************************************* //
#endif
