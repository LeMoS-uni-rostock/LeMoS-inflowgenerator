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
 */

#ifndef FOAM_SYNTHETICTURBULENCECASEELEMENTS_H
#define FOAM_SYNTHETICTURBULENCECASEELEMENTS_H

#include "openfoam/boundaryconditioncaseelements.h"

namespace insight
{
    

class SyntheticTurbulenceVelocityInletBC
: public VelocityInletBC
{
 
public:

#include "syntheticturbulencecaseelements__VelocityInletBC__Parameters.h"
/*
PARAMETERSET>>> SyntheticTurbulenceVelocityInletBC Parameters

turbulence_synthesis=selectablesubset {{

 uniform 
 set { 
   values=array
   [
    set {
     time=double 0 "time description"
     value=vector (1 0 0) "value description"
    } "desc"
   ] * 1  "values1 description"
 }


}} uniform "desc"
<<<PARAMETERSET
*/


  
protected:
  Parameters p_;
  
public:
  SyntheticTurbulenceVelocityInletBC
  (
    OpenFOAMCase& c,
    const std::string& patchName, 
    const OFDictData::dict& boundaryDict, 
    Parameters const& p = Parameters()
  );

//   virtual void setField_U(OFDictData::dict& BC) const;
//   virtual void setField_p(OFDictData::dict& BC) const;
//   virtual void addIntoFieldDictionaries(OFdicts& dictionaries) const;
};


}

#endif
