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
#ifndef TESTINFLOWGENERATORRS_H
#define TESTINFLOWGENERATORRS_H


#include <openfoam/openfoamanalysis.h>


namespace insight
{
  
  
class TestInflowGeneratorRS
    : insight::OpenFOAMAnalysis
{
public:
  
#include "testinflowgeneratorrs__TestInflowGeneratorRS__Parameters.h"
  
/*
PARAMETERSET>>> TestInflowGeneratorRS Parameters
inherits OpenFOAMAnalysis::Parameters

geometry=set
{
  H 		= double 	2 		"total height of channel inlet face"
  B 		= double 	3 		"total width of channel inlet face"
} "Geometrical properties of the jet mixer"

mesh=set
{
nh = int 64 "# cells in vertical direction"
fixbuf = bool false "fix cell layer size inside buffer layer"
dzplus = double 15 "Dimensionless grid spacing in spanwise direction"
ypluswall = double 2 "yPlus at the wall grid layer"
} "Properties of the computational mesh"
operation=set
{
 spottype=selection (
  hatSpot
  gaussianSpot
  decayingTurbulenceSpot
  decayingTurbulenceVorton
  anisotropicVorton
  modalTurbulence
  ) anisotropicVorton "Type of inflow generator"

} "Definition of the operation point under consideration"
fluid=set
{
  rho		= double 	998.0 		"[kg/m^3] Density of the fluid"
  nu		= double 	1e-6 		"[m^2/s] Viscosity of the fluid"
} "Parameters of the fluid"
evaluation=set
{
  meantime	= double 	10.0 		"[T] length of time period for averaging of velocity and RMS (as multiple of flow-through time)"
} "Parameters of the result evaluation"
<<<PARAMETERSET
*/

protected:
    // derived data
  
public:
    declareType("TestInflowGeneratorRS");
    
    TestInflowGeneratorRS(const ParameterSet& ps, const boost::filesystem::path& exepath);
    
    static std::string category() { return "Validation Cases/Inflow Generator"; }
    
    virtual void calcDerivedInputData();
    
    virtual void createCase(insight::OpenFOAMCase& cm);
    virtual void createMesh(insight::OpenFOAMCase& cm);
    
    virtual ResultSetPtr evaluateResults(OpenFOAMCase& cmp);
};
}
#endif // TESTINFLOWGENERATORRS_H
