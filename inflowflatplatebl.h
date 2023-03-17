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
#ifndef INFLOWFLATPLATEBL_H
#define INFLOWFLATPLATEBL_H

#include "flatplatebl.h"
#include "openfoam/openfoamanalysis.h"

namespace insight
{
  
class InflowFlatPlateBL
: public FlatPlateBL
{
  
public:
#include "inflowflatplatebl__InflowFlatPlateBL__Parameters.h"
  
/*
PARAMETERSET>>> InflowFlatPlateBL Parameters
inherits insight::FlatPlateBL::Parameters

operation=set
{
  Ushape 	= selectablesubset
  {{
  
    block
    set {
    }
    
    powerlaw
    set {
    }
    
    tabulated
    set {
      uplus_vs_yplus  = matrix 1x4		"Profile of normalized mean velocity (u+,v+,w+) vs normalized wall distance y+"
    }
    
  }} powerlaw "Shape of the velocity"
  
  turb 	= selectablesubset
  {{
    uniformisotropic
    set {
      intensity 	= double 	0.05		"uniform and isotropic turbulence intensity"
    }
    
    universal
    set {
    }
    
    tabulated
    set {
      Rplus_vs_yplus  = matrix 1x7		"Profile of normalized Reynolds stresses (Rxx+, Rxy+, Rxz+, Ryy+, Ryz+, Rzz+) vs normalized wall distance y+"
    }
    
  }} universal "Turbulence specification in jet"
  
  lengthscales 	= selectablesubset
  {{
    uniform
    set {
      uniformLbydelta99 	= vector (1 1 1)		"uniform integral length scales"
    }
    
    tabulated
    set {
      Lbydelta99_vs_yplus  = matrix 1x4		"Profile of integral length scale, normalized by BL thickness delta_99 (Lxx, Lyy, Lzz) vs normalized wall distance y+"
    }
    
  }} uniform "Turbulence specification boundary layer"

  spottype=selection (
    hatSpot
    gaussianSpot
    decayingTurbulenceSpot
    decayingTurbulenceVorton
    anisotropicVorton_Analytic
    anisotropicVorton_PseudoInv
    anisotropicVorton_NumOpt
    anisotropicVorton2
    combinedVorton
    modalTurbulence
  ) anisotropicVorton_PseudoInv "Type of inflow generator"

} "Definition of the operation point under consideration"

<<<PARAMETERSET
*/

protected:
    struct supplementedInputData
        : public supplementedInputDataDerived<Parameters, FlatPlateBL::supplementedInputData>
    {

      supplementedInputData(
          std::unique_ptr<Parameters> pPtr,
          const boost::filesystem::path& workDir,
          ProgressDisplayer& progress = consoleProgressDisplayer
          );
    };

#ifndef SWIG
  defineDerivedClassWithSupplementedInputData(Parameters, supplementedInputData)
#endif

public:
    declareType("InflowFlatPlateBL");
    
    InflowFlatPlateBL(
            const ParameterSet& ps,
            const boost::filesystem::path& exepath,
            ProgressDisplayer& progress );

    static std::string category() { return "Validation Cases/Inflow Generator"; }

    virtual void createInflowBC(insight::OpenFOAMCase& cm, const OFDictData::dict& boundaryDict) const;
};

}
#endif // INFLOWFLATPLATE_H
