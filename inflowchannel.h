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
#ifndef INSIGHT_INFLOWCHANNEL_H
#define INSIGHT_INFLOWCHANNEL_H

#include <channel.h>

namespace insight
{

class ChannelInflow
: public ChannelBase
{
  
#ifndef SWIG
  /**
   * locations (x/H) of section evaluations 
   */
  const std::vector<double> sec_locs_ = {
      0, 0.1, 0.25, 0.5, 0.75, 1,
      1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16
    };

  const static int ntpc_ = 4;
  const static char* tpc_names_[ntpc_]; 
  const static double tpc_xlocs_[ntpc_];
#endif
  
public:
  static void modifyDefaults(ParameterSet& p);

#include "inflowchannel__ChannelInflow__Parameters.h"
/*
PARAMETERSET>>> ChannelInflow Parameters
inherits ChannelBase::Parameters
addTo_makeDefault { modifyDefaults(p); }

inflow_scalar = selectablesubset
{{

 uniform
 set {
 }
 
 inflow_generator
 set {
  theta_mean = path "theta_mean.xy" "Path to profile theta_mean(y)"
  theta_variance = path "theta_variance.xy" "Path to profile theta_mean(y)"
  theta_turbulentFlux = path "theta_turbulentFlux.xy" "Path to profile theta_mean(y)"
 }

}} uniform "Scalar inflow generation method"

<<<PARAMETERSET
*/

protected:
#ifndef SWIG
  defineDerivedClassWithSupplementedInputData(Parameters, supplementedInputData)
#endif

public:
  declareType("Channel Flow Test Case (Inflow Generator)");
  
  ChannelInflow(
          const ParameterSet& ps,
          const boost::filesystem::path& exepath,
          ProgressDisplayer& progress );

  static std::string category() { return "Validation Cases/Inflow Generator"; }

  virtual void createMesh
  (
    OpenFOAMCase& cm, ProgressDisplayer& progress
  );  
  
  virtual void createCase
  (
    OpenFOAMCase& cm, ProgressDisplayer& progress
  );

  void insertBCIntoResults(ResultSetPtr results) const;
  ResultSetPtr evaluateResults(OpenFOAMCase& cm, ProgressDisplayer& progress);

  virtual void applyCustomOptions(OpenFOAMCase& cm, std::shared_ptr<OFdicts>& dicts);
  virtual void applyCustomPreprocessing(OpenFOAMCase& cm, ProgressDisplayer& progress);
  
};


}
#endif // INSIGHT_INFLOWCHANNEL_H
