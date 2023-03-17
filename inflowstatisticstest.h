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
#ifndef INSIGHT_INFLOWSTATISTICSTEST_H
#define INSIGHT_INFLOWSTATISTICSTEST_H

#include "inflowchannel.h"

namespace insight
{
class InflowStatisticsTest
: public ChannelInflow
{
public:
#include "inflowstatisticstest__InflowStatisticsTest__Parameters.h"
/*
PARAMETERSET>>> InflowStatisticsTest Parameters

inherits ChannelInflow::Parameters

eval=set
{
  include_channelrefdata  = bool 		true		"Whether to include channel reference data in result plots"
  include_allRcmp  = bool 		false		"Whether to show all Reynolds stress components in result plots"
} "Parameters of the result evaluation"

<<<PARAMETERSET
*/
  
  double dt_;
  std::string startTimeName_;
  
protected:
    // derived data
  
public:
    declareType("InflowStatisticsTest");
    InflowStatisticsTest(const ParameterSet& ps, const boost::filesystem::path& exepath);

    static std::string category() { return "Validation Cases/Inflow Generator"; }
    
    virtual void calcDerivedInputData(ProgressDisplayer& progress);
    virtual void createCase(insight::OpenFOAMCase& cm, ProgressDisplayer& progress);
    virtual void applyCustomOptions(OpenFOAMCase& cm, std::shared_ptr<OFdicts>& dicts);
    virtual void applyCustomPreprocessing(OpenFOAMCase& cm, ProgressDisplayer& progress);
//     virtual void createMesh(insight::OpenFOAMCase& cm);
    virtual ResultSetPtr evaluateResults(OpenFOAMCase& cmp, ProgressDisplayer& progress);
};
}
#endif // INSIGHT_INFLOWSTATISTICSTEST_H
