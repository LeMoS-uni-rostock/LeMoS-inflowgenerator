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

#include "inflowstatisticstest.h"

#include "base/factory.h"

#include "openfoam/openfoamcase.h"
#include "openfoam/openfoamtools.h"
#include "openfoam/blockmesh.h"

#include "openfoam/caseelements/analysiscaseelements.h"

#include "refdata.h"

using namespace std;
using namespace boost;
using namespace boost::assign;

namespace insight
{
  
addToAnalysisFactoryTable(InflowStatisticsTest);

InflowStatisticsTest::InflowStatisticsTest(const ParameterSet& ps, const boost::filesystem::path& exepath)
: ChannelInflow(ps, exepath)
{
}


void InflowStatisticsTest::createCase(OpenFOAMCase& cm, ProgressDisplayer& progress)
{
    insight::ChannelInflow::createCase(cm, progress);
    cm.remove<extendedForces>(".*");
}

void InflowStatisticsTest::calcDerivedInputData(ProgressDisplayer& progress)
{
  Parameters p(parameters_);
  insight::ChannelInflow::calcDerivedInputData(progress);
  dt_=double(p.ChannelBase::Parameters::geometry.L/nax_)/Ubulk_;
  nax_=1;
}

void InflowStatisticsTest::applyCustomOptions(OpenFOAMCase& cm, std::shared_ptr< OFdicts >& dicts)
{
    insight::ChannelInflow::applyCustomOptions(cm, dicts);
    
    startTimeName_=str(format("%g") % avgStart_);
    OFDictData::dict& controlDict=dicts->lookupDict("system/controlDict");
    controlDict["application"]="checkInflowGeneratorStatistics";
    controlDict["startFrom"]="startTime";
    controlDict["startTime"]=startTimeName_;
    controlDict["endTime"]=avg2Start_;
    controlDict["deltaT"]=dt_; // value in ChannelBase will be computed from axial resolution nax=1...
    
    controlDict["writeInterval"]=1000;
    controlDict["writeControl"]="timeStep";
}

void InflowStatisticsTest::applyCustomPreprocessing(OpenFOAMCase& cm, ProgressDisplayer& progress)
{
    insight::ChannelInflow::applyCustomPreprocessing(cm, progress);
    
    boost::filesystem::remove_all(executionPath()/startTimeName_);
    copyFields
    (
      executionPath()/"0", 
      executionPath()/startTimeName_
    );
    
    OFDictData::dictFile checkInflowGeneratorStatisticsDict;
    
    OFDictData::list patchlist;
    patchlist.push_back(cycl_in_);
    checkInflowGeneratorStatisticsDict["patches"]=patchlist;
    
    int nsteps=std::max( int((avg2Start_-avgStart_)/dt_), 1 );
    checkInflowGeneratorStatisticsDict["nSteps"]=nsteps;

    checkInflowGeneratorStatisticsDict["writeInterval"]=std::min(100, nsteps);

    checkInflowGeneratorStatisticsDict["doFullTimeLoop"]=true;
    
    // then write to file
    boost::filesystem::path dictpath = executionPath() / "system" / "checkInflowGeneratorStatisticsDict";
//     if (!exists(dictpath.parent_path())) 
//     {
//       boost::filesystem::create_directories(dictpath.parent_path());
//     }    
    {
      std::ofstream f(dictpath.c_str());
      writeOpenFOAMDict(f, checkInflowGeneratorStatisticsDict, boost::filesystem::basename(dictpath));
    }    
}


ResultSetPtr InflowStatisticsTest::evaluateResults(OpenFOAMCase& cmp, ProgressDisplayer& progress)
{
    Parameters p(parameters_);
    
    ResultSetPtr results=insight::OpenFOAMAnalysis::evaluateResults(cmp, progress);
    Ordering o;

    insertBCIntoResults(results);
    
    RFieldName_="UPrime2Mean";
    UMeanName_="UMean";

    evaluateAtSection(cmp, results, (1e-6-0.5)*p.ChannelBase::Parameters::geometry.L, 0, o, p.eval.include_channelrefdata, p.eval.include_allRcmp);

    return results;
}

}
