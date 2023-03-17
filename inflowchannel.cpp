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
#include "inflowchannel.h"

#include "base/factory.h"

#include "openfoam/blockmesh.h"
#include "openfoam/openfoamtools.h"
#include "openfoam/setfields.h"
#include "openfoam/openfoamcase.h"
#include "openfoam/fielddata.h"

#include "openfoam/caseelements/boundaryconditions/turbulentvelocityinletbc.h"
#include "openfoam/caseelements/boundaryconditions/pressureoutletbc.h"
#include "openfoam/caseelements/analysiscaseelements.h"

#include "refdata.h"


using namespace arma;
using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace boost::filesystem;




namespace insight {




const char* ChannelInflow::tpc_names_[] = 
  {
    "tpc0_inlet",
    "tpc1_intermediate1",
    "tpc2_intermediate2",
    "tpc3_intermediate3"
  };




const double ChannelInflow::tpc_xlocs_[] = {0.0, 0.125, 0.25, 0.375};




ChannelInflow::ChannelInflow(
        const ParameterSet& ps,
        const boost::filesystem::path& exepath,
        ProgressDisplayer& progress )
: ChannelBase(ps, exepath, progress)
{}




void ChannelInflow::modifyDefaults(ParameterSet& p)
{
  std::unique_ptr<SubsetParameter> inflowparams(new SubsetParameter(TurbulentVelocityInletBC::defaultParameters(), "Inflow BC"));
  p.extend
  (
    boost::assign::list_of<ParameterSet::SingleEntry>
    ("inflow", inflowparams.release())
    .convert_to_container<ParameterSet::EntryList>()
  );
}




void ChannelInflow::createMesh
(
  OpenFOAMCase& cm, ProgressDisplayer& progress
)
{  
  ChannelBase::createMesh(cm, progress);
}




void ChannelInflow::createCase
(
  OpenFOAMCase& cm, ProgressDisplayer& progress
)
{
  // create local variables from ParameterSet
//  PSDBL(parameters(), "geometry", H);
//  PSDBL(parameters(), "geometry", B);
//  PSDBL(parameters(), "geometry", L);
//  PSDBL(parameters(), "operation", Re_tau);
  
  OFDictData::dict boundaryDict;
  cm.parseBoundaryDict(executionPath(), boundaryDict);
      
  cm.insert(new TurbulentVelocityInletBC( cm, sp().cycl_in_, boundaryDict, parameters().getSubset("inflow") ));
  
  cm.insert(new PressureOutletBC(cm, sp().cycl_out_, boundaryDict, PressureOutletBC::Parameters()
    .set_behaviour( PressureOutletBC::Parameters::behaviour_fixMeanValue_type(
                      0.0
                     ))
    .set_behaviour(PressureOutletBC::Parameters::behaviour_fixMeanValue_type())
  ));
  
  ChannelBase::createCase(cm, progress);
  
  for (double xH: sec_locs_)
  {
      double x = -0.5*p().geometry.L + xH*p().geometry.H;
      if (x<=(0.5*p().geometry.L))
      {
        std::vector<std::string> sample_fields = {"p", "U"};
        
        if (p().operation.wscalar)
        {
            sample_fields.push_back("theta");
        }
        
        std::vector<arma::mat> pl;
        std::transform(
            sp().probe_locations_.begin(), sp().probe_locations_.end(),
            std::back_inserter(pl),
            [&](const arma::mat& p)
            {
                return vec3(x, p(1), p(2));
            }
        );
        
        cm.insert(new probes(cm, probes::Parameters()
            .set_fields( sample_fields )
            .set_probeLocations(pl)
            .set_name( boost::str(boost::format("probes_xH_%1.2f") % xH) )
            .set_outputControl("timeStep")
            .set_outputInterval(10.0)
        ));
      }
  }
  
  if (p().ChannelBase::Parameters::run.eval2)
  {
    for (int i=0; i<ntpc_; i++)
    {
      cm.insert(new LinearTPCArray(cm, LinearTPCArray::Parameters()
        .set_R(0.5*p().geometry.H)
    // 	.set_x(-0.5*L+tpc_xlocs_[i]*H)
    // 	.set_z(-0.49*B)
        .set_p0(vec3(
                    -0.5*p().geometry.L + tpc_xlocs_[i]*p().geometry.H,
                    0,
                    -0.49*p().geometry.B) )
        .set_axSpan(0.5*p().geometry.L)
        .set_tanSpan(0.45*p().geometry.B)
        .set_name(tpc_names_[i])
        .set_timeStart( sp().avg2Start_ )
      ));
    }
  }
  
}




void ChannelInflow::insertBCIntoResults(ResultSetPtr results) const
{
  
    FieldData( parameters().getSubset("inflow/umean") ).insertGraphsToResultSet
    (
        results, executionPath(),
        "inflowUmean",
        "Prescribed inflow mean velocity",
        "\\langle U \\rangle"
    );
    
    FieldData(
        parameters().get<SelectableSubsetParameter>("inflow/turbulence")().getSubset("R")
    ).insertGraphsToResultSet
    (
        results, executionPath(),
        "inflowR",
        "Prescribed inflow Reynolds stresses",
        "\\langle \\mathbf R \\rangle"
    );

    FieldData(
        parameters().get<SelectableSubsetParameter>("inflow/turbulence")().getSubset("L")
    ).insertGraphsToResultSet
    (
        results, executionPath(),
        "inflowL",
        "Prescribed inflow length scales",
        "L"
    );
}



 
ResultSetPtr ChannelInflow::evaluateResults(OpenFOAMCase& cm, ProgressDisplayer& progress)
{
//  const ParameterSet& p=parameters_;
//  PSDBL(p, "geometry", H);
//  PSDBL(p, "geometry", B);
//  PSDBL(p, "geometry", L);
//  PSDBL(p, "operation", Re_tau);
  
  ResultSetPtr results = ChannelBase::evaluateResults(cm, progress);
  
  Ordering o(2.);
  
  insertBCIntoResults(results);
    
  {
    
   // Pressure fluctuations
   arma::mat pPrime=interiorPressureFluctuationProfile(
               cm, executionPath(), vec3(1,0,0), sp().nax_);
   pPrime.col(0)-=pPrime.col(0)(0); // ensure, that curve starts at x=0
    
   arma::mat pPrime2_vs_xp(join_rows(
      pPrime.col(0)*p().operation.Re_tau,
      pPrime.col(1)
    ));
    pPrime2_vs_xp.save( (executionPath()/"pPrime2_vs_xp.txt").c_str(), raw_ascii);
    
    addPlot
    (
      results, executionPath(), "chartMeanpPrime2",
      "$x^+$", "$\\langle p\\prime^2 \\rangle$",
      {
        PlotCurve(pPrime2_vs_xp, "pPrime2", "w l lt 1 lc -1 lw 2 not")
      },
      "Axial profile of pressure fluctuation",
      "set logscale y;"
    ) .setOrder(o.next());    
  }
  
  // ============= Longitudinal profile of Velocity an RMS ================
  Ordering lseco(20.);
  int nr=10;
  for (int i=0; i<nr; i++)
  {
    double r0=0.1, r1=0.997*0.5*p().geometry.H;
    double r=r0+(r1-r0)*double(i)/double(nr-1);
    double yByH=r/p().geometry.H;
    
    std::shared_ptr<ResultSection> section
    (
      new ResultSection
      (
        str(format("Longitudinal section at $y/H=%.2f$")%yByH), ""
      )
    );
    Ordering so;
    
    string title="longitudinal__yByH_"+str(format("%07.3f")%yByH);
    replace_all(title, ".", "_");

    boost::ptr_vector<sampleOps::set> sets;
    
    double delta_yp1=1./p().operation.Re_tau;
    sets.push_back(new sampleOps::linearAveragedUniformLine(sampleOps::linearAveragedUniformLine::Parameters()
      .set_start( vec3(-0.5*p().geometry.L, r*0.5*p().geometry.H, -0.49*p().geometry.B))
      .set_end(   vec3(0.5*p().geometry.L, r*0.5*p().geometry.H, -0.49*p().geometry.B))
      .set_dir1(vec3(1,0,0))
      .set_dir2(vec3(0,0,0.98*p().geometry.B))
      .set_nd1(1)
      .set_nd2(sp().n_hom_avg)
      .set_name("longitudinal"+lexical_cast<string>(i))
    ));
    
    sample(cm, executionPath(), 
      {"p", UMeanName_, RFieldName_},
      sets
    );
    
    sampleOps::ColumnDescription cd;
    arma::mat data = static_cast<sampleOps::linearAveragedUniformLine&>(*sets.begin())
      .readSamples(cm, executionPath(), &cd);

      
    // Mean velocity profiles
    {
      int c=cd[UMeanName_].col;
      
      double fac_yp=p().operation.Re_tau*2.0/p().geometry.H;
      double fac_Up=1./sp().utau_;
      
      addPlot
      (
	section, executionPath(), "chartMeanVelocity_"+title,
        "$x^+$", "$\\langle U^+ \\rangle$",
        {
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Up*data.col(c))), "Up", "w l lt 1 lc -1 t '$U^+$'"),
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Up*data.col(c+1))), "Vp", "w l lt 2 lc -1 t '$V^+$'" ),
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Up*data.col(c+2))), "Wp", "w l lt 3 lc -1 t '$W^+$'" )
        },
	"Longitudinal profiles of averaged velocities at y/H="+str(format("%g")%yByH)
      ) .setOrder(so.next());
    }
    
    // Mean reynolds stress profiles
    {
      double fac_yp=p().operation.Re_tau*2.0/p().geometry.H;
      double fac_Rp=1./pow(sp().utau_,2);
      int c=cd[RFieldName_].col;
      
      addPlot
      (
	section, executionPath(), "chartMeanRstress_"+title,
        "$x^+$", "$\\langle R^+ \\rangle$",
        {
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Rp*data.col(c))), "Ruup", "w l lt 1 lc -1 t '$R_{uu}^+$'"),
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Rp*data.col(c+3))), "Rvvp", "w l lt 2 lc -1 t '$R_{vv}^+$'" ),
          PlotCurve( arma::mat(join_rows(fac_yp*data.col(0), fac_Rp*data.col(c+5))), "Rwwp", "w l lt 3 lc -1 t '$R_{ww}^+$'" )
        },
	"Longitudinal profiles of averaged reynolds stresses at y/H="+str(format("%g")%yByH)
      ) .setOrder(so.next());    
    }
    
    results->insert(title, section) .setOrder(lseco.next());
  }
  
  int i=0;
  Ordering xseco(10.); xseco.next(); // first sec was created in base class
  for (double xH: sec_locs_)
  {
    double x=-0.5*p().geometry.L+xH*p().geometry.H;
    if (xH==0.0)
      x=-0.5*p().geometry.L + 1e-6;
    
    if (x<=(0.5*p().geometry.L))
    {
      if (0.5*p().geometry.L-x<1e-6) x=0.5*p().geometry.L-1e-6; // avoid placing the section directly on the exit plane
      evaluateAtSection( cm, results, x, i+1, xseco, true, false, boost::str(boost::format("probes_xH_%1.2f") % xH) );
    }
    
    i++;
  }
    
  for (int i=0; i<ntpc_; i++)
  {
    
    const LinearTPCArray* tpcs=cm.get<LinearTPCArray>( string(tpc_names_[i])+"TPCArray");
    
    if (!tpcs)
    {
      //throw insight::Exception("tpc FO array not found in case!");
    }
    else
    {
      tpcs->evaluate(cm, executionPath(), results,
	"two-point correlation of velocity at different radii at x/H="+str(format("%f")%tpc_xlocs_[i])
	    );
    }
    
  }
  
  return results;
}




void ChannelInflow::applyCustomPreprocessing(OpenFOAMCase& cm, ProgressDisplayer& progress)
{
  
  setFields(cm, executionPath(),
            {
              "volVectorFieldValue U "+OFDictData::to_OF(vec3(sp().Ubulk_, 0, 0))
            },
            {}
            );
  
//   cm.get<TurbulentVelocityInletBC>(cycl_in_+"BC")->initInflowBC(executionPath(), p.getSubset("inflow"));
  
  OpenFOAMAnalysis::applyCustomPreprocessing(cm, progress);
}




void ChannelInflow::applyCustomOptions(OpenFOAMCase& cm, std::shared_ptr<OFdicts>& dicts)
{ 
//  Parameters p(parameters_);
  
  ChannelBase::applyCustomOptions(cm, dicts);
  
  OFDictData::dictFile& controlDict=dicts->lookupDict("system/controlDict");
  controlDict["endTime"] = sp().end_;
  
  if (const auto* sip =
       boost::get<Parameters::inflow_scalar_inflow_generator_type>(
              &p().inflow_scalar) )
  {
      if (!p().operation.wscalar)
      {
          throw insight::Exception("Inconsistent configuration: scalar disabled, but scalar inflow specified!");
      }
      
      OFDictData::dictFile& theta = dicts->lookupDict("0/theta");
      OFDictData::dict& thetaBC = theta.subDict("boundaryField").subDict(sp().cycl_in_);
      
      thetaBC["type"]="turbulentScalar";
      thetaBC["mean"]="linearProfile (0 0 0) (0 1 0) \""
          +dicts->insertAdditionalInputFile(sip->theta_mean).string()+"\"";

      thetaBC["variance"]="linearProfile (0 0 0) (0 1 0) \""
          +dicts->insertAdditionalInputFile(sip->theta_variance).string()+"\"";

      thetaBC["turbulentFlux"]="linearProfile (0 0 0) (0 1 0) \""
          +dicts->insertAdditionalInputFile(sip->theta_turbulentFlux).string()+"\"";
  }
}

addToAnalysisFactoryTable(ChannelInflow);

}
