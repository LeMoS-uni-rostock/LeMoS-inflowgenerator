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
#include "inflowflatplatebl.h"

#include "base/factory.h"

#include "openfoam/blockmesh.h"
#include "openfoam/openfoamcase.h"
#include "openfoam/openfoamtools.h"

#include "openfoam/caseelements/boundaryconditions/turbulentvelocityinletbc.h"

#include "openfoam/caseelements/analysiscaseelements.h"

#include "refdata.h"

using namespace std;
using namespace boost;
using namespace boost::assign;

namespace insight 
{

addToAnalysisFactoryTable(InflowFlatPlateBL);

InflowFlatPlateBL::supplementedInputData::supplementedInputData(
    std::unique_ptr<Parameters> pPtr,
    const boost::filesystem::path& workDir,
    ProgressDisplayer& progress
    )
    : supplementedInputDataDerived<Parameters, FlatPlateBL::supplementedInputData>(
              std::move(pPtr),
              workDir,
              progress
              )
{
#warning Check initialization order and dependencies with base class

    if ( auto *utabulated
         = boost::get<Parameters::operation_type::Ushape_tabulated_type>(
             &p().operation.Ushape ) )
    {
      double uplus_inf=arma::as_scalar(arma::max(utabulated->uplus_vs_yplus.col(1)));
      cf_0_ = 2./pow(uplus_inf, 2);
      utau_0_=uinf_*sqrt(0.5*cf_0_);

      arma::mat uByUinf = join_rows
      (
        utabulated->uplus_vs_yplus.col(0) *p().fluid.nu/utau_0_,
        utabulated->uplus_vs_yplus.col(1) * utau_0_/uinf_
      );

      arma::mat delta123 = integrateDelta123( uByUinf );
      theta0_=delta123(1);

      delta99_0_=searchDelta99( uByUinf );
      Retheta0_= uinf_*theta0_/p().fluid.nu;

      struct Obj : public Objective1D
      {
        double Redelta2;
        virtual double operator()(double Rex) const
        {
      return FlatPlateBL::Redelta2(Rex)-Redelta2;
        }

      } o;
      o.Redelta2=Retheta0_;
      Rex_0_=nonlinearSolve1D(o, o.Redelta2, 1e6);
    }
    else
    {
      Retheta0_= p().FlatPlateBL::Parameters::operation.Retheta0;

      struct Obj : public Objective1D
      {
        double Redelta2;
        virtual double operator()(double Rex) const
        {
      return FlatPlateBL::Redelta2(Rex)-Redelta2;
        }

      } o;
      o.Redelta2=Retheta0_;
      Rex_0_=nonlinearSolve1D(o, o.Redelta2, 1e6);

      cf_0_=cf(Rex_0_);
      utau_0_=uinf_*sqrt(0.5*cf_0_);
      theta0_=Retheta0_*p().fluid.nu/uinf_;

      delta99_0_=Redelta99(Rex_0_)*p().fluid.nu/uinf_;
    }
}

InflowFlatPlateBL::InflowFlatPlateBL(
        const ParameterSet& ps,
        const boost::filesystem::path& exepath,
        ProgressDisplayer& progress)
: FlatPlateBL(
      std::make_unique<supplementedInputData>(
          std::make_unique<Parameters>(ps),
          exepath, progress
          ),
      exepath,
      "Inflow Generator Flat Plate",
      "")
{}



void InflowFlatPlateBL::createInflowBC(OpenFOAMCase& cm, const OFDictData::dict& boundaryDict) const
{

  TurbulentVelocityInletBC::Parameters::umean_type inflow_velocity;
  TurbulentVelocityInletBC::Parameters::turbulence_inflowGenerator_type::R_type inflow_R;
  TurbulentVelocityInletBC::Parameters::turbulence_inflowGenerator_type::L_type inflow_L;
  
  /*###################################################################################
    # Set velocities
    ################################################################################### */
  
  if ( auto *ublock
       = boost::get<Parameters::operation_type::Ushape_block_type>(
           &p().operation.Ushape) )
  {
      decltype(inflow_velocity)::fielddata_uniform_type umean_data;
      decltype(inflow_velocity)::fielddata_uniform_type::values_default_type inst;
      inst.time=0;
      inst.value=vec3(sp().uinf_, 0, 0);
      umean_data.values.push_back(inst);
      inflow_velocity.fielddata=umean_data;
  }
  else if ( auto *upl
       = boost::get<Parameters::operation_type::Ushape_powerlaw_type>(
                &p().operation.Ushape) )
  {
    boost::filesystem::path inlet_velocity_profile_tabfile(
                executionPath()/"inflow_velocity.dat");
    std::ofstream f(inlet_velocity_profile_tabfile.c_str());
    f<<" 0.0 0.0 0.0 0.0"<<endl;

    int n=20;
    for (int i=1; i<n; i++)
    {
        double eta=double(i)/double(n-1);
        
        double UByUinf = 
            //2.*eta - 2.*pow(eta,3) + pow(eta,4);
            pow(eta, 1./7.);
            
        f<<(sp().delta99_0_*eta)<<" "<<(sp().uinf_*UByUinf)<<" 0.0 0.0"<<endl;
    }

    if (sp().H_>sp().delta99_0_)
    {
        f<<sp().H_<<" "<<sp().uinf_<<" 0.0 0.0"<<endl;
    }
    decltype(inflow_velocity)::fielddata_linearProfile_type umean_data;

    // mean value profile
    umean_data.values.resize(1);
    umean_data.values[0].time=0;
    umean_data.values[0].profile->setOriginalFilePath(inlet_velocity_profile_tabfile);
    
    umean_data.p0=vec3(0,0,0);      
    umean_data.ep=vec3(0,1,0);    
      
    inflow_velocity.fielddata=umean_data;

  }
  else  if ( auto *utabulated
       = boost::get<Parameters::operation_type::Ushape_tabulated_type>(
                 &p().operation.Ushape) )
  {
    boost::filesystem::path inlet_velocity_profile_tabfile(
                executionPath()/"inflow_velocity.dat" );
    arma::mat u_vs_y = utabulated->uplus_vs_yplus;
    u_vs_y.col(0) 	*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    u_vs_y.cols(1,3) 	*= sp().utau_0_;
    u_vs_y.save(inlet_velocity_profile_tabfile.c_str(), arma::raw_ascii);
    
    decltype(inflow_velocity)::fielddata_linearProfile_type umean_data;

    // mean value profile
    umean_data.values.resize(1);
    umean_data.values[0].time=0;
    umean_data.values[0].profile->setOriginalFilePath(inlet_velocity_profile_tabfile);
    
    umean_data.p0=vec3(0,0,0);      
    umean_data.ep=vec3(0,1,0);
      
    inflow_velocity.fielddata=umean_data;    
  }

  
  
  
  /*###################################################################################
    # Set Reynolds stresses
    ################################################################################### */
  
  
  if ( auto *unifiso
       = boost::get<Parameters::operation_type::turb_uniformisotropic_type>(
           &p().operation.turb) )
  {
    double uprime2=pow( unifiso->intensity*sp().uinf_, 2);
    
    decltype(inflow_R)::fielddata_uniform_type::values_default_type inst;
    inst.time=0;
    inst.value 
      << uprime2 	<< arma::endr
      << 0 		<< arma::endr
      << 0 		<< arma::endr
      << uprime2 	<< arma::endr
      << 0 		<< arma::endr
      << uprime2 	<< arma::endr;

    decltype(inflow_R)::fielddata_uniform_type unifdata;
    unifdata.values.push_back(inst);
    
    inflow_R.fielddata=unifdata;
  }
  else  if ( auto *universalturb
       = boost::get<Parameters::operation_type::turb_universal_type>(
                 &p().operation.turb) )
  {
    arma::mat Ruu_vs_y=refdatalib.getProfile("Schlichting", "generic/RuuByUtauSqr_vs_yp");
    arma::mat Rvv_vs_y=refdatalib.getProfile("Schlichting", "generic/RvvByUtauSqr_vs_yp");
    arma::mat Rww_vs_y=refdatalib.getProfile("Schlichting", "generic/RwwByUtauSqr_vs_yp");
    Ruu_vs_y.col(0)*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    Rvv_vs_y.col(0)*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    Rww_vs_y.col(0)*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    Ruu_vs_y.col(1)*=pow(sp().utau_0_,2);
    Rvv_vs_y.col(1)*=pow(sp().utau_0_,2);
    Rww_vs_y.col(1)*=pow(sp().utau_0_,2);

    boost::filesystem::path inlet_R_profile_tabfile(executionPath()/"inflow_R.dat");
    arma::mat R_vs_y = arma::zeros(Ruu_vs_y.n_rows, 7);
    R_vs_y.col(0) = Ruu_vs_y.col(0);
    R_vs_y.col(1) = Ruu_vs_y.col(1);
    R_vs_y.col(4) = Interpolator(Rvv_vs_y.col(0), Rvv_vs_y.col(1))(R_vs_y.col(0));
    R_vs_y.col(6) = Interpolator(Rww_vs_y.col(0), Rww_vs_y.col(1))(R_vs_y.col(0));
    R_vs_y.save(inlet_R_profile_tabfile.c_str(), arma::raw_ascii);
    
    decltype(inflow_R)::fielddata_linearProfile_type data;
    // mean value profile
    data.values.resize(1);
    data.values[0].time=0;
    data.values[0].profile->setOriginalFilePath(inlet_R_profile_tabfile);

    
    data.p0=vec3(0,0,0);      
    data.ep=vec3(0,1,0);
      
    inflow_R.fielddata=data;    
  }
  else  if ( auto *tabturb
       = boost::get<Parameters::operation_type::turb_tabulated_type>(
                 &p().operation.turb) )
  {
    boost::filesystem::path inlet_R_profile_tabfile(
                executionPath()/"inflow_R.dat");
    arma::mat R_vs_y = tabturb->Rplus_vs_yplus;
    R_vs_y.col(0) 	*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    R_vs_y.cols(1,6) 	*= pow(sp().utau_0_, 2);
    R_vs_y.save(inlet_R_profile_tabfile.c_str(), arma::raw_ascii);
    
    decltype(inflow_R)::fielddata_linearProfile_type data;
    // mean value profile
    data.values.resize(1);
    data.values[0].time=0;
    data.values[0].profile->setOriginalFilePath(inlet_R_profile_tabfile);

    
    data.p0=vec3(0,0,0);      
    data.ep=vec3(0,1,0);
      
    inflow_R.fielddata=data;    
  }

  
  /*###################################################################################
    # Set length scales
    ################################################################################### */
  
  if ( auto *unifiso
       = boost::get<Parameters::operation_type::lengthscales_uniform_type>(
           &p().operation.lengthscales) )
  {    
    decltype(inflow_L)::fielddata_uniform_type::values_default_type inst;
    inst.time=0;
    inst.value 
      << unifiso->uniformLbydelta99(0)*sp().delta99_0_ 	<< arma::endr
      << unifiso->uniformLbydelta99(1)*sp().delta99_0_ 	<< arma::endr
      << unifiso->uniformLbydelta99(2)*sp().delta99_0_ 	<< arma::endr;

    decltype(inflow_L)::fielddata_uniform_type data;
    data.values.push_back(inst);
    
    inflow_L.fielddata=data;
  }
  else  if ( auto *tabL
       = boost::get<Parameters::operation_type::lengthscales_tabulated_type>(
                 &p().operation.lengthscales) )
  {
    boost::filesystem::path inlet_L_profile_tabfile(
                executionPath()/"inflow_L.dat" );
    arma::mat L_vs_y = tabL->Lbydelta99_vs_yplus;
    L_vs_y.col(0) 	*= p().FlatPlateBL::Parameters::fluid.nu / sp().utau_0_;
    L_vs_y.cols(1,3) 	*= sp().delta99_0_;
    L_vs_y.insert_cols(2,2);
    L_vs_y.insert_cols(5,1);
    L_vs_y.save(inlet_L_profile_tabfile.c_str(), arma::raw_ascii);
    
    decltype(inflow_L)::fielddata_linearProfile_type data;
    // mean value profile
    data.values.resize(1);
    data.values[0].time=0;
    data.values[0].profile->setOriginalFilePath(inlet_L_profile_tabfile);
    
    data.p0=vec3(0,0,0);      
    data.ep=vec3(0,1,0);
      
    inflow_L.fielddata=data;    
  }

  TurbulentVelocityInletBC::Parameters::turbulence_inflowGenerator_type inflturb;
  inflturb.R=inflow_R;
  inflturb.L=inflow_L;
  inflturb.volexcess=2;
  inflturb.uniformConvection=false;
  inflturb.type=int(p().operation.spottype);

  TurbulentVelocityInletBC::Parameters inflparams;
  inflparams.umean=inflow_velocity;
  inflparams.turbulence=inflturb;

  cm.insert(new TurbulentVelocityInletBC(cm, sp().in_, boundaryDict, inflparams));
  
}


}
