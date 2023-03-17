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

#include "modalTurbulenceFvPatchVectorField.H"

#include "transform.H"
#include "transformField.H"
#include "volFields.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "wallDist.H"

#include "base/vtktools.h"

#include <vector>
#include "boost/lexical_cast.hpp"
#include "boost/foreach.hpp"

using namespace std;
//using namespace boost;

namespace Foam 
{


void modalTurbulenceFvPatchVectorField::writeStateVisualization
(
  int i,
  const vectorField& u,
  const vectorField* uMean,
  const symmTensorField* uPrime2Mean,
  const scalarField *mult
) const
{
//   insight::vtk::vtkModel vtk_vortons;
//   
//   vtk_vortons.setPoints(
//     vortons_.size(),
//     vortons_.component(vector::X)().data(),
//     vortons_.component(vector::Y)().data(),
//     vortons_.component(vector::Z)().data()
//     );
//   
//   for (label k=0; k<3; k++)
//   {
//     std::vector<double> Lx, Ly, Lz;
//     forAll(vortons_, j)
//     {
//       const vector& L = vortons_[j].L(k);
//       Lx.push_back(L.x()); Ly.push_back(L.y()); Lz.push_back(L.z());
//     }
//     vtk_vortons.appendPointVectorField("L"+lexical_cast<string>(k), Lx.data(), Ly.data(), Lz.data());
//   }
//   
//   insight::vtk::vtkModel2d vtk_patch;
//   // set cells
//   const polyPatch& ppatch = patch().patch();
//   vtk_patch.setPoints
//   (
//     ppatch.localPoints().size(), 
//     ppatch.localPoints().component(vector::X)().data(),
//     ppatch.localPoints().component(vector::Y)().data(),
//     ppatch.localPoints().component(vector::Z)().data()
//   );
//   for(label fi=0; fi<ppatch.size(); fi++)
//   {
//     const face& f = ppatch.localFaces()[fi];
//     vtk_patch.appendPolygon(f.size(), f.cdata());
//   }
//   
//   vtk_patch.appendCellVectorField
//   (
//     "u", 
//     u.component(vector::X)().cdata(),
//     u.component(vector::Y)().cdata(),
//     u.component(vector::Z)().cdata()
//   );
//   if (uMean)
//   {
//     vtk_patch.appendCellVectorField
//     (
//       "uMean", 
//       uMean->component(vector::X)().cdata(),
//       uMean->component(vector::Y)().cdata(),
//       uMean->component(vector::Z)().cdata()
//     );
//   }
//   if (uPrime2Mean)
//   {
//     vtk_patch.appendCellTensorField
//     (
//       "R", 
//       uPrime2Mean->component(symmTensor::XX)().cdata(),
//       uPrime2Mean->component(symmTensor::XY)().cdata(),
//       uPrime2Mean->component(symmTensor::XZ)().cdata(),
//       uPrime2Mean->component(symmTensor::XY)().cdata(),
//       uPrime2Mean->component(symmTensor::YY)().cdata(),
//       uPrime2Mean->component(symmTensor::YZ)().cdata(),
//       uPrime2Mean->component(symmTensor::XZ)().cdata(),
//       uPrime2Mean->component(symmTensor::YZ)().cdata(),
//       uPrime2Mean->component(symmTensor::ZZ)().cdata()
//     );
//   }
//   {
//     
//     IOobject oo
//     (
//       "vortons_"+this->patch().name()+"_"+lexical_cast<string>(i)+".vtk",
//       this->db().time().timeName(),
//       this->db().time(),
//       IOobject::NO_READ,
//       IOobject::AUTO_WRITE
//     );
//     IOobject oop
//     (
//       "patch_"+this->patch().name()+"_"+lexical_cast<string>(i)+".vtk",
//       this->db().time().timeName(),
//       this->db().time(),
//       IOobject::NO_READ,
//       IOobject::AUTO_WRITE
//     );
//     mkDir(oo.path());
//     
//     Info<<"Writing "<<oo.objectPath()<<endl;
//     std::ofstream f(oo.objectPath().c_str());
//     vtk_vortons.writeLegacyFile(f);
//     f.close();
//     
//     Info<<"Writing "<<oop.objectPath()<<endl;
//     std::ofstream f2(oop.objectPath().c_str());
//     vtk_patch.writeLegacyFile(f2);
//     f2.close();
//   }
}

tmp<scalarField> modalTurbulenceFvPatchVectorField::calculateScalarFluctuations()
{
    FatalErrorIn("modalTurbulenceFvPatchVectorField::calculateScalarFluctuations()")
     << "Not implemented!"<<abort(FatalError);
    return tmp<scalarField>();
}

bool modalTurbulenceFvPatchVectorField::mode::operator!=(const mode& other) const
{
  return other.k!=k;
}

modalTurbulenceFvPatchVectorField::modalTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inflowGeneratorBaseFvPatchVectorField(p, iF),
    modes_(),
    tau_(-1)
{
}

modalTurbulenceFvPatchVectorField::modalTurbulenceFvPatchVectorField
(
    const modalTurbulenceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inflowGeneratorBaseFvPatchVectorField(ptf, p, iF, mapper),
    modes_(),
    tau_(-1)
{
}

modalTurbulenceFvPatchVectorField::modalTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inflowGeneratorBaseFvPatchVectorField(p, iF, dict),
    modes_(),
    tau_(-1)
{
    if (dict.found("modes"))
    {
      tau_=readScalar(dict.lookup("tau"));
      modes_=List<mode>(dict.lookup("modes"));
    }  
}

modalTurbulenceFvPatchVectorField::modalTurbulenceFvPatchVectorField
(
    const modalTurbulenceFvPatchVectorField& ptf
)
: inflowGeneratorBaseFvPatchVectorField(ptf),
  modes_(ptf.modes_),
  tau_(ptf.tau_)
{}

modalTurbulenceFvPatchVectorField::modalTurbulenceFvPatchVectorField
(
    const modalTurbulenceFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
: inflowGeneratorBaseFvPatchVectorField(ptf, iF),
  modes_(ptf.modes_),
  tau_(ptf.tau_)
{}


void modalTurbulenceFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
  inflowGeneratorBaseFvPatchVectorField::autoMap(m);
  forAll(modes_, i)
  {
    modes_[i].q.autoMap(m);
  }
}


void modalTurbulenceFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
  inflowGeneratorBaseFvPatchVectorField::rmap(ptf, addr);
  
  const modalTurbulenceFvPatchVectorField& tiptf = 
    refCast<const modalTurbulenceFvPatchVectorField >(ptf);
  forAll(modes_, i)
  {
    modes_[i].q.rmap(tiptf.modes_[i].q, addr);
  }
}

#define OFLD(f) \
 Pout<<#f<<": min="<<gMin(f)<<", max="<<gMax(f)<<", avg="<<gAverage(f)<<endl;

void modalTurbulenceFvPatchVectorField::createModes()
{
    const polyPatch& ppatch = patch().patch();

    wallDist ywall(patch().boundaryMesh().mesh());

    const scalarField& dw=UNIOF_WALLDIST_Y(ywall).boundaryField()[patch().index()];
//   OFLD(dw);
    scalarField delta_x( 2.*(patch().Cf() - ppatch.faceCellCentres()) & patch().nf() );
    scalarField delta_max_edge(edgeLengths(true));
    scalarField delta_max( max(delta_max_edge, delta_x) );
    scalarField lcut( 2.*min
                     (
                         max(delta_max_edge, 0.3*delta_max) + 0.1*dw,
                         delta_max
                     ));

//   OFLD(delta_max);
//   OFLD(lcut);

    scalar ce=3., ctau=2.;
    scalar alpha=0.025;
#warning viscosity is hard-coded!
    scalar nu=1e-5;

    scalarField lt( mag(L()) );
//   OFLD(lt);
    scalarField le( Foam::min(2.*dw, ce*lt) );
//   OFLD(le);
    scalarField ke(2.*M_PI/le);

    scalar lemax=gMax(le);
    scalar kmin = 2.*M_PI/lemax;
    scalarField kcut(2.*M_PI/lcut);
    scalar kcutmax=gMax(kcut);
    scalarField keta(2.*M_PI/ (pow(nu,3)*lt/pow(mag(Umean())+SMALL, 3)));

    tau_=ctau*lemax/ ( gSum(-Umean()&patch().Sf()) / gSum(patch().magSf()) );

//   OFLD(lt);
//   OFLD(le);
//   OFLD(ke);
//   OFLD(kcut);
//   OFLD(keta);
//   Info<<tau_<<endl;
//
    {
        std::vector<scalar> ks;
        scalar kn=0;

        const label nmax=10000;
        while ( (kn=kmin*::pow(1.+alpha, ks.size())) < 1.5*kcutmax )
        {
//       Info<<kn<<endl;
            ks.push_back(kn);
            if (ks.size()>=nmax)
            {
                FatalErrorIn("modalTurbulenceFvPatchVectorField")
                        << "maximum number of modes reached ("<<nmax<<"): kn="<<kn  <<", limit kn="<< 1.5*kcutmax
                        <<abort(FatalError);
            }
        }

        modes_.resize(ks.size());

        // var_nor: Generator for random numbers with normal distribution
        boost::mt19937 rng;
        boost::normal_distribution<> nd(2.*M_PI, 2.*M_PI);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

        scalarField qsum(size(), 0.0);

        for (size_t i=0; i<ks.size(); i++)
        {
            mode& m = modes_[i];
            m.k=ks[i];
            const scalar& k=m.k;

            // compute random numbers on master processor and distribute afterwards
            // => need to be consistent over all processors
            if (Pstream::master())
            {
                m.d = 2.*vector(ranGen_()-0.5, ranGen_()-0.5, ranGen_()-0.5);
                m.d/=SMALL+mag(m.d);

                m.phi=ranGen_()*2.*M_PI;

                m.s=var_nor();

                m.sigma=vector(ranGen_()-0.5, ranGen_()-0.5, 0);
                m.sigma.z() = - m.sigma&m.d / m.d.z(); // ensure orthogonality of sigma and d
                m.sigma /= SMALL+mag(m.sigma);
            }
            else
            {
                m.d=vector::zero;
                m.phi=0.0;
                m.s=0.0;
                m.sigma=vector::zero;
            }

            reduce(m.d, sumOp<vector>());
            reduce(m.phi, sumOp<scalar>());
            reduce(m.s, sumOp<scalar>());
            reduce(m.sigma, sumOp<vector>());

            scalarField feta(exp(-pow(12.*k/keta,2)));
            scalarField fcut(exp(-pow(4*max(k-0.9*kcut,0.0)/kcut, 3)));
            scalarField Ek( feta*fcut*
                            pow(k/ke, 4)
                            /
                            pow(1.+2.4*pow(k/ke, 2), 17./6.));

            scalar delta_k=k;
            if (i>0) delta_k-=ks[i-1];

            m.q=Ek*delta_k;
            qsum+=m.q;

//       Pout<<"mode #"<<label(i)<<" : "<<endl;
//       Pout<<m.d<<", "<<m.phi<<", "<<m.s<<", "<<m.sigma<<endl;
//       OFLD(m.q);
        }
        forAll(modes_, i)
        {
            modes_[i].q/=max(SMALL, qsum);
        }
        Info<<"Created n="<<modes_.size()<<" modes."<<endl;
        
        Rnorm_=symmTensorField(size(), symmTensor::zero);
        
        forAll(modes_, m)
        {
            vectorField cm( Foam::sqrt(modes_[m].q)*modes_[m].sigma );
            Rnorm_ += 0.5*symm(cm * cm);
        }
        Rnorm_*=6;
        Pout<<"Rnorm = "<<Rnorm_<<endl;
    }
}


tmp<vectorField> modalTurbulenceFvPatchVectorField::continueFluctuationProcess(scalar t, ProcessStepInfo *info, double multidensity)
{

  if (modes_.size()==0) createModes();
  
  tmp<vectorField> tfluctuations(new vectorField(size(), vector::zero));
  vectorField& fluctuations = UNIOF_TMP_NONCONST(tfluctuations);
  
  
  forAll(modes_, mi)
  {
    const mode& m = modes_[mi];
    fluctuations += sqrt(m.q) * (m.sigma * cos( ((m.k*m.d)&patch().Cf()) + m.phi + m.s*t/tau_));
  }
  fluctuations *= 2.*sqrt(3./2.);

  if (info) 
  {
//     info->n_removed=n_removed;
//     info->n_generated=n_generated;
//     info->n_total=n_total;
  }
  
  // cholesky transformation
  symmTensorField Rtarg(R());
  tensorField LT(size(), tensor::zero);
  
  scalarField R11(Rtarg.component(symmTensor::XX));
  scalarField R12(Rtarg.component(symmTensor::XY));
  scalarField R13(Rtarg.component(symmTensor::XZ));
  scalarField R22(Rtarg.component(symmTensor::YY));
  scalarField R23(Rtarg.component(symmTensor::YZ));
  scalarField R33(Rtarg.component(symmTensor::ZZ));
  
  scalarField rho11(Rnorm_.component(symmTensor::XX));
  scalarField rho12(Rnorm_.component(symmTensor::XY));
  scalarField rho13(Rnorm_.component(symmTensor::XZ));
  scalarField rho22(Rnorm_.component(symmTensor::YY));
  scalarField rho23(Rnorm_.component(symmTensor::YZ));
  scalarField rho33(Rnorm_.component(symmTensor::ZZ));
  
  scalarField a11( sqrt(R11/rho11) );
  scalarField a22( sqrt( (R22-sqr(R12)/R11)/(rho22-sqr(rho12)/rho11) ) );
  scalarField a21( R12/sqrt(R11*rho11) - a22*rho12/rho11 );
  
  scalarField denom( sqr(rho13)*rho22 -2.*rho12*rho13*rho23 +rho11*sqr(rho23) +sqr(rho12)*rho33 -rho11*rho22*rho33 );
  
  scalarField a33
  (  sqrt(
        (
            ( R33*sqr(a11*a22)*(sqr(rho12)-rho22*rho11) +sqr(R13)*(sqr(a21)*rho11 +sqr(a22)*rho22 +2*a21*a22*rho12) ) /denom              
        + (-2.*R13*R23*a11*(a22*rho12+a21*rho11) +sqr(R23)*sqr(a11)*rho11) /denom
        ) / sqr(a11*a22)
     )
  );
  
  scalarField a32( (a33*(rho13*rho12-rho23*rho11) +R23*rho11/a22 -(a21*rho11+a22*rho12)*R13/a11/a22)/(rho22*rho11-sqr(rho12)) );
  scalarField a31( (a33*(rho23*rho12-rho13*rho22) -R23*rho12/a22 +(a21*rho12+a22*rho22)*R13/a11/a22)/(rho22*rho11-sqr(rho12)) );
  
  LT.replace(tensor::XX, a11);
  LT.replace(tensor::YY, a22);
  LT.replace(tensor::ZZ, a33);
  LT.replace(tensor::YX, a21);
  LT.replace(tensor::ZX, a31);
  LT.replace(tensor::ZY, a32);
  
  if (debug>=3)
  {
   Pout<<" fluctuations: min/max/avg = "<<min(fluctuations)<<" / "<<max(fluctuations) << " / "<<average(fluctuations)<<endl;
   forAll(fluctuations, j)
    if (mag(fluctuations[j])>1e3) Pout<<j<<": "<<tfluctuations<<endl;
  }

  return ( LT & tfluctuations );
}

void modalTurbulenceFvPatchVectorField::restartFluctuationProcess()
{
}


void modalTurbulenceFvPatchVectorField::write(Ostream& os) const
{
  if (modes_.size()>0)
  {
    os.writeKeyword("tau") << tau_ << token::END_STATEMENT <<nl;
    os.writeKeyword("modes") << modes_ << token::END_STATEMENT <<nl;
  }
    
  inflowGeneratorBaseFvPatchVectorField::write(os);
}

Ostream& operator<<(Ostream& os, const modalTurbulenceFvPatchVectorField::mode& m)
{
  dictionary d;
  d.add("k", m.k);
  d.add("phi", m.phi);
  d.add("s", m.s);
  d.add("d", m.d);
  d.add("sigma", m.sigma);
  d.add("q", m.q);
  os << d;
  return os;
}

Istream& operator>>(Istream& os, modalTurbulenceFvPatchVectorField::mode& m)
{
  dictionary d(os);
  m.k=readScalar(d.lookup("k"));
  m.phi=readScalar(d.lookup("phi"));
  m.s=readScalar(d.lookup("s"));
  m.d=vector(d.lookup("d")); 
  m.sigma=vector(d.lookup("sigma"));
  m.q=scalarField(d.lookup("q"));
  return os;
}

makePatchTypeField
(
    fvPatchVectorField,
    modalTurbulenceFvPatchVectorField
);

}
