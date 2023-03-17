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


#include "turbulentStructure.H"
#include <armadillo>
#include "inflowGeneratorBaseFvPatchVectorField.H"

#include "base/exception.h"

namespace Foam
{

  
ESAnalyze::ESAnalyze(const symmTensor& L)
{
  arma::mat mL;
  mL 
    << L.xx() << L.xy() << L.xz() << arma::endr
    << L.xy() << L.yy() << L.yz() << arma::endr
    << L.xz() << L.yz() << L.zz() << arma::endr
    ;
  
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, mL);
  
  arma::uvec idx=arma::sort_index(eigval, /*"descend"*/ 1);
 
//    std::cout<<"eval="<<eigval<<"evec="<<eigvec<<"diag1="<< (eigvec.t()*mL*eigvec)<<"diag2="<< (eigvec*mL*eigvec.t()) <<std::endl; // only diag1 is right!
  
  e1_=vector(eigvec.col(idx(0))(0), eigvec.col(idx(0))(1), eigvec.col(idx(0))(2));
  e2_=vector(eigvec.col(idx(1))(0), eigvec.col(idx(1))(1), eigvec.col(idx(1))(2));
  e3_=vector(eigvec.col(idx(2))(0), eigvec.col(idx(2))(1), eigvec.col(idx(2))(2));
  
  e1_/=mag(e1_)+SMALL;
  e2_/=mag(e2_)+SMALL;
  e3_/=mag(e3_)+SMALL;
  
  es_=tensor
  (
   e1_ * eigval(idx(0)),
   e2_ * eigval(idx(1)),
   e3_ * eigval(idx(2))
  );
}



bool ESAnalyze::clip(scalar minL)
{
  vector L1=c1();
  vector L2=c2();
  vector L3=c3();

  scalar mL1=mag(L1);
  scalar mL2=mag(L2);
  scalar mL3=mag(L3);
  if (mL1<SMALL || mL2<SMALL || mL3<SMALL)
  {
   FatalErrorIn("ESAnalyze::clip")
   <<"zero principal length scale was prescribed! (L1="<<mL1<<", L2="<<mL2<<", L3="<<mL3<<")"
    <<abort(FatalError);
  }

  bool clipped=false;
  
  if (mL1<minL)
  {
    L1*=minL/mL1;
    clipped=true;
  }
  if (mL2<minL)
  {
    L2*=minL/mL2;
    clipped=true;
  }
  if (mL3<minL)
  {
    L3*=minL/mL3;
    clipped=true;
  }
  
  es_=tensor(L1, L2, L3);
  
  return clipped;
}

scalar ESAnalyze::Lalong(const vector& x) const
{
  return Lalong(x, c1(), c2(), c3());
}

scalar ESAnalyze::Lalong(const vector& x, const vector& L1, const vector& L2, const vector& L3)
{
  vector e=x/(mag(x)+SMALL);
  diagTensor Alpha(mag(L1), mag(L2), mag(L3));
  tensor Q( L1/Alpha.xx(), L2/Alpha.yy(), L3/Alpha.zz() );
  return (e & (Q.T()&Alpha&Q) & e);
}



void turbulentStructure::Parameters::initialize()
{
    e1_=vectorField(inputData().num(), vector::zero); 
    e2_=vectorField(inputData().num(), vector::zero); 
    e3_=vectorField(inputData().num(), vector::zero); 
    yw_=vectorField(inputData().num(), vector::zero); 
    Rp_=vectorField(inputData().num(), vector::zero);
    LalongN_=scalarField(inputData().num(), 0.0);
    
    vectorField flowDir( inputData().flowDirection() );
    
    for (label i=0; i<inputData().num(); i++)
    {
        ESAnalyze ea( inputData().R()[i] );
        Rp_[i][0]=mag(ea.c1());
        Rp_[i][1]=mag(ea.c2());
        Rp_[i][2]=mag(ea.c3());
        
        if (mag(Rp_[i][0])<SMALL || mag(Rp_[i][1])<SMALL || mag(Rp_[i][2])<SMALL)
        {
            e1_[i]=vector(1,0,0);
            e2_[i]=vector(0,1,0);
            e3_[i]=vector(0,0,1);
        }
        else
        {
            e1_[i]=ea.c1(); e1_[i]/=mag(e1_[i]);
            e2_[i]=ea.c2(); e2_[i]/=mag(e2_[i]);
            e3_[i]=ea.c3(); e3_[i]/=mag(e3_[i]);
        }  
        
        {
            vector l=inputData().L()[i];
            
            vector L1=l.x()*e1_[i],
            L2=l.y()*e2_[i], 
            L3=l.z()*e3_[i];
            
            vector x=flowDir[i];
            vector e=x/(mag(x)+SMALL);
            diagTensor Alpha(mag(L1), mag(L2), mag(L3));
            tensor Q( L1/Alpha.xx(), L2/Alpha.yy(), L3/Alpha.zz() );
            LalongN_[i] = (e & (Q.T()&Alpha&Q) & e);
        }
    }
}



turbulentStructure::Parameters::Parameters(const inflowInputDataField& ifp, const dictionary&)
: ifp_(ifp)
{
    initialize();
}


turbulentStructure::Parameters::~Parameters()
{}


void turbulentStructure::Parameters::autoMap
(
    const fvPatchFieldMapper&
)
{
}

//- Reverse map the given fvPatchField onto this fvPatchField
void turbulentStructure::Parameters::rmap
(
    const fvPatchField<vector>&,
    const labelList&
)
{
}

tensor turbulentStructure::Parameters::Lund(label iface) const
{
    symmTensor R=ifp_.R()[iface];
    
    tensor LT = tensor::zero;
    LT.xx()=sqrt(R.xx());
    LT.yx()=R.xy()/(SMALL+LT.xx());
    LT.yy()=sqrt(R.yy()-sqr(LT.yx()));
    LT.zx()=R.xz()/(SMALL+LT.xx());
    LT.zy()=(R.yz() - LT.yx()*LT.zx() )/(SMALL+LT.yy());
    LT.zz()=sqrt(R.zz() - sqr(LT.zx()) - sqr(LT.zy()));
    
    return LT;
}


void turbulentStructure::Parameters::write(Ostream&) const
{
}



int turbulentStructure::Parameters::nParamFields() const
{
  return 5;
}


string turbulentStructure::Parameters::paramFieldName(int i) const
{
  switch(i)
    {
    case 0: return "e1";
    case 1: return "e2";
    case 2: return "e3";
    case 3: return "Rp";
    case 4: return "yw";
    }
  return "";
}


turbulentStructure::Parameters::ParamField turbulentStructure::Parameters::paramField(int i) const
{
  switch(i)
    {
    case 0: return &e1_;
    case 1: return &e2_;
    case 2: return &e3_;
    case 3: return &Rp_;
    case 4: return &yw_;
    }

  throw insight::Exception("Invalid parameter field index!");
}








turbulentStructure::turbulentStructure(const Parameters& p, Istream& is)
: p_(p)
{
  is >> *this;
}



turbulentStructure::turbulentStructure
(
  const Parameters& p,
  BoostRandomGen&,
  const point& footPoint, 
  const vector& initialDelta,
  label creaface
)
: p_(p),
  creaFace_(creaface)
{  
  initialPositioning(footPoint, initialDelta);
}

turbulentStructure::turbulentStructure(const turbulentStructure& o)
: point(o),
  p_(o.p_),
  startPoint_(o.startPoint_),
  footPoint_(o.footPoint_),
  creaFace_(o.creaFace_)
{
}

turbulentStructure::~turbulentStructure()
{}


scalar turbulentStructure::travelledDistance() const
{
  return mag(*this - startPoint_);
}

vector turbulentStructure::motion() const
{
  return (*this - footPoint_);
}

scalar turbulentStructure::linearMotion() const
{
  vector dir = p_.inputData().flowDirection(creaFace_);
  return motion() & dir;
}
    
scalar turbulentStructure::scalar_fluctuation(const vector& x) const
{
    return fluctuation(x).x() / sqrt( p_.inputData().R()[creaFace_].xx() );
}


void turbulentStructure::randomize(BoostRandomGen& rand)
{
  epsilon_s_ = 2.0*(rand() - 0.5);
}

void turbulentStructure::operator=(const turbulentStructure& rhs)
{
  // Check for assignment to self
  if (this == &rhs)
  {
      FatalErrorIn("turbulentStructure::operator=(const turbulentStructure&)")
	  << "Attempted assignment to self"
	  << abort(FatalError);
  }

  point::operator=(rhs);
  startPoint_=rhs.startPoint_;
  footPoint_=rhs.footPoint_;
  creaFace_=rhs.creaFace_;

}

bool turbulentStructure::operator<(const turbulentStructure& o) const
{
  return this < &o;
}


Ostream& operator<<(Ostream& s, const turbulentStructure& ht)
{
  writeComponents(s, static_cast<const point&>(ht));
  writeComponents(s, ht.startPoint_);
  writeComponents(s, ht.footPoint_);
  s<<ht.creaFace_<<endl;
  return s;
}

Istream& operator>>(Istream& s, turbulentStructure& ht)
{
  readComponents(s, static_cast<point&>(ht));

  readComponents(s, ht.startPoint_);
  readComponents(s, ht.footPoint_);
  ht.creaFace_=readLabel(s);

  return s;
}


}
