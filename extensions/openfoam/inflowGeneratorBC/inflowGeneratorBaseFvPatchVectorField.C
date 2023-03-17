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

#include "inflowGeneratorBaseFvPatchVectorField.H" 

#include "transform.H"
#include "transformField.H"
#include "volFields.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "globalMeshData.H"
#include "globalIndex.H"

#include "base/vtktools.h"

#include <vector>
#include "boost/lexical_cast.hpp"

#include "uniof.h"

//using namespace std;
//using namespace boost;

namespace Foam 
{

defineTypeNameAndDebug(inflowGeneratorBaseFvPatchVectorField, 0);

vector arma_eigenValues(const symmTensor& r)
{
  arma::mat mr;
  mr 
    << r.xx() << r.xy() << r.xz() << arma::endr
    << r.xy() << r.yy() << r.yz() << arma::endr
    << r.xz() << r.yz() << r.zz() << arma::endr
    ;
  
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, mr);
  
  arma::uvec idx=arma::sort_index(eigval, /*"descend"*/ 1);

  return vector
  (
    eigval(idx(0)),
    eigval(idx(1)),
    eigval(idx(2))
  );
}



inflowInputDataField::inflowInputDataField(const pointField& x, const Time& time)
: x_(x), time_(time)
{}

inflowInputDataField::inflowInputDataField(const pointField& x, const Time& time, const dictionary& dict)
: x_(x), time_(time),
  Umean_(FieldDataProvider<vector>::New(dict.lookup("Umean"))),
  R_(FieldDataProvider<symmTensor>::New(dict.lookup("R"))),
  L_(FieldDataProvider<vector>::New(dict.lookup("L"))),
  c_(FieldDataProvider<scalar>::New(dict.lookup("c")))
{}

inflowInputDataField::inflowInputDataField(const inflowInputDataField& o)
: x_(o.x_), 
  time_(o.time_),
  Umean_(o.Umean_().clone()),
  R_(o.R_().clone()),
  L_(o.L_().clone()),
  c_(o.c_().clone())
{
}

const Field<point>& inflowInputDataField::x() const
{
  return x_;
}

const Field<vector>& inflowInputDataField::Umean() const 
{ 
    if ( (!Umeancur_.valid()) || (curUmeanAccessTimeIndex_ != time_.timeIndex()) )
    {
        curUmeanAccessTimeIndex_ = time_.timeIndex();
        Umeancur_ = Umean_()(time_.timeOutputValue(), x_); 
    }
    return Umeancur_();
}

tmp<Field<vector> > inflowInputDataField::flowDirection() const
{
    return Umean()/mag(Umean());
}

vector inflowInputDataField::flowDirection(label I) const
{
    return Umean()[I]/mag(Umean()[I]);
}

const Field<symmTensor>& inflowInputDataField::R() const 
{ 
    if ( (!Rcur_.valid()) || (curRAccessTimeIndex_ != time_.timeIndex()) )
    {
        curRAccessTimeIndex_ = time_.timeIndex();
        Rcur_ = R_()(time_.timeOutputValue(), x_); 
    }
    return Rcur_();
}

const Field<vector>& inflowInputDataField::L() const 
{ 
    if ( (!Lcur_.valid()) || (curLAccessTimeIndex_ != time_.timeIndex()) )
    {
        curLAccessTimeIndex_ = time_.timeIndex();
        Lcur_ = L_()(time_.timeOutputValue(), x_); 
    }
    return Lcur_();
}

const scalarField& inflowInputDataField::c() const 
{ 
    if ( (!ccur_.valid()) || (curcAccessTimeIndex_ != time_.timeIndex()) )
    {
        curcAccessTimeIndex_ = time_.timeIndex();
        ccur_ = c_()(time_.timeOutputValue(), x_); 
    }
    return ccur_();
}

void inflowInputDataField::write(Ostream& os) const
{
    Umean_().writeEntry("Umean", os);
    R_().writeEntry("R", os);
    L_().writeEntry("L", os);
    c_().writeEntry("c", os);
}

void inflowGeneratorBaseFvPatchVectorField::computeConditioningFactor(int writeInterval, int nsteps)
{
//   scalar relax=0.5;
      
  symmTensorField r_prescribed(R());
  
  // collect statistics
  vectorField u(size(), vector::zero);
  vectorField uMean(size(), vector::zero);
  symmTensorField uPrime2Mean(size(), symmTensor::zero);

  restartFluctuationProcess();

  scalarField ratio(size(), 1.0);

  for (int i=1; i<nsteps; i++)
  {
    ProcessStepInfo info;
    u = continueFluctuationProcess(0.0, &info, 1.0);

    scalar alpha = scalar(i - 1)/scalar(i);
    scalar beta = 1.0/scalar(i);

//     uPrime2Mean += sqr(uMean);
    uMean = alpha*uMean + beta*u;
//     N = alpha*N + beta*scalar(info.n_induced);
    uPrime2Mean = alpha*uPrime2Mean + beta*sqr(u) /*- sqr(uMean)*/; //uMean shoudl be zero

    Info<<"i="<<i<<": Averages: uMean="
        <<gSum(uMean*patch().magSf())/gSum(patch().magSf())
        <<" \t R^2="
        <<gSum(uPrime2Mean*patch().magSf())/gSum(patch().magSf())
        << endl;

    forAll(uPrime2Mean, i)
    {
      const scalar& r_cur = tr(uPrime2Mean[i]);

      ratio[i]=
        tr(r_prescribed[i])
          /
        r_cur;
  //       ratio[i].x()=min(max( ( r_prescribed[i].xx() / (mag(R1)+SMALL) ), mir), mar);
  //       ratio[i].y()=min(max( ( r_prescribed[i].yy() / (mag(R2)+SMALL) ), mir), mar);
  //       ratio[i].z()=min(max( ( r_prescribed[i].zz() / (mag(R3)+SMALL) ), mir), mar);
    }

    if (i%writeInterval==0)
      writeStateVisualization(i, u, &uMean, &uPrime2Mean, &ratio);
  }

//  writeStateVisualization(0, u, &uMean, &uPrime2Mean, &ratio);

}




inflowGeneratorBaseFvPatchVectorField::inflowGeneratorBaseFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    inflowInputDataField(patch().Cf(), this->db().time()),
    ranGen_(Pstream::myProcNo()),
    curTimeIndex_(-1),
    scaleToMassflow_(true)
{
}




inflowGeneratorBaseFvPatchVectorField::inflowGeneratorBaseFvPatchVectorField
(
    const inflowGeneratorBaseFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    inflowInputDataField(ptf),
    ranGen_(Pstream::myProcNo()),
    curTimeIndex_(ptf.curTimeIndex_),
    scaleToMassflow_(ptf.scaleToMassflow_)
{
}




inflowGeneratorBaseFvPatchVectorField::inflowGeneratorBaseFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    inflowInputDataField(patch().Cf(), this->db().time(), dict),
    ranGen_(Pstream::myProcNo()),
    curTimeIndex_(-1),
    scaleToMassflow_(dict.lookupOrDefault<Switch>("scaleToMassflow", true))
{  
}




inflowGeneratorBaseFvPatchVectorField::inflowGeneratorBaseFvPatchVectorField
(
    const inflowGeneratorBaseFvPatchVectorField& ptf
)
: fixedValueFvPatchField<vector>(ptf),
  inflowInputDataField(ptf),
  ranGen_(Pstream::myProcNo()),
  curTimeIndex_(ptf.curTimeIndex_),
  scaleToMassflow_(ptf.scaleToMassflow_)
{  
}




inflowGeneratorBaseFvPatchVectorField::inflowGeneratorBaseFvPatchVectorField
(
    const inflowGeneratorBaseFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
: fixedValueFvPatchField<vector>(ptf, iF),
  inflowInputDataField(ptf),
  ranGen_(Pstream::myProcNo()),
  curTimeIndex_(ptf.curTimeIndex_),
  scaleToMassflow_(ptf.scaleToMassflow_)
{
}


 
void inflowGeneratorBaseFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
//     rescale_.autoMap(m);
//     rescaleMean_.autoMap(m);
//     RplainMean_.autoMap(m);
}


void inflowGeneratorBaseFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);
    
//     const inflowGeneratorBaseFvPatchVectorField& tiptf = 
//       refCast<const inflowGeneratorBaseFvPatchVectorField >(ptf);

//     rescale_.rmap(tiptf.rescale_, addr);
//     rescaleMean_.rmap(tiptf.rescaleMean_, addr);
//     RplainMean_.rmap(tiptf.RplainMean_, addr);
}

const Field<vector>& inflowGeneratorBaseFvPatchVectorField::fluctuations() const
{ 
    if ( (!curfluctuations_.valid()) || (curflucAccessTimeIndex_ != this->db().time().timeIndex()) )
    {
        curflucAccessTimeIndex_ = this->db().time().timeIndex();
        curfluctuations_ = const_cast<inflowGeneratorBaseFvPatchVectorField&>(*this).continueFluctuationProcess(this->db().time().value()); 
    }
    return curfluctuations_();
}

const scalarField& inflowGeneratorBaseFvPatchVectorField::normalizedScalarFluctuations() const
{ 
    if ( (!cursfluctuations_.valid()) || (cursflucAccessTimeIndex_ != this->db().time().timeIndex()) )
    {
        cursflucAccessTimeIndex_ = this->db().time().timeIndex();
        fluctuations();
        cursfluctuations_ = const_cast<inflowGeneratorBaseFvPatchVectorField&>(*this).calculateScalarFluctuations(); 
    }
    return cursfluctuations_();
}

tmp<scalarField> inflowGeneratorBaseFvPatchVectorField::edgeLengths(bool maxL) const
{
  tmp<scalarField> res(new scalarField(size(), maxL?0.0:GREAT));
  scalarField& delta_edge = UNIOF_TMP_NONCONST(res);

  const polyPatch& ppatch = patch().patch();

  forAll(ppatch.edges(), ei)
  {
    const edge& e = ppatch.edges()[ei];
    scalar this_edge_len=e.mag(ppatch.localPoints());
    forAll(ppatch.edgeFaces()[ei], j)
    {
      label fi = ppatch.edgeFaces()[ei][j];
      if (maxL)
        delta_edge[fi] = max(delta_edge[fi], this_edge_len);
      else
        delta_edge[fi] = min(delta_edge[fi], this_edge_len);
    }
  }
  
  return res;
}

vector inflowGeneratorBaseFvPatchVectorField::randomTangentialDeflection(label fi)
{
  vector n=patch().Sf()[fi]; n/=mag(n);
  vector e1=n^vector(1,1,1);
  if (mag(e1)<SMALL) e1=n^vector(0,1,0);
  vector e2=n^e1;
  
  scalar dist=Foam::sqrt(patch().magSf()[fi]);

  return (0.5-ranGen_())*dist*e1 + (0.5-ranGen_())*dist*e2 ;
}


void inflowGeneratorBaseFvPatchVectorField::updateCoeffs()
{
  if (this->updated())
  {
      return;
  }

  if (curTimeIndex_ != this->db().time().timeIndex())
  {
    if (this->db().time().outputTime())
    {
      writeStateVisualization(0, fluctuations());
    }
    
    vectorField turbField( Umean() + fluctuations() );
   
    scalar rescale=1.0;

    if (scaleToMassflow_)
    { 
     scalar meanflux = gSum(Umean() & patch().Sf());
     scalar turbflux = gSum(turbField & patch().Sf());
     rescale = meanflux/turbflux;
     Info
        <<" Inflow generator ["<<patch().name()<<"]:"
        <<" mean="<<meanflux<<"/turb="<<turbflux
        <<", scaling turbulent fluctuations by "
        << rescale 
        << " to ensure constant flux across boundary."<<endl;
    }

    fixedValueFvPatchField<vector>::operator==( turbField * rescale );
    curTimeIndex_ = this->db().time().timeIndex();
    
  }

  fixedValueFvPatchField<vector>::updateCoeffs();
}


void inflowGeneratorBaseFvPatchVectorField::write(Ostream& os) const
{    
    inflowInputDataField::write(os);

    os.writeKeyword("scaleToMassflow") << scaleToMassflow_ << token::END_STATEMENT << nl;
        
    fixedValueFvPatchField<vector>::write(os);
}

}
