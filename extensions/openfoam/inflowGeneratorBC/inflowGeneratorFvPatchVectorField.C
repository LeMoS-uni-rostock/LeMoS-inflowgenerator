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



#include "inflowGeneratorFvPatchVectorField.H"
#include "transform.H"
#include "transformField.H"
#include "volFields.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "PstreamCombineReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "vtkSurfaceWriter.H"

#include "base/vtktools.h"

#include <vector>
#include "boost/lexical_cast.hpp"
#include "boost/foreach.hpp"
#include "boost/iterator/counting_iterator.hpp"

#include "uniof.h"

//using namespace std;
//using namespace boost;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{




//template<class TurbulentStructure>
//void inflowGeneratorFvPatchVectorField<TurbulentStructure>::writeStateVisualization
//(
//  int i,
//  const vectorField& u,
//  const vectorField* uMean,
//  const symmTensorField* uPrime2Mean,
//  const scalarField* mult
//) const
//{
//  // #############################################################################
//  // #############################################################################
//  // #############################################################################
//  //
//  //                          VTK representation of vorton cloud
//  //
//  // #############################################################################
//  // #############################################################################
//  // #############################################################################

//  if ( time().outputTime() && lagrangianOutput_ )
//    {
//      {
//        IOobject po
//        (
//          "positions",
//          time().timeName(),
//          *this,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//        );
//        mkDir ( po.path() );
//        OFstream os ( po.objectPath() );
////       po.writeHeader(os);

//        // manually write header to ensure class name "Cloud"
//        os
//            << "FoamFile\n{\n"
//            << "    version     " << os.version() << ";\n"
//            << "    format      " << os.format() << ";\n"
//            << "    class       " << "Cloud" << ";\n"
//            << "    location    " << po.instance() /po.db().dbDir() /po.local() << ";\n"
//            << "    object      " << po.name() << ";\n"
//            << "}" << nl;


//        os<<label ( nvorton() ) <<nl<<token::BEGIN_LIST<<nl;
//        for (const auto& vl: vortons_)
//            for (const auto& vorton: vl)
//            {
//              os << point ( *vorton ) << " "<<label ( 0 ) <<nl;
//            }
//        os<<token::END_LIST<<nl;
//      }

//      {
//        IOobject ioo
//        (
//          "faceID",
//          time().timeName(),
//          *this,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//        );
//        IOField<Foam::scalar> faceid ( ioo, nvorton() );

//        label i=0;
//        for (const auto& vl: vortons_)
//            for (const auto& vorton: vl)
//            {
//              faceid[i++] = vorton->creaFace();
//            }
//        faceid.write();
//      }


//       {
//        IOobject ioo
//        (
//          "L",
//          time().timeName(),
//          *this,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//        );
//        IOField<Foam::tensor> L ( ioo, nvorton() );

//        label i=0;
//        for (const auto& vl: vortons_)
//            for (const auto& vorton: vl)
//            {
//              label j=vorton->creaFace();
//              vector l=this->L()[j];
//              const vector& L1 = l[0]*p_->e1()[j];
//              const vector& L2 = l[1]*p_->e2()[j];
//              const vector& L3 = l[2]*p_->e3()[j];
//              tensor tt ( L1/ ( l[0]+SMALL ), L2/ ( l[1]+SMALL ), L3/ ( l[2]+SMALL ) ); // vectors are rows!

//              L[i++]=tt.T() &diagTensor ( l[0],l[1],l[2] ) &tt;
//            }
//        L.write();
//      }

//      {
//        IOobject ioo
//        (
//          "Lp",
//          time().timeName(),
//          *this,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//        );
//        IOField<Foam::vector> Lp ( ioo, nvorton() );

//        label i=0;
//        for (const auto& vl: vortons_)
//            for (const auto& vorton: vl)
//            {
//              Lp[i++]=p_->inputData().L()[vorton->creaFace()];
//            }
//        Lp.write();
//      }

//      {
//        IOobject ioo
//        (
//          "Rp",
//          time().timeName(),
//          *this,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//        );
//        IOField<Foam::vector> Rp ( ioo, nvorton() );

//        label i=0;
//        for (const auto& vl: vortons_)
//            for (const auto& vorton: vl)
//            {
//              Rp[i++]=p_->Rp()[vorton->creaFace()];
//            }
//        Rp.write();
//      }

//    }

//  // #############################################################################
//  // #############################################################################
//  // #############################################################################
//  //
//  //                          VTK representation of patch
//  //
//  // #############################################################################
//  // #############################################################################
//  // #############################################################################
//  insight::vtk::vtkModel2d vtk_patch;
//  // set cells
//  const polyPatch& ppatch = patch().patch();
//  vtk_patch.setPoints
//  (
//    ppatch.localPoints().size(),
//    UNIOF_TMP_NONCONST(ppatch.localPoints().component(vector::X)).data(),
//    UNIOF_TMP_NONCONST(ppatch.localPoints().component(vector::Y)).data(),
//    UNIOF_TMP_NONCONST(ppatch.localPoints().component(vector::Z)).data()
//  );
//  for(label fi=0; fi<ppatch.size(); fi++)
//  {
//    const face& f = ppatch.localFaces()[fi];
//    vtk_patch.appendPolygon(f.size(), f.cdata());
//  }

//  for (int i=0; i<p_->nParamFields(); i++)
//      {
//          string name = p_->paramFieldName(i);
//          typename TurbulentStructure::Parameters::ParamField pf = p_->paramField(i);
//          if (const scalarField** sf = boost::get<const scalarField*>(&pf))
//              {
//                  vtk_patch.appendCellScalarField
//                  (
//                    name, (*sf)->cdata()
//                  );
//              }
//          else if (const vectorField** vf = boost::get<const vectorField*>(&pf))
//              {
//                  vtk_patch.appendCellVectorField
//                  (
//                    name,
//                    (*vf)->component(vector::X)().cdata(),
//                    (*vf)->component(vector::Y)().cdata(),
//                    (*vf)->component(vector::Z)().cdata()
//                  );
//              }
//      }
  
//  vtk_patch.appendCellVectorField
//  (
//    "u",
//    u.component(vector::X)().cdata(),
//    u.component(vector::Y)().cdata(),
//    u.component(vector::Z)().cdata()
//  );
//  if (uMean)
//  {
//    vtk_patch.appendCellVectorField
//    (
//      "uMean",
//      uMean->component(vector::X)().cdata(),
//      uMean->component(vector::Y)().cdata(),
//      uMean->component(vector::Z)().cdata()
//    );
//  }
//  if (uPrime2Mean)
//  {
//    vtk_patch.appendCellTensorField
//    (
//      "R",
//      uPrime2Mean->component(symmTensor::XX)().cdata(),
//      uPrime2Mean->component(symmTensor::XY)().cdata(),
//      uPrime2Mean->component(symmTensor::XZ)().cdata(),
//      uPrime2Mean->component(symmTensor::XY)().cdata(),
//      uPrime2Mean->component(symmTensor::YY)().cdata(),
//      uPrime2Mean->component(symmTensor::YZ)().cdata(),
//      uPrime2Mean->component(symmTensor::XZ)().cdata(),
//      uPrime2Mean->component(symmTensor::YZ)().cdata(),
//      uPrime2Mean->component(symmTensor::ZZ)().cdata()
//    );
//  }

//  {
//    vtk_patch.appendCellScalarField
//    (
//      "c",
//      ccorr()().cdata()
//    );
//  }
//  if (mult)
//  {
//    vtk_patch.appendCellScalarField
//    (
//      "mult",
//      (*mult).cdata()
//    );
//  }
//  {
//    symmTensorField Rp=R();
//    vtk_patch.appendCellTensorField
//    (
//      "UPrime2Mean",
//      Rp.component(symmTensor::XX)().cdata(),
//      Rp.component(symmTensor::XY)().cdata(),
//      Rp.component(symmTensor::XZ)().cdata(),
//      Rp.component(symmTensor::XY)().cdata(),
//      Rp.component(symmTensor::YY)().cdata(),
//      Rp.component(symmTensor::YZ)().cdata(),
//      Rp.component(symmTensor::XZ)().cdata(),
//      Rp.component(symmTensor::YZ)().cdata(),
//      Rp.component(symmTensor::ZZ)().cdata()
//    );
//  }
////   {
////     symmTensorField Rp=RplainMean_;
////     vtk_patch.appendCellTensorField
////     (
////       "Rmeas",
////       Rp.component(symmTensor::XX)().cdata(),
////       Rp.component(symmTensor::XY)().cdata(),
////       Rp.component(symmTensor::XZ)().cdata(),
////       Rp.component(symmTensor::XY)().cdata(),
////       Rp.component(symmTensor::YY)().cdata(),
////       Rp.component(symmTensor::YZ)().cdata(),
////       Rp.component(symmTensor::XZ)().cdata(),
////       Rp.component(symmTensor::YZ)().cdata(),
////       Rp.component(symmTensor::ZZ)().cdata()
////     );
////   }
//  {
    
//    IOobject oop
//    (
//      "patch_"+this->patch().name()+"_"+lexical_cast<string>(i)+".vtk",
//      this->db().time().timeName(),
//      this->db().time(),
//      IOobject::NO_READ,
//      IOobject::AUTO_WRITE
//    );
        
//    Info<<"Writing "<<oop.objectPath()<<endl;
//    std::ofstream f2(oop.objectPath().c_str());
//    vtk_patch.writeLegacyFile(f2);
//    f2.close();
//  }
//}







//template<class TurbulentStructure>
//typename inflowGeneratorFvPatchVectorField<TurbulentStructure>::VortonList
//inflowGeneratorFvPatchVectorField<TurbulentStructure>::filterVortons
//(
//    const inflowGeneratorFvPatchVectorField<TurbulentStructure>& ptf,
//    const fvPatchFieldMapper& mapper,
//    const VortonList& vlist
//) const
//{
//    VortonList newlist(size());

//    if (mapper.direct() && &mapper.directAddressing() && mapper.directAddressing().size())
//    {

//        const UNIOF_LABELULIST& addr=mapper.directAddressing();

//        HashTable<label, label, Foam::Hash<label> > inverseAddr;
//        forAll(addr, j)
//        {
//            inverseAddr.insert(addr[j], j);
//        }

//        if (debug)
//        {
//            Info
//                    <<"mapinfo: "<<mapper.direct()<<", "
//                    <<size()<<", "
//                    <<mapper.directAddressing().size()<<", "
//                    <<inverseAddr.size()
//                    <<endl;
//        }

//        label n_initial=0, n_kept=0;
//        for (const auto& vl: vlist)
//        {
//            for (const auto& ov: vl)
//            {
//                n_initial++;

//                if (inverseAddr.found(ov->creaFace()))
//                {
//                    n_kept++;
//                    label myf=inverseAddr[ov->creaFace()];

//                    TurbulentStructurePtr nv(new TurbulentStructure(*ov));
//                    nv->setCreaFace(myf);
//                    newlist[myf].push_back(nv);
//                }
//            }
//        }

//        if (debug)
//        {
//            Pout<<"Kept " << n_kept << " out of " << n_initial << " vortons." << endl;
//        }
//    }
//    else
//    {
//        // create new vortons, if meshes don't match
//    }

//    return newlist;
//}






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




template<class TurbulentStructure>
inflowGeneratorFvPatchVectorField<TurbulentStructure>::inflowGeneratorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),

    inflowGenerator_(new structureBasedInflowGenerator(
        this->db().time().value(),
        uniformField<vector>(vector(1,0,0)),
        uniformField<symmTensor>(symmTensor(1,0,0,1,0,1)),
        uniformField<vector>(vector(1,1,1)),
        uniformField<scalar>(1.),
        uniformField<scalar>(1.),
        &this->patch().boundaryMesh().mesh(),
        this->patch().patch()
        ) ),

    cloud(this->db(), this->patch().name()+"VortonCloud"),

    lagrangianOutput_(false),
    curTimeIndex_(-1),
    scaleToMassflow_(true)
{}




// this is called during decomposePar
template<class TurbulentStructure>
inflowGeneratorFvPatchVectorField<TurbulentStructure>::inflowGeneratorFvPatchVectorField
(
    const inflowGeneratorFvPatchVectorField<TurbulentStructure>& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),

    inflowGenerator_( new structureBasedInflowGenerator(
        this->db().time().value(),
        ptf.inflowGenerator_->UMeanInput(),
        ptf.inflowGenerator_->RInput(),
        ptf.inflowGenerator_->LInput(),
        ptf.inflowGenerator_->calibrationFactorInput(),
        ptf.inflowGenerator_->vortonDensityInput(),
        &this->patch().boundaryMesh().mesh(),
        this->patch().patch()
        )
    ),

    cloud(this->db(), this->patch().name()+"VortonCloud"),

    lagrangianOutput_(ptf.lagrangianOutput_),
    curTimeIndex_(-1),
    scaleToMassflow_(ptf.scaleToMassflow_)
{
    // upon decomposePar:
    // pick out vortons on current processor mesh:
    if (mapper.direct() && (&mapper.directAddressing()) && mapper.directAddressing().size())
    {

        const UNIOF_LABELULIST& addr = mapper.directAddressing();

        HashTable<label, label, Foam::Hash<label> > inverseAddr;
        forAll(addr, j)
        {
            inverseAddr.insert(addr[j], j);
        }

        label nInitial=0, nKept=0;
        for (label otherFaceI=0; otherFaceI<ptf.fvPatchField::size(); ++otherFaceI)
        {
            if (inverseAddr.found(otherFaceI))
            {
                nKept++;
                label myFaceI=inverseAddr[otherFaceI];
                this->inflowGenerator_->queues_.set(
                            myFaceI,
                            new vortonQueue<TurbulentStructure>(
                                ptf.inflowGenerator_->queues_[otherFaceI]
                                ));
            }
        }

        reduce(nInitial, sumOp<label>());
        reduce(nKept, sumOp<label>());

        if (debug)
        {
            Pout<<"Kept " << nKept << " out of " << nInitial << " vortons." << endl;
        }
    }
    else
    {
        // create new vortons, if meshes don't match
    }
}




template<class TurbulentStructure>
inflowGeneratorFvPatchVectorField<TurbulentStructure>::inflowGeneratorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedValueFvPatchVectorField(p, iF, dict),

    inflowGenerator_(
        dict.found("ySource") ?
        new structureBasedInflowGenerator(
            this->db().time().value(),
            FieldDataProvider<vector>::New(dict.lookup("UMeanSource")),
            FieldDataProvider<symmTensor>::New(dict.lookup("RSource")),
            typename structureBasedInflowGenerator::WallDistanceAndLengthScaleInput::input_type(
             FieldDataProvider<scalar>::New(dict.lookup("LSource")),
             FieldDataProvider<scalar>::New(dict.lookup("ySource"))
            ),
            FieldDataProvider<scalar>::New(dict.lookup("calibrationFactorSource")),
            FieldDataProvider<scalar>::New(dict.lookup("vortonDensitySource")),
            &this->patch().boundaryMesh().mesh(),
            this->patch().patch()
        ) :
        new structureBasedInflowGenerator(
            this->db().time().value(),
            FieldDataProvider<vector>::New(dict.lookup("UMeanSource")),
            FieldDataProvider<symmTensor>::New(dict.lookup("RSource")),
            FieldDataProvider<vector>::New(dict.lookup("LSource")),
            FieldDataProvider<scalar>::New(dict.lookup("calibrationFactorSource")),
            FieldDataProvider<scalar>::New(dict.lookup("vortonDensitySource")),
            &this->patch().boundaryMesh().mesh(),
            this->patch().patch()
        )

    ),

    cloud(this->db(), patch().name()+"VortonCloud"),

    lagrangianOutput_(dict.lookupOrDefault("lagrangianOutput", false)),
    curTimeIndex_(-1),
    scaleToMassflow_(dict.lookupOrDefault("scaleToMassflow", true))
{
    if (dict.found("vortons"))
    {
        this->inflowGenerator_->readVortons(dict.lookup("vortons"));
    }
}




template<class TurbulentStructure>
inflowGeneratorFvPatchVectorField<TurbulentStructure>::inflowGeneratorFvPatchVectorField
(
    const inflowGeneratorFvPatchVectorField& ptf
)
: fixedValueFvPatchVectorField(ptf),

  inflowGenerator_(new structureBasedInflowGenerator(
      static_cast<const structureBasedInflowGenerator&>(ptf.inflowGenerator_()) // typecast to enforce call of copy constructor
      ) ),

  cloud(this->db(), patch().name()+"VortonCloud"),

  lagrangianOutput_(ptf.lagrangianOutput_),
  curTimeIndex_(ptf.curTimeIndex_),
  scaleToMassflow_(ptf.scaleToMassflow_)
{}




template<class TurbulentStructure>
inflowGeneratorFvPatchVectorField<TurbulentStructure>::inflowGeneratorFvPatchVectorField
(
    const inflowGeneratorFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
: fixedValueFvPatchVectorField(ptf, iF),

  inflowGenerator_(new structureBasedInflowGenerator(
      static_cast<const structureBasedInflowGenerator&>(ptf.inflowGenerator_()) // typecast to enforce call of copy constructor
      ) ),

  cloud(this->db(), patch().name()+"VortonCloud"),

  lagrangianOutput_(ptf.lagrangianOutput_),
  curTimeIndex_(ptf.curTimeIndex_),
  scaleToMassflow_(ptf.scaleToMassflow_)
{}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class TurbulentStructure>
void inflowGeneratorFvPatchVectorField<TurbulentStructure>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


template<class TurbulentStructure>
void inflowGeneratorFvPatchVectorField<TurbulentStructure>::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
    
    const inflowGeneratorFvPatchVectorField<TurbulentStructure>& tiptf =
      refCast<const inflowGeneratorFvPatchVectorField<TurbulentStructure> >(ptf);

    // insert vortons
    for (label otherFaceI=0; otherFaceI<tiptf.fvPatchField::size(); ++otherFaceI)
    {
        label myFaceI = addr[otherFaceI];

        inflowGenerator_->queues_.set(
                myFaceI,
                new vortonQueue<TurbulentStructure>(
                    tiptf.inflowGenerator_->queues_[otherFaceI]
                )
        );
    } 
}







template<class TurbulentStructure>
void inflowGeneratorFvPatchVectorField<TurbulentStructure>::updateCoeffs()
{
  if (this->updated())
  {
      return;
  }

  if (curTimeIndex_ != this->db().time().timeIndex())
  {

    this->inflowGenerator_->advance( this->db().time().value() );

//    if (this->db().time().outputTime())
//    {
//      writeStateVisualization(0, fluctuations());
//    }

    vectorField turbField( this->inflowGenerator_->turbulentField() );

    scalar rescale=1.0;

    if (scaleToMassflow_)
    {
     scalar meanflux = gSum(this->inflowGenerator_->UMean(this->db().time().value()) & this->patch().Sf());
     scalar turbflux = gSum(turbField & this->patch().Sf());
     rescale = meanflux/turbflux;
     Info
        <<" Inflow generator ["<<this->patch().name()<<"]:"
        <<" mean="<<meanflux<<"/turb="<<turbflux
        <<", scaling turbulent fluctuations by "
        << rescale
        << " to ensure constant flux across boundary."<<endl;
    }

    this->fixedValueFvPatchField<vector>::operator==( turbField * rescale );
    curTimeIndex_ = this->db().time().timeIndex();

  }

  this->fixedValueFvPatchField<vector>::updateCoeffs();
}


template<class TurbulentStructure>
void inflowGeneratorFvPatchVectorField<TurbulentStructure>::write(Ostream& os) const
{
    os.writeKeyword("lagrangianOutput") << lagrangianOutput_ <<token::END_STATEMENT<<nl;
    os.writeKeyword("scaleToMassflow") << scaleToMassflow_ <<token::END_STATEMENT<<nl;

    this->inflowGenerator_->UMeanInput().writeEntry("UMeanSource", os);
    this->inflowGenerator_->RInput().writeEntry("RSource", os);

    auto Lspec = this->inflowGenerator_->LInput();
    if (const auto * Ldirect =
            boost::get<
            const typename InflowGeneratorBase<face,SubList,const pointField&>
             ::DirectLengthScaleInput::input_type&>(&Lspec) )
    {
        Ldirect->writeEntry("LSource", os);
    }
    else
    {
        FatalErrorIn("void inflowGeneratorFvPatchVectorField<TurbulentStructure>::write(Ostream& os) const")
                << "Unsupported length scale specification!" << abort(FatalError);
    }
    this->inflowGenerator_->calibrationFactorInput().writeEntry("calibrationFactorSource", os);
    this->inflowGenerator_->vortonDensityInput().writeEntry("vortonDensitySource", os);

    this->inflowGenerator_->writeVortons( os.writeKeyword("vortons") << token::SPACE ) << token::END_STATEMENT << nl;

    this->fixedValueFvPatchVectorField::write(os);
}




} // End namespace Foam

// ************************************************************************* //
