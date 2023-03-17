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

#ifndef turbulentStructure_H
#define turbulentStructure_H

#include "fvCFD.H"
#include "boostRandomGen.H"
#include "IOdictionary.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <boost/graph/graph_concepts.hpp>
#include <boost/variant.hpp>

namespace Foam
{

class inflowGeneratorBaseFvPatchVectorField;
class inflowInputDataField;
class turbulentStructure;

Ostream& operator<<(Ostream&, const turbulentStructure&);
Istream& operator>>(Istream&, turbulentStructure&);


class ESAnalyze
{
  tensor es_;
  vector e1_, e2_, e3_;
public:
  ESAnalyze(const symmTensor& t);
  bool clip(scalar minL);
  inline vector e1() const { return e1_; }
  inline vector e2() const { return e2_; }
  inline vector e3() const { return e3_; }
  inline vector c1() const { return es_.x(); }
  inline vector c2() const { return es_.y(); }
  inline vector c3() const { return es_.z(); }
  inline tensor Leig() const { return es_; }
  scalar Lalong(const vector& x) const;

  static scalar Lalong(const vector& x, const vector& L1, const vector& L2, const vector& L3);  
//   static tensor eigenSystem(const symmTensor& L);
};



/**
 * base class for turbulent structures.
 * derived from point => represents its location
 */
class turbulentStructure
: public point
{
public:
    
    /**
     * stores the local parameters of the vorton
     */
    class Parameters
    {
    protected:        
        const inflowInputDataField& ifp_;
        
        /**
         * principal directions of Reynolds stresses
         */
        vectorField e1_, e2_, e3_;
        
        /**
         * wall distance and direction, for possible usage in scalaing factor determination
         */
        vectorField yw_;
        
        /**
        * principal components of prescribed Reynolds Stresses
        */
        vectorField Rp_;
        
        /**
         * spot extent in face normal direction
         */
        scalarField LalongN_;
        
        void initialize();
        
    public:
        Parameters(const inflowInputDataField& ifp, const dictionary&);
        virtual ~Parameters();

        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<vector>&,
            const labelList&
        );
    
        tensor Lund(label iface) const;
        
        inline const inflowInputDataField& inputData() const { return ifp_; }
        inline const vectorField& e1() const { return e1_; }
        inline const vectorField& e2() const { return e2_; }
        inline const vectorField& e3() const { return e3_; }
        inline const vectorField& Rp() const { return Rp_; }
        inline const scalarField& LalongN() const { return LalongN_; }
        
        virtual void write(Ostream&) const;

        /**
         * @brief nParamFields
         * @return number of parameter fields
         */
        virtual int nParamFields() const;

        /**
         * @brief paramFieldName
         * @param i
         * @return name of the i-th parameter field
         */
        virtual string paramFieldName(int i) const;


        typedef boost::variant<const scalarField*,const vectorField*> ParamField;
        /**
         * @brief paramField
         * @param i
         * @return values of the i-th parameter field
         */
        virtual ParamField paramField(int i) const;
    };
    
protected:
  /**
   * spot parameters
   */
  const Parameters& p_;

  /**
   * start point of TS motion
   */
  point startPoint_; 

  /**
   * intersection point of motion path and (global) patch
   */
  point footPoint_;
  
  /**
   * label of the face on which the spot was created
   */
  label creaFace_;
  
  /**
   * random parameter for scalar fluctuation
   * (needs to be different for statistical independence)
   */
  scalar epsilon_s_;
      
  
public:    
    
    /**
     * Create from stream (during restart)
     */
    turbulentStructure(const Parameters& p, Istream&);
    
    /**
     * create at some location (during run)
     */
    turbulentStructure
    (
      const Parameters& p,
      BoostRandomGen&, 
      const point& footPoint, 
      const vector& initialDelta,
      label creaface
    );

    /**
     * copy constructor
     */
    turbulentStructure(const turbulentStructure& p);
    
    virtual ~turbulentStructure();

    inline void setLocation(const point& x) { point::operator=(x); }
    
    inline const point& footPoint() const { return footPoint_; }
    
    inline void initialPositioning(const point& x, const vector& initialDelta) 
    { 
      startPoint_=x+initialDelta; 
      point::operator=(startPoint_); 
      footPoint_=x; 
    }
    
    inline const point& location() const { return *this; }

    
    scalar travelledDistance() const;
    vector motion() const;

    /**
     * @brief linearMotion
     * distance to footPoint, only component parallel to motion direction
     * @return
     */
    scalar linearMotion() const;
    scalar passedThrough() const;
    inline bool noFluctuation() const
    {
      return (mag(p_.Rp()[creaFace_])<SMALL);
    }
    
    virtual vector fluctuation(const vector&) const =0;
    virtual scalar scalar_fluctuation(const vector&) const;
    
    inline label creaFace() const { return creaFace_; }
    inline void setCreaFace(label cf) { creaFace_=cf; }
    
    // Edit
    void randomize(BoostRandomGen&);
    
    void moveForward(scalar dt);

    void operator=(const turbulentStructure&);
    
    /**
     * compare structures
     * required for storage in ptr_map
     */
    bool operator<(const turbulentStructure& o) const;
    
    friend Ostream& operator<<(Ostream&, const turbulentStructure&);
    friend Istream& operator>>(Istream&, turbulentStructure&);
        

};


template<class V>
void writeComponents(Ostream& os, const V& v)
{
  for (int i=0; i<pTraits<V>::nComponents; i++)
  {
    if (i>0) os << token::SPACE;
    os << v[i];
  }
  os << endl;
}

template<class V>
void readComponents(Istream& is, V& v)
{
  for (int i=0; i<pTraits<V>::nComponents; i++)
  {
    v[i]=readScalar(is);
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
