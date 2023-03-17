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

#include "combinedVorton.H"
#include "transformField.H"
#include "gsl/gsl_multimin.h"
#include <armadillo>

#include "base/linearalgebra.h"

namespace Foam
{
  
std::unique_ptr<insight::Interpolator> readTable(fileName fn)
{
  arma::mat xy;
  fn.expand();
  xy.load(fn.c_str(), arma::raw_ascii);

  return std::unique_ptr<insight::Interpolator>(new insight::Interpolator(xy, true));
}


combinedVorton::Parameters::Parameters(const inflowInputDataField& ifp, const dictionary& dict)
  : turbulentStructure::Parameters(ifp, dict),
    rx1_(ifp.num()), ry1_(ifp.num()), rz1_(ifp.num()),
    rx2_(ifp.num()), ry2_(ifp.num()), rz2_(ifp.num()),
    sx1_(ifp.num()), sy1_(ifp.num()), sz1_(ifp.num()),
    sx2_(ifp.num()), sy2_(ifp.num()), sz2_(ifp.num()),
    residual_sigma_(ifp.num()),
    error_gamma_1_(ifp.num()), error_gamma_2_(ifp.num())
{
  if (dict.found("sigmaTableFile") && dict.found("gammaTableFile"))
    {
      sigmaFileName_=fileName(dict.lookup("sigmaTableFile"));
      gammaFileName_=fileName(dict.lookup("gammaTableFile"));
      coordinateSystem_.read(dict.lookup("coordinateSystem"));

      std::unique_ptr<insight::Interpolator> sigma_vals = readTable(sigmaFileName_);
      std::unique_ptr<insight::Interpolator> gamma_vals = readTable(gammaFileName_);

      for(label fi = 0; fi < inputData().num(); fi++)
      {
        double t = coordinateSystem_.t( inputData().x()[fi] );

        arma::mat s = (*sigma_vals)(t);
        sx1_[fi]=s(0);
        sy1_[fi]=s(1);
        sz1_[fi]=s(2);
        sx2_[fi]=s(3);
        sy2_[fi]=s(4);
        sz2_[fi]=s(5);

        arma::mat g = (*gamma_vals)(t);
        rx1_[fi]=g(0);
        ry1_[fi]=g(1);
        rz1_[fi]=g(2);
        rx2_[fi]=g(3);
        ry2_[fi]=g(4);
        rz2_[fi]=g(5);

      }
   }
  else
    {
      struct Obj : public insight::ObjectiveND
      {
        vector L_, R_;

        const int X=0;
        const int Y=1;
        const int Z=2;

        const double lolim=0.001;
        const double hilim=0.999;

        double operator()(const arma::mat& xy) const
        {
          return residual(xy, false);
        }

        double residual(const arma::mat& xy, bool show=false) const
        {
//          std::cout << xy << std::flush;

          arma::mat c = coeff(xy);

          double penalty =0.;

          {
            arma::mat x = xy.rows(0,2);
            arma::mat y = xy.rows(3,5);
            penalty =
                pow( 1.*(lolim-::min(lolim, xy.min())), 2)
                +
                pow( 1.*(::max(hilim, xy.max())-hilim), 2)

  //              +
  //              pow( 0.1*arma::norm(arma::ones(6)-xy, 2), 2)
              /*  +
                pow( 0.1*arma::norm(arma::ones(3)*arma::mean(x)-x, 2), 2)
                +
                pow( 0.1*arma::norm(arma::ones(3)*arma::mean(y)-y, 2), 2)*/
                ;
          }

          arma::mat x = arma::clamp(xy.rows(0,2), SMALL, 1.-SMALL), mx = arma::ones(3)-x;
          arma::mat y = arma::clamp(xy.rows(3,5), SMALL, 1.-SMALL), my = arma::ones(3)-y;

          double K=4.; ///::pow(M_PI, 1.5);

          double rh11= K * R_[X] *  x(X) *  y(Y)*L_[Y] *  y(Z)*L_[Z] / stabilise( y(X)*L_[X], SMALL);
          double rh21= K * R_[X] * mx(X) * my(Y)*L_[Y] * my(Z)*L_[Z] / stabilise(my(X)*L_[X], SMALL);

          double residual =
              sqr( rh11 - sqr( c(Y)*sqr(y(Y)*L_[Y]) - c(Z)*sqr(y(Z)*L_[Z]) ) )
              +
              sqr( rh21 - sqr( c(3+Y)*sqr(my(Y)*L_[Y]) - c(3+Z)*sqr(my(Z)*L_[Z]) ) )
              ;

          for(int i=0; i<3; i++)
            {
              residual += sqr (R_[i]*L_[i] - x(i)*R_[i]*y(i)*L_[i] - mx(i)*R_[i]*my(i)*L_[i] );
            }

//          std::cout << "/" << residual <<std::endl;

          if (show)
            {
              for (int i=0; i<3; i++)
                {
                  Info<<i<<": "<<x(i)<<" "<<(x(i)*R_[i])<<" "<<(mx(i)*R_[i])<<" "<<(R_[i]-x(i)*R_[i]-mx(i)*R_[i])<<endl;
                  Info<<i<<": "<<y(i)<<" "<<(x(i)*R_[i]*y(i)*L_[i])<<" "<<(mx(i)*R_[i]*my(i)*L_[i])<<" "<<(R_[i]*L_[i] - x(i)*R_[i]*y(i)*L_[i] - mx(i)*R_[i]*my(i)*L_[i])<<endl;
                }
            }

          return residual + penalty;
        }

        arma::mat coeff(const arma::mat& xy) const
        {
          arma::mat x = arma::clamp(xy.rows(0,2), SMALL, 1.-SMALL), mx = arma::ones(3)-x;
          arma::mat y = arma::clamp(xy.rows(3,5), SMALL, 1.-SMALL), my = arma::ones(3)-y;

          double K=4.; ///::pow(M_PI, 1.5);


          double rh12= K * R_[Y] *  x(Y) *  y(X)*L_[X] *  y(Z)*L_[Z] / stabilise( y(Y)*L_[Y], SMALL);
          double rh13= K * R_[Z] *  x(Z) *  y(X)*L_[X] *  y(Y)*L_[Y] / stabilise( y(Z)*L_[Z], SMALL);

          double rh22= K * R_[Y] * mx(Y) * my(X)*L_[X] * my(Z)*L_[Z] / stabilise(my(Y)*L_[Y], SMALL);
          double rh23= K * R_[Z] * mx(Z) * my(X)*L_[X] * my(Y)*L_[Y] / stabilise(my(Z)*L_[Z], SMALL);

          double g11 = 0.0;
          double g12 = ( g11 * sqr(y(Y)*L_[Y]) - ::sqrt(rh13) ) / stabilise(sqr( y(Y)*L_[Y] ), SMALL);
          double g13 = ( g11 * sqr(y(X)*L_[X]) - ::sqrt(rh12) ) / stabilise(sqr( y(Z)*L_[Z] ), SMALL);

          double g21 = 0.0;
          double g22 = ( g21 * sqr(my(Y)*L_[Y]) - ::sqrt(rh23) ) / stabilise(sqr( my(Y)*L_[Y] ), SMALL);
          double g23 = ( g21 * sqr(my(X)*L_[X]) - ::sqrt(rh22) ) / stabilise(sqr( my(Z)*L_[Z] ), SMALL);

          const double K2=1./::sqrt(M_PI);
          arma::mat c;
          c
              << K2* y(X)*L_[X] << K2* y(Y)*L_[Y] << K2* y(Z)*L_[Z]
              << K2*my(X)*L_[X] << K2*my(Y)*L_[Y] << K2*my(Z)*L_[Z]
              << g11 << g12 << g13
              << g21 << g22 << g23
                 ;
          return c;
        }

        int numP() const { return 6; }

      } o;


      double maxeps=-GREAT;
      arma::mat x0;
      for(label fi = 0; fi < inputData().num(); fi++)
      {
          o.L_=inputData().L()[fi];
          o.R_=this->Rp()[fi];

          if (x0.n_rows==0)
            {
              x0 = arma::ones(6)*0.5;
            }
          arma::mat steps = arma::ones(6)*0.2;
          arma::mat x = insight::nonlinearMinimizeND(o, x0, 1e-5, steps);
          x0=x;
          arma::mat c = o.coeff(x);

          double eps=o.residual(x, true);
          maxeps=::max(maxeps,eps);
          std::cout<<x<< " / eps="<<eps<<std::endl;

          residual_sigma_[fi]=eps;

          sx1_[fi]=c(0);
          sy1_[fi]=c(1);
          sz1_[fi]=c(2);
          sx2_[fi]=c(3);
          sy2_[fi]=c(4);
          sz2_[fi]=c(5);

          rx1_[fi]=c(6);
          ry1_[fi]=c(7);
          rz1_[fi]=c(8);
          rx2_[fi]=c(9);
          ry2_[fi]=c(10);
          rz2_[fi]=c(11);

          error_gamma_1_[fi] = 0;
          error_gamma_2_[fi] = 0;
      }

      Info << "maximum residual = " << maxeps <<endl;
    }
}

void combinedVorton::Parameters::write(Ostream& os) const
{
  if (!sigmaFileName_.empty() && ! gammaFileName_.empty())
    {
      os.writeKeyword("sigmaTableFile") << token::SPACE << sigmaFileName_ << token::END_STATEMENT << nl;
      os.writeKeyword("gammaTableFile") << token::SPACE << gammaFileName_ << token::END_STATEMENT << nl;
      os.writeKeyword("coordinateSystem") << token::SPACE;
    }
  coordinateSystem_.writeSup(os);
  os << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


combinedVorton::combinedVorton
(
    const Parameters& p,
    Istream& s
)
: turbulentStructure(p, s),
  p_(p)
{
  s>>(*this);
}


combinedVorton::combinedVorton
(
    const Parameters& p,
    BoostRandomGen& rg,
    const point& footPoint,
    const vector& initialDelta,
    label creaface
)
: turbulentStructure
  (
    p,
    rg,
    footPoint,
    initialDelta,
    creaface
  ),
  p_(p),
  epsilon1_(1.0),
  epsilon2_(1.0)
{  
}

combinedVorton::combinedVorton(const combinedVorton& o)
: turbulentStructure(o),
  p_(o.p_),
  epsilon1_(o.epsilon1_),
  epsilon2_(o.epsilon2_)
{ 
}

vector combinedVorton::fluctuation(const vector& x) const
{
  vector delta_x = x - location();


  vector
      e1=p_.e1()[creaFace_],
      e2=p_.e2()[creaFace_],
      e3=p_.e3()[creaFace_];
  
  double Xx=delta_x&e1;
  double Yy=delta_x&e2;
  double Zz=delta_x&e3;


  tensor trs = tensor(e1, e2, e3).T();
  
  vector up=vector::zero;
  
  double xx2=Xx*Xx;
  double yy2=Yy*Yy;
  double zz2=Zz*Zz;

  const double low=SMALL;

  vector L=p_.inputData().L()[creaFace_];
  if
  (
      (mag(Xx)  < 4.*L[0]) &&
      (mag(Yy)  < 4.*L[1]) &&
      (mag(Zz)  < 4.*L[2])
  )
    {
      {
        double rx=p_.rx1_[creaFace_];
        double ry=p_.ry1_[creaFace_];
        double rz=p_.rz1_[creaFace_];
        double sx=p_.sx1_[creaFace_];
        double sy=p_.sy1_[creaFace_];
        double sz=p_.sz1_[creaFace_];

        double sx2=stabilise( sx*sx, low);
        double sy2=stabilise( sy*sy, low);
        double sz2=stabilise( sz*sz, low);

        double e=exp( -0.5 * (xx2/sx2 + yy2/sy2 + zz2/sz2) );

        up +=  epsilon1_ * e * vector
            (
                Yy*Zz*( ry/sz2 - rz/sy2 ),
                Xx*Zz*( rz/sx2 - rx/sz2 ),
                Xx*Yy*( rx/sy2 - ry/sx2 )
            );
      }

      {
        double rx=p_.rx2_[creaFace_];
        double ry=p_.ry2_[creaFace_];
        double rz=p_.rz2_[creaFace_];
        double sx=p_.sx2_[creaFace_];
        double sy=p_.sy2_[creaFace_];
        double sz=p_.sz2_[creaFace_];

        double sx2=stabilise( sx*sx, low);
        double sy2=stabilise( sy*sy, low);
        double sz2=stabilise( sz*sz, low);

        double e=exp( -0.5 * (xx2/sx2 + yy2/sy2 + zz2/sz2) );

        up +=  epsilon2_ * e * vector
            (
                Yy*Zz*( ry/sz2 - rz/sy2 ),
                Xx*Zz*( rz/sx2 - rx/sz2 ),
                Xx*Yy*( rx/sy2 - ry/sx2 )
            );
      }
    }

//   Info<<er1_<<"\t"<<er2_<<"\t"<<er3_<<"\t"<<trs<<"\t"<<Rp_<<endl;
//   Info<<(trs&trs.T())<<endl;

  return vector
  (
    transform
    ( 
      trs, 
      up
    )
  );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //



void combinedVorton::randomize(BoostRandomGen& rand)
{
  epsilon1_ = 2.0*(rand() - 0.5);
  epsilon2_ = 2.0*(rand() - 0.5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

combinedVorton::~combinedVorton()
{}


autoPtr<combinedVorton> combinedVorton::clone() const
{
  return autoPtr<combinedVorton>
  (
    new combinedVorton(*this)
  );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void combinedVorton::operator=(const combinedVorton& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("combinedVorton::operator=(const combinedVorton&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    turbulentStructure::operator=(rhs);
    epsilon1_=rhs.epsilon1_;
    epsilon2_=rhs.epsilon2_;
}

// bool combinedVorton::operator!=(const combinedVorton& o) const
// {
//     return 
// //     turbulentStructure::operator!=(o)
// // ||
//         (location()!=o.location())
//         ||
//         (epsilon_!=o.epsilon_)
// 	||
//         (rx_!=o.rx_)
// 	||
//         (ry_!=o.ry_)
// 	||
//         (rz_!=o.rz_)
// 	||
//         (sx_!=o.sx_)
// 	||
//         (sy_!=o.sy_)
// 	||
//         (sz_!=o.sz_);
// }


int combinedVorton::Parameters::nParamFields() const
{
  return turbulentStructure::Parameters::nParamFields()+15;
}

string combinedVorton::Parameters::paramFieldName(int i) const
{
  int i0 = i - turbulentStructure::Parameters::nParamFields();
  switch(i0)
    {
    case 0: return "rx1"; break;
    case 1: return "ry1"; break;
    case 2: return "rz1"; break;
    case 3: return "sx1"; break;
    case 4: return "sy1"; break;
    case 5: return "sz1"; break;

    case 6: return "rx2"; break;
    case 7: return "ry2"; break;
    case 8: return "rz2"; break;
    case 9: return "sx2"; break;
    case 10: return "sy2"; break;
    case 11: return "sz2"; break;

    case 12: return "residual_sigma"; break;
    case 13: return "error_gamma_1"; break;
    case 14: return "error_gamma_2"; break;
    }
  return turbulentStructure::Parameters::paramFieldName(i);
}

turbulentStructure::Parameters::ParamField combinedVorton::Parameters::paramField(int i) const
{
  int i0 = i - turbulentStructure::Parameters::nParamFields();
  switch(i0)
    {
    case 0: return &rx1_; break;
    case 1: return &ry1_; break;
    case 2: return &rz1_; break;
    case 3: return &sx1_; break;
    case 4: return &sy1_; break;
    case 5: return &sz1_; break;

    case 6: return &rx2_; break;
    case 7: return &ry2_; break;
    case 8: return &rz2_; break;
    case 9: return &sx2_; break;
    case 10: return &sy2_; break;
    case 11: return &sz2_; break;

    case 12: return &residual_sigma_; break;
    case 13: return &error_gamma_1_; break;
    case 14: return &error_gamma_2_; break;
    }
  return turbulentStructure::Parameters::paramField(i);
}


Ostream& operator<<(Ostream& s, const combinedVorton& ht)
{
    s << static_cast<const turbulentStructure&>(ht);
    s<<ht.epsilon1_<<endl;
    s<<ht.epsilon2_<<endl;
    return s;
}

Istream& operator>>(Istream& s, combinedVorton& ht)
{
    s >> static_cast<turbulentStructure&>(ht);
    ht.epsilon1_=readScalar(s);
    ht.epsilon2_=readScalar(s);
    return s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
