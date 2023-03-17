#ifndef ANISOTROPICVORTON_H
#define ANISOTROPICVORTON_H

#include "turbulentstructure.h"
#include "vortonqueue.h"

#include "base/linearalgebra.h"

namespace Foam
{


struct Eigensystem
        : public vector /* eigenvalues */
{
    vector e_[3];

    Eigensystem();
    Eigensystem(const tensor& t);
};


template<class ResolutionCriterion = MaxLengthResolvedByMaxEdge>
struct StaticAnisotropicParametersBase
        : public turbulentStructure::StaticVortonParameters<ResolutionCriterion>
{
    vector Rprincip;
    vector e[3];
    scalar sx, sy, sz;
    scalar rx, ry, rz;

    StaticAnisotropicParametersBase()
        : turbulentStructure::StaticVortonParameters<ResolutionCriterion>()
    {}

    StaticAnisotropicParametersBase(
            const vector& umean,
            const symmTensor& r,
            const vector& l,
            scalar c,
            scalar faceArea,
            const vector& faceEdgeLength )
        : turbulentStructure::StaticVortonParameters<ResolutionCriterion>(
              umean, r, l, c,
              faceArea, faceEdgeLength )
    {
        Eigensystem es(r);
        Rprincip=es;
        for (int i=0; i<3; ++i)
            e[i]=es.e_[i];
    }

    scalar maxL() const
    {
        return 2.*sqrt(M_PI)*max(max(sx, sy), sz);
    }

    scalar influenceLength() const
    {
        return maxL();
    }
};



template<class AnisotropicVortonParametersComputer>
class anisotropicVorton
        : public turbulentStructure
{
public:

    typedef
        AnisotropicVortonParametersComputer
        StaticVortonParameters;

private:
    const StaticVortonParameters& prop_;
    scalar epsilon_;

public:
    anisotropicVorton(const StaticVortonParameters& prop, const point& p)
        : turbulentStructure(p),
          prop_(prop),
          epsilon_( 2.0*(rand() - 0.5) ) // -1 .. +1
    {}

    anisotropicVorton(const StaticVortonParameters& prop, Istream& is)
        : turbulentStructure(is),
          prop_(prop)
    {
        is >> epsilon_;
    }

    // compute fluctuation induced at x
    vector fluctuation(const point& p) const
    {
        vector delta_x = p - *this;

        scalar x = delta_x & prop_.e[0];
        scalar y = delta_x & prop_.e[1];
        scalar z = delta_x & prop_.e[2];

        scalar SPI=sqrt(M_PI);
        if  (
             (mag(x) < prop_.influenceLength()) &&
             (mag(y) < prop_.influenceLength()) &&
             (mag(z) < prop_.influenceLength())
             )
        {

            double sx2=prop_.sx*prop_.sx;
            double sy2=prop_.sy*prop_.sy;
            double sz2=prop_.sz*prop_.sz;
            double x2=x*x;
            double y2=y*y;
            double z2=z*z;

            return  transform
            (
                tensor(prop_.e[0], prop_.e[1], prop_.e[2]).T(),
                epsilon_
                    * exp( -0.5 * (x2/sx2 + y2/sy2 + z2/sz2) )
                    * vector
                    (
                        y*z*( prop_.ry/sz2 - prop_.rz/sy2 ),
                        x*z*( prop_.rz/sx2 - prop_.rx/sz2 ),
                        x*y*( prop_.rx/sy2 - prop_.ry/sx2 )
                    )
            );
        }
        else
        {
            return pTraits<vector>::zero;
        }
    }

    void write(Ostream& os) const override
    {
        turbulentStructure::write(os);
        os
                << token::SPACE << epsilon_;
    }
};




template<class ResolutionCriterion = MaxLengthResolvedByMaxEdge>
struct PseudoInvsereAnisotropicVortonParametersComputer
: public StaticAnisotropicParametersBase<ResolutionCriterion>
{
    PseudoInvsereAnisotropicVortonParametersComputer()
        : StaticAnisotropicParametersBase<ResolutionCriterion>()
    {}

    PseudoInvsereAnisotropicVortonParametersComputer(
            const vector& umean,
            const symmTensor& r,
            const vector& l,
            scalar c,
            scalar faceArea,
            const vector& faceEdgeLength )
        : StaticAnisotropicParametersBase<ResolutionCriterion>(
              umean, r, l, c,
              faceArea, faceEdgeLength )
    {
            double
                L1 = this->L[0],
                L2 = this->L[1],
                L3 = this->L[2]
                ;

            double
                R1 = std::max(this->Rprincip[0], 1e-3),
                R2 = std::max(1e-3*R1, std::min(this->Rprincip[1], 0.999*R1)),
                R3 = std::max(1e-3*R1, std::min(this->Rprincip[2], 0.999*R1));

            this->sx = L1/sqrt(M_PI)/2.;
            this->sy = L2/sqrt(M_PI)/2.;
            this->sz = L3/sqrt(M_PI)/2.;

            /*
            * determine ri via pseudo-inverse
            */
            arma::mat cM, rhs;

            Info<<"L1/L2/L3="<<L1<<"/"<<L2<<"/"<<L3<<endl;

            cM << 0.    << L2/L3	<< -L3/L2 << arma::endr
               << L1/L3	<< 0.       << -L3/L1 << arma::endr
               << L1/L2	<< -L2/L1 	<< 0.     << arma::endr;
            cM *= sqrt( c / 96. );


            rhs << sqrt(R1) << arma::endr
                << sqrt(R2) << arma::endr
                << sqrt(R3) << arma::endr;

            arma::mat rxryrz = arma::pinv(cM)*rhs;

            this->rx = rxryrz(0);
            this->ry = rxryrz(1);
            this->rz = rxryrz(2);

            Info<<this->Rprincip<<this->e[0]<<this->e[1]<<this->e[2]
               <<this->sx<<" "<<this->sy<<" "<<this->sz<<" "
              <<this->rx<<" "<<this->ry<<" "<<this->rz<<endl;
            std::cout<<"R="<<pow(cM*rxryrz,2)<<std::endl;

    }
};



template<class ResolutionCriterion = MaxLengthResolvedByMaxEdge>
struct AnalyticAnisotropicVortonParametersComputer
: public StaticAnisotropicParametersBase<ResolutionCriterion>
{
    AnalyticAnisotropicVortonParametersComputer()
        : StaticAnisotropicParametersBase<ResolutionCriterion>()
    {}

    AnalyticAnisotropicVortonParametersComputer(
            const vector& umean,
            const symmTensor& r,
            const vector& l,
            scalar c,
            scalar faceArea,
            const vector& faceEdgeLength )
        : StaticAnisotropicParametersBase<ResolutionCriterion>(
              umean, r, l, c,
              faceArea, faceEdgeLength )
    {

        double
            L1 = this->L[0],
            L2 = this->L[1],
            L3 = this->L[2];

        double
            R1 = std::max(this->Rprincip[0], 1e-3),
            R2 = std::max(1e-3*R1, std::min(this->Rprincip[1], 0.999*R1)),
            R3 = std::max(1e-3*R1, std::min(this->Rprincip[2], 0.999*R1));

        // 1.) modify L3
        //L3 = L1*L2*sqrt(R3) / ( L2*sqrt(R1) + L1*sqrt(R2) );
        // 1.) modify L2
        L2 = L1*L3*sqrt(R2) / ( L3*sqrt(R1) + L1*sqrt(R3) );

        this->rx=-4.*sqrt(6.)*L3*sqrt(R2)/L1/sqrt(c);
        this->ry=-4.*sqrt(6.)*(L1*L3*sqrt(R2)-L1*L2*sqrt(R3))/pow(L2,2)/sqrt(c);
        this->rz=0.0;



        this->sx=L1/sqrt(M_PI)/2.;
        this->sy=L2/sqrt(M_PI)/2.;
        this->sz=L3/sqrt(M_PI)/2.;

        /*
            Info<<this->Rprincip<<this->e[0]<<this->e[1]<<this->e[2]
               <<this->sx<<" "<<this->sy<<" "<<this->sz<<" "
              <<this->rx<<" "<<this->ry<<" "<<this->rz<<endl;
            Info
                    << "R1="<<(1./96.*c*pow(this->sy*this->sz*(this->rz/this->sy/this->sy - this->ry/this->sz/this->sz),2))<<" ("<<R1<<")\n"
                    << "R2="<<(1./96.*c*pow(this->sx*this->sz*(this->rz/this->sx/this->sx - this->rx/this->sz/this->sz),2))<<" ("<<R2<<")\n"
                    << "R3="<<(1./96.*c*pow(this->sx*this->sy*(this->ry/this->sx/this->sx - this->rx/this->sy/this->sy),2))<<" ("<<R3<<")\n"
                    << "L1="<<(2.*this->sx*sqrt(M_PI))<<" ("<<L1<<")\n"
                    << "L2="<<(2.*this->sy*sqrt(M_PI))<<" ("<<L2<<")\n"
                    << "L3="<<(2.*this->sz*sqrt(M_PI))<<" ("<<L3<<")\n"
                       ;
        */
        
        // recalc delta with modified length scales
        if (this->enabled)
        {
            this->L = 2.*sqrt(M_PI)*vector(this->sx, this->sy, this->sz);

            this->delta = pow(2.*sqrt(M_PI),3)*this->sx*this->sy*this->sz
                    /
                    ( c * faceArea );
        }
    }
};



template<class ResolutionCriterion = MaxLengthResolvedByMaxEdge>
struct NumericalAnisotropicVortonParametersComputer
: public StaticAnisotropicParametersBase<ResolutionCriterion>
{
    NumericalAnisotropicVortonParametersComputer()
        : StaticAnisotropicParametersBase<ResolutionCriterion>()
    {}

    NumericalAnisotropicVortonParametersComputer(
            const vector& umean,
            const symmTensor& r,
            const vector& l,
            scalar c,
            scalar faceArea,
            const vector& faceEdgeLength )
        : StaticAnisotropicParametersBase<ResolutionCriterion>(
              umean, r, l, c,
              faceArea, faceEdgeLength )
    {

        double
            L1 = this->L[0],
            L2 = this->L[1],
            L3 = this->L[2];

        double
            R1 = std::max(this->Rprincip[0], 1e-3),
            R2 = std::max(1e-3*R1, std::min(this->Rprincip[1], 0.999*R1)),
            R3 = std::max(1e-3*R1, std::min(this->Rprincip[2], 0.999*R1));


        arma::mat rxryrzsxsysz0;
        rxryrzsxsysz0
                << 1.
                << 1.
                << 1.
                << (L1/sqrt(M_PI)/2.)
                << (L2/sqrt(M_PI)/2.)
                << (L3/sqrt(M_PI)/2.) ;

        arma::mat rxryrzsxsysz = insight::nonlinearMinimizeND(
                    [&](const arma::mat& rxryrzsxsysz)
                    {
                        double
                                rx=rxryrzsxsysz(0),
                                ry=rxryrzsxsysz(1),
                                rz=rxryrzsxsysz(2),
                                sx=rxryrzsxsysz(3),
                                sy=rxryrzsxsysz(4),
                                sz=rxryrzsxsysz(5)
                                ;
                        return
                                 pow( (R1 - 1./96.*c*pow(sy*sz*(rz/sy/sy - ry/sz/sz),2))/R1, 2)
                                +pow( (R2 - 1./96.*c*pow(sx*sz*(rz/sx/sx - rx/sz/sz),2))/R2, 2)
                                +pow( (R3 - 1./96.*c*pow(sx*sy*(ry/sx/sx - rx/sy/sy),2))/R3, 2)
                                +pow( (L1 - 2.*sqrt(M_PI)*sx)/L1, 2)
                                +pow( (L2 - 2.*sqrt(M_PI)*sy)/L2, 2)
                                +pow( (L3 - 2.*sqrt(M_PI)*sz)/L3, 2)
                                ;
                    },
                    rxryrzsxsysz0
        );

        this->rx=rxryrzsxsysz(0);
        this->ry=rxryrzsxsysz(1);
        this->rz=rxryrzsxsysz(2);
        this->sx=rxryrzsxsysz(3);
        this->sy=rxryrzsxsysz(4);
        this->sz=rxryrzsxsysz(5);
        
        /*
        Info<<this->Rprincip<<this->e[0]<<this->e[1]<<this->e[2]
           <<this->sx<<" "<<this->sy<<" "<<this->sz<<" "
          <<this->rx<<" "<<this->ry<<" "<<this->rz<<endl;

        Info
                << "R1="<<(1./96.*c*pow(this->sy*this->sz*(this->rz/this->sy/this->sy - this->ry/this->sz/this->sz),2))<<" ("<<R1<<")\n"
                << "R2="<<(1./96.*c*pow(this->sx*this->sz*(this->rz/this->sx/this->sx - this->rx/this->sz/this->sz),2))<<" ("<<R2<<")\n"
                << "R3="<<(1./96.*c*pow(this->sx*this->sy*(this->ry/this->sx/this->sx - this->rx/this->sy/this->sy),2))<<" ("<<R3<<")\n"
                << "L1="<<(2.*this->sx*sqrt(M_PI))<<" ("<<L1<<")\n"
                << "L2="<<(2.*this->sy*sqrt(M_PI))<<" ("<<L2<<")\n"
                << "L3="<<(2.*this->sz*sqrt(M_PI))<<" ("<<L3<<")\n"
                   ;
        */
        // recalc delta with modified length scales
        if (this->enabled)
        {
            this->L = 2.*sqrt(M_PI)*vector(this->sx, this->sy, this->sz);

            this->delta = pow(2.*sqrt(M_PI),3)*this->sx*this->sy*this->sz
                    /
                    ( c * faceArea );
        }

    }


};

}


#endif // ANISOTROPICVORTON_H
