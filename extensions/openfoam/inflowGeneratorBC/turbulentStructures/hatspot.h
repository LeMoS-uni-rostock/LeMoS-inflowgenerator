#ifndef HATSPOT_H
#define HATSPOT_H

#include "turbulentstructure.h"
#include "vortonqueue.h"

namespace Foam
{



class hatSpot
        : public turbulentStructure
{
public:
    struct StaticVortonParameters
            : public turbulentStructure::StaticVortonParameters<MaxLengthResolvedByMaxEdge>
    {
        tensor Lund;
        vector e[3];

        StaticVortonParameters()
            : turbulentStructure::StaticVortonParameters<MaxLengthResolvedByMaxEdge>()
        {}

        StaticVortonParameters(
                const vector& umean,
                const symmTensor& r,
                const vector& l,
                scalar c,
                scalar faceArea,
                const vector& faceEdgeLength )
            : turbulentStructure::StaticVortonParameters<MaxLengthResolvedByMaxEdge>(
                  umean, r, l, c,
                  faceArea, faceEdgeLength )
        {
            Lund = LundTransform(r)
                    / sqrt( c / 81. );
            e[0]=vector(1,0,0);
            e[1]=vector(0,1,0);
            e[2]=vector(0,0,1);
        }

        scalar maxL() const
        {
            return cmptMax(L);
        }

        scalar influenceLength() const
        {
            return maxL();
        }
    };

private:
    const StaticVortonParameters& prop_;
    vector epsilon_;

public:
    hatSpot(const StaticVortonParameters& prop, const point& p)
        : turbulentStructure(p),
          prop_(prop),
          epsilon_(
              2.0*(rand() - 0.5),
              2.0*(rand() - 0.5),
              2.0*(rand() - 0.5)
              )
    {}

    hatSpot(const StaticVortonParameters& prop, Istream& is)
        : turbulentStructure(is),
          prop_(prop)
    {
        is >> epsilon_.x() >> epsilon_.y() >> epsilon_.z();
    }

    // compute fluctuation induced at x
    vector fluctuation(const point& p) const
    {
        vector delta_x = p - *this;
        scalar x = mag(delta_x & prop_.e[0]);
        scalar y = mag(delta_x & prop_.e[1]);
        scalar z = mag(delta_x & prop_.e[2]);

        if  (
             (x  < (prop_.L[0] / 2.0)) &&
             (y  < (prop_.L[1] / 2.0)) &&
             (z  < (prop_.L[2] / 2.0))
             )
        {

            vector uprime = cmptMultiply(epsilon_,

                     (1.0 - 2.0*x/prop_.L[0] )
                    *(1.0 - 2.0*y/prop_.L[1] )
                    *(1.0 - 2.0*z/prop_.L[2] )
                    * pTraits<vector>::one

                    );

            return prop_.Lund & uprime;
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
                << token::SPACE << epsilon_.x()
                << token::SPACE << epsilon_.y()
                << token::SPACE << epsilon_.z();
    }
};

}

#endif // HATSPOT_H
