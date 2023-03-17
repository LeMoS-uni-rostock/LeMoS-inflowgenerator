#ifndef TURBULENTSTRUCTURE_H
#define TURBULENTSTRUCTURE_H

#include "boostRandomGen.H"

#include "point.H"
#include "symmTensor.H"


namespace Foam
{


class MaxTransversalResolvedByAverageEdge
{
public:
    bool operator()(const vector& Lprincip, const vector& edgeLength) const;
};

class MaxTransversalResolvedByMaxEdge
{
public:
    bool operator()(const vector& Lprincip, const vector& edgeLength) const;
};

class MaxLengthResolvedByMaxEdge
{
public:
    bool operator()(const vector& Lprincip, const vector& edgeLength) const;
};



class turbulentStructure
        : public point
{
public:

    template<class ResolutionCriterion>
    struct StaticVortonParameters
    {
        bool enabled;
        int disableReason;

        vector UMean;
        vector motionDirection;
        vector L;
        scalar c;

        scalar delta; // average distance between subsequent vortons

        StaticVortonParameters()
            : enabled(false),
              disableReason(-1)
        {}

        StaticVortonParameters(
                const vector& umean,
                const symmTensor& r,
                const vector& l,
                scalar _c,
                scalar faceArea,
                const vector& faceEdgeLength
                )
            : UMean(umean),
              L(l),
              c(_c),
              enabled(true),
              disableReason(0)
        {
            if ( faceArea<SMALL )
            {
                disableReason += 1;
                enabled=false;
            }
            if ( !ResolutionCriterion()(L, faceEdgeLength) /*cmptMin(L)*/ /*Foam::max(L.y(), L.z())*/ /*cmptMax(L) < faceEdgeLength*/ )
            {
                disableReason += 2;
                enabled=false;
            }
            if ( 0.5*tr(r) < 0.001 * 0.5*(umean&umean) )
            {
                disableReason += 4;
                enabled=false;
            }
            if ( mag(UMean) < 0.001 )
            {
                disableReason += 8;
                enabled=false;
            }

            if (enabled)
            {
                delta = L.x()*L.y()*L.z()
                    /
                    ( c * faceArea );

                motionDirection = UMean/mag(UMean);
            }
        }

        operator bool() const
        {
            return enabled;
        }

    };

public:
    turbulentStructure(const point& p);
    turbulentStructure(Istream& is);
    virtual ~turbulentStructure();

    virtual void write(Ostream& os) const;
};


}

#endif // TURBULENTSTRUCTURE_H
