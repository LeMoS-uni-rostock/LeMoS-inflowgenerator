#ifndef FOAM_VORTONQUEUE_H
#define FOAM_VORTONQUEUE_H

#include <deque>
#include <functional>
#include <memory>

#include "fvCFD.H"

#include "boostRandomGen.H"

namespace Foam {


tensor LundTransform(const symmTensor& R);


template<class TurbulentStructure>
class vortonQueue
        : public std::deque<std::shared_ptr<TurbulentStructure> >
{
public:
    typedef std::shared_ptr<TurbulentStructure> TurbulentStructurePtr;

    vector edgeLength() const
    {
        scalar minL=GREAT;
        scalar avgL=0.0;
        scalar maxL=-GREAT;
        const auto& el = face_.edges();
        forAll(el, ei)
        {
            scalar l=el[ei].mag(points_);
            minL = min(minL, l);
            avgL += l;
            maxL = max(maxL, l);
        }
        return vector(
                    minL,
                    avgL / scalar(el.size()),
                    maxL
                    );
    }

private:
    const face& face_;
    const pointField& points_;
    point basePoint_;

    typename TurbulentStructure::StaticVortonParameters vprops_;


    TurbulentStructurePtr createVorton( const point& p)
    {
        return std::make_shared<TurbulentStructure>(vprops_, p);
    }

    void moveVortons(scalar deltaT)
    {
        for (auto& v: *this)
        {
            (*v) += vprops_.UMean * deltaT;
        }
    }

    vector randomLateralDisplacement() const
    {

        label triI;
        scalar u = rand(),
               v = rand();

        faceList tris(face_.nTriangles());
        {
            label tril=0;
            face_.triangles(points_, tril, tris); // split face into triangles

            scalar A_tot=0.0;
            std::vector<double> A_f;
            forAll(tris, ti)
            {
              scalar A=tris[ti].mag(points_);
              A_tot+=A;
              A_f.push_back(A);
            }
            forAll(tris, ti)
            {
              A_f[ti]/=A_tot;
            }

            // setup discrete random generator with area-weighted probability

            boost::random::discrete_distribution<> tri_selector(A_f.begin(), A_f.end());
            triI = tri_selector(rand);
        }
        const face& t = tris[triI];

        vector a = points_[t[1]] - points_[t[0]];
        vector b = points_[t[2]] - points_[t[0]];
        if ( (u+v) > 1. )
        {
            u = 1.-u;
            v = 1.-v;
        }
        return
                points_[t[0]] + u*a + v*b
                - basePoint_;
    }


    /**
     * @brief fillDistance
     * @param dist
     * negative value: fill upstream
     * @param nmin
     * insert at least nmin new vortons
     */
    void fillDistance(scalar dist, label nmin=0)
    {
        if (this->size()==0)
        {
            this->reset();
        }


        std::function<typename vortonQueue::const_reference(void)> lastVorton;
        std::function<void(typename vortonQueue::value_type&&)> appendVorton;
        std::function<bool()> lastVortonTooClose;

        if (dist<0.)
        {
            lastVorton = std::bind(
                        static_cast<typename vortonQueue::const_reference (vortonQueue::*)() const>(
                            &vortonQueue::front ),
                        this );
            appendVorton = std::bind(
                        static_cast<void (vortonQueue::*)(const typename vortonQueue::value_type&)>(
                            &vortonQueue::push_front ),
                        this,
                        std::placeholders::_1);
            lastVortonTooClose = [&]() {
                return ( ((*lastVorton() - basePoint_)&vprops_.motionDirection) < -dist );
            };
        }
        else
        {
            lastVorton = std::bind(
                        static_cast<typename vortonQueue::const_reference (vortonQueue::*)() const>(
                            &vortonQueue::back),
                        this );
            appendVorton = std::bind(
                        static_cast<void (vortonQueue::*)(const typename vortonQueue::value_type&)>(
                            &vortonQueue::push_back ),
                        this,
                        std::placeholders::_1 );
            lastVortonTooClose = [&]() {
                return ( -((*lastVorton() - basePoint_)&vprops_.motionDirection) < dist );
            };
        }

        for (label nAdded=0; lastVortonTooClose()||(nAdded < nmin); ++nAdded )
        {
            vector lastVortonDistance = *lastVorton() - basePoint_;

            lastVortonDistance = (lastVortonDistance & vprops_.motionDirection) * vprops_.motionDirection;

            appendVorton( createVorton(
                        basePoint_
                        + lastVortonDistance
                        + 2.*rand()*std::copysign( vprops_.delta, -dist )*vprops_.motionDirection
                        + randomLateralDisplacement() ) );
        }
    }


    void removeDistantVortons()
    {
        if (this->size()>0)
        {
            while ( ( (*this->front() - basePoint_)&vprops_.motionDirection ) > vprops_.maxL() )
            {
                if (this->size()==1)
                {
                    // about to remove last vorton,
                    // insert one upstream to maintain chain
                    fillDistance( vprops_.maxL(), 1 );
                    this->pop_front();
                    break;
                }
                this->pop_front();
            }
        }
    }

public:

    vortonQueue(
            const face& face,
            const pointField& points,
            const vector& UMean,
            const symmTensor& R,
            const vector& L,
            scalar c
            )
        : face_(face),
          points_(points),
          basePoint_( face_.centre(points_) ), // moving towards inside, normal pointing outside
          vprops_( UMean, R, L, c,
                   face_.mag(points_), edgeLength() )
    {}



    vortonQueue(
            const vortonQueue& o
            )
        : face_(o.face_),
          points_(o.points_),
          basePoint_( o.basePoint_ ), // moving towards inside, normal pointing outside
          vprops_( o.vprops_ )
    {
        for (const auto& v: o)
        {
            this->push_back( std::make_shared<TurbulentStructure>(*v) );
        }
    }

    void reset()
    {
        // remove all
        std::deque<TurbulentStructurePtr> empty;
        std::swap( *this, empty );

        if (vprops_)
        {
            // insert first vorton
            this->push_back(
                        createVorton(
                            basePoint_
                            + 2.*(rand()-0.5)*vprops_.delta*vprops_.motionDirection
                            )
                        );

            fillDistance( -vprops_.influenceLength() ); // upstream
            fillDistance( vprops_.influenceLength() ); // downstream
        }
    }


    void setProperties(
            const vector& UMean,
            const symmTensor& R,
            const vector& L,
            scalar c
            )
    {
        vprops_ = TurbulentStructure::StaticVortonParameters(
                    UMean, R, L, c,
                    face_.mag(points_), edgeLength() );
    }


    const typename TurbulentStructure::StaticVortonParameters& properties() const
    {
        return vprops_;
    }



    void advance(scalar deltaT)
    {
        moveVortons(deltaT);

        if (vprops_)
        {
            fillDistance( vprops_.influenceLength() );
        }

        removeDistantVortons();
    }

    autoPtr<vortonQueue> clone() const
    {
        return autoPtr<vortonQueue>(new vortonQueue(*this) );
    }


    void read(Istream& is)
    {
        label n;
        is >> n;
        for (label j=0; j<n; ++j)
        {
            this->push_back( std::make_shared<TurbulentStructure>(vprops_, is) );
        }
    }


    void write(Ostream& os) const
    {
        label size = this->std::deque<std::shared_ptr<TurbulentStructure> >::size();
        os << size;
        for (auto i=this->begin(); i!=this->end(); ++i)
        {
            os << token::SPACE;
            (*i)->write(os);
        }
    }
};

} // namespace Foam

#endif // FOAM_VORTONQUEUE_H
