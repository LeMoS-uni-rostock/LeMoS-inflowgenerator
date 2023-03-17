
#include "turbulentstructure.h"

namespace Foam
{

bool MaxTransversalResolvedByAverageEdge::operator()
(const vector& Lprincip, const vector& edgeLength) const
{
    return Lprincip.y() > edgeLength.y();
}

bool MaxTransversalResolvedByMaxEdge::operator()
(const vector& Lprincip, const vector& edgeLength) const
{
    return max(Lprincip.y(),Lprincip.z()) > edgeLength.z();
}

bool MaxLengthResolvedByMaxEdge::operator()
(const vector& Lprincip, const vector& edgeLength) const
{
    return cmptMax(Lprincip) > edgeLength.z();
}


turbulentStructure::turbulentStructure(const point& p)
    : point(p)
{}

turbulentStructure::turbulentStructure(Istream &is)
{
    is >> point::x() >> point::y() >> point::z();
}

turbulentStructure::~turbulentStructure()
{}

void turbulentStructure::write(Ostream &os) const
{
    os
            << point::x() << token::SPACE
            << point::y() << token::SPACE
            << point::z();
}


}
