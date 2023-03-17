#include "vortonqueue.h"

namespace Foam {


tensor LundTransform(const symmTensor& R)
{
    if (R.xx()<0. || R.yy()<0. || R.zz()<0.)
    {
        FatalErrorIn("tensor LundTransform(const symmTensor& R)")
                <<"all diagonal elements of the Reynolds stress tensor are expected to be positve."<<endl
                <<" Got: R="<<R<<endl
               <<abort(FatalError);
    }
    tensor LT = tensor::zero;
    LT.xx()=sqrt(R.xx());
    LT.yx()=R.xy()/(SMALL+LT.xx());

    scalar rad1=R.yy()-sqr(LT.yx());
    if (rad1<0.)
    {
        FatalErrorIn("tensor LundTransform(const symmTensor& R)")
                <<"error: radicand in Lyy expression negative"<<endl
                <<" Input: R="<<R<<endl
               <<abort(FatalError);
    }
    LT.yy()=sqrt(rad1);

    LT.zx()=R.xz()/(SMALL+LT.xx());
    LT.zy()=(R.yz() - LT.yx()*LT.zx() )/(SMALL+LT.yy());

    scalar rad2=R.zz() - sqr(LT.zx()) - sqr(LT.zy());
    if (rad2<0.)
    {
        FatalErrorIn("tensor LundTransform(const symmTensor& R)")
                <<"error: radicand in Lzz expression negative"<<endl
                <<" Input: R="<<R<<endl
               <<abort(FatalError);
    }
    LT.zz()=sqrt(rad2);

    return LT;
}

} // namespace Foam
