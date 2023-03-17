
#include "anisotropicvorton.h"


namespace Foam
{


Eigensystem::Eigensystem()
{}

Eigensystem::Eigensystem(const tensor& t)
{
    arma::mat at;
    at
      << t.xx() << t.xy() << t.xz() << arma::endr
      << t.xy() << t.yy() << t.yz() << arma::endr
      << t.xz() << t.yz() << t.zz() << arma::endr
      ;

    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, at);

    arma::uvec idx=arma::sort_index(eigval,
#if ARMA_VERSION_MAJOR>=9
                                    "descend"
#else
                                    1
#endif
                                    );

    for (int i=0; i<3; ++i)
    {
        (*this)[i]=eigval(idx(i));
        e_[i]=vector(eigvec.col(idx(i))(0), eigvec.col(idx(i))(1), eigvec.col(idx(i))(2));
        e_[i]/=mag(e_[i])+SMALL;
    }
}





};
