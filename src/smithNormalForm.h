#ifndef SNF_SMITH_NORMAL_FORM_H
#define SNF_SMITH_NORMAL_FORM_H

#define ARMA_64BIT_WORD
#include <armadillo>


class SNF {
public:
    arma::imat calculate_naive(arma::imat m);
    arma::imat calculate_storjohann(arma::imat m);
    void calculate_probabilistic(arma::imat & m);
};

#endif // SNF_SMITH_NORMAL_FORM_H