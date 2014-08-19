#ifndef SNF_SMITH_NORMAL_FORM_H
#define SNF_SMITH_NORMAL_FORM_H

#include <armadillo>


class SNF {
public:
    void calculate_naive(arma::imat & m);
    void calculate_probabilistic(arma::imat & m);
};

#endif // SNF_SMITH_NORMAL_FORM_H