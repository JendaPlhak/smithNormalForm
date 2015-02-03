#ifndef SNF_NAIVE_H
#define SNF_NAIVE_H

#define ARMA_64BIT_WORD
#include <armadillo>

namespace naive
{
void rowReduce(arma::imat & m);
void diagonalize(arma::imat & m);
void ensure_divisibility(arma::imat & m);
void make_gcd(arma::subview<arma::sword> m);
void qsort_diagonal(arma::imat & m);
}

#endif // SNF_NAIVE_H