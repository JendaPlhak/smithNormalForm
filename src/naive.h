#ifndef SNF_NAIVE_H
#define SNF_NAIVE_H

#include <armadillo>

void diagonalize(arma::imat & m);
void ensure_divisibility(arma::imat & m);
void make_gcd(arma::subview<arma::sword> m);
void qsort_diagonal(arma::imat & m);

#endif // SNF_NAIVE_H