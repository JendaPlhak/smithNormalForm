#include "probabilistic.h"

#define ARMA_64BIT_WORD
#include <vector>
#include <cmath>
#include <iostream>
#include <armadillo>

/*
    This code is altogether based on paper [1] "On Efficient Sparse Integer
    Matrix Smith Normal Form Computations"
 */

double ovalsCassiniBound(arma::imat & A);

int
valence(arma::imat & A)
{
    arma::imat B;
    if (A.n_rows < A.n_cols) {
        B = A * A.t();
    } else {
        B = A.t() * A;
    }
    double ocb = ovalsCassiniBound(A);
    return (int) ocb;
}

/*
    Ovals-of-Cassini-Bound
    Output: β ∈ R, such that for every eigenvalue λ of AA', |λ| ≤ β.
    See section 3.4 in [1]
 */
double
ovalsCassiniBound(arma::imat & A)
{
    // First calculate centers
    arma::ivec centers(A.n_rows);
    for (unsigned int i = 0; i < A.n_rows; ++i) {
        centers(i) = arma::dot(A.row(i), A.row(i));
    }

    A = arma::abs(A);

    arma::ivec ones = arma::ivec(A.n_rows, arma::fill::ones);
    arma::ivec v    = A * A.t() * ones;
    arma::ivec r    = v - centers;

    int max_center = arma::max(centers);
    int max1 = 0;
    int max2 = 0;
    for (int x : r) {
        if ( x > max1 ) {
            max1 = x;
        } else if ( x > max2 ) {
            max2 = x;
        }
    }
    return max_center + std::sqrt(max1 * max2);
}