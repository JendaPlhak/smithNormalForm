#pragma once

#define ARMA_64BIT_WORD
#include <armadillo>


/*
    Input:  general integer matrix
    Action: Calculates triangular Hermite normal form.
 */
void triangularize(arma::imat & A);
