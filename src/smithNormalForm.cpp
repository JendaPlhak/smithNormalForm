#include "smithNormalForm.h"
#include "naive.h"
#include "probabilistic.h"
#include "storjohannTriangular.h"
#include <iostream>
#include <armadillo>
#include <algorithm>


arma::imat
SNF::calculate_naive(arma::imat m)
{
    std::cout << "Performing diagonalization...\n";
    diagonalize(m);
    std::cout << "Sorting diagonal...\n";
    m.diag() = arma::sort(m.diag());
    std::cout << "Ensuring divisibility...\n";
    ensure_divisibility(m);
    std::cout << "Complete!\n";
    return m;
}

arma::imat
SNF::calculate_storjohann(arma::imat matrix)
{
    makeHermiteNormalForm(matrix);
    hermiteTriangToSNF(matrix);
    return matrix;
}

void
SNF::calculate_probabilistic(arma::imat & m)
{
    int val = valence(m);
    std::cout << val;
}