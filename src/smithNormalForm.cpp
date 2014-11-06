#include "smithNormalForm.h"
#include "naive.h"
#include "probabilistic.h"
#include <iostream>
#include <armadillo>
#include <algorithm>


void
SNF::calculate_naive(arma::imat & m)
{
    std::cout << "Performing diagonalization...\n";
    diagonalize(m);
    std::cout << "Sorting diagonal...\n";
    m.diag() = arma::sort(m.diag());
    std::cout << "Ensuring divisibility...\n";
    ensure_divisibility(m);
    std::cout << "Complete!\n";
}

void
SNF::calculate_probabilistic(arma::imat & m)
{
    int val = valence(m);
    std::cout << val;
}