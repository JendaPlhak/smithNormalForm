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
    naive::diagonalize(m);
    std::cout << "Sorting diagonal...\n";
    m.diag() = arma::sort(m.diag());
    std::cout << "Ensuring divisibility...\n";
    naive::ensure_divisibility(m);
    std::cout << "Complete!\n";
    return m;
}

arma::imat
SNF::calculate_storjohann(arma::imat matrix)
{
    std::cout << "Input matrix: \n" << matrix << std::endl;
    std::cout << "Performing naive row reduction...";
    naive::rowReduce(matrix);
    std::cout << "Result: \n" << matrix << std::endl;
    std::cout << "Making Hermite normal form...\n";
    makeHermiteNormalForm(matrix);
    std::cout << "Result: \n" << matrix << std::endl;
    std::cout << "Converting Hermite matrix to SNF...\n";
    hermiteTriangToSNF(matrix.submat(0, 0, matrix.n_rows-1, matrix.n_cols-1));
    std::cout << "Result: \n" << matrix << std::endl;
    std::cout << "Eliminating extra columns...\n";
    eliminateExtraColumns(matrix);
    std::cout << "Result: \n" << matrix << std::endl;
    std::cout << "Reducing final non-trivial square matrix to SNF...\n";
    reduceResultingSquareToSNF(matrix);
    std::cout << "Result: \n" << matrix << std::endl;
    return matrix;
}

void
SNF::calculate_probabilistic(arma::imat & m)
{
    int val = valence(m);
    std::cout << val;
}