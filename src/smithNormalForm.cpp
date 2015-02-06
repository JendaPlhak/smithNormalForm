#include "smithNormalForm.h"
#include "probabilistic.h"
#include "storjohannTriangular.h"
#include "util.h"

#define ARMA_64BIT_WORD
#include <iostream>
#include <armadillo>
#include <algorithm>

arma::imat
SNF::calculate_storjohann(arma::imat matrix)
{
    I_ std::cout << "Input matrix: \n" << matrix << std::endl;
    // std::cout << "Making Hermite normal form...\n";
    // makeHermiteNormalForm(matrix);
    // std::cout << "Result: \n" << matrix << std::endl;
    I_ std::cout << "Stripping zero rows...\n";
    uint n_stripped = stripZeroRows(matrix);
    I_ std::cout << "   Done!...\n";
    I_ std::cout << "Converting Hermite matrix to SNF...\n";
    hermiteTriangToSNF(matrix);
    I_ std::cout << "Result: \n" << matrix << std::endl;
    I_ std::cout << "Eliminating extra columns...\n";
    eliminateExtraColumns(matrix);
    I_ std::cout << "Result: \n" << matrix << std::endl;
    I_ std::cout << "Reducing final non-trivial square matrix to SNF...\n";
    D_ std::cout << "Input matrix: \n"
                 << matrix << std::endl;
    reduceResultingSquareToSNF(matrix);
    I_ std::cout << "Result: \n" << matrix << std::endl;

    // give matrix its original size
    matrix.resize(matrix.n_rows + n_stripped, matrix.n_cols);
    return matrix;
}

void
SNF::calculate_probabilistic(arma::imat & m)
{
    int val = valence(m);
    std::cout << val;
}