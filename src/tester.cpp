#include "smithNormalForm.h"
#include "storjohannTriangular.h"

#include <iostream>
#include <armadillo>
#include <chrono>
#include <stdlib.h>
#include <time.h>

typedef std::chrono::duration<int,std::micro> mu_t;

void print_wolfram_matrix(arma::imat matrix, uint size);

int main(int argc, char const *argv[])
{
    // unsigned int size = atoi(argv[1]);
    uint size = 4;

    srand(time(NULL));
    // arma::imat matrix = arma::randi<arma::imat>(size, size);

    // for (int & c : matrix) {
    //     c = c % 10;
    //     c += 1;
    // }
    std::vector<int> matrix_array = {8, 11286,  4555,  46515,
                                     0,     1, 66359, 153094,
                                     0,     0,     9,  43651,
                                     0,     0,     0,     77};
    arma::imat matrix(matrix_array.data(), size, size);
    matrix = matrix.t();

    print_wolfram_matrix(matrix, size);

    makeHermiteNormalForm(matrix);
    std::cout << matrix << std::endl;
    hermiteTriangToSNF(matrix);
    std::cout << matrix << std::endl;

    SNF snf;

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    snf.calculate_naive(matrix);

    end = std::chrono::high_resolution_clock::now();
    mu_t duration(std::chrono::duration_cast <mu_t>(end - start));

    std::cout << matrix << std::endl;
    printf("\nCalculation took %f s\n", duration.count() / 1000000.);

    return 0;
}

void print_wolfram_matrix(arma::imat matrix, uint size) {
    std::cout << "";
    for (uint i = 0; i < size; ++i) {
        std::cout << "";
        for (uint j = 0; j < size; ++j) {
            std::cout << matrix(i, j);
            if (j != size - 1) {
                std::cout << " ";
            }
        }
        std::cout << "";
        if (i != size - 1) {
            std::cout << " \n";
        }
    }
    std::cout << "\n";
}