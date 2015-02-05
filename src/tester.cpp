#include "smithNormalForm.h"
#include "storjohannTriangular.h"
#include "storjohannNumeric.h"
#include "triangularization.h"

#define ARMA_64BIT_WORD
#include <iostream>
#include <armadillo>
#include <chrono>
#include <stdlib.h>
#include <time.h>

typedef std::chrono::duration<int,std::micro> mu_t;

void print_wolfram_matrix(arma::imat matrix, uint size);

int main(int argc, char const *argv[])
{
    if (argc < 2) {
        std::cout << "./tester <matrix_size>\n";
        return 0;
    }
    uint size = atoi(argv[1]);

    srand(time(NULL));
    // for (int i = 0; i < 1; ++i) {
    while (true) {
    arma::imat matrix = arma::randi<arma::imat>(size, size + 2);


    // for (int & c : matrix) {
    //     c = c % 10;
    //     c += 1;
    // }
    // std::vector<int> matrix_array = {8, 11286,  4555,  46515,
    //                                  0,     1, 66359, 153094,
    //                                  0,     0,     9,  43651,
    //                                  0,     0,     0,     77};
    // std::vector<int> matrix_array = {2,     -6,     0,
    //                                  0,    2,     -6,
    //                                  -6,   0,     2};
    // std::vector<int> matrix_array = {2,     1, 7,
    //                                  0,     6, 12};
    // std::vector<int> matrix_array = {36,113344, 95472,  42884, 12373,12503,12303,
    //                                  0,     41,  1576,  98594, 1172,1172,11872,
    //                                  0,     0,     13,  99206, 952,94662,94192,
    //                                  0,     0,     0,   94,  770, 7780, 7030};



    // uint size = std::sqrt(matrix_array.size());
    // arma::imat matrix(matrix_array.data(), size, size);
    // matrix = matrix.t();
    matrix.transform(PositiveModulo(5));
    // float det = std::abs(arma::det(arma::conv_to<arma::mat>::from(matrix)));
    // if (0.01f > det) {
    //     continue;
    // }

    // printf("Determinant: %f\n", det);
    // std::cout << matrix << std::endl;
    // print_wolfram_matrix(matrix, size);
    triangularize(matrix);
    // float new_det = std::abs(arma::det(arma::conv_to<arma::mat>::from(matrix)));
    // if (0.01f < std::abs(det - new_det) ) {
    //     printf("Determinants differ!!! New det: %f\n", new_det);
    //     break;
    // }

    SNF snf;

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    // arma::imat naive_result      = snf.calculate_naive(matrix);
    arma::imat storjohann_result = snf.calculate_storjohann(matrix);

    end = std::chrono::high_resolution_clock::now();
    mu_t duration(std::chrono::duration_cast <mu_t>(end - start));

    // std::cout << "Naive method result:" << std::endl;
    // std::cout << naive_result << std::endl;
    std::cout << "Storjohann method: " << std::endl;
    std::cout << storjohann_result << std::endl;
    printf("\nCalculation took %f s\n", duration.count() / 1000000.);
    std::cout << "###################################################################\n";

}
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