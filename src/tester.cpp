#include "smithNormalForm.h"

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <chrono>
#include <stdlib.h>
#include <time.h>

typedef std::chrono::duration<int,std::micro> mu_t;



int main(int argc, char const *argv[])
{   
    unsigned int size = atoi(argv[1]);
    Eigen::MatrixXi matrix(size, size);
    // Eigen::MatrixXi matrix(3, 3);
    srand(time(NULL));
    matrix.setRandom();
    // matrix <<  2,  4, 4,
    //           -6,  6, 12,
    //            10,-4,-16;
           

    for (unsigned int i = 0; i < size; ++i) {
        for (unsigned int j = 0; j < size; ++j)
        {
            matrix(i, j) %= 4;
        }
    }
    // std::cout << matrix << std::endl;
    SNF snf;

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    snf.calculate_naive(matrix);

    end = std::chrono::high_resolution_clock::now();
    mu_t duration(std::chrono::duration_cast <mu_t>(end - start));

    // std::cout << matrix << std::endl;
    printf("\nCalculation took %f s\n", duration.count() / 1000000.);

    return 0;
}