#include "smithNormalForm.h"

#include <iostream>
#include <armadillo>
#include <chrono>
#include <stdlib.h>
#include <time.h>

typedef std::chrono::duration<int,std::micro> mu_t;



int main(int argc, char const *argv[])
{   
    unsigned int size = atoi(argv[1]);

    srand(time(NULL));
    arma::imat matrix = arma::randi<arma::imat>(size, size);    

    for (int & c : matrix) {
        c %= 4; 
    }  

    SNF snf;

    std::cout << matrix << std::endl;
    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    snf.calculate_probabilistic(matrix);

    end = std::chrono::high_resolution_clock::now();
    mu_t duration(std::chrono::duration_cast <mu_t>(end - start));

    // std::cout << matrix << std::endl;
    printf("\nCalculation took %f s\n", duration.count() / 1000000.);

    return 0;
}