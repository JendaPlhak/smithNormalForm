#include "smithNormalForm.h"
#include "numeric.h"
#include "util.h"

#define ARMA_64BIT_WORD
#include <iostream>
#include <armadillo>
#include <chrono>
#include <exception>
#include <stdlib.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <NTL/mat_ZZ.h>

typedef std::chrono::duration<int,std::micro> mu_t;

void print_wolfram_matrix(arma::imat matrix, uint size);
arma::imat loadMatrixFromFile(const std::string file_path);
void saveMatrixToFile(const std::string file_path, const arma::imat & mat);


// Define simple matrix size structure
struct Size {
    uint n_rows;
    uint n_cols;

    Size(const uint rows, const uint cols) : n_rows(rows), n_cols(cols) {}
    // operator= will be used to assign to the vector
    Size& operator=(const std::string& str)
    {
    std::istringstream iss(str);
    if (!(iss >> n_rows >> n_cols))
        throw TCLAP::ArgParseException(str + " doesn't define size!");

    return *this;
    }
};

// Create an ArgTraits for the Size type that declares it to be
// of string like type
namespace TCLAP {
template<>
struct ArgTraits<Size> {
    typedef StringLike ValueCategory;
};
}

int main(int argc, char const *argv[])
{
    try {
    // srand(10);
    srand(time(NULL));

    TCLAP::CmdLine cmd("SNF library tester", ' ', "1.0");

    // Define a value argument and add it to the command line.
    TCLAP::SwitchArg randArg("r","random",
                                "Generates random matrix and calculates its SNF",
                                false);
    TCLAP::ValueArg<Size> sizeArg("s", "size",
                                    "Specifies size of random matrix",
                                    false, Size(5, 5), "rows columns", cmd);
    TCLAP::ValueArg<std::string> fileArg("f", "file",
                                            "loads matrix from file",
                                            false, "", "file_path");
    TCLAP::ValueArg<std::string> saveArg("", "save",
                                            "Saves matrix to specified file",
                                            false, "", "file_path", cmd);

    TCLAP::MultiSwitchArg verboseArg("v", "",
                                        "Specifies verbosity level. -v will"
                                        " print only basic info, -vv will produce"
                                        " step by step partial results",
                                        cmd, 0);

    // We want either generate random matrix or load matrix from file.
    cmd.xorAdd(randArg, fileArg);
    // Parse the argv array.
    cmd.parse(argc, argv);

    arma::imat matrix;
    int count = 0;

    if (randArg.isSet()) {
CYCLE:;
        // srand(count);
        printf("Round: %d\n", count);
        count += 1;
        Size size = sizeArg.getValue();

        matrix = arma::randi<arma::imat>(size.n_rows, size.n_cols);
        matrix.transform(PositiveModulo(10));
    } else if (fileArg.isSet()) {
        matrix = loadMatrixFromFile(fileArg.getValue());
    }

    SNF snf;

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();

    // float det = arma::det(arma::conv_to<arma::mat>::from(matrix));

    // arma::imat naive_result      = snf.calculate_naive(matrix);
    NTL::mat_ZZ storjohann_result = snf.calculate(matrix);

    end = std::chrono::high_resolution_clock::now();
    mu_t duration(std::chrono::duration_cast <mu_t>(end - start));

    // std::cout << "Naive method result:" << std::endl;
    // std::cout << naive_result << std::endl;
    P_ std::cout << "Storjohann method: " << std::endl;
    P_ std::cout << storjohann_result << std::endl;
    // printf("Matrix determinant is: %f\n", det);
    printf("\nCalculation took %f s\n", duration.count() / 1000000.);
    P_ std::cout << "###################################################################\n";

    // save file to file
    if (saveArg.isSet()) {
        saveMatrixToFile(saveArg.getValue(), matrix);
    }
    if (randArg.isSet()) {
        // goto CYCLE;
    }

    } catch (TCLAP::ArgException &e) {  // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
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

arma::imat
loadMatrixFromFile(const std::string file_path)
{
    std::ifstream file(file_path);
    if (not file.is_open()) {
        throw std::runtime_error("Failed to open specified file!");
    }

    arma::imat mat;
    mat.load(file, arma::raw_ascii);
    return mat;
}

void
saveMatrixToFile(const std::string file_path, const arma::imat & mat)
{
    std::ofstream file(file_path);
    if (not file.is_open()) {
        throw std::runtime_error("Failed to open specified file!");
    }
    mat.save(file, arma::raw_ascii);
}