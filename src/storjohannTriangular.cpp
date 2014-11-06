#include <armadillo>
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/math/common_factor.hpp>

#include "storjohannNumeric.h"
#include "storjohannTriangular.h"

void triangularReduction(arma::subview<arma::sword> A);

class IncorrectForm : public std::exception {
    std::string error_message;
    virtual const char* what() const throw() {
        return error_message.c_str();
    }
public:
    IncorrectForm(const std::string & message) : error_message(message) {}
};

//! check if A is upper diagonal. If not, IncorrectForm exception is thrown
void
checkUpperDiagonal(const arma::imat & A)
{
    for (uint i = 0; i < A.n_rows; ++i) {
        for (uint j = 0; j < i; ++j) {
            if (A(i, j) != 0) {
                throw IncorrectForm("Input matrix is not upper diagonal!");
            }
        }
    }
}

//! calculates multiple product of elements on diagonal
int
diagonalMultiple(const arma::diagview<int> & diag)
{
    int det = 1;
    for (uint i = 0; i < diag.n_rows; ++i) {
        det *= diag[i];
    }
    return det;
}

//! checks whether matrix A is in correct form to perform conversion to Hermite
//! normal form.
void
checkCorrectForm(const arma::imat & A)
{
    std::cout << gcdCombination(3, 5, 18) << std::endl;
    checkUpperDiagonal(A);
    // check if matrix is regular
    int det = diagonalMultiple(A.diag());
    if (0 == det) {
        throw IncorrectForm("Input matrix is singular!");
    }
    // check if off-diagonal entries are bounded by determinant
    for (uint i = 0; i < A.n_rows; ++i) {
        for (uint j = i + 1; j < A.n_cols; ++j) {
            if (A(i, j) > det) {
                throw IncorrectForm("Off-diagonal entries of input matrix are"
                                    " not bounded by its determinant!");
            }
        }
    }
}

/** \brief transform upper triangular matrix A to Hermite normal form
* \param A input upper triangular matrix
*
* Transform input matrix A to Hermite normal form. Input matrix is supposed to
* be upper triangular and all off-diagonal elements must be bounded by det(A).
* If conditions are not satisfied, exception is raised
*/
void
makeHermiteNormalForm(arma::imat & A)
{
    checkUpperDiagonal(A);
    for (int k = A.n_rows - 2; k >= 0; --k) {
        triangularReduction(A.submat(k, k, A.n_rows - 1, A.n_cols - 1));
    }
    checkCorrectForm(A);
}

//! Functor calculating positive modulo for given number e.
class PositiveModulo {
    int m_d;
public:
    PositiveModulo(int d) : m_d(d) {}
    int operator()(const int e) { return (e % m_d + m_d) % m_d; }
};

//! perform triangular reduction of input matrix T. For details please refer to
//! chapter 2, Lemma 2.
void
triangularReduction(arma::subview<arma::sword> T)
{
    // initialization
    int d = diagonalMultiple(T.submat(1, 1, T.n_rows - 1, T.n_cols - 1).diag());
    if (T(0, 0) < 0) {
        T.row(0) = (-1) * T.row(0);
    }
    T.row(0).subvec(1, T.n_cols - 1).transform(PositiveModulo(d));

    // Reduce off-diagonal entries in row 0
    for (uint j = 1; j < T.n_cols; ++j) {
        std::div_t div_result = std::div(T(0, j), T(j, j));
        T(0, j) = div_result.rem;

        if (j + 1 < T.n_cols) {
            arma::subview_row<int> sub_row = T.row(0).subvec(j+1, T.n_cols - 1);
            sub_row -= div_result.quot * T.row(j).subvec(j+1, T.n_cols - 1);
            sub_row.transform(PositiveModulo(d));
        }
    }
}
