#include <armadillo>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <boost/math/common_factor.hpp>

#include "storjohannNumeric.h"
#include "storjohannTriangular.h"

void triangularReduction(arma::subview<arma::sword> A);
void eliminateTrailingCol(arma::subview<arma::sword> T);
void ensureDivisibility(arma::subview<arma::sword> T);
void submatrixToSNF(arma::subview<arma::sword> T);
void processRow(arma::subview<arma::sword> T);
void eliminateExtraColumnsInFirstRow(arma::subview<arma::sword> T);

void
checkSNF(const arma::subview<arma::sword> T)
{
    for (uint i = 0; i < T.n_rows; ++i) {
        for (uint j = 0; j < T.n_rows; ++j) {
            if (i == j) {
                if (i >= 1 && T(i, i) % T(i-1, i-1) != 0) {
                    throw IncorrectForm("Diagonal entries don't form dividing "
                                        "sequence!");
                }
            } else {
                if (0 != T(i, j)) {
                    throw IncorrectForm("Matrix has non-zero off diagonal "
                                        "entries!");
                }
            }
        }
    }
}

void
checkConditionsTheorem6(const arma::subview<arma::sword> T)
{
    uint k = T.n_rows;
    arma::mat T_mat = arma::conv_to<arma::mat>::from(
                        T.submat(0, 0, T.n_rows - 1, T.n_rows - 1));

    if (k != arma::rank(T_mat)) {
        throw IncorrectForm("First k columns of matrix T don't have rank k!");
    }
    // first k-1 columns of T must be in Smith normal form
    try { checkSNF(T.submat(0, 0, k - 2, k - 2)); }
    catch (IncorrectForm & e) {
        std::stringstream s;
        s << "Theorem 6 requires principal (k-1)th matrix to be in SNF but "
          << e.str();
        throw IncorrectForm(s.str());
    }
    // off-diagonal entries in rows 0, 1,..., k-2 are reduced modulo the
    // diagonal entry in the same row.
    for (uint row = 0; row < k - 2; ++row) {
        for (uint col = row + 1; col < T.n_cols; ++col) {
            std::stringstream s;
            s << "Theorem 6 requires off-diagonal entries in rows 0,.., k-2 "
              << "to be reduced modulo the diagonal entry in the same row but ";
            if (T(row, col) < 0) {
                s << "T(" << row <<", " <<col<< ") = " << T(row, col) << " < 0";
                throw IncorrectForm(s.str());
            } else if (T(row, row) <= T(row, col)) {
                std::cout << T << std::endl;
                s << "T(" << row <<", " << col << ") = " << T(row, col)
                  << " >= " << T(row, col);
                throw IncorrectForm(s.str());
            }
        }
    }
    // off-diagonal entries in row k-1 must be bounded in magnitude by D, a
    // positive multiple of the determinant of the principal kth submatrix of T
    int determinant = std::abs(arma::det(T_mat));
    for (uint col = k; col < T.n_cols; ++col) {
        if (determinant < T(k - 1, col)) {
            throw IncorrectForm("off-diagonal entries in row k-1 must be"
                "bounded in magnitude by D, a positive multiple of the "
                "determinant of the principal kth submatrix of T");
        }
    }
}

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

//! checks whether matrix A is in correct form to perform conversion to Hermite
//! normal form.
void
checkCorrectFormHermiteTransform(const arma::imat & A)
{
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
    checkCorrectFormHermiteTransform(A);
}

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
    for (uint j = 1; j < T.n_rows; ++j) {
        std::div_t div_result = std::div(T(0, j), T(j, j));
        T(0, j) = div_result.rem;

        if (j + 1 < T.n_cols) {
            arma::subview_row<int> sub_row = T.row(0).subvec(j+1, T.n_cols - 1);
            sub_row -= div_result.quot * T.row(j).subvec(j+1, T.n_cols - 1);
            sub_row.transform(PositiveModulo(d));
        }
    }
}

void
hermiteTriangToSNF(arma::subview<arma::sword> A)
{
    for (uint i = 1; i < A.n_rows; ++i) {
        submatrixToSNF(A.submat(0, 0, i, A.n_cols - 1));
    }
}

void
submatrixToSNF(arma::subview<arma::sword> T)
{
    ensureDivisibility(T);
    eliminateTrailingCol(T);

    // check invariant
    try { checkSNF(T); }
    catch (IncorrectForm & e) {
        std::cerr << T << std::endl;
        throw IncorrectForm("Function: submatrixToSNF, Error: " + e.str());
    }
}

//! Implements Lemma 7 - ensures that
//! gcd(a_i, t_i) = gcd(a_i, t_i, t_i+1, ..., t_k) for 0 <= i <= k-1
void
ensureDivisibility(arma::subview<arma::sword> T)
{
    int k = T.n_rows;
    arma::subview_col<int> t = T.col(k - 1);

    for (int i = k - 2; i >= 0; --i) {
        int a = T.diag()[i];
        int c = gcdCombination(t[i], t[i + 1], a);

        arma::subview_row<int> sub_row = T.row(i).subvec(k - 1, T.n_cols - 1);
        sub_row += c * T.row(i+1).subvec(k - 1, T.n_cols - 1);
        sub_row.transform(PositiveModulo(a));
    }

    // check invariant
    int gcd_t = t[k - 1];
    for (int i = k - 2; i >= 0; --i) {
        using boost::math::gcd;
        gcd_t = gcd(gcd_t, t[i]);
        if (gcd(t[i], T.diag()[i]) != gcd(gcd_t, T.diag()[i])) {
            throw std::logic_error("Invariant of Lemma 7 doesn't hold!");
        }
     }
}

//! implements Lemma 9
void
eliminateTrailingCol(arma::subview<arma::sword> T)
{
    // std::cout << "Eliminating trailing column\n";
    for (uint i = 0; i < T.n_rows - 1; ++i) {
        // std::cout << "Processing row " << i << std::endl;
        // std::cout << T.submat(i, i, T.n_rows - 1, T.n_cols - 1) << std::endl;
        processRow(T.submat(i, i, T.n_rows - 1, T.n_cols - 1));
    }
    // std::cout << "Resulting sub-matrix is: \n";
    // std::cout << T << std::endl;
}

//! process first row of given view and convert it to form required by Lemma 9
void
processRow(arma::subview<arma::sword> T)
{
    const uint k = T.n_rows;
    arma::subview_col<int> t_col = T.col(k - 1);
    arma::diagview<int>    diag  = T.diag();
    int s = 0, t = 0, s1 = 0;
    extendedGCD(s, t, s1, diag[0], t_col[0]);

    t_col[0] = 0;
    for (uint i = 1; i < k; ++i) {
        int q = (t * t_col[i] / s1) % diag[i];
        // std::cout << "    q = " << q << std::endl;
        // printf("Before mod:\n");
        // std::cout << T << std::endl;
        if (k < T.n_cols) {
            arma::subview_row<int> sub_row = T.row(i).subvec(k, T.n_cols - 1);
            sub_row -= q * T.row(0).subvec(k, T.n_cols - 1);
        }
        // std::cout << "    t[i] = " << t_col[i] << ", diag[i] = " << diag[i]
                  // << ", s1 = " << s1 << std::endl;
        // printf("After mod:\n");
        // std::cout << T << std::endl;
        t_col[i] = t_col[i] * diag[0] / s1;
        if (k <= T.n_cols) {
            if (k - 1 == i) {
                if (k < T.n_cols) {
                    arma::subview_row<int> sub_row = T.row(i).subvec(k, T.n_cols - 1);
                    sub_row.transform(PositiveModulo(diag[i]));
                }
            } else {
                arma::subview_row<int> sub_row = T.row(i).subvec(k - 1, T.n_cols - 1);
                sub_row.transform(PositiveModulo(diag[i]));
            }

        }
    }
    diag[0] = s1;

    if (k <= T.n_cols) {
        T.row(0).subvec(k - 1, T.n_cols - 1).transform(PositiveModulo(s1));
    }
    // check invariant of Lemma 9
    // Resulting matrix has to satisfy all 4 conditions of Theorem 6
    checkConditionsTheorem6(T);

    // T(0,0) has to divide all entries in the principal kth submatrix of T
    for (uint row = 0; row < k; ++row) {
        for (uint col = 0; col < k; ++col) {
            if (T(row, col) % T(0, 0) != 0) {
                std::stringstream e;
                e << "T(" << row << "," << col << ") = " << T(row, col)
                  << " is not divisible by T(0,0) = " << T(0, 0) << "!";
                throw IncorrectForm(e.str());
            }
        }
    }
    // T(0, k-1) has to be zero - thats what this function is supposed to do.
    if (0 != T(0, k-1)) {
        throw IncorrectForm("T(0, k-1) != 0!");
    }
}

//! Implements Theorem 10 - eliminate columns k,.., m
void
eliminateExtraColumns(arma::imat & T)
{
    uint k = T.n_rows;
    uint m = T.n_cols;
    bool identity = true;
    for (uint i = 0; i < k; ++i) {
        // For first identity block of matrix T no need for elimination.
        if (identity && 1 == T.diag()[i]) {
            continue;
        } else {
            identity = false;
            eliminateExtraColumnsInFirstRow(T.submat(i, i, k - 1, m - 1));
        }
    }
}

//! Process first row and eliminate all extra entries on indexes k,.., m
void
eliminateExtraColumnsInFirstRow(arma::subview<arma::sword> T)
{
    uint k = T.n_rows;
    uint m = T.n_cols;
    for (uint j = k; j < m; ++j) {
        int s = 0; int t = 0; int s1 = 0;
        extendedGCD(s, t, s1, T(0, 0), T(0, j));

        // initialize
        int a   = (-1) * T(0, j) / s1;
        int b   = T(0, 0) / s1;
        T(0, 0) = s1;
        T(0, j) = 0;

        // TODO implement using column-wise operations!
        for (uint i = 1; i < k; ++i) {
            int C   = (s * T(i, 0) + t * T(i, j)) % T.diag()[i];
            T(i, j) = (a * T(i, 0) + b * T(i, j)) % T.diag()[i];
            T(i, 0) = C;
        }
    }
}

void
reduceResultingSquareToSNF(arma::imat & T)
{
    for (uint i = 0; i < T.n_rows; ++i) {
        if (1 != T.diag()[i]) {
            // transpose the nontrivial matrix
            T.submat(i, i, T.n_rows - 1, T.n_rows - 1) =
                T.submat(i, i, T.n_rows - 1, T.n_rows - 1).t();
            hermiteTriangToSNF(T.submat(i, i, T.n_rows - 1, T.n_rows - 1));
            break;
        }
    }
    checkSNF(T.submat(0, 0, T.n_rows - 1, T.n_cols - 1));
}
