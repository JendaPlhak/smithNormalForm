#include <boost/math/common_factor.hpp>
#include <vector>

#define ARMA_64BIT_WORD
#include "triangularization.h"
#include "storjohannTriangular.h"
#include "storjohannNumeric.h"

void conditioningRoutine(arma::subview<arma::sword>);
void columnReduction(arma::subview<arma::sword> B, const int k, const int col2);

template <typename Matrix_Type>
int
imat_rank(Matrix_Type A)
{
    return arma::rank(arma::conv_to<arma::mat>::from(A));
}

template <typename Matrix_Type>
float
imat_det(Matrix_Type A)
{
    return arma::det(arma::conv_to<arma::mat>::from(A));
}

uint
get_next_rank(arma::subview<arma::sword> A)
{
    const int_t N = A(0,0);
    uint j = 1;
    for (; j < A.n_cols; ++j) {
        const int_t N_ = A(0,j);
        for (uint i = 1; i < A.n_rows; ++i) {
            if (N * A(i, j) - N_ * A(i, 0) != 0) {
                goto PROFILE_FOUND;
            }
        }
    }
    return j - 1;
PROFILE_FOUND:;
    return j;
}

void
triangularize(arma::imat & A)
{
    float det_orig = std::abs(imat_det(A));

    // First reshape the matrix in correspondence with requirements of RowReducedEchelonForm
    uint n_rows = A.n_rows + 2;
    uint n_cols = A.n_cols + 2;
    gcdCombination(12, 13, 405);

    arma::imat A_tmp(n_rows, n_cols);
    A_tmp.zeros();
    A_tmp(0,0) = 1;
    A_tmp(n_rows - 1, n_cols - 1) = 1;
    A_tmp.submat(1, 1, n_rows - 2, n_cols - 2) = A;

    uint prev_rank = 0;
    for (uint k = 0; k < n_rows - 2; ++k) {
        // std::cout << "-------------------------------------------------------------------------\n";
        // std::cout << "Round " << k << " input:" << std::endl;
        // std::cout << A_tmp << std::endl;

        auto sub_A = A_tmp.submat(0, prev_rank, n_rows - 1, n_cols - 1);

        uint new_rank_offset = get_next_rank(sub_A.submat(k, 0,
                                                            sub_A.n_rows - 1,
                                                            sub_A.n_cols - 1));
        // std::cout << "New rank offset: " << new_rank_offset << std::endl;
        columnReduction(sub_A, n_rows - k - 2, new_rank_offset);
        prev_rank += new_rank_offset;

        // printf("Round %d result: (det = %d)\n", k, imat_det(A_tmp));
        // std::cout << A_tmp << std::endl;
    }
    A = A_tmp.submat(1, 1, n_rows - 2, n_cols - 2);

    if (0.1f < std::abs(det_orig - std::abs(imat_det(A)))) {
        fprintf(stderr, "det_orig = %f, new_det = %f\n", det_orig, imat_det(A));
        throw IncorrectForm("Determinant mustn't change during triangularization!");
    }
}

template <typename Matrix_Type>
void
conditioningRoutineInputCheck(Matrix_Type B)
{
    if (2 != B.n_cols) {
        throw IncorrectForm("Matrix must be (k+2)x2 shaped!");
    }
    if (not (B(0,0) > 0)) {
        throw IncorrectForm("Element B[0,0] must be positive!");
    }
    if (2 != imat_rank(B)) {
        std::cout << B << std::endl;
        throw IncorrectForm("Matrix must have rank 2!");
    }
    if (0 == B(0,0) * B(1,1) - B(0,1) * B(1,0)) {
        throw IncorrectForm("Principal matrix must have rank 2!");
    }
}

template <typename Matrix_Type>
void
conditioningRoutineOutputCheck(Matrix_Type B)
{
    using boost::math::gcd;

    if (0 == imat_det(B.submat(0, 0, 1, 1))) {
        std::cout << B << std::endl;
        throw IncorrectForm("Principal matrix must have rank 2!");
    }
    int_t x = gcd(B(1,0), B(0,0));
    for (uint i = 2; i < B.n_rows; i++) {
        int_t y = gcd(x, B(i,0));
        if (y != x) {
            throw IncorrectForm("gcd(a_k,N) != gcd(a_k,N, b_1,..., b_k)");
        }
    }

}

/** \brief Implements Theorem 5 - The conditioning routine
* \param (k+2)x2 matrix subview
*
*/
void
conditioningRoutine(arma::subview<arma::sword> B, const uint col2)
{
    // Element B(0,0) has to be positive
    if (not (B(0,0) > 0)) {
        B.row(0) *= -1;
    }

    // std::cout << "Ensuring principal matrix to has full rank...\n";
    // std::cout << B << std::endl;
    // ensure that principal 2x2 matrix has full rank
    if (0 == B(0,0) * B(1,col2) - B(0,col2) * B(1,0)) {
        for (uint i = 2; i < B.n_rows; i++) {
            if (0 != B(0,0) * B(i, col2) - B(0,col2) * B(i,0)) {
                B.row(1) += B.row(i); // TODO this might be necessary to take into account later
                break;
            }
        }
    }
    // std::cout << "Result:\n";
    // std::cout << B << std::endl;
    conditioningRoutineInputCheck(arma::imat(arma::join_rows(B.col(0), B.col(col2))));

    if (B.n_rows == 2) { // There is no work to be done
        return;
    }

    const int_t & N  = B(0,0);
    const int_t & N_ = B(0,col2);
    int_t & a  = B(1,0);
    int_t & a_ = B(1,col2);
    for (uint i = 2; i < B.n_rows; i++) {
        using boost::math::gcd;

        int_t g = gcd(a, B(i,0));
        if (g == 0) { // a and B(i,0) must be 0, therefore no action needed.
            continue;
        }
        int_t a_tmp = (a / g)      % N;
        int_t b_tmp = (B(i,0) / g) % N;

        int_t t = gcdCombination(a_tmp, b_tmp, N);
        if (0 == N * (a_ + t * B(i,col2)) - N_ * (a + t * B(i,0))) {
            t = -gcdCombination(a_tmp, -b_tmp, N);
        }
        B.row(1) += t * B.row(i);
    }

    conditioningRoutineOutputCheck(arma::imat(arma::join_rows(B.col(0), B.col(col2))));
}

void
columnReductionOutputCheck(arma::subview<arma::sword> B,
                            const int k,
                            const int col2)
{
    const int_t t1 = B(B.n_rows - k - 2, 0);
    const int_t t2 = B(B.n_rows - k - 1, col2);
    if (t1 <= 0 || t2 <= 0) {
        std::cout << B << std::endl;
        std::cerr << "t1 = " << t1 << ", "
                  << "t2 = " << t2 << std::endl;
        throw IncorrectForm("t1 and t2 must be positive!");
    }
    for (uint i = 0; i < B.n_rows - k - 2; i++) {
        if (B(i, 0) > t1 - 1) {
            throw IncorrectForm("Elements in first column above t1 must be bounded by t1!");
        } else if (B(i, 0) < 0) {
            throw IncorrectForm("Elements in first column above t1 must be nonnegative!");
        }
    }
    for (uint i = 0; i < B.n_rows; i++) {
        if (i > B.n_rows - k - 2 && 0 != B(i, 0) ) {
            throw IncorrectForm("Elements in first column under t1 must be zero!");
        } else if (i == B.n_rows - k - 1) {
            continue;
        } else if (B(i, col2) > t2 - 1) {
            throw IncorrectForm("Elements in second column above must be bounded by t2!");
        } else if (B(i, col2) < 0) {
            throw IncorrectForm("Elements in second column must be nonnegative!");
        }
    }
}

void
columnReduction(arma::subview<arma::sword> B, const int k, const int col2)
{
    // printf("Before conditioningRoutine: \n");
    // std::cout << B << std::endl;
    conditioningRoutine(B.submat(B.n_rows - k - 2, 0, B.n_rows - 1, B.n_cols - 1), col2);
    // printf("After conditioningRoutine: \n");
    // std::cout << A_tmp << std::endl;

    uint offset = B.n_rows - k - 2;
    int_t m1, m2, t1;

    const int_t N = B(offset, 0);
    const int_t a = B(offset + 1, 0);

    extendedGCD(m1, m2, t1, N, a, true); // gcd(N, a_k)

    arma::irowvec vec = B.row(offset);
    B.row(offset)     = m1 * B.row(offset) + m2       * B.row(offset + 1);
    B.row(offset + 1) = (-a / t1) * vec    + (N / t1) * B.row(offset + 1);


    if (B(offset, 0) < 0) { // ensure t1 to be positive
        B.row(offset) *= -1;
        t1            *= -1;
    }
    if (B(offset + 1, col2) < 0) { // ensure t2 to be positive
        B.row(offset + 1) *= -1;
    }
    // std::cout << "Before elimination:\n";
    // std::cout << B << std::endl;

    const int_t t2 = B(offset + 1, col2);
    // Reduce first and second column
    for (uint i = 0; i < B.n_rows; i++) {
        if (i == offset + 1) {
            continue;
        } if (i == offset) {
            if (t2 != 0) {
                B.row(i) -= floored_factor(B(i, col2), t2) * B.row(offset + 1);
            }
        } else {
            B.row(i) -= floored_factor(B(i, 0), t1) * B.row(offset);
            if (t2 != 0) {
                B.row(i) -= floored_factor(B(i, col2), t2) * B.row(offset + 1);
            }
        }
    }
    // std::cout << "After elimination:\n";
    // std::cout << B << std::endl;
    columnReductionOutputCheck(B, k, col2);
}