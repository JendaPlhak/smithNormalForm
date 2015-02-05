#include <boost/math/common_factor.hpp>
#include <vector>

#define ARMA_64BIT_WORD
#include "triangularization.h"
#include "storjohannTriangular.h"
#include "storjohannNumeric.h"
#include "util.h"

template <typename Matrix_Type>
void conditioningRoutineInputCheck(Matrix_Type B);

template <typename Matrix_Type>
void conditioningRoutineOutputCheck(Matrix_Type B);

template <typename Matrix_Type>
void
columnReductionInputCheck(Matrix_Type B, const int k);

template <typename Matrix_Type>
void
columnReductionOutputCheck(Matrix_Type B, const int k);

template <typename Matrix_Type>
void hermiteTriangularFormCheck(Matrix_Type A);

void conditioningRoutine(arma::subview<arma::sword>);
void columnReduction(arma::subview<arma::sword> B, const int k, const int col2);
void reduceBetweenProfiles(arma::imat & A,
                            const std::vector<uint> & rank_profile,
                            const uint new_prof);
void reshape(arma::imat & A);

template <typename Matrix_Type>
int
mat_rank(Matrix_Type A)
{
    return arma::rank(arma::conv_to<arma::mat>::from(A));
}

template <typename Matrix_Type>
float
mat_det(Matrix_Type A)
{
    return arma::det(arma::conv_to<arma::mat>::from(A));
}

template <>
float
mat_det<>(arma::imat A)
{
    if (A.is_square()) {
        return arma::det(arma::conv_to<arma::mat>::from(A));
    } else {
        return 0;
    }
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
    float det_orig = std::abs(mat_det(A));

    // First reshape the matrix in correspondence with requirements of RowReducedEchelonForm
    uint n_rows = A.n_rows + 2;
    uint n_cols = A.n_cols + 2;

    // add additional columns and rows.
    arma::imat A_tmp(n_rows, n_cols);
    A_tmp.zeros();
    A_tmp(0,0) = 1;
    A_tmp(n_rows - 1, n_cols - 1) = 1;
    A_tmp.submat(1, 1, n_rows - 2, n_cols - 2) = A;

    std::vector<uint> rank_profile = { 0 };

    uint prev_rank = 0;
    for (uint k = 0; k < n_rows - 1; ++k) {
        D_ std::cout << "-------------------------------------------------------------------------\n";
        D_ std::cout << "Round " << k << " input:" << std::endl;
        D_ std::cout << A_tmp << std::endl;

        auto sub_A = A_tmp.submat(0, prev_rank, n_rows - 1, n_cols - 1);
        uint new_rank_offset = get_next_rank(sub_A.submat(k, 0,
                                                            sub_A.n_rows - 1,
                                                            sub_A.n_cols - 1));

        D_ std::cout << "New rank offset: " << new_rank_offset << std::endl;

        // if new_rank_offset is zero, we have reached last column of matrix.
        if (0 == new_rank_offset) {
            break;
        }

        columnReduction(sub_A, n_rows - k - 2, new_rank_offset);

        if (1 != new_rank_offset) {
            reduceBetweenProfiles(A_tmp, rank_profile,
                                    prev_rank + new_rank_offset);
        }

        prev_rank += new_rank_offset;
        rank_profile.push_back(prev_rank);

        D_ printf("Round %d result: (det = %d)\n", k, (int) mat_det(A_tmp));
        D_ std::cout << A_tmp << std::endl;
    }
    A = A_tmp.submat(1, 1, n_rows - 2, n_cols - 2);
    reshape(A);

    D_ std::cout << "Resulting matrix: " << std::endl;
    D_ std::cout << A << std::endl;
    P_ hermiteTriangularFormCheck(A);

    if (0.1f < std::abs(det_orig - std::abs(mat_det(A)))) {
        fprintf(stderr, "det_orig = %f, new_det = %f\n", det_orig, mat_det(A));
        throw HermiteTriangularError("Determinant mustn't change during triangularization!");
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

    D_ std::cout << "Ensuring principal matrix has full rank...\n";
    D_ std::cout << B << std::endl;
    // ensure that principal 2x2 matrix has full rank
    if (0 == B(0,0) * B(1,col2) - B(0,col2) * B(1,0)) {
        for (uint i = 2; i < B.n_rows; i++) {
            if (0 != B(0,0) * B(i, col2) - B(0,col2) * B(i,0)) {
                B.row(1) += B.row(i); // TODO this might be necessary to take into account later
                break;
            }
        }
    }
    D_ std::cout << "Result:\n";
    D_ std::cout << B << std::endl;
    P_ conditioningRoutineInputCheck(arma::imat(arma::join_rows(B.col(0), B.col(col2))));

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
        int_t a_tmp = PositiveModulo::mod(a / g, N);
        int_t b_tmp = PositiveModulo::mod(B(i,0) / g, N);

        int_t t = gcdCombination(a_tmp, b_tmp, N);
        if (0 == N * (a_ + t * B(i,col2)) - N_ * (a + t * B(i,0))) {
            t = -gcdCombination(a_tmp, -b_tmp, N);
        }
        B.row(1) += t * B.row(i);
    }

    P_ conditioningRoutineOutputCheck(
            arma::imat(arma::join_rows(B.col(0), B.col(col2))));
}

void
columnReduction(arma::subview<arma::sword> B, const int k, const int col2)
{
    D_ printf("columnReduction scope:\n");
    D_ std::cout << B.submat(B.n_rows - k - 2, 0, B.n_rows - 1, B.n_cols - 1)
                 << std::endl;
    D_ printf("Before conditioningRoutine: \n");
    D_ std::cout << arma::join_rows(B.col(0), B.col(col2)) << std::endl;
    conditioningRoutine(B.submat(B.n_rows - k - 2, 0,
                                 B.n_rows - 1, B.n_cols - 1),
                        col2);
    D_ printf("After conditioningRoutine: \n");
    D_ std::cout << arma::join_rows(B.col(0), B.col(col2)) << std::endl;

    P_ columnReductionInputCheck(
            arma::imat(arma::join_rows(B.col(0), B.col(col2))), k);

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
    D_ std::cout << "Before elimination:\n";
    D_ std::cout << B << std::endl;

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
    D_ std::cout << "After elimination:\n";
    D_ std::cout << B << std::endl;
    P_ columnReductionOutputCheck(
            arma::imat(arma::join_rows(B.col(0), B.col(col2))), k);
}
/**
 * Reduce columns between matrix profiles modulo profile entries
 * @param A matrix to edit
 * @param rank_profile previous rank profiles
 * @param new_prof current rank profile
 */
void
reduceBetweenProfiles(arma::imat & A,
                        const std::vector<uint> & rank_profile,
                        const uint new_prof)
{
    I_ std::cout << "Reducing elements between profiles...\n";
    uint n_rows    = rank_profile.size();
    uint prev_prof = rank_profile.back();
    // first row can be skipped since it contains always only leading 1
    for (int i = (int) n_rows - 1; i > 0; --i) {
        uint prof   = rank_profile[i];
        D_ printf("   Rank profile: %d\n", prof);

        int_t pivot = A(i, prof);
        for (uint j = prev_prof + 1; j < new_prof; ++j) {
            int_t q = floored_factor(A(i, j), pivot);
            A.col(j).subvec(0, i) -= q * A.col(prof).subvec(0, i);
        }
    }
    I_ std::cout << "   Complete!\n";
    D_ std::cout << "Result:\n";
    D_ std::cout << A << std::endl;
}
/**
 * Reorganizes columns so the resulting matrix has principal rxr of rank r,
 * where r is also rank of the original matrix and maintenances triangular shape
 * @param A matrix to reshape
 */
void
reshape(arma::imat & A)
{
    for (uint i = 0; i < std::min(A.n_rows, A.n_cols); ++i) {
        if (0 != A(i, i)) { // Diagonal entry is nonzero so continue.
            continue;
        }
        for (uint j = i + 1; j < A.n_cols; ++j) {
            if (0 != A(i, j)) {
                A.swap_cols(i, j);
                break;
            }
        }
    }
}

template <typename Matrix_Type>
void
conditioningRoutineInputCheck(Matrix_Type B)
{
    if (2 != B.n_cols) {
        throw ConditioningRoutineError("Matrix must be (k+2)x2 shaped!");
    }
    if (not (B(0,0) > 0)) {
        throw ConditioningRoutineError("Element B[0,0] must be positive!");
    }
    if (2 != mat_rank(B)) {
        std::cout << B << std::endl;
        throw ConditioningRoutineError("Matrix must have rank 2!");
    }
    if (0 == B(0,0) * B(1,1) - B(0,1) * B(1,0)) {
        throw ConditioningRoutineError("Principal matrix must have rank 2!");
    }
}

template <typename Matrix_Type>
void
conditioningRoutineOutputCheck(Matrix_Type B)
{
    using boost::math::gcd;

    if (0 == mat_det(B.submat(0, 0, 1, 1))) {
        std::cout << B.submat(0, 0, 1, 1) << std::endl;
        std::cout << B << std::endl;
        throw ConditioningRoutineError("Principal matrix must have rank 2!");
    }
    int_t x = gcd(B(1,0), B(0,0));
    for (uint i = 2; i < B.n_rows; i++) {
        int_t y = gcd(x, B(i,0));
        if (y != x) {
            throw ConditioningRoutineError("gcd(a_k,N) != gcd(a_k,N, b_1,..., b_k)");
        }
    }

}

template <typename Matrix_Type>
void
columnReductionInputCheck(Matrix_Type B, const int k)
{
    if (B(B.n_rows - k - 2, 0) <= 0) {
        throw ColumnReductionError("N = B(B.n_rows - k - 2, 0) has to be positive!");
    } else if (2 != mat_rank(B.submat(B.n_rows - k - 2, 0, B.n_rows - k - 1, 1))) {
        std::cerr << B << std::endl;
        std::cerr << B.submat(B.n_rows - k - 2, 0, B.n_rows - k - 1, 1) << std::endl;
        throw ColumnReductionError("trailing (k+2)x2 submatrix must have rank 2!");
    }
}

template <typename Matrix_Type>
void
columnReductionOutputCheck(Matrix_Type B, const int k)
{
    const int_t t1 = B(B.n_rows - k - 2, 0);
    const int_t t2 = B(B.n_rows - k - 1, 1);
    if (t1 <= 0 || t2 <= 0) {
        std::cout << B << std::endl;
        std::cerr << "t1 = " << t1 << ", "
                  << "t2 = " << t2 << std::endl;
        throw ColumnReductionError("t1 and t2 must be positive!");
    }
    for (uint i = 0; i < B.n_rows - k - 2; i++) {
        if (B(i, 0) > t1 - 1) {
            throw ColumnReductionError("Elements in first column above t1 must be bounded by t1!");
        } else if (B(i, 0) < 0) {
            throw ColumnReductionError("Elements in first column above t1 must be nonnegative!");
        }
    }
    for (uint i = 0; i < B.n_rows; i++) {
        if (i > B.n_rows - k - 2 && 0 != B(i, 0) ) {
            throw ColumnReductionError("Elements in first column under t1 must be zero!");
        } else if (i == B.n_rows - k - 1) {
            continue;
        } else if (B(i, 1) > t2 - 1) {
            throw ColumnReductionError("Elements in second column above must be bounded by t2!");
        } else if (B(i, 1) < 0) {
            throw ColumnReductionError("Elements in second column must be nonnegative!");
        }
    }
}

template <typename Matrix_Type>
void hermiteTriangularFormCheck(Matrix_Type A)
{
    for (uint j = 0; j < A.n_cols; ++j) {
        for (uint i = 0; i < A.n_rows; ++i) {
            if (A(i, j) < 0) {
                throw HermiteTriangularError("Hermite form requires elements to be positive!");
            }
            else if (i > j) {
                if (A(i, j) != 0) {
                    throw HermiteTriangularError("Hermite form requires matrix to be triangular!");
                }
            } else if (i < j && j < A.n_rows && 0 != A(j, j)) {
                if (A(j, j) < A(i, j)) {
                    throw HermiteTriangularError("Hermite form requires off diagonal entries to be bounded by diagonal entry!");
                }
            }
        }
    }
}