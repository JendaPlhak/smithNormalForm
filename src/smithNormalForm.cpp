#include "smithNormalForm.h"
#include "probabilistic.h"
#include "storjohannTriangular.h"
#include "triangularization.h"
#include "util.h"
#include "numeric.h"


#include <iostream>
#include <armadillo>
#include <algorithm>
#include <math.h>
#include "matrix.h"
#include <vector>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

arma::imat
SNF::calculate(arma::imat matrix)
{
    // generate primes for moduli system
    std::vector<NTL::ZZ> primes = get_primes(matrix);

    std::vector<uint> rank_profile = get_rank_profile(matrix, primes);
    I_ std::cout << "Rank profile: " << std::endl;
    for (uint i = 0; i < rank_profile.size(); ++i) {
        I_ std::cout << rank_profile[i] << " ";
    }
    I_ std::cout << std::endl;
    for (const NTL::ZZ & p : primes) {
        IMat a(arma::imat(matrix), ZZ_to_int_t(p));
        calculate_storjohann(a, rank_profile, ZZ_to_int_t(p));
    }
    return matrix;
}


// before mod implemented  0.020107 for 10x10
void
SNF::calculate_storjohann(IMat & matrix,
                            std::vector<uint> & rank_profile,
                            int_t p)
{
    I_ std::cout << "Input matrix: \n" << matrix << std::endl;
    std::cout << "Making Hermite triangular normal form...\n";
    // makeHermiteNormalForm(matrix);
    triangularize(matrix, rank_profile, p);
    std::cout << "Result: \n" << matrix << std::endl;

    I_ std::cout << "Stripping zero rows...\n";
    uint n_stripped = stripZeroRows(matrix);
    I_ std::cout << "Result: \n" << matrix << std::endl;

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
}

void
SNF::calculate_probabilistic(arma::imat & m)
{
    int val = valence(m);
    std::cout << val;
}

std::vector<NTL::ZZ>
SNF::get_primes(const arma::imat & matrix)
{
    long A_max = std::max(matrix.max(), std::abs(matrix.min()));
    long m = matrix.n_cols;
    long l = 6 + positive_log(1 + m*std::floor(positive_log(std::sqrt(m) * A_max)));

    NTL::ZZ b;
    NTL::power(b, std::ceil(std::sqrt(m)) * A_max, m);
    long s;
    {
        using namespace std;
        long denom = ceil(positive_log(2) + m * positive_log(sqrt(m) * A_max));
        s = 2 * ceil(denom / (double) (l - 1));
    }

    NTL::ZZ total(NTL::INIT_VAL, 1);
    // we will start the search from 2^(l-1)
    NTL::ZZ prev_prime(NTL::INIT_VAL, 1 << (l - 2));

    std::vector<NTL::ZZ> primes;
    if (b <= total) {
        primes.push_back(NextPrime(prev_prime + 1));
    } else {
        while (total < b) {
            prev_prime = NextPrime(prev_prime + 1);
            total *= prev_prime;
            primes.push_back(prev_prime);
        }
    }

    I_ std::cout << "Number of primes: "      << primes.size() << std::endl;
    I_ std::cout << "Max moduli bit length: " << l << std::endl;
    return primes;
}

std::vector<uint>
SNF::get_rank_profile(const arma::imat & matrix,
                        const std::vector<NTL::ZZ> & primes)
{
    // force NTL::ZZ_p to use mod
    NTL::mat_ZZ_p ntl_matrix;
    std::vector<long> ranks;
    std::vector<std::vector<uint>> rank_profiles;

    long A_max = std::max(matrix.max(), std::abs(matrix.min()));
    NTL::ZZ p_multi(NTL::INIT_VAL, 1);
    NTL::ZZ i_bound;

    long rank = 0;
    int s_i   = 0;
    for (long i = 1; i <= (long) primes.size(); ++i) {
        NTL::power(i_bound, std::ceil(std::sqrt(i)) * A_max, i);

        while (p_multi <= i_bound) {
            NTL::ZZ_p::init(primes[s_i]);
            matrix_convert(matrix, ntl_matrix);

            rank_profiles.push_back(get_rank_profile_mod(ntl_matrix, rank));

            ranks.push_back(rank);
            p_multi *= primes[s_i++];
        }
        for (uint j = 0; j < ranks.size(); ++j) {
            if (i <= ranks[j]) {
                break;
            } else {
                goto FINAL_STAGE;
            }
        }
    }
    I_ std::cout << "primes used: " << s_i << std::endl;
FINAL_STAGE:;
    // find minimal rank profile
    uint min_prof_idx = 0;

    for (uint i = 1; i < rank_profiles.size(); ++i) {
        bool is_smaller = false;
        if (ranks[min_prof_idx] < ranks[i]) {
            is_smaller = true;
        } else if (ranks[min_prof_idx] == ranks[i])
        {
            std::vector<uint> & a = rank_profiles[min_prof_idx];
            std::vector<uint> & b = rank_profiles[i];
            for (uint j = 0; j < a.size(); ++j) {
                if (a[j] == b[j]) {
                    continue;
                } else if (b[j] < a[j]) {
                    is_smaller = true;
                }
                break;
            }
        }
        if (is_smaller) {
            min_prof_idx = i;
        }
    }
    if (rank_profiles.empty()) {
        return std::vector<uint>();
    } else {
        return rank_profiles[min_prof_idx];
    }
}

std::vector<uint>
SNF::get_rank_profile_mod(NTL::mat_ZZ_p & ntl_matrix, long & rank)
{
    rank = NTL::gauss(ntl_matrix);
    std::vector<uint> rank_profile;

    uint j = 0;
    for (uint i = 0; i < rank; ++i) {
        while (j < ntl_matrix.NumCols()) {
            if (0 != ntl_matrix(i + 1, j + 1)) {
                rank_profile.push_back(j++);
                goto NEXT_ROW;
            } else {
                ++j;
            }
        }
    NEXT_ROW:;
    }

    return rank_profile;
}