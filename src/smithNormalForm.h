#ifndef SNF_SMITH_NORMAL_FORM_H
#define SNF_SMITH_NORMAL_FORM_H

#include <armadillo>
#include <vector>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

#include "matrix.h"
#include "numeric.h"

class SNF {
public:
    /**
     * Computes Smith Normal Form
     * @param matrix Input matrix whose SNF should be computed
     * @return Returns matrix in Smith Normal Form
     */
    arma::imat calculate(arma::imat matrix);
private:
    void calculate_storjohann(IMat & m,
                                std::vector<uint> & rank_profile,
                                int_t p);
    void calculate_probabilistic(arma::imat & m);
    /**
     * Computes moduli primes whose multiple bounds matrix determinant
     * @param matrix Input matrix for which primes will be computed
     * @return Returns vector containing primes bounding matrix determinant.
     */
    std::vector<NTL::ZZ> get_primes(const arma::imat & matrix);
    /**
     * Computes rank profile of given matrix using primes bounding determinant
     * @param matrix Input matrix whose rank profile should be computed
     * @param primes Vector of primes bounding determinant
     * @return Returns vector containing rank profile.
     */
    std::vector<uint> get_rank_profile(const arma::imat & matrix,
                                        const std::vector<NTL::ZZ> & primes);
    /**
     * Computes rank profile of given matrix
     * @param matrix Input matrix whose rank profile should be computed
     * @param p      Compute the rank profile in Z_p
     * @return Returns vector containing rank profile.
     */
    std::vector<uint> get_rank_profile_mod(NTL::mat_ZZ_p & ntl_matrix,
                                            long & rank);
};

#endif // SNF_SMITH_NORMAL_FORM_H