#ifndef SNF_SMITH_NORMAL_FORM_H
#define SNF_SMITH_NORMAL_FORM_H

#include <armadillo>
#include <mutex>
#include <thread>
#include <vector>
#include <queue>
#include <memory>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ.h>
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
    NTL::mat_ZZ calculate(arma::imat matrix);
private:
    static void calculate_storjohann(IMat & m,
                                const std::vector<uint> & rank_profile,
                                int_t p);
    void calculate_probabilistic(arma::imat & m);
    static void
    worker(const std::shared_ptr<std::mutex> q_mut,
            const std::shared_ptr<std::mutex> v_mut,
            std::queue<NTL::ZZ> & primes,
            arma::imat matrix,
            std::vector<std::pair<NTL::ZZ, arma::imat>> & moduli_results,
            std::vector<uint> rank_profile);

    /**
     * Computes moduli primes whose multiple bounds matrix determinant
     * @param matrix Input matrix for which primes will be computed
     * @return Returns vector containing primes bounding matrix determinant.
     */
    static std::vector<NTL::ZZ> get_primes(const arma::imat & matrix);
    /**
     * Computes rank profile of given matrix using primes bounding determinant
     * @param matrix Input matrix whose rank profile should be computed
     * @param primes Vector of primes bounding determinant
     * @return Returns vector containing rank profile.
     */
    static std::vector<uint>
    get_rank_profile(const arma::imat & matrix,
                        const std::vector<NTL::ZZ> & primes);
    /**
     * Computes rank profile of given matrix
     * @param matrix Input matrix whose rank profile should be computed
     * @param p      Compute the rank profile in Z_p
     * @return Returns vector containing rank profile.
     */
    static std::vector<uint>
    get_rank_profile_mod(NTL::mat_ZZ_p & ntl_matrix,
                            long & rank);
    /**
     * Computes resulting matrix from moduli results using Chinese reminder
     * theorem.
     * @param moduli_results Partial SNF moduli results
     * @param primes Primes modulo which the results where obtained
     * @return Returns resulting matrix.
     */
    static NTL::mat_ZZ
    apply_iso(const std::vector<std::pair<NTL::ZZ, arma::imat>> & moduli_results,
                const std::vector<NTL::ZZ> & primes);
};

#endif // SNF_SMITH_NORMAL_FORM_H