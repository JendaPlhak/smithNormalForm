#ifndef SNF_STORJOHANN_NUMERIC_H
#define SNF_STORJOHANN_NUMERIC_H

typedef long long int_t;

#include <NTL/mat_ZZ_p.h>

//! Functor calculating positive modulo for given number e.
class PositiveModulo {
    int_t m_d;
public:
    PositiveModulo(int_t d) : m_d(d) {}
    static int_t mod(const int_t e, const int_t d)
    {
        return (e % d + d) % d;
    }
    int_t operator()(const int_t e) { return mod(e, m_d); }
};

//! Calculate log in correspondence with chapter 4.4 in Storjohann
long positive_log(const long number);

int_t gcdCombination(int_t a, int_t b, int_t N);

void extendedGCD(int_t & s_out, int_t & t_out, int_t & gcd_out,
                    const int_t a, const int_t b, bool first_nonzero = false);

int_t CRT(const int_t x, const int_t y);

int_t diagonalMultiple(const arma::diagview<int_t> & diag);

//! implements valid x / y over Z
int_t floored_factor(const int_t x, const int_t y);

template <typename T>
int_t ZZ_to_int_t(const T & big_num);

void matrix_convert(const NTL::mat_ZZ_p & A_from, arma::imat & A_to);
void matrix_convert(const arma::imat & A_from, NTL::mat_ZZ_p & A_to);

#include "numeric.tpp"

#endif // SNF_STORJOHANN_NUMERIC_H