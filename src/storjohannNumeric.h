#ifndef SNF_STORJOHANN_NUMERIC_H
#define SNF_STORJOHANN_NUMERIC_H

typedef long long int_t;

//! Functor calculating positive modulo for given number e.
class PositiveModulo {
    int_t m_d;
public:
    PositiveModulo(int_t d) : m_d(d) {}
    int_t operator()(const int_t e) { return (e % m_d + m_d) % m_d; }
};

int_t gcdCombination(int_t a, int_t b, int_t N);

void extendedGCD(int_t & s_out, int_t & t_out, int_t & gcd_out,
                    const int_t a, const int_t b, bool first_nonzero = false);

int_t CRT(const int_t x, const int_t y);

int_t diagonalMultiple(const arma::diagview<int_t> & diag);

//! implements valid x / y over Z
int_t floored_factor(const int_t x, const int_t y);

#endif // SNF_STORJOHANN_NUMERIC_H