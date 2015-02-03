#ifndef SNF_STORJOHANN_NUMERIC_H
#define SNF_STORJOHANN_NUMERIC_H

//! Functor calculating positive modulo for given number e.
class PositiveModulo {
    int m_d;
public:
    PositiveModulo(int d) : m_d(d) {}
    int operator()(const int e) { return (e % m_d + m_d) % m_d; }
};

int gcdCombination(int a, int b, int N);
void extendedGCD(int & s_out, int & t_out, int & gcd_out,
                    const int a, const int b, bool first_nonzero = false);
int CRT(const int x, const int y);
int diagonalMultiple(const arma::diagview<int> & diag);

int floored_factor(const int x, const int y);

#endif // SNF_STORJOHANN_NUMERIC_H