#ifndef SNF_STORJOHANN_NUMERIC_H
#define SNF_STORJOHANN_NUMERIC_H

int gcdCombination(int a, int b, int N);
void extendedGCD(int & s_out, int & t_out, int & gcd_out,
                    const int a, const int b);
int CRT(const int x, const int y);
int diagonalMultiple(const arma::diagview<int> & diag);

#endif // SNF_STORJOHANN_NUMERIC_H