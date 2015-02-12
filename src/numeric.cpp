#include <armadillo>
#include <cmath>
#include <exception>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <boost/math/common_factor.hpp>

#include "numeric.h"

long
positive_log(const long x)
{
    if (x == 0) {
        return 1;
    } else {
        return 1 + std::floor(std::log2(std::abs(x)));
    }
}

//! For given integers a,b,N, N positive and gcd(a,b,N) = d calculates integer c
//! 0 <= c < N and such that gcd(a + cb, N) = d
int_t
gcdCombination(int_t a, int_t b, int_t N)
{
    using boost::math::gcd;
    // input arguments check
    if (N <= 0) {
        throw std::logic_error("Argument N has to be positive!");
    }
    // if d != 1, we divide a, b, N by d which reduces the problem to Theorem 5
    int_t d = gcd(a, gcd(b, N));
    if (1 != d) {
        a /= d; b /= d; N /= d;
    }

    int_t c = -1;
    if (1 == gcd(a, N)) {
        c = 0;
    } else if (1 == gcd(a + b, N)) {
        c = 1;
    } else {
        int_t g = gcd(a, N);
        int_t x = g;
        int_t y = N / g;      // ensuring that gcd(g, y) = 1

        // ensure that gcd(x, y) = 1
        while (1 != gcd(x, y)) {
            int_t gcd_x_y = gcd(x, y);
            y = y / gcd_x_y;
            x = x * gcd_x_y;
        }
        // check invariant
        if (1 != gcd(x, y)) {
            throw std::logic_error("gcd(x, y) == 1 doesn't hold!");
        } else if (1 != gcd(g, y)) {
            throw std::logic_error("gcd(g, y) == 1 doesn't hold!");
        }
        // printf("Factorization N = %d * %d\n", x, y);
        // apply Chinese remaindering algorithm to obtain c that satisfies
        // (c = 1) mod x, (c = 0) mod y
        c = CRT(x, y);
    }

    // check invariant
    if (1 != gcd(a + c * b, N)) {
        throw std::logic_error("gcd(a + c * b, N) == 1 doesn't hold!");
    }
    return c;
}

//! Solves equations of the format (c = 1) mod x, (c = 0) mod y
int_t
CRT(const int_t x, const int_t y)
{
    std::vector<int_t> r    = {1, 0};
    std::vector<int_t> mods = {x, y};
    int_t M = 1;
    for(uint i = 0; i < mods.size(); i++) {
        M *= mods[i];
    }
    std::vector<int_t> m, s;
    for(uint i = 0; i < mods.size(); i++){
        m.push_back(M / mods[i]);
        int_t temp = m[i] % mods[i];
        int_t k    = 0;
        /* if there is a possibility of k being very big, then prime factorize m[i],
         * find modular inverse of 'temp' of each of the factors
         * 'k' equals to the multiplication ( modular mods[i] ) of modular inverses
         */
        while((k * temp) % mods[i] != 1) {
            k++;
        }
        s.push_back(k);
    }
    int_t c = 0;
    for(uint i = 0; i < s.size(); i++) {
        c += ((m[i]*s[i]) % M * r[i]) % M;
        if(c >= M) {
            c -= M;
        }
    }
    return c;
}

// Given a and b, calculate s, t, and gcd(a, b) such that
// s * a + t * b = gcd(a, b)
// Also, if a = +/- gcd, then choose t = 0.
void
extendedGCD(int_t & s_out, int_t & t_out, int_t & gcd_out,
                const int_t a, const int_t b, bool first_nonzero)
{
    static int_t one ( 1 );
    static int_t zero ( 0 );
    static int_t neg_one ( -1 );
    // Relies upon Euclidean Division being available, and also comparison.
    // The user doesn't need to sort the input.
    bool reversed = false;
    if ( a < b ) reversed = true;

    int_t x, y;
    if ( reversed ) {
        x = a;
        y = b;
    } else {
        x = b;
        y = a;
    }

   // For extended euclidean algorithm
    int_t s0 = one; int_t s1 = zero;
    int_t t0 = zero; int_t t1 = one;
    while ( x != zero ) {
        int_t q = std::div(y, x).quot;
        int_t r = y - x * q;
        int_t s = s0 - q * s1;
        int_t t = t0 - q * t1;

        s0 = s1; s1 = s;
        t0 = t1; t1 = t;
        y = x;
        x = r;
    }
    // Set output
    // The Bezout coefficients s and t are the second to last ones
    // to be calculated (the last ones give 0 = s*a + t*b)
    if ( not reversed ) {
        s_out = s0;
        t_out = t0;
    } else {
        t_out = s0;
        s_out = t0;
    }
    gcd_out = y;
    // TODO: generalize this to all unit multiples, somehow, for general Rs
    // TODO: should this really be here? perhaps the caller should worry about this
    if ( gcd_out == a ) {
        s_out = one;
        t_out = zero;
    }
    if ( gcd_out == -a ) {
        s_out = neg_one;
        t_out = zero;
    }
    if ( first_nonzero ) {
        if (s_out == 0) {
            if (t_out > 0) {
                s_out += b / gcd_out;
                t_out -= a / gcd_out;
            } else {
                s_out -= b / gcd_out;
                t_out += a / gcd_out;
            }
        }
    }
}

//! calculates multiple product of elements on diagonal
int_t
diagonalMultiple(const arma::diagview<int_t> & diag)
{
    int_t det = 1;
    for (uint i = 0; i < diag.n_rows; ++i) {
        det *= diag[i];
    }
    return det;
}


int_t
floored_factor(const int_t x, const int_t y)
{
    auto r = std::div(x, y);
    if (x >= 0 || r.rem == 0) {
        return r.quot;
    } else {
        return r.quot - 1;
    }
}

void
matrix_convert(const NTL::mat_ZZ_p & A_from, arma::imat & A_to)
{
    A_to = arma::imat(A_from.NumRows(), A_from.NumCols());
    for (uint i = 0; i < A_from.NumRows(); ++i) {
        for (uint j = 0; j < A_from.NumCols(); ++j) {
            A_to(i, j) = ZZ_to_int_t(A_from(i + 1, j + 1));
        }
    }
}

void
matrix_convert(const arma::imat & A_from, NTL::mat_ZZ_p & A_to)
{
    A_to.SetDims(A_from.n_rows, A_from.n_cols);

    for (uint i = 0; i < A_from.n_rows; ++i) {
        for (uint j = 0; j < A_from.n_cols; ++j) {
            A_to(i + 1, j + 1) = A_from(i, j);
        }
    }
}
