#include <cmath>
#include <exception>
#include <iostream>
#include <vector>
#include <boost/math/common_factor.hpp>

#include "storjohannNumeric.h"

//! For given integers a,b,N, N positive and gcd(a,b,N) = d calculates integer c
//! 0 <= c < N and such that gcd(a + cb, N) = d
int
gcdCombination(int a, int b, int N)
{
    using boost::math::gcd;
    // input arguments check
    if (N <= 0) {
        throw ("Argument N has to be positive!");
    }
    // if d != 1, we divide a, b, N by d which reduces the problem to Theorem 5
    int d = gcd(a, gcd(b, N));
    if (1 != d) {
        a /= d; b /= d; N /= d;
    }

    int c = -1;
    if (1 == gcd(a, N)) {
        c = 0;
    } else if (1 == gcd(a + b, N)) {
        c = 1;
    } else {
        int g = gcd(a, N);
        int x = g;
        int y = N / g;      // ensuring that gcd(g, y) = 1

        // ensure that gcd(x, y) = 1
        int gcd_x_y = gcd(x, y);
        y = y / gcd_x_y;
        x = x * gcd_x_y;
        // check invariant
        if (1 != gcd(x, y)) {
            throw ("gcd(x, y) == 1 doesn't hold!");
        } else if (1 != gcd(g, y)) {
            throw ("gcd(g, y) == 1 doesn't hold!");
        }
        // printf("Factorization N = %d * %d\n", x, y);
        // apply Chinese remaindering algorithm to obtain c that satisfies
        // (c = 1) mod x, (c = 0) mod y
        c = CRT(x, y);
    }

    // check invariant
    if (1 != gcd(a + c * b, N)) {
        throw "gcd(a + c * b, N) == 1 doesn't hold!";
    }
    return c;
}

//! Solves equations of the format (c = 1) mod x, (c = 0) mod y
int
CRT(const int x, const int y)
{
    std::vector<int> r    = {1, 0};
    std::vector<int> mods = {x, y};
    int M = 1;
    for(uint i = 0; i < mods.size(); i++) {
        M *= mods[i];
    }
    std::vector<int> m, s;
    for(uint i = 0; i < mods.size(); i++){
        m.push_back(M / mods[i]);
        int temp = m[i] % mods[i];
        int k    = 0;
        /* if there is a possibility of k being very big, then prime factorize m[i],
         * find modular inverse of 'temp' of each of the factors
         * 'k' equals to the multiplication ( modular mods[i] ) of modular inverses
         */
        while((k * temp) % mods[i] != 1) {
            k++;
        }
        s.push_back(k);
    }
    int c = 0;
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
void extendedGCD(int & s_out, int & t_out, int & gcd_out,
                    const int a, const int b)
{
    static int one ( 1 );
    static int zero ( 0 );
    static int neg_one ( -1 );
    // Relies upon Euclidean Division being available, and also comparison.
    // The user doesn't need to sort the input.
    bool reversed = false;
    if ( a < b ) reversed = true;

    int x, y;
    if ( reversed ) {
        x = a;
        y = b;
    } else {
        x = b;
        y = a;
    }

   // For extended euclidean algorithm
    int s0 = one; int s1 = zero;
    int t0 = zero; int t1 = one;
    while ( x != zero ) {
        int q = std::div(y, x).quot;
        int r = y - x * q;
        int s = s0 - q * s1;
        int t = t0 - q * t1;

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
}