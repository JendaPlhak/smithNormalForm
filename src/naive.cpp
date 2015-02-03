#include "naive.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <armadillo>
#include <boost/math/common_factor.hpp>

void
naive::rowReduce(arma::imat & m)
{
    for (unsigned int i = 0; i < m.n_rows; ++i) {
        naive::make_gcd(m.submat(i, i, m.n_rows - 1, m.n_cols - 1));
        if (m.diag()[i] < 0) {
            m.row(i) *= (-1);
        }
    }
}

void
naive::diagonalize(arma::imat & m)
{
    for (unsigned int i = 0; i < m.n_rows; ++i) {
        naive::make_gcd(m.submat(i, i, m.n_rows - 1, m.n_cols - 1));
        // Since -1 is invertible element in Z we can take absolute value of
        // diagonal element.
        m(0, 0) = std::abs(m(0, 0));
        for (unsigned int j = 1; j < m.n_cols; ++j) {
            m(0, j) = 0;
        }
    }
}

/*
    Makes GCD at position r, c using row and column-wise
    operations. Expects input in format
     ..c..
    .x0000
    .0y000
    r00z..
    .00...
    .00...
 */
void
naive::make_gcd(arma::subview<arma::sword> m)
{
    for (unsigned int i = 1; i < m.n_rows; ++i) {
        // Check Euclidean algorithm invariant |a| >= |b|
        while (std::abs(m(0, 0)) < std::abs(m(i, 0)) ) {
            m.swap_rows(0, i);
        }
        // Apply euclidean algorithm
        while (0 != m(i,0)) {
            m.row(0) -= (m(0,0) / m(i,0)) * m.row(i);
            m.swap_rows(0, i);
        }
    }
}

/*
    Input:  diagonal matrix m
    Action: Ensures that zero diagonal elements are shifted to very end of
    diagonal.
 */
uint
eliminate_zero_elements(arma::imat & m)
{
    arma::imat diag = m.diag();
    int zero_n = 0;
    for (int entry : diag) {
        if (entry == 0) {
            ++zero_n;
        } else {
            break;
        }
    }

    for (uint i = 0; i < m.n_rows; ++i) {
        if (i < m.n_rows - zero_n) {
            m(i, i) = diag[i + zero_n];
        } else {
            m(i, i) = 0;
        }
    }
    return zero_n;
}

/*
    Input:  diagonal matrix m
    Action: Ensures that elements on diagonal are devising other elements in
            ascending order. Zero entries are placed at very end of diagonal.
 */
void
naive::ensure_divisibility(arma::imat & m)
{
    uint zero_n  = eliminate_zero_elements(m);

    for (bool updated = true; true == updated; ) {
        updated = false;
        for (uint i = 0; i < m.n_rows - 1 - zero_n; ++i) {
            int a = m.diag()[i];
            int b = m.diag()[i + 1];
            m.diag()[i]     = boost::math::gcd(a, b);
            m.diag()[i + 1] = boost::math::lcm(a, b);
            if (not updated && (a != m.diag()[i] || b != m.diag()[i + 1])) {
                updated = true;
            }
        }
    }
}

void
naive::qsort_diagonal(arma::imat & m)
{
    #define  MAX_LEVELS  1000
    int piv;
    int beg[MAX_LEVELS];
    int end[MAX_LEVELS];
    int i = 0;
    int L;
    int R;

    arma::imat diag = m.diag();

    beg[0] = 0;
    end[0] = diag.size();
    while (i >= 0) {
        L = beg[i];
        R = end[i] - 1;
        if (L < R) {
            piv = diag[L];
            if (i == MAX_LEVELS - 1)
                return;
            while (L < R) {
                while (diag[R] >= piv && L < R)
                    R--;
                if (L < R)
                    diag[L++] = diag[R];
                while (diag[L] <= piv && L < R)
                    L++;
                if (L < R)
                    diag[R--] = diag[L];
            }
            diag[L]    = piv;
            beg[i + 1] = L + 1;
            end[i + 1] = end[i];
            end[i++]   = L;
        } else {
            i--;
        }
    }
}
