#include "smithNormalForm.h"
#include <eigen3/Eigen/Dense>
#include <boost/math/common_factor.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace Eigen;

void
SNF::calculate_naive(MatrixXi & m)
{   
    std::cout << "Performing diagonalization...\n";
    diagonalize(m);
    std::cout << "Sorting diagonal...\n";
    qsort_diagonal(m);
    std::cout << "Ensuring divisibility...\n";
    // std::cout << m << std::endl;
    ensure_divisibility(m);
    std::cout << "Complete!\n";

    std::cout << m(m.cols() - 1, m.cols() - 1) << std::endl;
}

void
SNF::diagonalize(MatrixXi & m)
{
    for (unsigned int i = m.rows(); 0 < i; --i) {
        make_gcd(m.bottomRightCorner(i, i));
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
SNF::make_gcd(MatrixXi & m)
{   
    for (unsigned int i = 1; i < m.cols(); ++i) {
        // Check Euclidean algorithm invariant |a| >= |b|
        if ( std::abs(m(0,0)) < std::abs(m(i,0)) ) {
            m.row(0).swap(m.row(i));
        }
        // Apply euclidean algorithm
        while (0 != m(i,0)) {
            m.row(0) -= (m(0,0) / m(i,0)) * m.row(i);
            m.row(0).swap(m.row(i));
        }
    }
    // Since -1 is invertible element in Z we can take absolute value of 
    // diagonal element.
    m(0, 0) = std::abs(m(0, 0));
    for (unsigned int j = 1; j < m.cols(); ++j) {
        m(0, j) = 0;
    }
    // cout << m << endl << endl;
}

/*
    Input:  diagonal matrix m
    Action: Ensures that elements on diagonal are devising other elements in 
            ascending order.
 */
void
SNF::ensure_divisibility(MatrixXi & m) 
{
    auto diag = m.diagonal();
    for (unsigned int i = 0; i < m.rows() - 1; ++i) {
        unsigned int gcd = boost::math::gcd(diag[i], diag[i + 1]);
        diag[i + 1]      = boost::math::lcm(diag[i], diag[i + 1]);
        diag[i]          = gcd;
    }
}

void
SNF::qsort_diagonal(MatrixXi & m) 
{
    #define  MAX_LEVELS  1000
    int piv;
    int beg[MAX_LEVELS];
    int end[MAX_LEVELS];
    int i = 0;
    int L;
    int R;

    auto diag = m.diagonal();

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