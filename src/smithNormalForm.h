#ifndef SNF_SMITH_NORMAL_FORM_H
#define SNF_SMITH_NORMAL_FORM_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

class SNF {
public:
    void calculate_naive(MatrixXi & m);
private:
    void diagonalize(MatrixXi & m);
    void ensure_divisibility(MatrixXi & m);
    void make_gcd(Block<MatrixXi> m);
    void qsort_diagonal(MatrixXi & m);
};

#endif // SNF_SMITH_NORMAL_FORM_H