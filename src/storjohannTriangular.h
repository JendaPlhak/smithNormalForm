#ifndef SNF_STORJOHANN_TRIANGULAR_H
#define SNF_STORJOHANN_TRIANGULAR_H

#include <armadillo>

void makeHermiteNormalForm(arma::imat & A);
void hermiteTriangToSNF(arma::subview<arma::sword> A);
void eliminateExtraColumns(arma::imat & T);
void reduceResultingSquareToSNF(arma::imat & T);


#endif // SNF_STORJOHANN_TRIANGULAR_H