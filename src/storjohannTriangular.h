#ifndef SNF_STORJOHANN_TRIANGULAR_H
#define SNF_STORJOHANN_TRIANGULAR_H

void makeHermiteNormalForm(arma::imat & A);
void hermiteTriangToSNF(arma::imat & A);
void eliminateExtraColumns(arma::imat & T);


#endif // SNF_STORJOHANN_TRIANGULAR_H