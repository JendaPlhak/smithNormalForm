#ifndef SNF_STORJOHANN_TRIANGULAR_H
#define SNF_STORJOHANN_TRIANGULAR_H

#include <armadillo>

class IncorrectForm : public std::logic_error {
public:
    IncorrectForm(const std::string & message)
     : std::logic_error(message) { }

    virtual ~IncorrectForm() { }
    std::string str() const { return what(); }
};

void makeHermiteNormalForm(arma::imat & A);
void hermiteTriangToSNF(arma::imat & A);
void eliminateExtraColumns(arma::imat & T);
void reduceResultingSquareToSNF(arma::imat & T);

/**
 * Strips all trailing zero columns off the matrix T
 * @param T input matrix
 * @return Number of stripped rows
 */
uint stripZeroRows(arma::imat & T);


#endif // SNF_STORJOHANN_TRIANGULAR_H