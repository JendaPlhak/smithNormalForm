#ifndef SNF_STORJOHANN_TRIANGULAR_H
#define SNF_STORJOHANN_TRIANGULAR_H

#define ARMA_64BIT_WORD
#include <armadillo>

class IncorrectForm : public std::exception {
    std::string error_message;
    virtual const char* what() const throw() {
        return error_message.c_str();
    }
public:
    IncorrectForm(const std::string & message) : error_message(message) {}
    std::string str() const { return error_message; }
};

void makeHermiteNormalForm(arma::imat & A);
void hermiteTriangToSNF(arma::subview<arma::sword> A);
void eliminateExtraColumns(arma::imat & T);
void reduceResultingSquareToSNF(arma::imat & T);


#endif // SNF_STORJOHANN_TRIANGULAR_H