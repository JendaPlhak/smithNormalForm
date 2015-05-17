#ifndef SNF_STORJOHANN_TRIANGULAR_H
#define SNF_STORJOHANN_TRIANGULAR_H

#include <armadillo>
#include "matrix.h"

class IncorrectForm : public std::logic_error {
public:
    IncorrectForm(const std::string & message)
     : std::logic_error(message) { }

    virtual ~IncorrectForm() { }
    std::string str() const { return what(); }
};

void makeHermiteNormalForm(IMat& A);
void hermiteTriangToSNF(IMat & A, int_t p);
void eliminateExtraColumns(IMat & T);
void reduceResultingSquareToSNF(IMat & T, int_t p);

/**
 * Strips all trailing zero columns off the matrix T
 * @param T input matrix
 * @return Number of stripped rows
 */
uint stripZeroRows(IMat & T);


#endif // SNF_STORJOHANN_TRIANGULAR_H