#pragma once

#include <armadillo>

#include "numeric.h"
#include "storjohannTriangular.h"

class HermiteTriangularError : public IncorrectForm {
public:
    HermiteTriangularError(const std::string & message)
     : IncorrectForm(message) { }
};

class ColumnReductionError : public HermiteTriangularError {
public:
    ColumnReductionError(const std::string & message)
     : HermiteTriangularError(message) { }
};

class ConditioningRoutineError : public HermiteTriangularError {
public:
    ConditioningRoutineError(const std::string & message)
     : HermiteTriangularError(message) { }
};

/*
    Input:  general integer matrix
    Action: Calculates triangular Hermite normal form.
 */
void triangularize(IMat & A, std::vector<uint> rank_profile, int_t p);
