#pragma once

#include <armadillo>
#include "numeric.h"


/**
 * General matrix interface
 */
class MatrixInterface {

};

class SubRow;
class RowVec;


class SubCol {
public:
    SubCol(arma::subview_row<int_t> data, int_t p_)
     : m_data(data), p(p_)
    {}

    inline int_t
    get_mod() const {
        return p;
    }

    inline const arma::subview_row<int_t> &
    get_base() const {
        return m_data;
    }

    inline bool
    any() const {
        return arma::any(m_data);
    }

    inline SubCol
    subvec(const uint from, const uint to) {
        return SubCol(m_data.subvec(from, to), p);
    }

    inline arma::eOp<arma::subview_row<long long>, arma::eop_scalar_times>
    operator *(const int_t n) const {
        return m_data * n;
    }

    inline void
    operator +=(const SubCol & rhs) {
        m_data += rhs.m_data;
    }

    template <typename T>
    inline void
    operator +=(const T & rhs) {
        m_data += rhs;
    }

    inline void
    operator -=(const SubCol & rhs) {
        m_data -= rhs.m_data;
    }

    template <typename T>
    inline void
    operator -=(const T & rhs) {
        m_data -= rhs;
    }

    inline int_t &
    operator [](const uint i) {
        return m_data[i];
    }

    inline int_t
    operator [](const uint i) const {
        return m_data[i];
    }

private:
    arma::subview_row<int_t> m_data;
    int_t p;
};

inline arma::eOp<arma::subview_row<long long>, arma::eop_scalar_times>
operator *(const int_t n, const SubCol & rhs) {
    return rhs * n;
}

class SubRow {
public:
    SubRow(arma::subview_col<int_t> data, int_t p_)
     : m_data(data), p(p_)
    {}

    inline int_t
    get_mod() const {
        return p;
    }

    inline const arma::subview_col<int_t> &
    get_base() const {
        return m_data;
    }

    inline bool
    any() const {
        return arma::any(m_data);
    }

    inline SubRow
    subvec(const uint begin, const uint end) const {
        return SubRow(m_data.subvec(begin, end), p);
    }

    template <typename OperatorT>
    inline void
    transform(const OperatorT & o) {
        m_data.transform(o);
    }

    inline void
    mulAdd(const SubRow & rhs, const int_t mul) {
        m_data += mul * rhs.m_data;
    }

    template <typename T>
    inline void
    operator =(const T & rhs) {
        m_data = rhs;
    }

    inline void
    operator =(const RowVec & rhs);

    inline void
    operator *=(const int_t n) {
        m_data *= n;
    }

    inline arma::eOp<arma::subview_col<long long>, arma::eop_scalar_times>
    operator *(const int_t n) const {
        return n * m_data;
    }

    inline void
    operator +=(const SubRow & rhs) {
        m_data += rhs.m_data;
    }

    void
    operator +=(const RowVec & rhs);

    template <typename T>
    inline void
    operator +=(const T & rhs) {
        m_data += rhs;
    }

    inline void
    operator -=(const SubRow & rhs) {
        m_data -= rhs.m_data;
    }

    template <typename T>
    inline void
    operator -=(const T & rhs) {
        m_data -= rhs;
    }

private:
    arma::subview_col<int_t> m_data;
    int_t p;

    friend class RowVec;
};

inline arma::eOp<arma::subview_col<long long>, arma::eop_scalar_times>
operator *(const int_t n, const SubRow & rhs) {
    return rhs * n;
}

class RowVec {
public:
    RowVec(arma::icolvec data, int_t p_)
     : m_vec(data), p(p_)
    {}

    RowVec(const SubRow & r)
     : m_vec(r.m_data), p(r.get_mod())
    {}

    inline int_t
    get_mod() const {
        return p;
    }

    inline const arma::icolvec &
    get_base() const {
        return m_vec;
    }

    inline bool
    any() const {
        return arma::any(m_vec);
    }

    inline void
    operator =(const SubRow & rhs) {
        *this = RowVec(rhs);
    }

    inline void
    operator *=(const int_t n) {
        m_vec *= n;
    }

    inline RowVec
    operator *(const int_t n) const {
        return RowVec(n * m_vec, p);
    }

    template <typename T>
    inline T
    operator +(const T & rhs) const {
        return m_vec + rhs;
    }

    inline void
    operator +=(const SubRow & rhs) {
        m_vec += rhs.m_data;
    }

    template <typename T>
    inline void
    operator +=(const T & rhs) {
        m_vec += rhs;
    }

    inline void
    operator -=(const SubRow & rhs) {
        m_vec -= rhs.m_data;
    }

    template <typename T>
    inline void
    operator -=(const T & rhs) {
        m_vec -= rhs;
    }

private:
    arma::icolvec m_vec;
    int_t p;

    friend class SubRow;
};

inline RowVec
operator *(const int_t n, const RowVec & rhs) {
    return rhs * n;
}


class SubMat {
public:
    uint n_rows;
    uint n_cols;

    inline const arma::subview<int_t> &
    get_base() const {
        return m_submat;
    }

    inline float
    det() const {
        return arma::det(arma::conv_to<arma::mat>::from(m_submat));
    }

    inline int_t
    rank() const {
        return arma::rank(arma::conv_to<arma::mat>::from(m_submat));
    }

    inline int_t
    diagMultiple() const {
        return diagonalMultiple(m_submat.diag());
    }

    inline SubRow
    row(uint i) {
        return SubRow(m_submat.col(i), p);
    }
    inline SubCol
    col(uint i) {
        return SubCol(m_submat.row(i), p);
    }

    inline arma::diagview<int_t>
    diag() {
        return m_submat.diag();
    }

    inline void
    transpose() {
        m_submat = m_submat.t();
    }

    template <typename OperatorT>
    inline void
    transform(const OperatorT & o) {
        m_submat.transform(o);
    }

    inline SubMat
    submat(uint row_begin, uint col_begin, uint row_end, uint col_end) {
        return SubMat(m_submat.submat(col_begin, row_begin, col_end, row_end), p);
    }

    inline SubMat
    submat(uint row_begin, uint col_begin, uint row_end, uint col_end) const {
        return SubMat(m_submat.submat(col_begin, row_begin, col_end, row_end), p);
    }

    inline int_t &
    operator ()(const uint row, const uint col) {
        return m_submat(col, row);
    }

    inline int_t
    operator ()(const uint row, const uint col) const {
        return m_submat(col, row);
    }

    template <typename T>
    inline void
    operator =(const T & rhs) {
        m_submat = rhs.get_base();
    }

private:
    SubMat(arma::subview<int_t> submat, const int_t p_)
     : n_rows(submat.n_cols), n_cols(submat.n_rows), m_submat(submat), p(p_)
    {}

    arma::subview<int_t> m_submat;
    int_t p;

    friend class IMat;
    friend std::ostream& operator<<(std::ostream& os, const SubMat& mat);
};

/**
 * Matrix from arma::imat with columns and rows reverted.
 */
class IMat : public MatrixInterface {
public:
    uint n_rows;
    uint n_cols;

    IMat(const arma::imat & matrix, int_t p_)
     : n_rows(matrix.n_rows), n_cols(matrix.n_cols), M(matrix.t()), p(p_)
    {
        // M.transform(Modulo(p));
    }
    IMat(uint rows, uint cols, int_t p_)
     : n_rows(rows), n_cols(cols), M(n_cols, n_rows), p(p_)
    {
        // M.transform(Modulo(p));
    }

    inline int_t
    get_mod() const {
        return p;
    }

    inline const arma::imat &
    get_base() const {
        return M;
    }

    inline float
    is_square() const {
        return M.is_square();
    }

    inline arma::diagview<int_t>
    diag() {
        return M.diag();
    }

    inline int_t
    diagMultiple() const {
        return diagonalMultiple(M.diag());
    }

    template <typename OperatorT>
    inline void
    transform(const OperatorT & o) {
        M.transform(o);
    }

    inline void
    resize(uint rows_, uint cols_) {
        M.resize(cols_, rows_);
        n_rows = rows_;
        n_cols = cols_;
    }

    inline void
    zeros() {
        M.zeros();
    }

    inline SubMat
    submat(uint row_begin, uint col_begin, uint row_end, uint col_end) {
        return SubMat(M.submat(col_begin, row_begin, col_end, row_end), p);
    }

    inline SubRow
    row(uint i) {
        return SubRow(M.col(i), p);
    }
    inline SubCol
    col(uint i) {
        return SubCol(M.row(i), p);
    }

    inline void
    swap_cols(uint i, uint j) {
        M.swap_rows(i, j);
    }

    inline int
    rank() const {
        return arma::rank(arma::conv_to<arma::mat>::from(M));
    }

    inline float
    det() const {
        return arma::det(arma::conv_to<arma::mat>::from(M));
    }

    inline int_t &
    operator ()(uint row, uint col) {
        return M(col, row);
    }

    inline const int_t &
    operator ()(uint row, uint col) const {
        return M(col, row);
    }

    template <typename T>
    inline void
    operator =(const T & rhs) {
        M = rhs.get_base();
    }

private:
    arma::imat M;
    const int_t p;
    friend std::ostream& operator<<(std::ostream& os, const IMat& mat);
};



inline std::ostream& operator<<(std::ostream& os, const IMat& mat)
{
  os << mat.M.t();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const SubMat& mat)
{
  os << mat.m_submat.t();
  return os;
}

namespace mat {
    // TODO probably will create wrong object
    inline IMat join_rows(const SubCol & r1, const SubCol & r2) {
        return IMat(arma::imat(arma::join_cols(r1.get_base(), r2.get_base())).t(),
                    r1.get_mod());
    }
}

#include "matrix.tpp"