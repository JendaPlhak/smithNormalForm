
inline void
SubRow::operator +=(const RowVec & rhs) {
    m_data += rhs.m_vec;
}

inline void
SubRow::operator =(const RowVec & rhs) {
    m_data = rhs.m_vec;
}