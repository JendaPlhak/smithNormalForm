template <typename T>
inline int_t
ZZ_to_int_t(const T & big_num) {
    int_t output;

    std::stringstream big_num_str;
    big_num_str << big_num;
    big_num_str >> output;

    return output;
}