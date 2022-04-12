#ifndef KAORI_UTILS_HPP
#define KAORI_UTILS_HPP

#include <bitset>

namespace {

template<size_t N>
constexpr std::bitset<N> A_(0);

template<size_t N>
constexpr std::bitset<N> C_(1);

template<size_t N>
constexpr std::bitset<N> G_(2);

template<size_t N>
constexpr std::bitset<N> T_(3);

inline char reverse_complement(char b) {
    char output;
    switch (b) {
        case 'A': case 'a':
            output = 'T';
            break;
        case 'C': case 'c':
            output = 'G';
            break;
        case 'G': case 'g':
            output = 'C';
            break;
        case 'T': case 't':
            output = 'A';
            break;
        default:
            throw std::runtime_error("unknown base '" + std::string(1, b) + "'");
            break;
    }
    return output;
}

template<size_t N>
void add_base(std::bitset<N>& x, char b) {
    x <<= 2;
    switch (b) {
        case 'A': case 'a':
            x |= A_<N>;
            break;
        case 'C': case 'c':
            x |= C_<N>;
            break;
        case 'G': case 'g':
            x |= G_<N>;
            break;
        case 'T': case 't':
            x |= T_<N>;
            break;
        default:
            throw std::runtime_error("unknown base '" + std::string(1, b) + "'");
            break;
    }
    return;
}

}

#endif
