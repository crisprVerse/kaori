#ifndef KAORI_EXACT_LIBRARY_HPP
#define KAORI_EXACT_LIBRARY_HPP

#include <unordered_map>
#include <bitset>

namespace kaori {

template<size_t N>
struct ExactLibrary {
public:
    ExactLibrary(int seq_len, size_t expected_num_seq, bool forward, bool reverse) : length(seq_len) {
        if (seq_len * 2 > N) {
            throw std::runtime_error("sequence length must be no greater than the hard limit of " + std::to_string(N/2));
        }
        if (forward || reverse) {
            if (forward && reverse) {
                expected_num_seq *= 2;
            }
            mappings.reserve(expected_num_seq);
        }
    }

    void add(const char* seq) {
        if (forward) {
            std::bitset<N> current;
            for (int i = 0; i < length; ++i) {
                add_base<N>(current, seq[i]);
            }
            add(current, false);
        }

        if (reverse) {
            std::bitset<N> current;
            for (int i = 0; i < length; ++i) {
                add_base<N>(current, reverse_complement(seq[length - i - 1]));
            }
            add(current, true);
        }

        ++counter;
        return;
    }

    void finish() {
        return;
    }

    std::pair<int, bool> query(const std::bitset<N>& hash) const {
        auto it = mappings.find(hash);
        if (it == mappings.end()) {
            return std::make_pair(-1, false);
        } else {
            return it->second;
        }
    }

    int hash_length() const {
        return length;
    }

private:
    std::unordered_map<std::bitset<N>, std::pair<int, bool> > mappings;
    int counter = 0;
    int length;

    void add(const std::bitset<N>& x, bool reverse) {
        auto it = mappings.find(x);
        if (it != mappings.end()) {
            (it->second).first = -1;
        } else {
            mappings[x] = std::make_pair(counter, reverse);
        }
    }
};

}

#endif
