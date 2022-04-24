#ifndef KAORI_SEQUENCE_SET_HPP
#define KAORI_SEQUENCE_SET_HPP

#include <vector>
#include <string>

/**
 * @file SequenceSet.hpp
 *
 * @brief Represent sequences for a variable region.
 */

namespace kaori {

/**
 * @brief Choice of sequences for a variable region.
 *
 * The construct template can be realized into a construct sequence by replacing each variable region with a known sequence.
 * The `SequenceSet` class defines the set of possible sequences that may be used for a given variable region.
 * All sequences in this set are assumed to have the same length.
 */
struct SequenceSet {
    /**
     * Default constructor.
     */
    SequenceSet() {}
    
    /**
     * @param c Vector of pointers to sequences of length `l`.
     * @param l Length of each sequence.
     */
    SequenceSet(std::vector<const char*> c, size_t l) : choices(std::move(c)), length(l) {}

    /**
     * @param s Vector of sequences of the same length.
     *
     * It is assumed that the lifetime of `s` (and its contents) exceeds that of the constructed `SequenceSet`. 
     */
    SequenceSet(const std::vector<std::string>& s) {
        if (s.size()) {
            length = s.front().size();
            choices.reserve(s.size());
            for (const auto& x : s) {
                if (x.size() != length) {
                    throw std::runtime_error("sequences for a given variable region should be of a constant length");
                }
                choices.push_back(x.c_str());
            }
        }
    }

    /**
     * Vector containing pointers to sequences of length `len`.
     */
    std::vector<const char*> choices;

    /**
     * Length of each sequence in `choices`.
     */
    size_t length = 0;

    /**
     * @return Number of sequences in the set.
     */
    size_t size() const {
        return choices.size();
    }

    /**
     * @param i Index of the sequence of interest.
     * @return Pointer to the `i`-th sequence in the set.
     */
    const char* operator[](size_t i) const {
        return choices[i];
    }
};

}

#endif
