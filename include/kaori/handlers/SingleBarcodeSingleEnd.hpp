#ifndef KAORI_SINGLE_BARCODE_SINGLE_END_HPP
#define KAORI_SINGLE_BARCODE_SINGLE_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

/**
 * @file SingleBarcodeSingleEnd.hpp
 *
 * @brief Process single-end single barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end single barcodes.
 *
 * A construct contains the barcode sequence and is subjected to single-end sequencing.
 * This handler will search the read for the barcode - possibly on both strands - and count the frequency.
 *
 * @tparam N Size of the bitset to use for each constant template.
 * The maximum size of the template is defined as `N / 4`, see `ConstantTemplate` for details.
 */
template<size_t N>
class SingleBarcodeSingleEnd {
public:
    /**
     * @param[in] constant Template sequence for the first barcode.
     * This should contain one variable regions.
     * @param size Length of the template.
     * @param strand Strand to use for searching the read sequence - forward (0), reverse (1) or both (2).
     * @param variable Known sequences for the variable region.
     * @param mismatches Maximum number of mismatches allowed across the target sequence.
     */
    SingleBarcodeSingleEnd(const char* constant, size_t size, int strand, const SequenceSet& variable, int mismatches = 0) : 
        matcher(constant, size, strand != 1, strand != 0, variable, mismatches), counts(variable.size()) {}
        
    /**
     * @param t Whether to search only for the first match.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `SingleBarcodeSingleEnd` instance.
     */
    SingleBarcodeSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    /**
     * @cond
     */
    struct State {
        State() {}

        State(typename SimpleSingleMatch<N>::SearchState s, size_t nvar) : search(std::move(s)), counts(nvar) {}

        typename SimpleSingleMatch<N>::SearchState search;
        std::vector<int> counts;
        int total = 0;
    };

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        bool found = false;
        if (use_first) {
            found = matcher.search_first(x.first, x.second - x.first, state.search);
        } else {
            found = matcher.search_best(x.first, x.second - x.first, state.search);
        }
        if (found) {
            ++(state.counts[state.search.index]);
        }
        ++state.total;
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    State initialize() const {
        return State(matcher.initialize(), counts.size());
    }

    void reduce(State& s) {
        matcher.reduce(s.search);
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] += s.counts[i];
        }
        total += s.total;
    }
    /**
     * @endcond
     */

private:
    SimpleSingleMatch<N> matcher;
    std::vector<int> counts;
    int total = 0;
    bool use_first = true;

public:
    /**
     * @return Vector containing the frequency of each barcode.
     * This has length equal to the number of valid barcodes (i.e., the length of `variable` in the constructor).
     */
    const std::vector<int>& get_counts() const {
        return counts;        
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    int get_total() const {
        return total;
    }
};

}

#endif
