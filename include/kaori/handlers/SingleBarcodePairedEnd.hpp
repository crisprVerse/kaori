#ifndef KAORI_SINGLE_BARCODE_PAIRED_END_HPP
#define KAORI_SINGLE_BARCODE_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include <vector>

/**
 * @file SingleBarcodePairedEnd.hpp
 *
 * @brief Process paired-end single barcodes.
 */

namespace kaori {

/**
 * @brief Handler for paired-end single barcodes.
 *
 * A construct contains the barcode sequence on one side and is subjected to paired-end sequencing.
 * However, either read could contain the barcode sequence.
 * This handler will search both reads for the barcode and count the frequency.
 *
 * @tparam N Size of the bitset to use for each constant template.
 * The maximum size of the template is defined as `N / 4`, see `ConstantTemplate` for details.
 */
template<size_t N>
class SingleBarcodePairedEnd {
public:
    /**
     * @param[in] constant Template sequence for the first barcode.
     * This should contain one variable regions.
     * @param size Length of the template.
     * @param reverse Whether to search the reverse strand of the read sequence.
     * @param variable Known sequences for the variable region.
     * @param mismatches Maximum number of mismatches allowed across the target sequence.
     */
    SingleBarcodePairedEnd(const char* constant, size_t size, bool reverse, const SequenceSet& variable, int mismatches = 0) : 
        matcher(constant, size, !reverse, reverse, variable, mismatches), counts(variable.size()) {}
        
    /**
     * @param t Whether to search only for the first match.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `SingleBarcodePairedEnd` instance.
     */
    SingleBarcodePairedEnd& set_first(bool t = true) {
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

    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        if (use_first) {
            if (matcher.search_first(r1.first, r1.second - r1.first, state.search)) {
                ++state.counts[state.search.index];
            } else if (matcher.search_first(r2.first, r2.second - r2.first, state.search)) {
                ++state.counts[state.search.index];
            }
        } else {
            bool found1 = matcher.search_best(r1.first, r1.second - r1.first, state.search);
            auto id1 = state.search.index;
            auto mm1 = state.search.mismatches;

            bool found2 = matcher.search_best(r2.first, r2.second - r2.first, state.search);
            auto id2 = state.search.index;
            auto mm2 = state.search.mismatches;

            if (found1 && !found2) {
                ++state.counts[id1];
            } else if (!found1 && found2) {
                ++state.counts[id2];
            } else if (found1 && found2) {
                if (mm1 < mm2) {
                    ++state.counts[id1];
                } else if (mm1 > mm2) {
                    ++state.counts[id2];
                } else if (id1 == id2) {
                    ++state.counts[id1];
                }
            }
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
