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
 * In this design, the target sequence is created from a template with a single variable region drawn from a pool of barcode sequences.
 * The construct containing the target sequence is then subjected to paired-end sequencing, where either end could contain the target sequence.
 * This handler will search both reads for the target sequence and count the frequency of each barcode.
 *
 * @tparam max_size Maximum length of the template sequence.
 */
template<size_t N>
class SingleBarcodePairedEnd {
public:
    /**
     * @param[in] template_seq Template sequence for the first barcode.
     * This should contain exactly one variable region.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param reverse Whether to search the reverse strand of each read.
     * @param barcode_pool Known barcode sequences for the variable region.
     * @param max_mismatches Maximum number of mismatches allowed across the target sequence.
     */
    SingleBarcodePairedEnd(const char* template_seq, size_t template_length, bool reverse, const SequenceSet& barcode_pool, int max_mismatches = 0) : 
        matcher(template_seq, template_length, !reverse, reverse, barcode_pool, max_mismatches), counts(barcode_pool.size()) {}
        
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
     * This has length equal to the number of valid barcodes (i.e., the length of `barcode_pool` in the constructor).
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
