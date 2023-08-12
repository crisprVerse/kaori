#ifndef KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP
#define KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include "../utils.hpp"

#include <array>
#include <vector>

/**
 * @file CombinatorialBarcodesPairedEnd.hpp
 *
 * @brief Process paired-end combinatorial barcodes.
 */

namespace kaori {

/**
 * @brief Handler for paired-end combinatorial barcodes.
 *
 * In this design, each read contains a target sequence created from a template with a single variable region.
 * For one read, the barcode is drawn from one pool of options, while the other read contains a barcode from another pool.
 * Combinations are assembled randomly by library construction, where the large number of combinations provide many unique identifiers for cell-tracing applications.
 * This handler will capture the frequencies of each barcode combination. 
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class CombinatorialBarcodesPairedEnd { 
public:
    /**
     * @param[in] template_seq1 Pointer to a character array containing the first template sequence. 
     * This should contain exactly one variable region.
     * @param template_length1 Length of the first template.
     * This should be less than or equal to `max_size`.
     * @param reverse1 Whether to search the reverse strand of the read for the first template.
     * @param barcode_pool1 Pool of known barcode sequences for the variable region in the first template.
     * @param max_mismatches1 Maximum number of mismatches across the target sequence corresponding to the first template.
     * @param[in] template_seq2 Pointer to a character array containing the second template sequence. 
     * This should contain exactly one variable region.
     * @param template_length2 Length of the second template.
     * This should be less than or equal to `max_size`.
     * @param reverse2 Whether to search the reverse strand of the read for the second template.
     * @param barcode_pool2 Pool of known barcode sequences for the variable region in the second template.
     * @param max_mismatches2 Maximum number of mismatches across the target sequence corresponding to the second template.
     * @param random Whether the reads are randomized with respect to the first/second target sequences.
     * If `false`, the first read is searched for the first template only, and the second read is searched for the second template only.
     * If `true`, an additional search will be performed in the opposite orientation.
     * @param duplicates How duplicates in `barcode_pool1` and `barcode_pool2` should be handled.
     */
    CombinatorialBarcodesPairedEnd(
        const char* template_seq1, size_t template_length1, bool reverse1, const BarcodePool& barcode_pool1, int max_mismatches1, 
        const char* template_seq2, size_t template_length2, bool reverse2, const BarcodePool& barcode_pool2, int max_mismatches2,
        bool random = false,
        DuplicateAction duplicates = DuplicateAction::ERROR
    ) :
        matcher1(template_seq1, template_length1, !reverse1, reverse1, barcode_pool1, max_mismatches1, duplicates),
        matcher2(template_seq2, template_length2, !reverse2, reverse2, barcode_pool2, max_mismatches2, duplicates),
        randomized(random)
    {
        num_options[0] = barcode_pool1.size();
        num_options[1] = barcode_pool2.size();
    }

    /**
     * @param t Whether to search only for the first match to the target sequence in each read.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `CombinatorialBarcodesPairedEnd` instance.
     */
    CombinatorialBarcodesPairedEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    /**
     *@cond
     */
    struct State {
        State() {}

        State(typename SimpleSingleMatch<max_size>::State s1, typename SimpleSingleMatch<max_size>::State s2) : search1(std::move(s1)), search2(std::move(s2)) {}

        std::vector<std::array<int, 2> >collected;
        int barcode1_only = 0;
        int barcode2_only = 0;
        int total = 0;

        /**
         * @cond
         */
        typename SimpleSingleMatch<max_size>::State search1;
        typename SimpleSingleMatch<max_size>::State search2;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(matcher1.initialize(), matcher2.initialize());
    }

    void reduce(State& s) {
        matcher1.reduce(s.search1);
        matcher2.reduce(s.search2);
        combinations.insert(combinations.end(), s.collected.begin(), s.collected.end());
        total += s.total;
        barcode1_only += s.barcode1_only;
        barcode2_only += s.barcode2_only;
    }

    constexpr static bool use_names = false;
    /**
     *@endcond
     */

public:
    /**
     *@cond
     */
    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        if (use_first) {
            bool m1 = matcher1.search_first(r1.first, r1.second - r1.first, state.search1);
            bool m2 = matcher2.search_first(r2.first, r2.second - r2.first, state.search2);

            if (m1 && m2) {
                state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
            } else if (randomized) {
                bool n1 = matcher1.search_first(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_first(r1.first, r1.second - r1.first, state.search2);
                if (n1 && n2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else {
                    if (m1 || n1) {
                        ++state.barcode1_only;
                    } else if (m2 || n2) {
                        ++state.barcode2_only;
                    }
                }
            } else {
                if (m1) {
                    ++state.barcode1_only;
                } else if (m2) {
                    ++state.barcode2_only;
                }
            }

        } else {
            bool m1 = matcher1.search_best(r1.first, r1.second - r1.first, state.search1);
            bool m2 = matcher2.search_best(r2.first, r2.second - r2.first, state.search2);

            if (!randomized) {
                if (m1 && m2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else if (m1) {
                    ++state.barcode1_only;
                } else if (m2) {
                    ++state.barcode2_only;
                }
            } else if (m1 && m2) {
                std::array<int, 2> candidate{ state.search1.index, state.search2.index };
                int mismatches = state.search1.mismatches + state.search2.mismatches;

                bool n1 = matcher1.search_best(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_best(r1.first, r1.second - r1.first, state.search2);

                if (n1 && n2) {
                    int rmismatches = state.search1.mismatches + state.search2.mismatches;
                    if (mismatches > rmismatches) {
                        state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                    } else if (mismatches < rmismatches) {
                        state.collected.emplace_back(candidate);
                    } else if (candidate[0] == state.search1.index && candidate[1] == state.search2.index) {
                        // If the mismatches are the same, it may not be ambiguous
                        // if the indices would be the same anyway.
                        state.collected.emplace_back(candidate);
                    }
                } else {
                    state.collected.emplace_back(candidate);
                }
            } else {
                bool n1 = matcher1.search_best(r2.first, r2.second - r2.first, state.search1);
                bool n2 = matcher2.search_best(r1.first, r1.second - r1.first, state.search2);

                if (n1 && n2) {
                    state.collected.emplace_back(std::array<int, 2>{ state.search1.index, state.search2.index });
                } else if (m1 || n1) {
                    ++state.barcode1_only;
                } else if (m2 || n2) {
                    ++state.barcode2_only;
                }
            }
        }

        ++state.total;
    }
    /**
     *@endcond
     */

public:
    /**
     * @return Sort the combinations for easier frequency counting.
     * Combinations are sorted by the first index, and then the second index.
     */
    void sort() {
        sort_combinations(combinations, num_options);
    }

    /**
     * @return All combinations encountered by the handler.
     * In each array, the first and second element contains the indices of known barcodes in the first and second pools, respectively.
     */
    const std::vector<std::array<int, 2> >& get_combinations() const {
        return combinations;
    }

    /**
     * @return Total number of read pairs processed by the handler.
     */
    int get_total() const {
        return total;
    }

    /**
     * @return Number of read pairs with a valid match to the first barcode but no valid match to the second barcode.
     */
    int get_barcode1_only() const {
        return barcode1_only;
    }

    /**
     * @return Number of read pairs with a valid match to the second barcode but no valid match to the first barcode.
     */
    int get_barcode2_only() const {
        return barcode2_only;
    }
private:
    SimpleSingleMatch<max_size> matcher1, matcher2;
    std::array<size_t, 2> num_options;

    bool randomized;
    bool use_first = true;

    std::vector<std::array<int, 2> > combinations;
    int total = 0;
    int barcode1_only = 0;
    int barcode2_only = 0;
};

}

#endif
