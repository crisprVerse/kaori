#ifndef KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP
#define KAORI_COMBINATORIAL_BARCODES_PAIRED_END_HPP

#include "../SimpleSingleMatch.hpp"
#include "../utils.hpp"

/**
 * @file CombinatorialBarcodesPairedEnd.hpp
 *
 * @brief Process paired-end combinatorial barcodes.
 */

namespace kaori {

/**
 * @brief Handler for paired-end combinatorial barcodes.
 *
 * One of the paired reads contains a barcode from one pool of options, while the other read contains a barcode from another pool.
 * Typically, these combinations are assembled randomly by library construction.
 * This handler will capture the frequencies of each barcode combination. 
 *
 * @tparam N Size of the bitset to use for each constant template.
 * The maximum size of the template is defined as `N / 4`, see `ConstantTemplate` for details.
 */
template<size_t N>
class CombinatorialBarcodesPairedEnd { 
public:
    /**
     * @param[in] con1 Template sequence for the first barcode.
     * This should contain one variable region.
     * @param n1 Length of the first barcode template.
     * @param rev1 Whether to search the reverse strand for the first barcode template.
     * @param var1 Set of known sequences for the variable region in the first barcode.
     * @param mm1 Maximum number of mismatches for the first barcode.
     * @param[in] con2 Template sequence for the second barcode.
     * This should contain one variable region.
     * @param n2 Length of the second barcode template.
     * @param rev2 Whether to search the reverse strand for the second barcode template.
     * @param var2 Set of known sequences for the variable region in the second barcode.
     * @param mm2 Maximum number of mismatches for the second barcode.
     * @param random Whether the reads are randomized with respect to the first/second barcode.
     * If `false`, the first read is searched for the first barcode only, and the second read is searched for the second barcode only.
     * If `true`, an additional search will be performed in the opposite orientation.
     * @param duplicates Whether to allow duplicates in `var1` and `var2`, see `MismatchTrie` for details.
     */
    CombinatorialBarcodesPairedEnd(
        const char* con1, size_t n1, bool rev1, const SequenceSet& var1, int mm1, 
        const char* con2, size_t n2, bool rev2, const SequenceSet& var2, int mm2,
        bool random = false,
        bool duplicates = false
    ) :
        matcher1(con1, n1, !rev1, rev1, var1, mm1, duplicates),
        matcher2(con2, n2, !rev2, rev2, var2, mm2, duplicates),
        randomized(random)
    {
        num_options[0] = var1.size();
        num_options[1] = var2.size();
    }

    /**
     * @param t Whether to search only for the first match in each read.
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

        State(typename SimpleSingleMatch<N>::SearchState s1, typename SimpleSingleMatch<N>::SearchState s2) : search1(std::move(s1)), search2(std::move(s2)) {}

        std::vector<std::array<int, 2> >collected;
        int barcode1_only = 0;
        int barcode2_only = 0;
        int total = 0;

        /**
         * @cond
         */
        typename SimpleSingleMatch<N>::SearchState search1;
        typename SimpleSingleMatch<N>::SearchState search2;
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
     */
    void sort() {
        sort_combinations(combinations, num_options);
    }

    /**
     * @return All combinations encountered by the handler.
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
    SimpleSingleMatch<N> matcher1, matcher2;
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
