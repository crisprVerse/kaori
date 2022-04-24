#ifndef KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP
#define KAORI_DUAL_BARCODES_WITH_DIAGNOSTICS_HPP

#include "DualBarcodes.hpp"
#include "CombinatorialBarcodesPairedEnd.hpp"
#include "../utils.hpp"

/**
 * @file DualBarcodesWithDiagnostics.hpp
 *
 * @brief Process dual barcodes with extra diagnostics.
 */

namespace kaori {

/**
 * @brief Handler for dual barcodes with extra diagnostics.
 *
 * This provides the same information as `DualBarcodes` but also captures the frequency of the invalid combinations.
 * These frequences can be helpful for diagnosing problems with library construction.
 * The handler also counts the number of reads where only one barcode construct matches to a read.
 *
 * @tparam N Size of the bitset to use for each constant template.
 * The maximum size of the template is defined as `N / 4`, see `ConstantTemplate` for details.
 */
template<size_t N>
class DualBarcodesWithDiagnostics { 
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
     *
     * `var2` and `var1` are assumed to have the same number of barcodes.
     * Corresponding values across `var1` and `var2` define a particular combination. 
     */
    DualBarcodesWithDiagnostics(
        const char* con1, size_t n1, bool rev1, const SequenceSet& var1, int mm1, 
        const char* con2, size_t n2, bool rev2, const SequenceSet& var2, int mm2,
        bool random = false
    ) :
        dual_handler(con1, n1, rev1, var1, mm1, con2, n2, rev2, var2, mm2, random),
        combo_handler(con1, n1, rev1, var1, mm1, con2, n2, rev2, var2, mm2, random, true) // we allow duplicates in the trie.
    {}

    /**
     * @param t Whether to search only for the first match across reads (for valid combinations) or in each read (for invalid combinations).
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `DualBarcodesWithDiagnostics` instance.
     */
    DualBarcodesWithDiagnostics& set_first(bool t = true) {
        dual_handler.set_first(t);
        combo_handler.set_first(t);
        return *this;
    }

private:
    DualBarcodes<N> dual_handler;
    CombinatorialBarcodesPairedEnd<N> combo_handler;

public:
    /**
     *@cond
     */
    struct State {
        State() {}
        State(typename DualBarcodes<N>::State ds, typename CombinatorialBarcodesPairedEnd<N>::State cs) : dual_state(std::move(ds)), combo_state(std::move(cs)) {}

        /**
         * @cond
         */
        typename DualBarcodes<N>::State dual_state;
        typename CombinatorialBarcodesPairedEnd<N>::State combo_state;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(dual_handler.initialize(), combo_handler.initialize());
    }

    void reduce(State& s) {
        dual_handler.reduce(s.dual_state);
        combo_handler.reduce(s.combo_state);
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
        // Only searching for combinations if we couldn't find a proper dual barcode match.
        if (!dual_handler.process(state.dual_state, r1, r2)) {
            combo_handler.process(state.combo_state, r1, r2);
        }
    }
    /**
     *@endcond
     */

public:
    /**
     * @return Sort the invalid combinations for easier frequency counting.
     */
    void sort() {
        combo_handler.sort();
    }

    /**
     * @return Vector containing the frequency of each combination.
     * This has length equal to the number of valid combinations (i.e., the length of `var1` and `var2` in the constructor).
     */
    const std::vector<int>& get_counts() const {
        return dual_handler.get_counts();
    }

    /**
     * @return All invalid combinations encountered by the handler.
     */
    const std::vector<std::array<int, 2> >& get_combinations() const {
        return combo_handler.get_combinations();
    }

    /**
     * @return Total number of read pairs processed by the handler.
     */
    int get_total() const {
        return dual_handler.get_total();
    }

    /**
     * @return Number of read pairs with a valid match to the first barcode but no valid match to the second barcode.
     */
    int get_barcode1_only() const {
        return combo_handler.get_barcode1_only();
    }

    /**
     * @return Number of read pairs with a valid match to the second barcode but no valid match to the first barcode.
     */
    int get_barcode2_only() const {
        return combo_handler.get_barcode2_only();
    }
};

}

#endif
