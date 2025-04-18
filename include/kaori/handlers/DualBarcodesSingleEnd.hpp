#ifndef KAORI_DUAL_BARCODES_SINGLE_END_HPP
#define KAORI_DUAL_BARCODES_SINGLE_END_HPP

#include "../ScanTemplate.hpp"
#include "../BarcodeSearch.hpp"
#include "../utils.hpp"

#include <array>
#include <vector>

/**
 * @file DualBarcodesSingleEnd.hpp
 *
 * @brief Process single-end dual barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end dual barcodes.
 *
 * In this design, the barcoding element is created from a template with multiple variable regions.
 * Each region contains a barcode from a different pool of options, where the valid combinations of barcodes across variable regions are known beforehand.
 * This differs from `CombinatorialBarcodesSingleEnd` where the combinations are assembled randomly.
 * Despite its name, this handler can actually handle any number (>= 2) of variable regions in the combination.
 * It will count the frequency of each barcode combination, along with the total number of reads. 
 *
 * @tparam max_size_ Maximum length of the template sequence.
 */
template<SeqLength max_size_>
class DualBarcodesSingleEnd {
public:
    /**
     * @brief Optional parameters for `DualBarcodeSingleEnd`.
     */
    struct Options {
        /**
         * Maximum number of mismatches allowed across the barcoding element.
         */
        SeqLength max_mismatches = 0;

        /** @param Whether to search only for the first match.
         * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
         */
        bool use_first = true;

        /**
         * Strand(s) of the read sequence to search for the barcoding element.
         */
        SearchStrand strand = SearchStrand::FORWARD;

        /**
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;
    };

public:
    /**
     * @param[in] template_seq Template sequence containing any number (usually 2 or more) of variable regions.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size_`.
     * @param barcode_pools Array containing the known barcode sequences for each of the variable regions, in the order of their appearance in the template sequence.
     * Each pool should have the same length, and corresponding values across pools define a specific combination of barcodes. 
     * @param options Optional parameters.
     */
    DualBarcodesSingleEnd(const char* template_seq, SeqLength template_length, const std::vector<BarcodePool>& barcode_pools, const Options& options) :
        my_forward(search_forward(options.strand)),
        my_reverse(search_reverse(options.strand)),
        my_max_mm(options.max_mismatches),
        my_use_first(options.use_first),
        my_constant_matcher(template_seq, template_length, options.strand)
    {
        const auto& regions = my_constant_matcher.forward_variable_regions();
        my_num_variable = regions.size();
        if (barcode_pools.size() != my_num_variable) {
            throw std::runtime_error("length of 'barcode_pools' should equal the number of variable regions");
        }

        for (decltype(my_num_variable) i = 0; i < my_num_variable; ++i) {
            SeqLength rlen = regions[i].second - regions[i].first;
            SeqLength vlen = barcode_pools[i].length();
            if (vlen != rlen) {
                throw std::runtime_error("length of variable region " + std::to_string(i + 1) + " (" + std::to_string(rlen) + 
                    ") should be the same as its sequences (" + std::to_string(vlen) + ")");
            }
        }

        BarcodeIndex num_choices = 0;
        if (my_num_variable) {
            num_choices = barcode_pools[0].size();
            for (decltype(my_num_variable) i = 1; i < my_num_variable; ++i) {
                if (num_choices != barcode_pools[i].size()) {
                    throw std::runtime_error("all entries of 'barcode_pools' should have the same length");
                }
            }
        }
        my_counts.resize(num_choices);

        // Constructing the combined varlib.
        std::vector<std::string> combined(num_choices); 
        for (decltype(my_num_variable) v = 0; v < my_num_variable; ++v) {
            const auto& curpool = barcode_pools[v];
            SeqLength n = curpool.length();
            for (BarcodeIndex c = 0; c < num_choices; ++c) {
                auto ptr = curpool[c];
                combined[c].insert(combined[c].end(), ptr, ptr + n);
            }
        }

        SimpleBarcodeSearch::Options bopt;
        bopt.max_mismatches = options.max_mismatches;
        bopt.duplicates = options.duplicates;

        if (my_forward) {
            bopt.reverse = false;
            my_forward_lib = SimpleBarcodeSearch(combined, bopt);
        }

        if (my_reverse) {
            bopt.reverse = true;
            my_reverse_lib = SimpleBarcodeSearch(combined, bopt);
        }
    }

private:
    bool my_forward;
    bool my_reverse;
    SeqLength my_max_mm;
    bool my_use_first;

    ScanTemplate<max_size_> my_constant_matcher;
    int my_num_variable;

    SimpleBarcodeSearch my_forward_lib, my_reverse_lib;
    std::vector<Count > my_counts;
    Count my_total = 0;

public:
    /**
     * @cond
     */
    struct State {
        State(decltype(counts) n) : counts(n) {}

        std::vector<Count> counts;
        Count total = 0;

        std::string buffer;

        // Default constructors should be called in this case, so it should be fine.
        typename SimpleBarcodeSearch::State forward_details, reverse_details;
    };
    /**
     * @endcond
     */

private:
    std::pair<Barcodeindex, SeqLength> find_match(
        bool reverse,
        const char* seq, 
        SeqLength position, 
        SeqLength obs_mismatches, 
        const SimpleBarcodeSearch& lib, 
        typename SimpleBarcodeSearch::State state, 
        std::string& buffer
    ) const {
        const auto& regions = my_constant_matcher.variable_regions(reverse);
        buffer.clear();

        for (decltype(my_num_variable) r = 0; r < my_num_variable; ++r) {
            auto start = seq + position;
            buffer.insert(buffer.end(), start + regions[r].first, start + regions[r].second);
        }

        lib.search(buffer, state, my_max_mm - obs_mismatches);
        return std::make_pair(state.index, obs_mismatches + state.mismatches);
    }

    std::pair<Barcodeindex, SeqLength> forward_match(const char* seq, const typename ScanTemplate<max_size_>::State& deets, State& state) const {
        return find_match(false, seq, deets.position, deets.forward_mismatches, my_forward_lib, state.forward_details, state.buffer);
    }

    std::pair<Barcodeindex, SeqLength> reverse_match(const char* seq, const typename ScanTemplate<max_size_>::State& deets, State& state) const {
        return find_match(true, seq, deets.position, deets.reverse_mismatches, my_reverse_lib, state.reverse_details, state.buffer);
    }

private:
    bool process_first(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = my_constant_matcher.initialize(x.first, x.second - x.first);

        while (!deets.finished) {
            my_constant_matcher.next(deets);

            if (my_forward && deets.forward_mismatches <= my_max_mm) {
                auto id = forward_match(x.first, deets, state).first;
                if (id >= 0) {
                    ++state.counts[id];
                    return true;
                }
            }

            if (my_reverse && deets.reverse_mismatches <= my_max_mm) {
                auto id = reverse_match(x.first, deets, state).first;
                if (id >= 0) {
                    ++state.counts[id];
                    return true;
                }
            }
        }
        return false;
    }

    bool process_best(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = my_constant_matcher.initialize(x.first, x.second - x.first);
        bool found = false;
        SeqLength best_mismatches = my_max_mm + 1;
        BarcodeIndex best_id = UNMATCHED;

        auto update = [&](std::pair<BarcodeIndex, SeqLength> match) -> void {
            if (match.first < 0){ 
                return;
            }
            if (match.second == best_mismatches) {
                if (best_id != match.first) { // ambiguous.
                    found = false;
                }
            } else if (match.second < best_mismatches) { 
                // A further optimization at this point would be to narrow
                // max_mm to the current 'best_mismatches'. But
                // this probably isn't worth it.

                found = true;
                best_mismatches = match.second;
                best_id = match.first;
            }
        };

        while (!deets.finished) {
            my_constant_matcher.next(deets);

            if (my_forward && deets.forward_mismatches <= my_max_mm) {
                update(forward_match(x.first, deets, state));
            }

            if (my_reverse && deets.reverse_mismatches <= my_max_mm) {
                update(reverse_match(x.first, deets, state));
            }
        }

        if (found) {
            ++state.counts[best_id];
        }
        return found;
    }

public:
    /**
     * @cond
     */
    State initialize() const {
        return State(my_counts.size());
    }

    void reduce(State& s) {
        if (my_forward) {
            my_forward_lib.reduce(s.forward_details);
        }
        if (my_reverse) {
            my_reverse_lib.reduce(s.reverse_details);
        }

        for (decltype(my_counts.size()) i = 0, end = my_counts.size(); i < end; ++i) {
            my_counts[i] += s.counts[i];
        }
        my_total += s.total;
        return;
    }

    bool process(State& state, const std::pair<const char*, const char*>& x) const {
        ++state.total;
        if (my_use_first) {
            return process_first(state, x);
        } else {
            return process_best(state, x);
        }
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * @return Counts for each combination.
     */
    const std::vector<Count>& get_counts() const {
        return my_counts;
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    Count get_total() const {
        return my_total;
    }
};

}

#endif
