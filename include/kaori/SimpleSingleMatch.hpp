#ifndef KAORI_SIMPLE_SINGLE_MATCH_HPP
#define KAORI_SIMPLE_SINGLE_MATCH_HPP

#include "ScanTemplate.hpp"
#include "BarcodePool.hpp"
#include "BarcodeSearch.hpp"
#include "utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file SimpleSingleMatch.hpp
 *
 * @brief Defines the `SimpleSingleMatch` class.
 */

namespace kaori {

/**
 * @brief Search for a template with a single variable region.
 *
 * This class implements the most common use case for barcode matching, where the template sequence has a single variable region.
 * It will find a match to any valid barcoding element, i.e., the realization of the template where the variable region is replaced with one sequence from a pool of known barcodes.
 * No restrictions are placed on the distribution of mismatches throughout the barcoding element.
 * 
 * @tparam max_size_ Maximum length of the template sequence.
 */
template<SeqLength max_size_>
class SimpleSingleMatch {
public:
    /**
     * @brief Optional parameters for `SimpleSingleMatch`.
     */
    struct Options {
        /**
         * Maximum number of mismatches for any search performed by `SimpleBarcodeSearch::search`.
         */
        SeqLength max_mismatches = 0;

        /** 
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;

        /** 
         * Strand(s) of the read sequence to search.
         */
        SearchStrand strand = SearchStrand::FORWARD;
    };

public:
    /**
     * @param[in] template_seq Pointer to a character array containing the template sequence, see `ScanTemplate`.
     * @param template_length Length of the array pointed to by `barcode_length`.
     * This should be less than or equal to `max_size`.
     * @param barcode_pool Known sequences for the single variable region in `template_seq`.
     * @param options Optional parameters.
     */
    SimpleSingleMatch(
        const char* template_seq, 
        SeqLength template_length, 
        const BarcodePool& barcode_pool, 
        const Options& options
    ) : 
        my_forward(search_forward(options.strand)), 
        my_reverse(search_reverse(options.strand)),
        my_max_mm(options.max_mismatches),
        my_constant(template_seq, template_length, options.strand)
    {
        // Exact strandedness doesn't matter here, just need the number and length.
        const auto& regions = my_constant.forward_variable_regions();
        if (regions.size() != 1) {
            throw std::runtime_error("expected one variable region in the constant template");
        }

        SeqLength var_length = regions[0].second - regions[0].first;
        if (var_length != barcode_pool.length()) {
            throw std::runtime_error("length of barcode_pool sequences (" + std::to_string(barcode_pool.length()) + 
                ") should be the same as the barcode_pool region (" + std::to_string(var_length) + ")");
        }

        SimpleBarcodeSearch::Options bs_opt;
        bs_opt.duplicates = options.duplicates;
        bs_opt.max_mismatches = my_max_mm;

        if (my_forward) {
            bs_opt.reverse = false;
            my_forward_lib = SimpleBarcodeSearch(barcode_pool, bs_opt);
        }
        if (my_reverse) {
            bs_opt.reverse = true;
            my_reverse_lib = SimpleBarcodeSearch(barcode_pool, bs_opt);
        }
    }

private:
    bool my_forward, my_reverse;
    SeqLength my_max_mm;
    ScanTemplate<max_size_> my_constant;
    SimpleBarcodeSearch my_forward_lib, my_reverse_lib;

public:
    /**
     * @brief State of the search on a read sequence.
     */
    struct State {
        /**
         * Index of the known barcode sequence that matches the variable region in the read sequence, after each call to `search()`.
         * This will be set to `UNMATCHED` if no match was found.
         */
        SeqIndex index = UNMATCHED;

        /**
         * Position of the match to the template after each call to `search()`.
         * This is reported as the position on the read at the start of the template.
         * This should only be used if `index != UNMATCHED`.
         */
        SeqLength position = 0;

        /**
         * Total number of mismatches after each call to `search()`.
         * This include both the constant and variable regions.
         * This should only be used if `index != UNMATCHED`.
         */
        SeqLength mismatches = 0;

        /**
         * Total number of mismatches in the variable region, after each call to `search()`.
         * This should only be used if `index != UNMATCHED`.
         */
        SeqLength variable_mismatches = 0;

        /**
         * Whether the match was found on the reverse strand of the read sequence.
         */
        bool reverse = false;

        /**
         * @cond
         */
        std::string buffer;
        typename SimpleBarcodeSearch::State forward_details, reverse_details;
        /**
         * @endcond
         */
    };

    /**
     * Initialize the search state for thread-safe execution.
     *
     * @return A new `State` object.
     */
    State initialize() const {
        return State();
    }

    /**
     * Incorporate search cache optimizations from `state`, see `SimpleBarcodeSearch::reduce()` for details.
     * This allows regular consolidation and sharing of optimizations across threads.
     * On return, optimizations from `state` are incorporated into this `SimpleSingleMatch` instance.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has been used in `search_first()` or `search_best()` at least once.
     */
    void reduce(State& state) {
        if (my_forward) {
            my_forward_lib.reduce(state.forward_details);
        }
        if (my_reverse) {
            my_reverse_lib.reduce(state.reverse_details);
        }
    }

private:
    bool has_match(SeqIndex obs_index, SeqLength obs_mismatches) const {
        return (obs_index != UNMATCHED && obs_mismatches <= my_max_mm);
    }

    void forward_match(const char* seq, const typename ScanTemplate<max_size_>::State& details, State& state) const {
        auto start = seq + details.position;
        const auto& range = my_constant.forward_variable_regions()[0];
        state.buffer.clear();
        state.buffer.insert(state.buffer.end(), start + range.first, start + range.second);
        my_forward_lib.search(state.buffer, state.forward_details, my_max_mm - details.forward_mismatches);
    }

    void reverse_match(const char* seq, const typename ScanTemplate<max_size_>::State& details, State& state) const {
        auto start = seq + details.position;
        const auto& range = my_constant.reverse_variable_regions()[0];
        state.buffer.clear();
        state.buffer.insert(state.buffer.end(), start + range.first, start + range.second);
        my_reverse_lib.search(state.buffer, state.reverse_details, my_max_mm - details.reverse_mismatches);
    }

public:
    /**
     * Search a read for the first match to a valid barcoding element.
     * A match is only reported if the number of mismatches of the entire barcoding element's sequence to the read is no greater than `max_mismatches` (see the constructor)
     * and there is exactly one barcode sequence with the fewest mismatches to the read sequence at the variable region.
     *
     * @param[in] read_seq Pointer to a character array containing the read sequence.
     * @param read_length Length of the read sequence.
     * @param state State object, used to store the search result.
     *
     * @return Whether an appropriate match was found.
     * If `true`, `state` is filled with the details of the first match.
     */
    bool search_first(const char* read_seq, SeqLength read_length, State& state) const {
        auto deets = my_constant.initialize(read_seq, read_length);
        bool found = false;
        state.index = UNMATCHED;
        state.mismatches = 0;
        state.variable_mismatches = 0;

        auto update = [&](bool rev, SeqLength const_mismatches, const typename SimpleBarcodeSearch::State& x) -> bool {
            if (x.index == UNMATCHED) {
                return false;
            }

            SeqLength total = const_mismatches + x.mismatches;
            if (total > my_max_mm) {
                return false;
            }

            found = true;
            state.position = deets.position;
            state.mismatches = total;
            state.reverse = rev;
            state.index = x.index;
            state.variable_mismatches = x.mismatches;
            return true;
        };

        while (!deets.finished) {
            my_constant.next(deets);

            if (my_forward && has_match(deets.position, deets.forward_mismatches)) {
                forward_match(read_seq, deets, state);
                if (update(false, deets.forward_mismatches, state.forward_details)) {
                    break;
                }
            }

            if (my_reverse && has_match(deets.position, deets.reverse_mismatches)) {
                reverse_match(read_seq, deets, state);
                if (update(true, deets.reverse_mismatches, state.reverse_details)) {
                    break;
                }
            }
        }

        return found;
    }

    /**
     * Search a read for the best match to a valid barcoding element. 
     * This is slower than `search_first()` but will find the matching position with the fewest mismatches.
     * If multiple positions are tied for the fewest mismatches, no match is reported.
     *
     * @param[in] read_seq Pointer to a character array containing the read sequence.
     * @param read_length Length of the read sequence.
     * @param state State object, used to store the search result.
     *
     * @return Whether a match was found.
     * If `true`, `state` is filled with the details of the best match.
     */
    bool search_best(const char* read_seq, SeqLength read_length, State& state) const {
        auto deets = my_constant.initialize(read_seq, read_length);
        state.index = UNMATCHED;
        bool found = false;
        SeqLength best = my_max_mm + 1;

        auto update = [&](bool rev, SeqLength const_mismatches, const typename SimpleBarcodeSearch::State& x) -> void {
            if (x.index == UNMATCHED) {
                return;
            }

            auto total = x.mismatches + const_mismatches;
            if (total == best) { 
                if (state.index != x.index) { // ambiguous, setting back to a mismatch.
                    found = false;
                    state.index = UNMATCHED;
                }
            } else if (total < best) {
                found = true;
                best = total; 
                // A further optimization at this point would be to narrow
                // max_mm to the current 'best'. But this probably
                // isn't worth it.

                state.index = x.index;
                state.mismatches = total;
                state.variable_mismatches = x.mismatches;
                state.position = deets.position;
                state.reverse = rev;
            }
        };

        while (!deets.finished) {
            my_constant.next(deets);

            if (my_forward && has_match(deets.position, deets.forward_mismatches)) {
                forward_match(read_seq, deets, state);
                update(false, deets.forward_mismatches, state.forward_details);
            }

            if (my_reverse && has_match(deets.position, deets.reverse_mismatches)) {
                reverse_match(read_seq, deets, state);
                update(true, deets.reverse_mismatches, state.reverse_details);
            }
        }

        return found;
    }
};

}

#endif
