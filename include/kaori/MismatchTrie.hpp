#ifndef KAORI_MISMATCH_TRIE_HPP
#define KAORI_MISMATCH_TRIE_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <cstddef>

#include "utils.hpp"

/**
 * @file MismatchTrie.hpp
 *
 * @brief Defines trie-based classes for mismatch-tolerant sequence matching.
 */

namespace kaori {

/**
 * Integer type for the trie indices.
 * This is semantically different from `BarcodeIndex` as it can have the special values `TRIE_STATUS_MISSING` and `TRIE_STATUS_AMBIGUOUS`.
 */
typedef std::size_t TrieIndex;

// We would prefer to use size_type but then the template would be recursive,
// i.e., std::vector<typename std::vector<[??recursion here??]>::size_type>.

/**
 * No match to a known barcode in the trie.
 */
inline constexpr TrieIndex TRIE_STATUS_MISSING = static_cast<TrieIndex>(-1);

/**
 * Ambiguous match to two or more known barcodes in the trie.
 */
inline constexpr TrieIndex TRIE_STATUS_AMBIGUOUS = static_cast<TrieIndex>(-2);

/**
 * @param index An index in the trie as returned by, e.g., `AnyMismatches::search()`. 
 * @return Whether `index` corresponds to a barcode.
 * If false, `index` represents one of the error codes, i.e., `TRIE_STATUS_MISSING` or `TRIE_STATUS_AMBIGUOUS`.
 */
inline bool is_trie_index_ok(TrieIndex index) {
    return index < TRIE_STATUS_AMBIGUOUS;
}

/** 
 * @brief Status of barcode sequence addition to the trie.
 */
struct TrieAddStatus {
    /**
     * Whether the newly added sequence contains ambiguous IUPAC codes.
     */
    bool has_ambiguous = false;

    /**
     * Whether the newly added sequence is a duplicate of an existing sequence in the trie.
     */
    bool is_duplicate = false;

    /**
     * Whether the newly added sequence replaced a duplicate in the trie.
     * Only set when `is_duplicate = true` and `duplicates` is set to `DuplicateAction::LAST` in `add()`.
     */
    bool duplicate_replaced = false;

    /**
     * Whether the newly added sequence caused an existing duplicate to be cleared from the trie.
     * Only set when `is_duplicate = true` and `duplicates` is set to `DuplicateAction::NONE` in `add()`.
     */
    bool duplicate_cleared = false;
};

/**
 * @cond
 */
template<char base_>
int trie_base_shift() {
    if constexpr(base_ == 'A') {
        return 0;
    } else if constexpr(base_ == 'C') {
        return 1;
    } else if constexpr(base_ == 'G') {
        return 2;
    } else { // i.e., base_ == 'T'
        return 3;
    }
}

class MismatchTrie {
public:
    MismatchTrie() = default;

    MismatchTrie(SeqLength barcode_length, DuplicateAction duplicates) : 
        my_length(barcode_length), 
        my_duplicates(duplicates),
        my_pointers(NUM_BASES, TRIE_STATUS_MISSING)
    {}

private:
    SeqLength my_length;
    DuplicateAction my_duplicates;
    std::vector<TrieIndex> my_pointers;
    TrieIndex my_counter = 0;

    TrieIndex next(TrieIndex node) {
        auto current = my_pointers[node]; // don't make this a reference as it gets invalidated by the resize.
        if (current == TRIE_STATUS_MISSING) {
            current = my_pointers.size();
            my_pointers[node] = current;
            my_pointers.resize(current + NUM_BASES, TRIE_STATUS_MISSING);
        }
        return current;
    }

    void end(TrieIndex node, TrieAddStatus& status) {
        auto& current = my_pointers[node];

        if (current == TRIE_STATUS_MISSING) {
            current = my_counter;
        } else if (current == TRIE_STATUS_AMBIGUOUS) {
            status.is_duplicate = true; 
        } else {
            status.is_duplicate = true;
            switch(my_duplicates) {
                case DuplicateAction::FIRST:
                    break;
                case DuplicateAction::LAST:
                    status.duplicate_replaced = true;
                    current = my_counter;
                    break;
                case DuplicateAction::NONE:
                    status.duplicate_cleared = true;
                    current = TRIE_STATUS_AMBIGUOUS;
                    break;
                case DuplicateAction::ERROR:
                    throw std::runtime_error("duplicate sequences detected (" + 
                        std::to_string(current + 1) + ", " + 
                        std::to_string(my_counter + 1) + ") when constructing the trie");
            }
        }
    }

    template<char base_>
    void process_ambiguous(SeqLength i, TrieIndex node, const char* barcode_seq, TrieAddStatus& status) {
        node += trie_base_shift<base_>();
        ++i;
        if (i == my_length) {
            end(node, status);
        } else {
            node = next(node);
            recursive_add(i, node, barcode_seq, status);
        }
    }

    void recursive_add(SeqLength i, TrieIndex node, const char* barcode_seq, TrieAddStatus& status) {
        // Processing a stretch of non-ambiguous codes, where possible.
        // This reduces the recursion depth among the (hopefully fewer) ambiguous codes.
        while (1) {
            switch (barcode_seq[i]) {
                case 'A': case 'a':
                    node += trie_base_shift<'A'>(); break;
                case 'C': case 'c':
                    node += trie_base_shift<'C'>(); break;
                case 'G': case 'g':
                    node += trie_base_shift<'G'>(); break;
                case 'T': case 't':
                    node += trie_base_shift<'T'>(); break;
                default:
                    goto ambiguous;
            }
            if ((++i) == my_length) {
                end(node, status);
                return;
            } else {
                node = next(node);
            }
        } 

ambiguous:
        // Processing the ambiguous codes.
        status.has_ambiguous = true;

        auto processA = [&]() -> void { process_ambiguous<'A'>(i, node, barcode_seq, status); };
        auto processC = [&]() -> void { process_ambiguous<'C'>(i, node, barcode_seq, status); };
        auto processG = [&]() -> void { process_ambiguous<'G'>(i, node, barcode_seq, status); };
        auto processT = [&]() -> void { process_ambiguous<'T'>(i, node, barcode_seq, status); };

        switch(barcode_seq[i]) {
            case 'R': case 'r':
                processA(); processG(); break;
            case 'Y': case 'y':
                processC(); processT(); break;
            case 'S': case 's':
                processC(); processG(); break;
            case 'W': case 'w':
                processA(); processT(); break;
            case 'K': case 'k':
                processG(); processT(); break;
            case 'M': case 'm':
                processA(); processC(); break;
            case 'B': case 'b':
                processC(); processG(); processT(); break;
            case 'D': case 'd':
                processA(); processG(); processT(); break;
            case 'H': case 'h':
                processA(); processC(); processT(); break;
            case 'V': case 'v':
                processA(); processC(); processG(); break;
            case 'N': case 'n':
                processA(); processC(); processG(); processT(); break;
            default:
                throw std::runtime_error("unknown base '" + std::string(1, barcode_seq[i]) + "' detected when constructing the trie");
        }
    }

public:
    TrieAddStatus add(const char* barcode_seq) {
        TrieAddStatus status;
        recursive_add(0, 0, barcode_seq, status);
        ++my_counter;
        return status;
    }

    SeqLength length() const {
        return my_length;
    }

    TrieIndex size() const {
        return my_counter;
    }

    const std::vector<TrieIndex>& pointers() const {
        return my_pointers;
    }

public:
    // To be called in the middle steps of the recursive search (i.e., for all but the last position).
    template<class SearchResult_>
    void replace_best_with_chosen(SearchResult_& best, TrieIndex& best_index, SeqLength best_score, const SearchResult_& chosen, TrieIndex chosen_index, SeqLength chosen_score) const {
        if (is_trie_index_ok(chosen_index)) {
            if (chosen_score < best_score) {
                best = chosen;
            } else if (chosen_score == best_score) { 
                if (chosen_index != best_index) { // protect against multiple occurrences of IUPAC code-containing barcodes.
                    if (my_duplicates == DuplicateAction::FIRST) {
                        if (chosen_index < best_index) {
                            best_index = chosen_index;
                        }
                    } else if (my_duplicates == DuplicateAction::LAST) {
                        if (chosen_index > best_index) {
                            best_index = chosen_index;
                        }
                    } else {
                        best_index = TRIE_STATUS_AMBIGUOUS; 
                    }
                }
            }

        } else if (chosen_index == TRIE_STATUS_AMBIGUOUS) {
            if (chosen_score < best_score) {
                best = chosen;
            } else if (chosen_score == best_score) {
                // Ambiguity is infectious. Each ambiguous status indicates that there
                // are already 2+ barcodes on this score, so it doesn't matter how
                // many other unambiguous barcodes are here; we're already ambiguous.
                best_index = TRIE_STATUS_AMBIGUOUS;
            }
        }
    }

    // To be called in the last step of the recursive search.
    void scan_final_position_with_mismatch(TrieIndex node, int refshift, TrieIndex& current_index, SeqLength current_mismatches, SeqLength& mismatch_cap) const {
        bool found = false;
        for (int s = 0; s < NUM_BASES; ++s) {
            if (s == refshift) { 
                continue;
            }

            auto candidate = my_pointers[node + s];
            if (is_trie_index_ok(candidate)) {
                if (found) { 
                    if (candidate != current_index) { // protect against multiple occurrences of IUPAC-containg barcodes.
                        if (my_duplicates == DuplicateAction::FIRST) {
                            if (current_index > candidate) {
                                current_index = candidate;
                            }
                        } else if (my_duplicates == DuplicateAction::LAST) {
                            if (current_index < candidate) {
                                current_index = candidate;
                            }
                        } else {
                            current_index = TRIE_STATUS_AMBIGUOUS; // ambiguous, so we quit early.
                            break;
                        }
                    }
                } else {
                    current_index = candidate;
                    mismatch_cap = current_mismatches;
                    found = true;
                }

            } else if (candidate == TRIE_STATUS_AMBIGUOUS) {
                // If an ambiguity is present on a base at the last position, 
                // and we're accepting a mismatch on the last position, then 
                // we already have at least two known barcodes that match the 
                // input sequence. The behavior of the other bases is irrelevant;
                // even if they refer to non-ambiguous barcodes, that just adds to the
                // set of 2+ barcodes that the input sequence already matches.
                // So, we have no choice but to fail the match due to ambiguity.
                current_index = TRIE_STATUS_AMBIGUOUS;
                mismatch_cap = current_mismatches;
                break;
            }
        }
    }

public:
    void optimize() {
        SeqLength maxed = 0;
        if (!is_optimal(0, 0, maxed)) {
            std::vector<TrieIndex> replacement;
            replacement.reserve(my_pointers.size());
            optimize(0, 0, replacement);
            my_pointers.swap(replacement);
        }
    }

private:
    // Optimization involves reorganizing the nodes so that the pointers are
    // always increasing. This promotes memory locality of similar sequences
    // in a depth-first search (which is what search() does anyway).
    bool is_optimal(SeqLength i, TrieIndex node, TrieIndex& maxed) const {
        ++i;
        if (i < my_length) {
            for (int s = 0; s < NUM_BASES; ++s) {
                auto v = my_pointers[node + s];
                if (!is_trie_index_ok(v)) {
                    continue;
                }

                if (v < maxed) {
                    return false;
                }

                maxed = v;
                if (!is_optimal(i, v, maxed)) {
                    return false;
                }
            }
        }
        return true;
    }

    void optimize(SeqLength i, TrieIndex node, std::vector<TrieIndex>& trie) const {
        auto it = my_pointers.begin() + node;
        TrieIndex new_node = trie.size();
        trie.insert(trie.end(), it, it + NUM_BASES);

        ++i;
        if (i < my_length) {
            for (int s = 0; s < NUM_BASES; ++s) {
                auto& v = trie[new_node + s];
                if (!is_trie_index_ok(v)) {
                    continue;
                }

                auto original = v;
                v = trie.size();
                optimize(i, original, trie);
            }
        }
    }
};
/**
 * @endcond
 */

/**
 * @brief Search for barcodes with mismatches anywhere.
 *
 * Given a (typically read-derived) sequence, this class performs a mismatch-aware search to a trie containing a pool of known barcode sequences.
 * It will then return the barcode with the fewest mismatches to the input sequence.
 * Any number of mismatches are supported, distributed anywhere throughout the sequence. 
 */
class AnyMismatches {
public:
    /**
     * Default constructor.
     * This is only provided to enable composition, the resulting object should not be used until it is copy-assigned to a properly constructed instance.
     */
    AnyMismatches() = default;

    /**
     * @param barcode_length Length of the barcode sequences.
     * @param duplicates How duplicate sequences across `add()` calls should be handled.
     */
    AnyMismatches(SeqLength barcode_length, DuplicateAction duplicates) : my_core(barcode_length, duplicates) {}

private:
    MismatchTrie my_core;

public:
    /**
     * @param[in] barcode_seq Pointer to a character array containing a barcode sequence.
     * The array should have length equal to `length()` and should only contain IUPAC nucleotides or their lower-case equivalents (excepting U or gap characters).
     *
     * @return The barcode sequence is added to the trie.
     * The index of the newly added sequence is defined as the number of sequences that were previously added. 
     * The status of the addition is returned.
     */
    TrieAddStatus add(const char* barcode_seq) {
        return my_core.add(barcode_seq);
    }

    /**
     * @return The length of the barcode sequences.
     */
    SeqLength length() const {
        return my_core.length();
    }

    /**
     * @return The number of barcode sequences added.
     */
    TrieIndex size() const {
        return my_core.size();
    }

    /**
     * Attempt to optimize the trie for more cache-friendly look-ups.
     * This is not necessary if sorted sequences are supplied in `add()`.
     */
    void optimize() {
        my_core.optimize();
    }

public:
    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches in the search.
     *
     * @return Pair containing:
     * 1. The index of the barcode sequence with the lowest number of mismatches to `search_seq`.
     *    This is only non-negative if a barcode was unambiguously matched with no more mismatches than `max_mismatches`.
     *    If multiple sequences have the same lowest number of mismatches, the match is ambiguous and `TRIE_STATUS_AMBIGUOUS` is returned.
     *    If all sequences have more mismatches than `max_mismatches`, `TRIE_STATUS_MISSING` is returned.
     * 2. The number of mismatches.
     */
    std::pair<TrieIndex, SeqLength> search(const char* search_seq, SeqLength max_mismatches) const {
        return search(search_seq, 0, 0, 0, max_mismatches);
    }

private:
    std::pair<TrieIndex, SeqLength> search(const char* seq, SeqLength i, TrieIndex node, SeqLength mismatches, SeqLength& max_mismatches) const {
        const auto& pointers = my_core.pointers();
        TrieIndex current;
        int shift;
        switch (seq[i]) {
            case 'A': case 'a':
                shift = trie_base_shift<'A'>(); current = pointers[node + shift]; break;
            case 'C': case 'c':
                shift = trie_base_shift<'C'>(); current = pointers[node + shift]; break;
            case 'G': case 'g':
                shift = trie_base_shift<'G'>(); current = pointers[node + shift]; break;
            case 'T': case 't':
                shift = trie_base_shift<'T'>(); current = pointers[node + shift]; break;
            default:
                shift = -1; current = TRIE_STATUS_MISSING; break;
        }

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        SeqLength length = my_core.length();
        ++i;

        if (i == length) {
            if (is_trie_index_ok(current) || current == TRIE_STATUS_AMBIGUOUS) {
                max_mismatches = mismatches; // this assignment should always decrease max_mismatches, otherwise the search would have terminated earlier.
                return std::make_pair(current, mismatches);
            }

            TrieIndex alt = TRIE_STATUS_MISSING;
            ++mismatches;
            if (mismatches <= max_mismatches) {
                my_core.scan_final_position_with_mismatch(node, shift, alt, mismatches, max_mismatches);
            }

            return std::make_pair(alt, mismatches);

        } else {
            std::pair<TrieIndex, SeqLength> best(TRIE_STATUS_MISSING, max_mismatches + 1);
            if (is_trie_index_ok(current)) {
                best = search(seq, i, current, mismatches, max_mismatches);
            }

            ++mismatches;
            if (mismatches <= max_mismatches) {
                for (int s = 0; s < NUM_BASES; ++s) {
                    if (shift == s) { 
                        continue;
                    } 

                    auto alt = pointers[node + s];
                    if (!is_trie_index_ok(alt)) {
                        continue;
                    }

                    if (mismatches <= max_mismatches) { // check again, just in case max_mismatches changed.
                        auto chosen = search(seq, i, alt, mismatches, max_mismatches);
                        my_core.replace_best_with_chosen(best, best.first, best.second, chosen, chosen.first, chosen.second);
                    }
                }
            }

            return best;
        }
    }
};

/**
 * @brief Search for barcodes with segmented mismatches.
 *
 * Given a (typically read-derived) sequence, this class will perform a mismatch-aware search to a pool of known barcode sequences.
 * However, the distribution of mismatches is restricted in different segments of the sequence, e.g., 1 mismatch in the first 4 bp, 3 mismatches for the next 10 bp, and so on.
 * The aim is to enable searching for concatenations of variable region sequences (and barcodes), where each segment is subject to a different number of mismatches.
 * The barcode with the fewest mismatches to the input sequence is then returned.
 *
 * @tparam num_segments_ Number of segments to consider.
 */
template<int num_segments_>
class SegmentedMismatches {
public:
    /**
     * Default constructor.
     * This is only provided to enable composition, the resulting object should not be used until it is copy-assigned to a properly constructed instance.
     */
    SegmentedMismatches() = default;

    /**
     * @param segments Length of each segment of the sequence.
     * Each entry should be positive and the sum should be equal to the total length of the barcode sequence.
     * @param duplicates How duplicate sequences across `add()` calls should be handled.
     */
    SegmentedMismatches(std::array<SeqLength, num_segments_> segments, DuplicateAction duplicates) : 
        my_core(std::accumulate(segments.begin(), segments.end(), 0), duplicates), 
        my_boundaries(segments)
    {
        for (int i = 1; i < num_segments_; ++i) {
            my_boundaries[i] += my_boundaries[i-1];
        }
    }

private:
    MismatchTrie my_core;
    std::array<SeqLength, num_segments_> my_boundaries;

public:
    /**
     * @param[in] barcode_seq Pointer to a character array containing a barcode sequence.
     * The array should have length equal to `length()` and should only contain IUPAC nucleotides or their lower-case equivalents (excepting U or gap characters).
     *
     * @return The barcode sequence is added to the trie.
     * The index of the newly added sequence is defined as the number of sequences that were previously added. 
     * The status of the addition is returned.
     */
    TrieAddStatus add(const char* barcode_seq) {
        return my_core.add(barcode_seq);
    }

    /**
     * @return The length of the barcode sequences.
     */
    SeqLength length() const {
        return my_core.length();
    }

    /**
     * @return The number of barcode sequences added.
     */
    TrieIndex size() const {
        return my_core.size();
    }

    /**
     * Attempt to optimize the trie for more cache-friendly look-ups.
     * This is not necessary if sorted sequences are supplied in `add()`.
     */
    void optimize() {
        my_core.optimize();
    }

public:
    /**
     * @brief Result of the segmented search.
     */
    struct Result {
        /**
         * @cond
         */
        Result() : per_segment() {}
        /**
         * @endcond
         */

        /**
         * Index of the known barcode sequence matching the input sequence in `search()`.
         * This is guaranteed to be non-negative only if an unambiguous match is found.
         */
        TrieIndex index = 0;

        /**
         * Total number of mismatches between the barcode sequence from `index` and the input sequence.
         */
        SeqLength total = 0;

        /**
         * Number of mismatches in each segment of the sequence.
         */
        std::array<SeqLength, num_segments_> per_segment;
    };

    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches for each segment.
     * Each entry should be non-negative.
     *
     * @return A `Result` containing the index of the barcode sequence where the number of mismatches in each segment is less than or equal to `max_mismatches`.
     * - If multiple barcode sequences satisfy this condition, the barcode sequence with the lowest total number of mismatches is reported in the form of a non-negative `Result::index`.
     * - If multiple barcode sequences share the same lowest total, the match is ambiguous and `TRIE_STATUS_AMBIGUOUS` is reported.
     * - If no barcode sequences satisfy the `max_mismatches` condition, `TRIE_STATUS_MISSING` is reported.
     */
    Result search(const char* search_seq, const std::array<SeqLength, num_segments_>& max_mismatches) const {
        SeqLength total_mismatches = std::accumulate(max_mismatches.begin(), max_mismatches.end(), static_cast<SeqLength>(0));
        return search(search_seq, 0, 0, Result(), max_mismatches, total_mismatches);
    }

private:
    Result search(const char* seq, SeqLength i, TrieIndex segment_id, Result state, const std::array<SeqLength, num_segments_>& segment_mismatches, SeqLength& total_mismatches) const {
        // Note that, during recursion, state.index does double duty 
        // as the index of the node on the trie.
        auto node = state.index;

        const auto& pointers = my_core.pointers();
        TrieIndex current;
        int shift;
        switch (seq[i]) {
            case 'A': case 'a':
                shift = trie_base_shift<'A'>(); current = pointers[node + shift]; break;
            case 'C': case 'c':
                shift = trie_base_shift<'C'>(); current = pointers[node + shift]; break;
            case 'G': case 'g':
                shift = trie_base_shift<'G'>(); current = pointers[node + shift]; break;
            case 'T': case 't':
                shift = trie_base_shift<'T'>(); current = pointers[node + shift]; break;
            default:
                shift = -1; current = TRIE_STATUS_MISSING; break;
        }

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        SeqLength length = my_core.length();
        ++i;

        if (i == length) {
            if (is_trie_index_ok(current) || current == TRIE_STATUS_AMBIGUOUS) {
                total_mismatches = state.total; // this assignment should always decrease total_mismatches, otherwise the search would have terminated earlier.
                state.index = current;
                return state;
            }

            state.index = TRIE_STATUS_MISSING;
            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
                my_core.scan_final_position_with_mismatch(node, shift, state.index, state.total, total_mismatches);
            }

            return state;

        } else {
            auto next_segment_id = segment_id;
            if (i == my_boundaries[segment_id]) {
                ++next_segment_id;
            }

            Result best;
            best.index = TRIE_STATUS_MISSING;
            best.total = total_mismatches + 1;

            if (is_trie_index_ok(current)) {
                state.index = current;
                best = search(seq, i, next_segment_id, state, segment_mismatches, total_mismatches);
            }

            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
                for (int s = 0; s < NUM_BASES; ++s) {
                    if (shift == s) { 
                        continue;
                    } 

                    auto alt = pointers[node + s];
                    if (!is_trie_index_ok(alt)) {
                        continue;
                    }

                    if (state.total <= total_mismatches) { // check again, just in case total_mismatches changed.
                        state.index = alt;
                        auto chosen = search(seq, i, next_segment_id, state, segment_mismatches, total_mismatches);
                        my_core.replace_best_with_chosen(best, best.index, best.total, chosen, chosen.index, chosen.total);
                    }
                }
            }

            return best;
        }
    }
};

}

#endif
