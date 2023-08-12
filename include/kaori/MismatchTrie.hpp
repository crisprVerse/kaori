#ifndef KAORI_MISMATCH_TRIE_HPP
#define KAORI_MISMATCH_TRIE_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <numeric>
#include "utils.hpp"
#include "BarcodePool.hpp"

/**
 * @file MismatchTrie.hpp
 *
 * @brief Defines the `MismatchTrie` class and its subclasses.
 */

namespace kaori {

/**
 * @brief Base class for the mismatch search.
 *
 * Given a (typically read-derived) sequence, this class will perform a mismatch-aware search to a pool of known barcode sequences.
 * It will then return the barcode with the fewest mismatches to the input sequence.
 * Any number of mismatches are supported; subclasses will decide how the mismatches can be distributed throughout the length of the sequence.
 */
class MismatchTrie {
protected:
    static constexpr int status_not_present = -1;

    static constexpr int status_ambiguous = -2;

public:
    /**
     * @param barcode_length Length of the barcodes in the pool.
     */
    MismatchTrie(size_t barcode_length = 0) : length(barcode_length), pointers(4, status_not_present), counter(0) {}

    /**
     * @param barcode_pool Pool of known barcode sequences.
     * @param duplicates How duplicated sequences in `barcode_pool` should be handled.
     */
    MismatchTrie(const BarcodePool& barcode_pool, DuplicateAction duplicates = DuplicateAction::ERROR) : MismatchTrie(barcode_pool.length) {
        for (auto s : barcode_pool.pool) {
            add(s, duplicates);
        }
    }

public:
    /** 
     * @brief Status of the barcode sequence addition.
     */
    struct AddStatus {
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

private:
    void next(int shift, int& position) {
        auto& current = pointers[position + shift];
        if (current < 0) {
            current = pointers.size();
            position = current;
            pointers.resize(position + 4, status_not_present);
        } else {
            position = current;
        }
    }

    void end(int shift, int position, DuplicateAction duplicates, AddStatus& status) {
        auto& current = pointers[position + shift];
        if (current >= 0) {
            status.is_duplicate = true;
            switch(duplicates) {
                case DuplicateAction::FIRST:
                    break;
                case DuplicateAction::LAST:
                    status.duplicate_replaced = true;
                    current = counter;
                    break;
                case DuplicateAction::NONE:
                    status.duplicate_cleared = true;
                    current = status_ambiguous;
                    break;
                case DuplicateAction::ERROR:
                    throw std::runtime_error("duplicate sequences detected (" + 
                        std::to_string(current + 1) + ", " + 
                        std::to_string(counter + 1) + ") when constructing the trie");
            }
        } else if (current == status_not_present) {
            current = counter;
        } else if (current == status_ambiguous) {
            status.is_duplicate = true;
        }
    }

    void recursive_add(size_t i, int position, const char* barcode_seq, DuplicateAction duplicates, AddStatus& status) {
        // Processing a stretch of non-ambiguous codes, where possible.
        // This reduces the recursion depth among the (hopefully fewer) ambiguous codes.
        while (1) {
            auto shift = base_shift<true>(barcode_seq[i]);
            if (shift == -1) {
               break;
            }

            if ((++i) == length) {
                end(shift, position, duplicates, status);
                return;
            } else {
                next(shift, position);
            }
        } 

        // Processing the ambiguous codes.
        status.has_ambiguous = true;

        auto process = [&](char base) -> void {
            auto shift = base_shift(base);
            if (i + 1 == length) {
                end(shift, position, duplicates, status);
            } else {
                auto curpos = position;
                next(shift, curpos);
                recursive_add(i + 1, curpos, barcode_seq, duplicates, status);
            }
        };

        switch(barcode_seq[i]) {
            case 'R': case 'r':
                process('A'); process('G'); break;
            case 'Y': case 'y':
                process('C'); process('T'); break;
            case 'S': case 's':
                process('C'); process('G'); break;
            case 'W': case 'w':
                process('A'); process('T'); break;
            case 'K': case 'k':
                process('G'); process('T'); break;
            case 'M': case 'm':
                process('A'); process('C'); break;
            case 'B': case 'b':
                process('C'); process('G'); process('T'); break;
            case 'D': case 'd':
                process('A'); process('G'); process('T'); break;
            case 'H': case 'h':
                process('A'); process('C'); process('T'); break;
            case 'V': case 'v':
                process('A'); process('C'); process('G'); break;
            case 'N': case 'n':
                process('A'); process('C'); process('G'); process('T'); break;
            default:
                throw std::runtime_error("unknown base '" + std::string(1, barcode_seq[i]) + "' detected when constructing the trie");
        }
    }

public:
    /**
     * @param[in] barcode_seq Pointer to a character array containing a barcode sequence.
     * The array should have length equal to `get_length()` and should only contain IUPAC nucleotides or their lower-case equivalents (excepting U or gap characters).
     * @param duplicates How duplicate sequences across `add()` calls should be handled.
     *
     * @return The barcode sequence is added to the trie.
     * The index of the newly added sequence is defined as the number of sequences that were previously added. 
     * The status of the addition is returned.
     */
    AddStatus add(const char* barcode_seq, DuplicateAction duplicates) {
        AddStatus status;
        recursive_add(0, 0, barcode_seq, duplicates, status);
        ++counter;
        return status;
    }

public:
    /**
     * @return The length of the barcode sequences.
     */
    size_t get_length() const {
        return length;
    }

    /**
     * @return The number of barcode sequences added.
     */
    int size() const {
        return counter;
    }

protected:
    /**
     * @cond
     */
    size_t length;
    std::vector<int> pointers;

    template<bool allow_unknown = false>
    static int base_shift(char base) {
        int shift = 0;
        switch (base) {
            case 'A': case 'a':
                break;
            case 'C': case 'c':
                shift = 1;
                break;
            case 'G': case 'g':
                shift = 2;
                break;
            case 'T': case 't':
                shift = 3;
                break;
            default:
                if constexpr(allow_unknown) {
                    shift = -1; 
                } else {
                    throw std::runtime_error("unknown base '" + std::string(1, base) + "' detected when constructing the trie");
                }
        }
        return shift;
    }
    /**
     * @endcond
     */

private:
    int counter;

public:
    /**
     * Attempt to optimize the trie for more cache-friendly look-ups.
     * This is not necessary if sorted sequences are supplied in `add()`.
     */
    void optimize() {
        int maxed = 0;
        if (!is_optimal(0, 0, maxed)) {
            std::vector<int> replacement;
            replacement.reserve(pointers.size());
            optimize(0, 0, replacement);
            pointers.swap(replacement);
        }
    }

private:
    // Optimization involves reorganizing the nodes so that the pointers are
    // always increasing. This promotes memory locality of similar sequences
    // in a depth-first search (which is what search() does anyway).
    bool is_optimal(int node, size_t pos, int& maxed) const {
        ++pos;
        if (pos < length) {
            for (int s = 0; s < 4; ++s) {
                auto v = pointers[node + s];
                if (v < 0) {
                    continue;
                }

                if (v < maxed) {
                    return false;
                }

                maxed = v;
                if (!is_optimal(v, pos, maxed)) {
                    return false;
                }
            }
        }
        return true;
    }

    void optimize(int node, size_t pos, std::vector<int>& trie) const {
        auto it = pointers.begin() + node;
        size_t new_node = trie.size();
        trie.insert(trie.end(), it, it + 4);

        ++pos;
        if (pos < length) {
            for (int s = 0; s < 4; ++s) {
                auto& v = trie[new_node + s];
                if (v < 0) {
                    continue;
                }

                auto original = v;
                v = trie.size();
                optimize(original, pos, trie);
            }
        }
    }
};

/**
 * @brief Search for barcodes with mismatches anywhere.
 *
 * This `MismatchTrie` subclass will search for the best match to known sequences in a barcode pool.
 * Any number of mismatches are supported, distributed anywhere throughout the sequence. 
 */
class AnyMismatches : public MismatchTrie {
public:
    /**
     * @param barcode_length Length of the barcode sequences.
     */
    AnyMismatches(size_t barcode_length = 0) : MismatchTrie(barcode_length) {}

    /**
     * @param barcode_pool Pool of known barcode sequences.
     * @param duplicates How duplicate sequences in `barcode_pool` should be handled.
     */
    AnyMismatches(const BarcodePool& barcode_pool, DuplicateAction duplicates = DuplicateAction::ERROR) : 
        MismatchTrie(barcode_pool, duplicates) {}

public:
    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `get_length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches in the search.
     *
     * @return Pair containing:
     * 1. The index of the barcode sequence with the lowest number of mismatches to `search_seq`.
     *    If multiple sequences have the same lowest number of mismatches, the match is ambiguous and -1 is returned.
     *    If all sequences have more mismatches than `max_mismatches`, -1 is returned.
     * 2. The number of mismatches.
     */
    std::pair<int, int> search(const char* search_seq, int max_mismatches) const {
        return search(search_seq, 0, 0, 0, max_mismatches);
    }

private:
    std::pair<int, int> search(const char* seq, size_t pos, int node, int mismatches, int& max_mismatches) const {
        int shift = base_shift<true>(seq[pos]);
        int current = (shift >= 0 ? pointers[node + shift] : -1);

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        if (pos + 1 == length) {
            if (current >= 0) {
                max_mismatches = mismatches;
                return std::make_pair(current, mismatches);
            }

            int alt = -1;
            ++mismatches;
            if (mismatches <= max_mismatches) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    }

                    int candidate = pointers[node + s];
                    if (candidate >= 0) {
                        if (found) { // ambiguous, so we quit early.
                            alt = -1;
                            break;
                        }
                        alt = candidate;
                        max_mismatches = mismatches;
                        found = true;
                    }
                }
            }
            return std::make_pair(alt, mismatches);

        } else {
            ++pos;

            std::pair<int, int> best(-1, max_mismatches + 1);
            if (current >= 0) {
                best = search(seq, pos, current, mismatches, max_mismatches);
            }

            ++mismatches;
            if (mismatches <= max_mismatches) {
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    } 
                    
                    int alt = pointers[node + s];
                    if (alt < 0) {
                        continue;
                    }

                    auto chosen = search(seq, pos, alt, mismatches, max_mismatches);
                    if (chosen.second < best.second) {
                        best = chosen;
                    } else if (chosen.second == best.second) {
                        best.first = -1;
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
 * This `MismatchTrie` subclass will search for the best match to known sequences in a barcode pool.
 * However, the distribution of mismatches is restricted in different segments of the sequence, e.g., 1 mismatch in the first 4 bp, 3 mismatches for the next 10 bp, and so on.
 * The intention is to enable searching for concatenations of variable region sequences (and barcodes), where each segment is subject to a different number of mismatches.
 *
 * @tparam num_segments Number of segments to consider.
 */
template<size_t num_segments>
class SegmentedMismatches : public MismatchTrie {
public:
    /**
     * Default constructor.
     */
    SegmentedMismatches() {}

    /**
     * @param segments Length of each segment of the sequence.
     * Each entry should be positive and the sum should be equal to the total length of the barcode sequence.
     */
    SegmentedMismatches(std::array<int, num_segments> segments) : MismatchTrie(std::accumulate(segments.begin(), segments.end(), 0)), boundaries(segments) {
        for (size_t i = 1; i < num_segments; ++i) {
            boundaries[i] += boundaries[i-1];
        }
    }

    /**
     * @param barcode_pool Possible set of known sequences for the variable region.
     * @param segments Length of each segment of the sequence.
     * Each entry should be positive and the sum should be equal to the total length of the barcode sequence.
     * @param duplicates How duplicated sequences in `barcode_pool` should be handled.
     */
    SegmentedMismatches(const BarcodePool& barcode_pool, std::array<int, num_segments> segments, DuplicateAction duplicates = DuplicateAction::ERROR) : 
        SegmentedMismatches(segments) 
    {
        if (length != barcode_pool.length) {
            throw std::runtime_error("length of barcode sequences should equal the sum of segment lengths");
        }
        for (auto s : barcode_pool.pool) {
            add(s, duplicates);
        }
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
         * If no unambiguous match is found, -1 is reported.
         */
        int index = 0;

        /**
         * Total number of mismatches between the barcode sequence from `index` and the input sequence.
         */
        int total = 0;

        /**
         * Number of mismatches in each segment of the sequence.
         */
        std::array<int, num_segments> per_segment;
    };

    /**
     * @param[in] search_seq Pointer to a character array containing a sequence to use for searching the barcode pool.
     * This is assumed to be of length equal to `get_length()` and is typically derived from a read.
     * @param max_mismatches Maximum number of mismatches for each segment.
     * Each entry should be non-negative.
     *
     * @return A `Result` containing the index of the barcode sequence where the number of mismatches in each segment is less than or equal to `max_mismatches`.
     * - If multiple barcode sequences satisfy this condition, the barcode sequence with the lowest total number of mismatches is reported.
     * - If multiple barcode sequences share the same lowest total, the match is ambiguous and -1 is reported.
     * - If no barcode sequences satisfy the `max_mismatches` condition, -1 is reported.
     */
    Result search(const char* search_seq, const std::array<int, num_segments>& max_mismatches) const {
        int total_mismatches = std::accumulate(max_mismatches.begin(), max_mismatches.end(), 0);
        return search(search_seq, 0, 0, Result(), max_mismatches, total_mismatches);
    }

private:
    Result search(
        const char* seq, 
        size_t pos, 
        size_t segment_id,
        Result state,
        const std::array<int, num_segments>& segment_mismatches, 
        int& total_mismatches
    ) const {
        // Note that, during recursion, state.index does double duty 
        // as the index of the node on the trie.
        int node = state.index;

        int shift = base_shift<true>(seq[pos]);
        int current = (shift >= 0 ? pointers[node + shift] : -1);

        // At the end: we prepare to return the actual values. We also refine
        // the max number of mismatches so that we don't search for things with
        // more mismatches than the best hit that was already encountered.
        if (pos + 1 == length) {
            if (current >= 0) {
                total_mismatches = state.total;
                state.index = current;
                return state;
            }

            state.index = -1;
            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
                bool found = false;
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    }

                    int candidate = pointers[node + s];
                    if (candidate >= 0) {
                        if (found) { // ambiguous, so we quit early.
                            state.index = -1;
                            break;
                        }
                        state.index = candidate;
                        total_mismatches = state.total;
                        found = true;
                    }
                }
            }
            return state;

        } else {
            auto next_pos = pos + 1;
            auto next_segment_id = segment_id;
            if (static_cast<int>(next_pos) == boundaries[segment_id]) { // TODO: boundaries should probably be size_t's, thus avoiding the need for this cast.
                ++next_segment_id;
            }

            Result best;
            best.index = -1;
            best.total = total_mismatches + 1;

            if (current >= 0) {
                state.index = current;
                best = search(seq, next_pos, next_segment_id, state, segment_mismatches, total_mismatches);
            }

            ++state.total;
            auto& current_segment_mm = state.per_segment[segment_id];
            ++current_segment_mm;

            if (state.total <= total_mismatches && current_segment_mm <= segment_mismatches[segment_id]) {
                for (int s = 0; s < 4; ++s) {
                    if (shift == s) { 
                        continue;
                    } 
                    
                    int alt = pointers[node + s];
                    if (alt < 0) {
                        continue;
                    }

                    state.index = alt;
                    auto chosen = search(seq, next_pos, next_segment_id, state, segment_mismatches, total_mismatches);
                    if (chosen.total < best.total) {
                        best = chosen;
                    } else if (chosen.total == best.total) { // ambiguous
                        best.index = -1;
                    }
                }
            }

            return best;
        }
    }
private:
    std::array<int, num_segments> boundaries;
};

}

#endif
