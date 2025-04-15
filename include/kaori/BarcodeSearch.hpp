#ifndef KAORI_VARIABLE_LIBRARY_HPP
#define KAORI_VARIABLE_LIBRARY_HPP

#include "BarcodePool.hpp"
#include "MismatchTrie.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <array>

/**
 * @file BarcodeSearch.hpp
 *
 * @brief Search for barcode sequences.
 */

namespace kaori {

/** 
 * @cond
 */
template<typename Trie_>
inline void fill_library(
    const std::vector<const char*>& options, 
    std::unordered_map<std::string, int>& exact,
    Trie_& trie,
    bool reverse
) {
    size_t len = trie.length();

    for (size_t i = 0; i < options.size(); ++i) {
        auto ptr = options[i];

        std::string current;
        if (!reverse) {
            current = std::string(ptr, ptr + len);
        } else {
            current.reserve(len);
            for (size_t j = 0; j < len; ++j) {
                current += complement_base<true, true>(ptr[len - j - 1]);
            }
        }

        // Note that this must be called, even if the sequence is duplicated;
        // otherwise the trie's internal counter will not be properly incremented.
        auto status = trie.add(current.c_str());

        if (!status.has_ambiguous) {
            if (!status.is_duplicate || status.duplicate_replaced) {
                exact[current] = i;
            } else if (status.duplicate_cleared) {
                exact[current] = -1;
            }
        }
    }

    trie.optimize();
    return;
}

template<class Methods_, class Cache_, class Trie_, class Result_, class Mismatch_>
void matcher_in_the_rye(const std::string& x, const Cache_& cache, const Trie_& trie, Result_& res, const Mismatch_& mismatches, const Mismatch_& max_mismatches) {
    // Seeing if it's any of the caches; otherwise searching the trie.
    auto cit = cache.find(x);
    if (cit == cache.end()) {
        auto lit = res.cache.find(x);
        if (lit != res.cache.end()) {
            Methods_::update(res, lit->second, mismatches);

        } else {
            auto missed = trie.search(x.c_str(), mismatches);

            // The trie search breaks early when it hits the mismatch cap,
            // but the cap might be different across calls. If we break
            // early and report a miss, the miss will be cached and
            // returned in cases where there is a higher cap (and thus
            // might actually be a hit). As such, we should only store a
            // miss in the cache when the requested number of mismatches is
            // equal to the maximum value specified in the constructor.
            if (Methods_::index(missed) >= 0 || mismatches == max_mismatches) {
                res.cache[x] = missed;
            }

            // No need to pass the requested number of mismatches,
            // as we explicitly searched for that in the trie.
            Methods_::update(res, missed);
        }
    } else {
        Methods_::update(res, cit->second, mismatches);
    }
    return;
}
/** 
 * @endcond
 */

/**
 * @brief Search for known barcode sequences.
 *
 * This supports exact and mismatch-aware searches for known sequences.
 * Mismatches may be distributed anywhere along the length of the sequence, see `AnyMismatches` for details.
 * Instances of this class use caching to avoid redundant work when a mismatching sequence has been previously encountered.
 */
class SimpleBarcodeSearch {
public:
    /**
     * @brief Optional parameters for `SimpleBarcodeSearch`.
     */
    struct Options {
        /**
         * Maximum number of mismatches for any search performed by `SimpleBarcodeSearch::search`.
         */
        int max_mismatches = 0;

        /** 
         * Whether to reverse-complement the barcode sequences before indexing them.
         */
        bool reverse = false;

        /** 
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;
    };

public:
    /**
     * Default constructor.
     * This is only provided for composition purposes; methods of this class should only be called on properly constructed instance.
     */
    SimpleBarcodeSearch() = default;

    /**
     * @param barcode_pool Pool of barcode sequences.
     * @param options Optional parameters for the search.
     */
    SimpleBarcodeSearch(const BarcodePool& barcode_pool, const Options& options) : 
        my_trie(barcode_pool.length(), options.duplicates), 
        my_max_mm(options.max_mismatches) 
    {
        fill_library(barcode_pool.pool(), my_exact, my_trie, options.reverse);
        return;
    }

public:
    /**
     * @brief State of the search.
     *
     * This contains both the results of the search for any given input sequence,
     * as well as cached mismatches to optimize searches for future inputs.
     */
    struct State {
        /**
         * Index of the known sequence that matches best to the input sequence in `search()` (i.e., fewest total mismatches).
         * If no match was found or if the best match is ambiguous, this will be set to -1.
         */
        int index = 0;

        /**
         * Number of mismatches with the matching known sequence.
         * This should only be used if `index != -1`.
         */
        int mismatches = 0;
        
        /**
         * @cond
         */
        std::unordered_map<std::string, std::pair<int, int> > cache;
        /**
         * @endcond
         */
    };

    /**
     * Initialize the search state for thread-safe execution.
     *
     * @return A new `SeachState()`.
     */
    State initialize() const {
        return State();
    }

    /**
     * Incorporate the mismatch cache from `state` into the cache for this `SimpleBarcodeSearch` instance.
     * This allows regular consolidation of optimizations across threads.
     * On return, the mismatch cache of `state` is combined with that of this instance.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has already been used in `search()` at least once.
     */
    void reduce(State& state) {
        my_cache.merge(state.cache);
        state.cache.clear();
    }

private:
    struct Methods {
        static int index(const std::pair<int, int>& val) {
            return val.first;
        }

        static void update(State& state, const std::pair<int, int>& val) {
            state.index = val.first;
            state.mismatches = val.second;
            return;
        }

        static void update(State& state, const std::pair<int, int>& val, int mismatches) {
            state.index = (val.second > mismatches ? -1 : val.first);
            state.mismatches = val.second;
            return;
        }
    };

public:
    /**
     * Search the known sequences in the barcode pool for an input sequence.
     * The number of allowed mismatches is equal to the maximum specified in the constructor.
     * 
     * @param search_seq The input sequence to use for searching.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * On return, `state` is filled with the details of the best-matching barcode sequence, if any exists.
     */
    void search(const std::string& search_seq, State& state) const {
        search(search_seq, state, my_max_mm);
        return;
    }

    /**
     * Search the known sequences in the barcode pool for an input sequence,
     * with potentially more stringent mismatch requirements.
     * This can improve efficiency in situations where some mismatches have already been consumed by matching the template sequence.
     *
     * @param search_seq The input sequence to use for searching.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * On return, `state` is filled with the details of the best-matching barcode sequence, if any exists.
     * @param allowed_mismatches Allowed number of mismatches.
     * This should not be greater than the maximum specified in the constructor.
     */
    void search(const std::string& search_seq, State& state, int allowed_mismatches) const {
        auto it = my_exact.find(search_seq);
        if (it != my_exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
        } else {
            matcher_in_the_rye<Methods>(search_seq, my_cache, my_trie, state, allowed_mismatches, my_max_mm);
        }
    }

private:
    std::unordered_map<std::string, int> my_exact;
    AnyMismatches my_trie;
    std::unordered_map<std::string, std::pair<int, int> > my_cache;
    int my_max_mm;
};

/**
 * @brief Search for known barcode sequences with segmented mismatches.
 *
 * This supports exact and mismatch-aware searches for known sequences.
 * Mismatches are restricted by segments along the sequence, see `SegmentedMismatches` for details. 
 * Instances of this class use caching to avoid redundant work when a mismatching sequence has been previously encountered.
 *
 * @tparam num_segments_ Number of segments to consider.
 */
template<size_t num_segments_>
class SegmentedBarcodeSearch {
public:
    /**
     * @brief Optional parameters for a `SegmentedBarcodeSearch`.
     */
    struct Options {
        /**
         * @param max_mismatch_per_segment Maximum number of mismatches per segment.
         * This is used to fill `max_mismatches`.
         */
        Options(int max_mismatch_per_segment = 0) {
            max_mismatches.fill(max_mismatch_per_segment);
        }
        
        /**
         * Maximum number of mismatches in each segment for `SegmentedBarcodeSearch::search()`.
         * All values should be non-negative.
         * Defaults to an all-zero array in the `Options()` constructor.
         */
        std::array<int, num_segments_> max_mismatches;

        /** 
         * Whether to reverse-complement the barcode sequences before indexing them.
         * Note that, even if `reverse = true`, the segment lengths in the `SegmentedBarcodeSearch()` constructor and `max_mismatches` are still reported in their order on the forward strand.
         */
        bool reverse = false;

        /** 
         * How duplicated barcode sequences should be handled.
         */
        DuplicateAction duplicates = DuplicateAction::ERROR;
    };

public:
    /**
     * Default constructor.
     */
    SegmentedBarcodeSearch() = default;

    /**
     * @param barcode_pool Pool of barcode sequences.
     * @param segments Size of each segment.
     * All values should be positive and their sum should be equal to the barcode length.
     * @param options Optional parameters.
     */
    SegmentedBarcodeSearch(
        const BarcodePool& barcode_pool, 
        std::array<int, num_segments_> segments, 
        const Options& options
    ) : 
        trie(
            (!options.reverse ? 
                segments :
                [&]{
                    auto copy = segments;
                    std::reverse(copy.begin(), copy.end());
                    return copy;
                }()
            ),
            options.duplicates
        ), 
        max_mm(
            (!options.reverse ?
                options.max_mismatches :
                [&]{
                    auto copy = options.max_mismatches;
                    std::reverse(copy.begin(), copy.end());
                    return copy;
                }()
            )
        )
    {
        if (barcode_pool.length() != trie.length()) {
            throw std::runtime_error("variable sequences should have the same length as the sum of segment lengths");
        }
        fill_library(barcode_pool.pool(), exact, trie, options.reverse);
        return;
    }

public:
    /**
     * @brief State of the search.
     *
     * This contains both the results of the search for any given input sequence,
     * as well as cached mismatches to optimize searches for future inputs.
     */
    struct State {
        /**
         * Index of the known sequence that matches best to the input sequence in `search()` (i.e., fewest total mismatches).
         * If no match was found or if the best match is ambiguous, this will be set to -1.
         */
        int index = 0;

        /**
         * Total number of mismatches with the matching known sequence, summed across all segments.
         * This should only be used if `index != -1`.
         */
        int mismatches = 0;

        /**
         * Number of mismatches in each segment.
         * This should only be used if `index != -1`.
         */
        std::array<int, num_segments_> per_segment;
        
        /**
         * @cond
         */
        State() : per_segment() {}

        std::unordered_map<std::string, typename SegmentedMismatches<num_segments_>::Result> cache;
        /**
         * @endcond
         */
    };

    /**
     * Initialize the search state for thread-safe execution.
     *
     * @return A new `SeachState()`.
     */
    State initialize() const {
        return State();
    }

    /**
     * Incorporate the mismatch cache from `state` into the cache for this `SimpleBarcodeSearch` instance.
     * This allows regular consolidation of optimizations across threads.
     * On return, the mismatch cache of `state` is combined with that of this instance.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has already been used in `search()`.
     */
    void reduce(State& state) {
        cache.merge(state.cache);
        state.cache.clear();
    }

private:
    typedef typename SegmentedMismatches<num_segments_>::Result SegmentedResult;

    struct Methods {
        static int index(const SegmentedResult& val) {
            return val.index;
        }

        static void update(State& state, const SegmentedResult& val) {
            state.index = val.index;
            state.mismatches = val.total;
            state.per_segment = val.per_segment;
            return;
        }

        static void update(State& state, const SegmentedResult& val, const std::array<int, num_segments_>& mismatches) {
            [&]{
                for (size_t s = 0; s < num_segments_; ++s) {
                    if (val.per_segment[s] > mismatches[s]) {
                        state.index = -1;
                        return;
                    }
                }
                state.index = val.index;
            }();
            state.mismatches = val.total;
            state.per_segment = val.per_segment;
            return;
        }
    };

public:
    /**
     * Search the known sequences in the barcode pool for an input sequence.
     * The number of allowed mismatches in each segment is equal to the maximum specified in the constructor.
     * 
     * @param search_seq The input sequence to use for searching.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * On return, `state` is filled with the details of the best-matching barcode sequence, if any exists.
     */
    void search(const std::string& search_seq, State& state) const {
        search(search_seq, state, max_mm);
        return;
    }

    /**
     * Search the known sequences in the barcode pool for an input sequence,
     * with potentially more stringent mismatch requirements for each segment.
     * This can improve efficiency in situations where some mismatches have already been consumed by matching the template sequence.
     *
     * @param search_seq The input sequence to use for searching.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * On return, `state` is filled with the details of the best-matching barcode sequence, if any exists.
     * @param allowed_mismatches Allowed number of mismatches in each segment.
     * Each value should not be greater than the corresponding maximum specified in the constructor.
     */
    void search(const std::string& search_seq, State& state, std::array<int, num_segments_> allowed_mismatches) const {
        auto it = exact.find(search_seq);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
            std::fill_n(state.per_segment.begin(), num_segments_, 0);
        } else {
            matcher_in_the_rye<Methods>(search_seq, cache, trie, state, allowed_mismatches, max_mm);
        }
    }

private:
    std::unordered_map<std::string, int> exact;
    SegmentedMismatches<num_segments_> trie;
    std::unordered_map<std::string, SegmentedResult> cache;
    std::array<int, num_segments_> max_mm;
};

}

#endif
