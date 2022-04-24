#ifndef KAORI_VARIABLE_LIBRARY_HPP
#define KAORI_VARIABLE_LIBRARY_HPP

#include "SequenceSet.hpp"
#include "MismatchTrie.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <array>

/**
 * @file VariableLibrary.hpp
 *
 * @brief Define libraries for the variable sequences.
 */

namespace kaori {

/** 
 * @cond
 */
template<class Trie>
void fill_library(
    const std::vector<const char*>& options, 
    std::unordered_map<std::string, int>& exact,
    Trie& trie,
    bool reverse,
    bool duplicates
) {
    size_t len = trie.get_length();

    for (size_t i = 0; i < options.size(); ++i) {
        auto ptr = options[i];

        std::string current;
        if (!reverse) {
            current = std::string(ptr, ptr + len);
        } else {
            for (int j = 0; j < len; ++j) {
                current += reverse_complement(ptr[len - j - 1]);
            }
        }

        auto it = exact.find(current);
        if (exact.find(current) != exact.end()) {
            if (!duplicates) {
                throw std::runtime_error("duplicate variable sequence '" + current + "'");
            }
        } else {
            exact[current] = i;
        }

        // Note that this must be called, even if the sequence is duplicated;
        // otherwise the trie's internal counter will not be properly incremented.
        trie.add(current.c_str(), duplicates);
    }
    return;
}

template<class Indexer, class Updater, class Cache, class Trie, class Result, class Mismatch>
void matcher_in_the_rye(const std::string& x, const Cache& cache, const Trie& trie, Result& res, const Mismatch& mismatches, const Mismatch& max_mismatches) {
    // Seeing if it's any of the caches; otherwise searching the trie.
    auto cit = cache.find(x);
    if (cit == cache.end()) {
        auto lit = res.cache.find(x);
        if (lit != res.cache.end()) {
            Updater::update(res, lit->second);
        } else {
            auto missed = trie.search(x.c_str(), mismatches);

            // The trie search breaks early when it hits the mismatch cap,
            // but the cap might be different across calls. If we break
            // early and report a miss, the miss will be cached and
            // returned in cases where there is a higher cap (and thus
            // might actually be a hit). As such, we should only store a
            // miss in the cache when the requested number of mismatches is
            // equal to the maximum value specified in the constructor.
            if (Indexer::index(missed) >= 0 || mismatches == max_mismatches) {
                res.cache[x] = missed;
            }

            Updater::update(res, missed);
        }
    } else {
        Updater::update(res, cit->second);
    }
    return;
}
/** 
 * @endcond
 */

/**
 * @brief Library of known sequences for variable regions.
 *
 * This supports exact and mismatch-aware searches for known sequences.
 * Mismatches may be distributed anywhere along the length of the sequence.
 * Instances of this class use caching to avoid searching the `MismatchTrie` when a mismatching sequence has been previously encountered.
 */
class SimpleVariableLibrary {
public:
    /**
     * Default constructor.
     */
    SimpleVariableLibrary() {}

    /**
     * @param sequences Set of known sequences for the variable region.
     * @param mismatch Maximum number of mismatches.
     * @param reverse Whether to reverse-complement the input `sequences` before adding them.
     * @param duplicates Whether duplicate `sequences` are supported, see `MismatchTrie`.
     */
    SimpleVariableLibrary(const SequenceSet& sequences, int mismatch = 0, bool reverse = false, bool duplicates = false) : 
        trie(sequences.length), 
        max_mismatches(mismatch) 
    {
        fill_library(sequences.choices, exact, trie, reverse, duplicates);
        return;
    }

public:
    /**
     * @brief State of the search.
     *
     * This contains both the results of the search for any given input sequence,
     * as well as cached mismatches to optimize searches for future inputs.
     */
    struct SearchState {
        /**
         * Index of the known sequence that matches to the input sequence in `SimpleVariableLibrary::match()`.
         * If no match was found, this will be -1.
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
    SearchState initialize() const {
        return SearchState();
    }

    /**
     * Incorporate the mismatch cache from `state` into the cache for this `SimpleVariableLibrary` instance.
     * This allows regular consolidation of optimizations across threads.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has already been used in `match()`.
     *
     * @return The mismatch cache of `state` is combined with that of this instance.
     */
    void reduce(SearchState& state) {
        cache.merge(state.cache);
        state.cache.clear();
    }

private:
    struct Index {
        static int index(const std::pair<int, int>& val) {
            return val.first;
        }
    };

    struct Updator {
        static void update(SearchState& state, const std::pair<int, int>& val) {
            state.index = val.first;
            state.mismatches = val.second;
            return;
        }
    };

public:
    /**
     * Match an input sequence to the known sequences in the library,
     * where the number of mismatches is equal to the maximum specified in the constructor.
     * 
     * @param x The input sequence.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     *
     * @return `state` is filled with the details of the matching known sequence, if any exists.
     */
    void match(const std::string& x, SearchState& state) const {
        match(x, state, max_mismatches);
        return;
    }

    /**
     * Match an input sequence to the known sequences in the library.
     *
     * @param x The input sequence.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * @param mismatches Maximum number of mismatches.
     *
     * @return `state` is filled with the details of the matching known sequence, if any exists.
     */
    void match(const std::string& x, SearchState& state, int mismatches) const {
        auto it = exact.find(x);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
        } else {
            matcher_in_the_rye<Index, Updator>(x, cache, trie, state, mismatches, max_mismatches);
        }
    }

private:
    std::unordered_map<std::string, int> exact;
    SimpleMismatchTrie trie;
    std::unordered_map<std::string, std::pair<int, int> > cache;
    int max_mismatches;
};

/**
 * @brief Library of known sequences for variable regions with segmented mismatches.
 *
 * This supports exact and mismatch-aware searches for known sequences.
 * Mismatches are segmented, see `SegmentedMismatchTrie` for details. 
 * Instances of this class use caching to avoid searching the trie when a mismatching sequence has been previously encountered.
 */
template<size_t num_segments>
class SegmentedVariableLibrary {
public:
    /**
     * Default constructor.
     */
    SegmentedVariableLibrary() {}

    /**
     * @param sequences Set of known sequences for the variable region.
     * @param segments Size of each segment.
     * All values should be positive and their sum should be equal to the sequence length.
     * @param mismatches Maximum number of mismatches in each segment.
     * All values should be non-negative.
     * @param reverse Whether to reverse-complement the input `sequences` before adding them.
     * @param duplicates Whether duplicate `sequences` are supported, see `MismatchTrie`.
     */
    SegmentedVariableLibrary(const SequenceSet& sequences, std::array<int, num_segments> segments, std::array<int, num_segments> mismatches, bool reverse = false, bool duplicates = false) : 
        trie(segments), 
        max_mismatches(mismatches) 
    {
        if (sequences.length != trie.get_length()) {
            throw std::runtime_error("variable sequences should have the same length as the sum of segment lengths");
        }
        fill_library(sequences.choices, exact, trie, reverse, duplicates);
        return;
    }

public:
    /**
     * @brief State of the search.
     *
     * This contains both the results of the search for any given input sequence,
     * as well as cached mismatches to optimize searches for future inputs.
     */
    struct SearchState {
        /**
         * Index of the known sequence that matches to the input sequence in `SegmentedVariableLibrary::match()`.
         * If no match was found, this will be -1.
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
        std::array<int, num_segments> per_segment;
        
        /**
         * @cond
         */
        SearchState() : per_segment() {}

        std::unordered_map<std::string, typename SegmentedMismatchTrie<num_segments>::SearchResult> cache;
        /**
         * @endcond
         */
    };

    /**
     * Initialize the search state for thread-safe execution.
     *
     * @return A new `SeachState()`.
     */
    SearchState initialize() const {
        return SearchState();
    }

    /**
     * Incorporate the mismatch cache from `state` into the cache for this `SimpleVariableLibrary` instance.
     * This allows regular consolidation of optimizations across threads.
     *
     * @param state A state object generated by `initialize()`.
     * Typically this has already been used in `match()`.
     *
     * @return The mismatch cache of `state` is combined with that of this instance.
     */
    void reduce(SearchState& state) {
        cache.merge(state.cache);
        state.cache.clear();
    }

private:
    typedef typename SegmentedMismatchTrie<num_segments>::SearchResult SegmentedResult;

    struct Index {
        static int index(const SegmentedResult& val) {
            return val.index;
        }
    };

    struct Updator {
        static void update(SearchState& state, const SegmentedResult& val) {
            state.index = val.index;
            state.mismatches = val.total;
            state.per_segment = val.per_segment;
            return;
        }
    };

public:
    /**
     * Match an input sequence to the known sequences in the library,
     * where the number of mismatches in each segment is equal to the maximum specified in the constructor.
     * 
     * @param x The input sequence.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     *
     * @return `state` is filled with the details of the matching known sequence, if any exists.
     */
    void match(const std::string& x, SearchState& state) const {
        match(x, state, max_mismatches);
        return;
    }

    /**
     * Match an input sequence to the known sequences in the library.
     *
     * @param x The input sequence.
     * This is expected to have the same length as the known sequences.
     * @param state A state object generated by `initialize()`.
     * @param mismatches Maximum number of mismatches in each segment.
     *
     * @return `state` is filled with the details of the matching known sequence, if any exists.
     */
    void match(const std::string& x, SearchState& state, std::array<int, num_segments> mismatches) const {
        auto it = exact.find(x);
        if (it != exact.end()) {
            state.index = it->second;
            state.mismatches = 0;
            std::fill_n(state.per_segment.begin(), num_segments, 0);
        } else {
            matcher_in_the_rye<Index, Updator>(x, cache, trie, state, mismatches, max_mismatches);
        }
    }

private:
    std::unordered_map<std::string, int> exact;
    SegmentedMismatchTrie<num_segments> trie;
    std::unordered_map<std::string, SegmentedResult> cache;
    std::array<int, num_segments> max_mismatches;
};

}

#endif
