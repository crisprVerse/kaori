#ifndef KAORI_MATCH_SEQUENCE_HPP
#define KAORI_MATCH_SEQUENCE_HPP

#include "ConstantTemplate.hpp"
#include "MismatchTrie.hpp"
#include "utils.hpp"

namespace kaori {

template<size_t N>
class MatchSequence {
public:
    MatchSequence(const char* s, size_t n, bool f, bool r, std::vector<int> cat, const std::vector<std::vector<const char*> >& opt) : 
        num_options(opt.size()),
        forward(f), 
        reverse(r),
        constant(s, n, f, r), 
        forward_categories(std::move(cat))
    {
        const auto& regions = constant.variable_regions();
        std::vector<int> var_lengths(num_options);
        size_t nchunks = regions.size();

        for (size_t f = 0; f < forward_categories.size(); ++f) {
            auto cat = forward_categories[f];
            if (cat < 0 || static_cast<size_t>(cat) >= num_options) {
                throw std::runtime_error("categories for variable regions are out of range");
            }
            reverse_categories.push_back(nchunks - cat - 1);
            var_lengths[cat] += regions[f].second - regions[f].first;
        }

        if (forward) {
            forward_variable.resize(num_options);
            forward_trie.resize(num_options);
        }
        if (reverse) {
            reverse_variable.resize(num_options);
            reverse_trie.resize(num_options);
        }

        for (size_t o = 0; o < num_options; ++o) {
            const auto& curopts = opt[o];
            auto len = var_lengths[o];

            if (forward) {
                auto& fseq = forward_variable[v];
                auto& ftrie = forward_trie[v];
                ftrie = MismatchTrie(len);

                for (size_t i = 0; i < curopts.size(); ++i) {
                    auto ptr = curopt[i];
                    std::string current(ptr, ptr + len);
                    if (forward_variable.find(current) != forward_variable.end()) {
                        throw std::runtime_error("already present");
                    }
                    forward_variable[current] = i;
                    ftrie.add(current.c_str(), len);
                }
            }

            if (reverse) {
                auto& rseq = reverse_variable[v];
                auto& rtrie = reverse_trie[v];
                rtrie = MismatchTrie(len);

                for (size_t i = 0; i < curopt.size(); ++i) {
                    auto ptr = curopt[i];
                    std::string current;
                    for (int j = 0; j < len; ++j) {
                        current += reverse_complement(ptr[len - j - 1]);
                    }
                    if (reverse_variable.find(current) != reverse_variable.end()) {
                        throw std::runtime_error("already present");
                    } 
                    reverse_variable[current] = i;
                    rtrie.add(current.c_str(), len);
                }
            }
        }
    }

public:
    struct SearchState {
        std::vector<std::pair<int, int> > results;

        /**
         * @cond
         */
        SearchState(size_t nopt) : results(nopt), seqbuffer(nopt), resbuffer(nopt) {}

        std::vector<std::string> seqbuffer;

        std::vector<std::pair<int, int> > resbuffer; // for use with search_all.
        /**
         * @endcond
         */
    };

    SearchState initialize() {
        return SearchState(num_options);
    }

private:
    int strand_match(
        const char* seq,
        size_t position,
        const std::vector<std::pair<int, int> >& ranges,
        const std::vector<int>& categories,
        const std::vector<std::unordered_map<std::string, int> >& variable,
        const std::vector<MismatchTrie>& trie,
        int obs_mismatches,
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer
    ) const {
        if (obs_mismatches < 0 || obs_mismatches > max_mismatches) {
            return -1;
        }

        // Assembling the variable regions.
        auto start = seq + position;
        for (size_t c = 0; c < categories.size(); ++c) {
            auto idx = categories[c];
            auto& target = seqbuffer[idx];

            const auto& range = ranges[c];
            auto begin = start + range.first;
            auto end = start + range.second;

            if (c == 0) {
                target= std::string(begin, end);
            } else {
                target.insert(target.end(), begin, end);
            }
        }

        for (size_t v = 0; v < varseq.size(); ++v) {
            // Searching for an exact match.
            auto it = variable.find(varseq[v]);
            if (it != variable.end()) {
                results[v].first = it->second;
                results[v].second = 0;
            } else if (obs_mismatches == max_mismatches) {
                // No hope to stay under the max in this case.
                return -1;
            } else {
                auto missed = trie.search(varseq[v].c_str(), max_mismatches - obs_mismatches);
                if (missed.first < 0) {
                    return -1;
                }
                results[v] = missed;
                obs_mismatches += missed.second;
            }
        }

        return obs_mismatches;
    }

    int forward_match(
        const char* seq,
        const ConstantTemplate<N>::MatchDetails& details
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer
    ) const {
        return strand_match(
            seq,
            details.position,
            constant.variable_regions(),
            forward_categories,
            forward_variable,
            forward_trie,
            details.forward_mismatches,
            max_mismatches,
            results,
            seqbuffer
        );
    }

    int reverse_match(
        const char* seq,
        const ConstantTemplate<N>::MatchDetails& details
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer
    ) const {
        return strand_match(
            seq,
            details.position,
            constant.variable_regions(true),
            reverse_categories,
            reverse_variable,
            reverse_trie,
            details.reverse_mismatches,
            max_mismatches,
            results,
            seqbuffer
        );
    }

public:
    bool search_first(const char* seq, size_t len, int max_mismatches, SearchState& state) {
        auto deets = constant.initialize(seq, len);

        while (!deets.finished) {
            if (forward) {
                auto out = forward_match(seq, deets, max_mismatches, state.results, state.seqbuffer);
                if (out >= 0) {
                    return true;
                }
            }

            if (reverse) {
                auto out = reverse_match(seq, deets, max_mismatches, state.results, state.seqbuffer);
                if (out >= 0) {
                    return true;
                }
            }

            constant.next(deets);
        }

        return false;
    }

    bool search_best(const char* seq, size_t len, int max_mismatches, bool fail_ambiguous, SearchState& state) {
        auto deets = constant.initialize(seq, len);
        int best = max_mismatches + 1;
        bool found = false;

        while (!deets.finished) {
            if (forward) {
                auto out = forward_match(seq, deets, max_mismatches, state.resbuffer, state.seqbuffer);
                if (out >= 0) {
                    if (out == best && fail_ambiguous) {
                        found = false;
                    } else if (out < best) {
                        best = out;
                        state.results.swap(state.resbuffer);
                        found = true;
                    } 
                }
            }

            if (reverse) {
                auto out = reverse_match(seq, deets, max_mismatches, state.resbuffer, state.seqbuffer);
                if (out >= 0) {
                    if (out == best && fail_ambiguous) {
                        found = false;
                    } else if (out < best) {
                        best = out;
                        state.results.swap(state.resbuffer);
                        found = true;
                    }
                }
            }

            constant.next(deets);
        }

        return found;
    }

private:
    size_t num_options;
    bool forward, reverse;

    ConstantTemplate<N> constant;
    std::vector<int> forward_categories, reverse_categories;
    std::vector<std::unordered_map<std::string, int> > forward_variable, reverse_variable;
    std::vector<MismatchTrie> forward_trie, reverse_trie;
};

}

#endif
