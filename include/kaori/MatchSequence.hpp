#ifndef KAORI_MATCH_SEQUENCE_HPP
#define KAORI_MATCH_SEQUENCE_HPP

#include "ConstantTemplate.hpp"
#include "MismatchTrie.hpp"
#include "utils.hpp"

#include <string>
#include <unordered_map>
#include <vector>

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
            forward_cache.resize(num_options);
        }
        if (reverse) {
            reverse_variable.resize(num_options);
            reverse_trie.resize(num_options);
            reverse_cache.resize(num_options);
        }

        for (size_t o = 0; o < num_options; ++o) {
            const auto& curopts = opt[o];
            auto len = var_lengths[o];

            if (forward) {
                auto& fseq = forward_variable[o];
                auto& ftrie = forward_trie[o];
                ftrie = MismatchTrie(len);

                for (size_t i = 0; i < curopts.size(); ++i) {
                    auto ptr = curopts[i];
                    std::string current(ptr, ptr + len);
                    if (fseq.find(current) != fseq.end()) {
                        throw std::runtime_error("already present");
                    }
                    fseq[current] = i;
                    ftrie.add(current.c_str());
                }
            }

            if (reverse) {
                auto& rseq = reverse_variable[o];
                auto& rtrie = reverse_trie[o];
                rtrie = MismatchTrie(len);

                for (size_t i = 0; i < curopts.size(); ++i) {
                    auto ptr = curopts[i];
                    std::string current;
                    for (int j = 0; j < len; ++j) {
                        current += reverse_complement(ptr[len - j - 1]);
                    }
                    if (rseq.find(current) != rseq.end()) {
                        throw std::runtime_error("already present");
                    } 
                    rseq[current] = i;
                    rtrie.add(current.c_str());
                }
            }
        }
    }

    MatchSequence(const char* s, size_t n, bool f, bool r, const std::vector<const char*>& opt) : 
        MatchSequence(s, n, f, r, std::vector<int>(1), std::vector<std::vector<const char*> >{ opt }) {}

public:
    struct SearchState {
        size_t position = 0;
        int mismatches = 0;
        bool reverse = false;
        std::vector<std::pair<int, int> > identity;

        /**
         * @cond
         */
        SearchState() {}

        SearchState(size_t nopt) : identity(nopt), seqbuffer(nopt), resbuffer(nopt), forward_cache(nopt), reverse_cache(nopt) {}

        std::vector<std::string> seqbuffer;

        std::vector<std::pair<int, int> > resbuffer; // for use with search_all.

        std::vector<std::unordered_map<std::string, std::pair<int, int> > > forward_cache, reverse_cache;
        /**
         * @endcond
         */
    };

    SearchState initialize() const {
        return SearchState(num_options);
    }

    void reduce(SearchState& state) {
        if (forward) {
            for (size_t v = 0; v < forward_cache.size(); ++v) {
                forward_cache[v].merge(state.forward_cache[v]);
                state.forward_cache[v].clear();
            }
        }
        if (reverse) {
            for (size_t v = 0; v < reverse_cache.size(); ++v) {
                reverse_cache[v].merge(state.reverse_cache[v]);
                state.reverse_cache[v].clear();
            }
        }
    }

private:
    int strand_match(
        const char* seq,
        size_t position,
        const std::vector<std::pair<int, int> >& ranges,
        const std::vector<int>& categories,
        const std::vector<std::unordered_map<std::string, int> >& variable,
        const std::vector<MismatchTrie>& trie,
        const std::vector<std::unordered_map<std::string, std::pair<int, int> > >& cache,
        int obs_mismatches,
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer,
        std::vector<std::unordered_map<std::string, std::pair<int, int> > >& local_cache
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

        for (size_t v = 0; v < seqbuffer.size(); ++v) {
            const auto& curseq = seqbuffer[v];
            const auto& curvar = variable[v];

            // Searching for an exact match.
            auto it = curvar.find(curseq);
            if (it != curvar.end()) {
                results[v].first = it->second;
                results[v].second = 0;
            } else if (obs_mismatches == max_mismatches) {
                // No hope to stay under the max in this case.
                return -1;
            } else {
                std::pair<int, int> missed;

                // Seeing if it's any of the caches; otherwise searching the trie.
                const auto& curcache = cache[v];
                auto cit = curcache.find(curseq);
                if (cit == curcache.end()) {
                    auto& curlocal = local_cache[v];
                    auto lit = curlocal.find(curseq);
                    if (lit != curlocal.end()) {
                        missed = lit->second;
                    } else {
                        const auto& curtrie = trie[v];
                        missed = curtrie.search(curseq.c_str(), max_mismatches - obs_mismatches);
                        curlocal[curseq] = missed;
                    }
                } else {
                    missed = cit->second;
                }

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
        const typename ConstantTemplate<N>::MatchDetails& details,
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer,
        std::vector<std::unordered_map<std::string, std::pair<int, int> > >& cache
    ) const {
        return strand_match(
            seq,
            details.position,
            constant.variable_regions(),
            forward_categories,
            forward_variable,
            forward_trie,
            forward_cache,
            details.forward_mismatches,
            max_mismatches,
            results,
            seqbuffer,
            cache
        );
    }

    int reverse_match(
        const char* seq,
        const typename ConstantTemplate<N>::MatchDetails& details,
        int max_mismatches,
        std::vector<std::pair<int, int> >& results,
        std::vector<std::string>& seqbuffer,
        std::vector<std::unordered_map<std::string, std::pair<int, int> > >& cache
    ) const {
        return strand_match(
            seq,
            details.position,
            constant.variable_regions(true),
            reverse_categories,
            reverse_variable,
            reverse_trie,
            reverse_cache,
            details.reverse_mismatches,
            max_mismatches,
            results,
            seqbuffer,
            cache
        );
    }

public:
    bool search_first(const char* seq, size_t len, int max_mismatches, SearchState& state) const {
        auto deets = constant.initialize(seq, len);
        bool found = false;

        auto update = [&](bool rev, int mismatches) -> void {
            found = true;
            state.position = deets.position;
            state.reverse = rev;
            state.mismatches = mismatches;
        };

        while (!deets.finished) {
            constant.next(deets);

            if (forward) {
                auto out = forward_match(seq, deets, max_mismatches, state.identity, state.seqbuffer, state.forward_cache);
                if (out >= 0) {
                    update(false, out);
                    break;
                }
            }

            if (reverse) {
                auto out = reverse_match(seq, deets, max_mismatches, state.identity, state.seqbuffer, state.reverse_cache);
                if (out >= 0) {
                    update(true, out);
                    break;
                }
            }
        }

        return found;
    }

    bool search_best(const char* seq, size_t len, int max_mismatches, SearchState& state) const {
        auto deets = constant.initialize(seq, len);
        int& best = state.mismatches;
        best = max_mismatches + 1;
        bool found = false;

        auto update = [&](bool rev, int mismatches) -> void {
            if (mismatches >= 0) {
                if (mismatches == best) {
                    found = false;
                } else if (mismatches < best) {
                    best = mismatches;
                    max_mismatches = best; // reducing it to truncate the search space for subsequent rounds.
                    state.identity.swap(state.resbuffer);
                    state.position = deets.position;
                    state.reverse = rev;
                    found = true;
                }
            }
        };

        while (!deets.finished) {
            constant.next(deets);

            if (forward) {
                auto out = forward_match(seq, deets, max_mismatches, state.resbuffer, state.seqbuffer, state.forward_cache);
                update(false, out);
            }

            if (reverse) {
                auto out = reverse_match(seq, deets, max_mismatches, state.resbuffer, state.seqbuffer, state.reverse_cache);
                update(true, out);
            }
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
    std::vector<std::unordered_map<std::string, std::pair<int, int> > > forward_cache, reverse_cache;
};

}

#endif
