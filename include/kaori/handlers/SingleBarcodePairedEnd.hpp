#ifndef KAORI_SINGLE_BARCODE_PAIRED_END_HPP
#define KAORI_SINGLE_BARCODE_PAIRED_END_HPP

#include "../MatchSequence.hpp"

namespace kaori {

template<size_t N>
class SingleBarcodePairedEnd {
public:
    SingleBarcodePairedEnd(const char* constant, size_t size, int strand, const std::vector<const char*>& variable) : 
        matcher(constant, size, strand != 1, strand != 0, variable), 
        counts(variable.size()) {}
        
    SingleBarcodePairedEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

    SingleBarcodePairedEnd& set_mismatches(int m = 0) {
        mismatches = m;
        return *this;
    }

public:
    struct State {
        State() {}

        State(typename MatchSequence<N>::SearchState s, size_t nvar) : search(std::move(s)), counts(nvar) {}

        typename MatchSequence<N>::SearchState search;
        std::vector<int> counts;
    };

    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) const {
        if (use_first) {
            if (matcher.search_first(r1.first, r1.second - r1.first, mismatches, state.search)) {
                ++state.counts[state.search.identity[0].first];
            } else if (matcher.search_first(r2.first, r2.second - r2.first, mismatches, state.search)) {
                ++state.counts[state.search.identity[0].first];
            }
        } else {
            bool found1 = matcher.search_best(r1.first, r1.second - r1.first, mismatches, state.search);
            auto res1 = state.search.identity[0];

            bool found2 = matcher.search_best(r2.first, r2.second - r2.first, mismatches, state.search);
            auto res2 = state.search.identity[0];

            if (found1 && !found2) {
                ++state.counts[res1.first];
            } else if (!found1 && found2) {
                ++state.counts[res2.first];
            } else if (found1 && found2) {
                auto diff = res1.second - res2.second;
                if (diff < 0) {
                    ++state.counts[res1.first];
                } else if (diff > 0) {
                    ++state.counts[res2.first];
                } else if (res2.first == res1.first) {
                    ++state.counts[res1.first];
                }
            }
        }
    }

    static constexpr bool use_names = false;

public:
    State initialize() const {
        return State(matcher.initialize(), counts.size());
    }

    void reduce(State& s) {
        matcher.reduce(s.search);
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] += s.counts[i];
        }
    }

private:
    MatchSequence<N> matcher;
    std::vector<int> counts;
    bool use_first = true;
    int mismatches = 0;

public:
    const std::vector<int>& results() const {
        return counts;        
    }
};

}

#endif
