#ifndef KAORI_SINGLE_BARCODE_SINGLE_END_HPP
#define KAORI_SINGLE_BARCODE_SINGLE_END_HPP

#include "../MatchSequence.hpp"

namespace kaori {

template<size_t N>
class SingleBarcodeSingleEnd {
public:
    SingleBarcodeSingleEnd(const char* constant, size_t size, int strand, const std::vector<const char*>& variable) : 
        matcher(constant, size, strand != 1, strand != 0, variable), counts(variable.size()) {}
        
    SingleBarcodeSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

    SingleBarcodeSingleEnd& set_mismatches(int m = 0) {
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

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        bool found = false;
        if (use_first) {
            found = matcher.search_first(x.first, x.second - x.first, mismatches, state.search);
        } else {
            found = matcher.search_best(x.first, x.second - x.first, mismatches, state.search);
        }
        if (found) {
            ++(state.counts[state.search.identity[0].first]);
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
