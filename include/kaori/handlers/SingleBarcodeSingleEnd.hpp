#ifndef KAORI_SINGLE_BARCODE_SINGLE_END_HPP
#define KAORI_SINGLE_BARCODE_SINGLE_END_HPP

#include "byteme/Reader.hpp"
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

        State(const MatchSequence<N>& m, bool f, int mm, size_t nvar) : 
            matcher(&m),
            state(m.initialize()),
            use_first(f),
            mismatches(mm),
            counts(nvar)
        {}

        void process(const std::pair<const char*, const char*>& x) {
            bool found = false;
            if (use_first) {
                found = matcher->search_first(x.first, x.second - x.first, mismatches, state);
            } else {
                found = matcher->search_best(x.first, x.second - x.first, mismatches, state);
            }
            if (found) {
                ++(counts[state.identity[0].first]);
            }
        }

        const MatchSequence<N>* matcher = NULL;
        typename MatchSequence<N>::SearchState state;
        bool use_first;
        int mismatches;

        std::vector<int> counts;
    };

    static constexpr bool use_names = false;

public:
    State initialize() const {
        return State(matcher, use_first, mismatches, counts.size());
    }

    void reduce(State& s) {
        matcher.reduce(s.state);
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
