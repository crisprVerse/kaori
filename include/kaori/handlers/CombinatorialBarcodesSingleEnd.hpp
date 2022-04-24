#ifndef KAORI_COMBINATORIAL_BARCODES_SINGLE_END_HPP
#define KAORI_COMBINATORIAL_BARCODES_SINGLE_END_HPP

#include "../ConstantTemplate.hpp"
#include "../VariableLibrary.hpp"
#include "../utils.hpp"

#include <array>
#include <unordered_map>

namespace kaori {

template<size_t N, size_t V>
class CombinatorialBarcodesSingleEnd {
public:
    CombinatorialBarcodesSingleEnd(const char* constant, size_t size, int strand, const std::array<SequenceSet, V>& variable, int mismatches = 0) : 
        forward(strand != 1),
        reverse(strand != 0),
        max_mismatches(mismatches),
        constant_matcher(constant, size, forward, reverse)
    {
        const auto& regions = constant_matcher.variable_regions();
        if (regions.size() != V) { 
            throw std::runtime_error("expected " + std::to_string(V) + " variable regions in the constant template");
        }
        for (size_t i = 0; i < V; ++i) {
            size_t rlen = regions[i].second - regions[i].first;
            size_t vlen = variable[i].length;
            if (vlen != rlen) {
                throw std::runtime_error("length of variable region " + std::to_string(i + 1) + " (" + std::to_string(rlen) + 
                    ") should be the same as its sequences (" + std::to_string(vlen) + ")");
            }
        }

        // We'll be using this later.
        for (size_t i = 0; i < V; ++i) {
            num_options[i] = variable[i].choices.size();
        }

        if (forward) {
            for (size_t i = 0; i < V; ++i) {
                const auto& current = regions[i];
                size_t len = current.second - current.first;
                forward_lib[i] = SimpleVariableLibrary(variable[i], mismatches);
            }
        }

        if (reverse) {
            const auto& rev_regions = constant_matcher.template variable_regions<true>();
            for (size_t i = 0; i < V; ++i) {
                const auto& current = rev_regions[i];
                size_t len = current.second - current.first;
                reverse_lib[i] = SimpleVariableLibrary(variable[V - i - 1], mismatches, true);
            }
        }
    }
        
    CombinatorialBarcodesSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

public:
    struct State {
        std::vector<std::array<int, V> >collected;
        int total = 0;

        /**
         * @cond
         */
        std::array<int, V> temp;

        // Default constructors should be called in this case, so it should be fine.
        std::array<typename SimpleVariableLibrary::SearchState, V> forward_details, reverse_details;
        /**
         * @endcond
         */
    };

private:
    template<bool reverse>
    std::pair<bool, int> find_match(
        const char* seq, 
        size_t position, 
        int obs_mismatches, 
        const std::array<SimpleVariableLibrary, V>& libs, 
        std::array<typename SimpleVariableLibrary::SearchState, V>& states, 
        std::array<int, V>& temp) 
    const {
        const auto& regions = constant_matcher.template variable_regions<reverse>();

        for (size_t r = 0; r < V; ++r) {
            auto range = regions[r];
            auto start = seq + position;
            std::string current(start + range.first, start + range.second);

            auto& curstate = states[r];
            libs[r].match(current, curstate, max_mismatches - obs_mismatches);
            if (curstate.index < 0) {
                return std::make_pair(false, 0);
            }
            
            obs_mismatches += curstate.mismatches;
            if (obs_mismatches > max_mismatches) {
                return std::make_pair(false, 0);
            }

            if constexpr(reverse) {
                temp[V - r - 1] = curstate.index;
            } else {
                temp[r] = curstate.index;
            }
        }

        return std::make_pair(true, obs_mismatches);
    }

    std::pair<bool, int> forward_match(const char* seq, const typename ConstantTemplate<N>::MatchDetails& deets, State& state) const {
        return find_match<false>(seq, deets.position, deets.forward_mismatches, forward_lib, state.forward_details, state.temp);
    }

    std::pair<bool, int> reverse_match(const char* seq, const typename ConstantTemplate<N>::MatchDetails& deets, State& state) const {
        return find_match<true>(seq, deets.position, deets.reverse_mismatches, reverse_lib, state.reverse_details, state.temp);
    }

private:
    void process_first(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);

        while (!deets.finished) {
            constant_matcher.next(deets);

            if (forward && deets.forward_mismatches <= max_mismatches) {
                if (forward_match(x.first, deets, state).first) {
                    state.collected.push_back(state.temp);
                    return;
                }
            }

            if (reverse && deets.reverse_mismatches <= max_mismatches) {
                if (reverse_match(x.first, deets, state).first) {
                    state.collected.push_back(state.temp);
                    return;
                }
            }
        }
    }

    void process_best(State& state, const std::pair<const char*, const char*>& x) const {
        auto deets = constant_matcher.initialize(x.first, x.second - x.first);
        bool found = false;
        int best_mismatches = max_mismatches + 1;
        std::array<int, V> best_id;

        auto update = [&](std::pair<int, int> match) -> void {
            if (match.first && match.second <= best_mismatches) {
                if (match.second == best_mismatches) {
                    if (best_id != state.temp) { // ambiguous.
                        found = false;
                    }
                } else { 
                    // A further optimization at this point would be to narrow
                    // max_mismatches to the current 'best_mismatches'. But
                    // this probably isn't worth it.

                    found = true;
                    best_mismatches = match.second;
                    best_id = state.temp;
                }
            }
        };

        while (!deets.finished) {
            constant_matcher.next(deets);

            if (forward && deets.forward_mismatches <= max_mismatches) {
                update(forward_match(x.first, deets, state));
            }

            if (reverse && deets.reverse_mismatches <= max_mismatches) {
                update(reverse_match(x.first, deets, state));
            }
        }

        if (found) {
            state.collected.push_back(best_id);
        }
    }

public:
    State initialize() const {
        return State();
    }

    void reduce(State& s) {
        if (forward) {
            for (size_t r = 0; r < V; ++r) {
                forward_lib[r].reduce(s.forward_details[r]);
            }
        }
        if (reverse) {
            for (size_t r = 0; r < V; ++r) {
                reverse_lib[r].reduce(s.reverse_details[r]);
            }
        }

        combinations.insert(combinations.end(), s.collected.begin(), s.collected.end());
        total += s.total;
        return;
    }

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        if (use_first) {
            process_first(state, x);
        } else {
            process_best(state, x);
        }
        ++state.total;
    }

    static constexpr bool use_names = false;

public:
    void sort() {
        sort_combinations(combinations, num_options);
    }

    const std::vector<std::array<int, V> >& get_combinations() const {
        return combinations;
    }

    int get_total() const {
        return total;
    }
private:
    bool forward;
    bool reverse;
    bool use_first = true;
    int max_mismatches;
    size_t nregions;

    ConstantTemplate<N> constant_matcher;
    std::array<SimpleVariableLibrary, V> forward_lib, reverse_lib;
    std::array<size_t, V> num_options;

    std::vector<std::array<int, V> > combinations;
    int total = 0;
};

}

#endif