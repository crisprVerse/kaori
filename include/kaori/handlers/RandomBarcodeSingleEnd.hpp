#ifndef KAORI_RANDOM_BARCODE_SINGLE_END_HPP
#define KAORI_RANDOM_BARCODE_SINGLE_END_HPP

#include "../ScanTemplate.hpp"
#include <vector>
#include <string>
#include <unordered_map>

/**
 * @file RandomBarcodeSingleEnd.hpp
 *
 * @brief Process single-end random barcodes.
 */

namespace kaori {

/**
 * @brief Handler for single-end random barcodes.
 *
 * In this design, the target sequence is created from a template with a single variable region containing a random barcode sequence.
 * The construct containing the target sequence is then subjected to single-end sequencing.
 * This handler will search the read for the target sequence and count the frequency of each random barcode.
 * Random barcodes containing N's are allowed and will be counted separately.
 *
 * @tparam max_size Maximum length of the template sequences on both reads.
 */
template<size_t max_size>
class RandomBarcodeSingleEnd {
public:
    /**
     * @param[in] template_seq Template sequence for the first barcode.
     * This should contain exactly one variable region.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size`.
     * @param strand Strand to use for searching the read sequence - forward (0), reverse (1) or both (2).
     * @param max_mismatches Maximum number of mismatches allowed across the target sequence.
     */
    RandomBarcodeSingleEnd(const char* template_seq, size_t template_length, int strand, int max_mismatches = 0) : 
        forward(strand != 1), 
        reverse(strand != 0),
        constant(template_seq, template_length, forward, reverse),
        max_mm(max_mismatches) {}

    /**
     * @param t Whether to search only for the first match.
     * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
     *
     * @return A reference to this `RandomBarcodeSingleEnd` instance.
     */
    RandomBarcodeSingleEnd& set_first(bool t = true) {
        use_first = t;
        return *this;
    }

private:
    bool use_first = true;

    std::unordered_map<std::string, int> counts;
    int total = 0;

    bool forward, reverse;
    ScanTemplate<max_size> constant;
    int max_mm;

    bool has_match(int obs_mismatches) const {
        return (obs_mismatches >= 0 && obs_mismatches <= max_mm);
    }

public:
    /**
     * @cond
     */
    struct State {
        State() {}
        State(size_t varsize) : buffer(varsize, ' ') {}

        std::unordered_map<std::string, int> counts;
        std::string buffer;
        int total = 0;
    };

    void forward_match(const char* seq, size_t position, State& state) const {
        auto start = seq + position;
        const auto& range = constant.variable_regions()[0];
        std::copy(start + range.first, start + range.second, state.buffer.data());

        auto it = state.counts.find(state.buffer);
        if (it != state.counts.end()) {
            ++(it->second);
        } else {
            state.counts[state.buffer] = 1;
        }
    }

    void reverse_match(const char* seq, size_t position, State& state) const {
        const auto& range = constant.variable_regions()[0];
        auto start = seq + position + range.first;
        size_t len = state.buffer.size();
        for (size_t j = 0; j < len; ++j) {
            state.buffer[j] = reverse_complement<true>(start[len - j - 1]);
        }

        auto it = state.counts.find(state.buffer);
        if (it != state.counts.end()) {
            ++(it->second);
        } else {
            state.counts[state.buffer] = 1;
        }
    }

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        auto read_seq = x.first;
        auto deets = constant.initialize(read_seq, x.second - x.first);

        if (use_first) {
            while (!deets.finished) {
                constant.next(deets);
                if (forward && has_match(deets.forward_mismatches)) {
                    forward_match(read_seq, deets.position, state);
                    break;
                }
                if (reverse && has_match(deets.reverse_mismatches)) {
                    reverse_match(read_seq, deets.position, state);
                    break;
                }
            }

        } else {
            int best = max_mm + 1;
            bool best_forward = true;
            size_t best_position = 0;
            bool best_tied = false;

            while (!deets.finished) {
                constant.next(deets);

                if (forward && has_match(deets.forward_mismatches)) {
                    if (deets.forward_mismatches < best) {
                        best = deets.forward_mismatches;
                        best_position = deets.position;
                        best_forward = true;
                        best_tied = false;
                    } else if (deets.forward_mismatches == best) {
                        best_tied = true;
                    }
                }

                if (reverse && has_match(deets.reverse_mismatches)) {
                    if (deets.reverse_mismatches < best) {
                        best = deets.reverse_mismatches;
                        best_position = deets.position;
                        best_forward = false;
                        best_tied = false;
                    } else if (deets.reverse_mismatches == best) {
                        best_tied = true;
                    }
                }
            }

            if (!best_tied && best <= max_mm) { 
                if (best_forward) {
                    forward_match(read_seq, best_position, state);
                } else {
                    reverse_match(read_seq, best_position, state);
                }
            }
        }

        ++state.total;
    }

    static constexpr bool use_names = false;
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    State initialize() const {
        const auto& range = constant.variable_regions()[0];
        return State(range.second - range.first);
    }

    void reduce(State& s) {
        for (const auto& pair : s.counts) {
            auto it = counts.find(pair.first);
            if (it != counts.end()) {
                it->second += pair.second;
            } else {
                counts[pair.first] = pair.second;
            }
        }
        total += s.total;
    }
    /**
     * @endcond
     */

public:
    /**
     * @return Unordered map containing the frequency of each random barcode.
     */
    const std::unordered_map<std::string, int>& get_counts() const {
        return counts;        
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    int get_total() const {
        return total;
    }
};

}

#endif
