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
 * In this design, the barcoding element is created from a template with a single variable region containing a random barcode sequence.
 * The construct containing the barcoding element is then subjected to single-end sequencing.
 * This handler will search the read for the barcoding element and count the frequency of each random barcode.
 * Random barcodes containing N's are allowed and will be counted separately.
 *
 * @tparam max_size_ Maximum length of the template sequences on both reads.
 */
template<size_t max_size_>
class RandomBarcodeSingleEnd {
public:
    /**
     * @brief Optional parameters for `SingleBarcodeSingleEnd`.
     */
    struct Options {
        /** 
         * Maximum number of mismatches allowed across the barcoding element.
         */
        int max_mismatches = 0;

        /** 
         * Whether to search only for the first match.
         * If `false`, the handler will search for the best match (i.e., fewest mismatches) instead.
         */
        bool use_first = true;

        /** 
         * Strand(s) of the read sequence to search.
         */
        SearchStrand strand = SearchStrand::FORWARD;
    };

public:
    /**
     * @param[in] template_seq Template sequence of the barcoding element.
     * This should contain exactly one variable region.
     * @param template_length Length of the template.
     * This should be less than or equal to `max_size_`.
     * @param options Optional parameters.
     */
    RandomBarcodeSingleEnd(const char* template_seq, size_t template_length, const Options& options) :
        my_forward(search_forward(options.strand)),
        my_reverse(search_reverse(options.strand)),
        my_constant(template_seq, template_length, options.strand),
        my_max_mm(options.max_mismatches),
        my_use_first(options.use_first)
    {}

private:
    std::unordered_map<std::string, int> my_counts;
    int my_total = 0;

    bool my_forward, my_reverse;
    ScanTemplate<max_size_> my_constant;
    int my_max_mm;
    bool my_use_first;

    bool has_match(int obs_mismatches) const {
        return (obs_mismatches >= 0 && obs_mismatches <= my_max_mm);
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
        const auto& range = my_constant.forward_variable_regions()[0];
        std::copy(start + range.first, start + range.second, state.buffer.data());

        auto it = state.counts.find(state.buffer);
        if (it != state.counts.end()) {
            ++(it->second);
        } else {
            state.counts[state.buffer] = 1;
        }
    }

    void reverse_match(const char* seq, size_t position, State& state) const {
        const auto& range = my_constant.forward_variable_regions()[0];
        auto start = seq + position + range.first;
        size_t len = state.buffer.size();
        for (size_t j = 0; j < len; ++j) {
            state.buffer[j] = complement_base<true>(start[len - j - 1]);
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
        auto deets = my_constant.initialize(read_seq, x.second - x.first);

        if (my_use_first) {
            while (!deets.finished) {
                my_constant.next(deets);
                if (my_forward && has_match(deets.forward_mismatches)) {
                    forward_match(read_seq, deets.position, state);
                    break;
                }
                if (my_reverse && has_match(deets.reverse_mismatches)) {
                    reverse_match(read_seq, deets.position, state);
                    break;
                }
            }

        } else {
            int best = my_max_mm + 1;
            bool best_forward = true;
            size_t best_position = 0;
            bool best_tied = false;

            while (!deets.finished) {
                my_constant.next(deets);

                if (my_forward && has_match(deets.forward_mismatches)) {
                    if (deets.forward_mismatches < best) {
                        best = deets.forward_mismatches;
                        best_position = deets.position;
                        best_forward = true;
                        best_tied = false;
                    } else if (deets.forward_mismatches == best) {
                        best_tied = true;
                    }
                }

                if (my_reverse && has_match(deets.reverse_mismatches)) {
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

            if (!best_tied && best <= my_max_mm) { 
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
        const auto& range = my_constant.forward_variable_regions()[0];
        return State(range.second - range.first);
    }

    void reduce(State& s) {
        for (const auto& pair : s.counts) {
            auto it = my_counts.find(pair.first);
            if (it != my_counts.end()) {
                it->second += pair.second;
            } else {
                my_counts[pair.first] = pair.second;
            }
        }
        my_total += s.total;
    }
    /**
     * @endcond
     */

public:
    /**
     * @return Unordered map containing the frequency of each random barcode.
     */
    const std::unordered_map<std::string, int>& get_counts() const {
        return my_counts;        
    }

    /**
     * @return Total number of reads processed by the handler.
     */
    int get_total() const {
        return my_total;
    }
};

}

#endif
