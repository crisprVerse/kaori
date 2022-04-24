#ifndef KAORI_CONSTANT_TEMPLATE_HPP
#define KAORI_CONSTANT_TEMPLATE_HPP

#include <bitset>
#include <deque>
#include "utils.hpp"

/**
 * @file ConstantTemplate.hpp
 *
 * @brief Defines the `ConstantTemplate` class.
 */

namespace kaori {

/**
 * @brief Search a read sequence for a constant template sequence.
 *
 * The template sequence is expected to contain constant regions interspersed with one or more variable regions.
 * This class will find the location on the read sequence that matches the constant regions of the template (give or take any number of substitutions).
 * Once a match is found, the variable regions can be extracted for look-up with, e.g., the `VariableLibrary` class.
 *
 * For efficiency, the search itself is done after converting all base sequences into a bit encoding.
 * The maximum bit size of the encoding window is determined at compile-time using `N`.
 * As four bits are used in the encoding, the maximum size of the constant template is `N / 4`.
 *
 * @tparam N Size of the bitset used to represent the constant template sequence.
 */
template<size_t N>
class ConstantTemplate { 
public:
    /**
     * Default constructor.
     * This is only provided to enable composition, the resulting object should not be used until it is copy-assigned to a properly constructed instance.
     */
    ConstantTemplate() {}

    /**
     * @param[in] s Pointer to a character array containing the constant template.
     * Constant sequences should only contain `A`, `C`, `G` or `T` (or their lower-case equivalents).
     * Variable regions should be marked with `-`.
     * @param n Length of the array pointed to by `s`.
     * @param f Should the search be performed on the forward strand of the read sequence?
     * @param r Should the search be performed on the reverse strand of the read sequence?
     */
    ConstantTemplate(const char* s, size_t n, bool f, bool r) : length(n), forward(f), reverse(r) {
        if (n * 4 > N) {
            throw std::runtime_error("maximum constant size should be " + std::to_string(N/4) + " nt");
        }

        if (forward) {
            for (size_t i = 0; i < n; ++i) {
                char b = s[i];
                if (b != '-') {
                    add_base(forward_ref, b);
                    add_mask(forward_mask, i);
                } else {
                    shift(forward_ref);
                    shift(forward_mask);
                    add_variable_base(forward_variables, i);
                }
            }
        } else {
            // Forward variable regions are always defined.
            for (size_t i = 0; i < n; ++i) {
                char b = s[i];
                if (b == '-') {
                    add_variable_base(forward_variables, i);
                }
            }
        }

        if (reverse) {
            for (size_t i = 0; i < n; ++i) {
                char b = s[n - i - 1];
                if (b != '-') {
                    add_base(reverse_ref, reverse_complement(b));
                    add_mask(reverse_mask, i);
                } else {
                    shift(reverse_ref);
                    shift(reverse_mask);
                    add_variable_base(reverse_variables, i);
                }
            }
        }
    }

public:
    /**
     * @brief Details on the current match to the read sequence.
     */
    struct MatchDetails {
        /**
         * Position of the match.
         */
        size_t position = static_cast<size_t>(-1); // overflow should be sane.

        /**
         * Number of mismatches on the forward strand.
         * If negative, no search was performed.
         */
        int forward_mismatches = -1;

        /**
         * Number of mismatches on the reverse strand.
         * If negative, no search was performed.
         */
        int reverse_mismatches = -1;

        /**
         * Whether the match is at the end of the read sequence.
         * If `true`, `next()` should not be called.
         */
        bool finished = false;

        /**
         * @cond
         */
        std::bitset<N> state, ambiguous;
        const char * seq;
        size_t len;
        std::deque<size_t> bad;
        /**
         * @endcond
         */
    };

    /**
     * Begin a new search for the template in a read sequence.
     *
     * @param[in] seq Pointer to an array containing the read sequence.
     * @param len Length of the read sequence.
     *
     * @return A `MatchDetails` object that can be passed to `next()`.
     */
    MatchDetails initialize(const char* seq, size_t len) const {
        MatchDetails out;
        out.seq = seq;
        out.len = len;

        if (length <= len) {
            for (size_t i = 0; i < length - 1; ++i) {
                char base = seq[i];

                if (is_good(base)) {
                    add_base(out.state, base);
                    if (!out.bad.empty()) {
                        shift(out.ambiguous);
                    }
                } else {
                    shift(out.state);
                    out.state |= other_<N>;
                    shift(out.ambiguous);
                    out.ambiguous |= other_<N>;
                    out.bad.push_back(i);
                }
            }
        } else {
            out.finished = true;
        }

        return out;
    }

    /**
     * Find the next match in the read sequence.
     * The first invocation will search for a match at position 0;
     * this can be repeatedly called until `match.finished` is `true`.
     *
     * @param match A `MatchDetails` object produced by `initialize()`.
     *
     * @return `match` is updated with the details of the current match at a particular position on the read sequence.
     */
    void next(MatchDetails& match) const {
        if (!match.bad.empty() && match.bad.front() == match.position) {
            match.bad.pop_front();
            if (match.bad.empty()) {
                // This should effectively clear the ambiguous bitset, allowing
                // us to skip its shifting if there are no more ambiguous
                // bases. We do it here because we won't get an opportunity to
                // do it later; as 'bad' is empty, the shift below is skipped.
                shift(match.ambiguous); 
            }
        }

        size_t right = match.position + length;
        char base = match.seq[right];
        if (is_good(base)) {
            add_base(match.state, base); // no need to trim off the end, the mask will handle that.
            if (!match.bad.empty()) {
                shift(match.ambiguous);
            }
        } else {
            shift(match.state);
            match.state |= other_<N>;
            shift(match.ambiguous);
            match.ambiguous |= other_<N>;
            match.bad.push_back(right);
        }

        ++match.position;
        full_match(match);
        if (right + 1 == match.len) {
            match.finished = true;
        }

        return;
    }

private:
    std::bitset<N> forward_ref, forward_mask;
    std::bitset<N> reverse_ref, reverse_mask;
    size_t length;
    int mismatches;
    bool forward, reverse;

    static void add_mask(std::bitset<N>& current, size_t pos) {
        shift(current);
        current |= other_<N>;
        for (int i = 0; i < 4; ++i) {
            current[i] = 1;
        }
    }

    static int strand_match(const MatchDetails& match, const std::bitset<N>& ref, const std::bitset<N>& mask) {
        // pop count here is equal to the number of non-ambiguous mismatches *
        // 2 + number of ambiguous mismatches * 3. This is because
        // non-ambiguous bases are encoded by 1 set bit per 4 bases (so 2 are
        // left after a XOR'd mismatch), while ambiguous mismatches are encoded
        // by all set bits per 4 bases (which means that 3 are left after XOR).
        int pcount = ((match.state & mask) ^ ref).count(); 

        // Counting the number of ambiguous bases after masking. Each ambiguous
        // base is represented by 4 set bits, so we divide by 4 to get the number
        // of bases; then we multiply by three to remove their contribution. The
        // difference is then divided by two to get the number of non-ambig mm's.
        if (!match.bad.empty()) {
            int acount = (match.ambiguous & mask).count();
            acount /= 4;
            return acount + (pcount - acount * 3) / 2;
        } else {
            return pcount / 2;
        }
    }

    void full_match(MatchDetails& match) const {
        if (forward) {
            match.forward_mismatches = strand_match(match, forward_ref, forward_mask);
        }
        if (reverse) {
            match.reverse_mismatches = strand_match(match, reverse_ref, reverse_mask);
        }
    }

private:
    std::vector<std::pair<int, int> > forward_variables, reverse_variables;

public:
    /**
     * Get the details about the variable regions in the constant template.
     *
     * @tparam reverse Should we return the coordinates of the variable regions when searching on the reverse strand?
     *
     * @return A vector of pairs where each pair specifies the start and one-past-the-end position of each variable region in the template.
     * Coordinates are reported relative to the start of the template.
     * If `reverse = true`, coordinates are reported after reverse-complementing the template sequence.
     */
    template<bool reverse = false>
    const std::vector<std::pair<int, int> >& variable_regions() const {
        if constexpr(reverse) { 
            return reverse_variables;
        } else {
            return forward_variables;
        }
    } 

private:
    static void add_variable_base(std::vector<std::pair<int, int> >& variables, int i) {
        if (!variables.empty()) {
            auto& last = variables.back().second;
            if (last == i) {
                ++last;
                return;
            }
        }
        variables.emplace_back(i, i + 1);
        return;
    }
};

}

#endif
