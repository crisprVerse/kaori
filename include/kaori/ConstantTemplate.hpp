#ifndef KAORI_CONSTANT_TEMPLATE_HPP
#define KAORI_CONSTANT_TEMPLATE_HPP

#include <bitset>
#include <deque>

namespace kaori {

template<class N>
class ConstantTemplate { 
public:
    ConstantTemplate(const char* s, size_t n, bool f, bool r) : length(n), forward(f), reverse(r) {
        if (n * 4 > N) {
            throw std::runtime_error("maximum constant size should be " + std::to_string(N/4) + " nt");
        }

        for (size_t i = 0; i < n; ++i) {
            char b = seq[i];
            if (b != "-") {
                add_base(forward_ref, b);
                forward_mask[i] = 1;
            } else {
                shift(forward_ref);
                if (forward_variables.empty()) {
                    forward_variables.emplace_back(i, i + 1);
                } else {
                    forward_variables.back().second = i + 1;
                }
            }
        }

        for (size_t i = 0; i < n; ++i) {
            char b = seq[n - i - 1];
            if (b != "-") {
                add_base(reverse_ref, reverse_complement(b));
                reverse_mask[i] = 1;
            } else {
                shift(reverse_ref);
                if (reverse_variables.empty()) {
                    reverse_variables.emplace_back(i, i + 1);
                } else {
                    reverse_variables.back().second = i + 1;
                }
            }
        }
    }

public:
    struct MatchDetails {
        std::bitset<N> state;
        const char * seq;
        size_t len;

        size_t position = 0;
        int forward_mismatches = -1;
        int reverse_mismatches = -1;
        bool finished = false;

        std::deque<size_t> bad_words;
    };

    MatchDetails initialize(const char* seq, size_t len) const {
        MatchDetails out;
        out.seq = seq;
        out.len = len;

        if (length <= len) {
            size_t limit = std::min(length, len);
            for (size_t i = 0; i < limit; ++i) {
                char base = seq[i];
                if (is_good(base)) {
                    add_base(out.state, base);
                } else {
                    shift(out.state);
                    out.bad.push_back(i);
                }
            }

            if (out.bad.empty()) {
                full_match(out);
            }
        } else {
            out.finished = true;
        }

        return out;
    }

    void next(MatchDetails& match) const {
        size_t right = match.position + length;
        while (right < match.len) {
            if (match.bad.empty() && match.bad.front() == match.position) {
                match.bad.pop_front();
            }

            char base = seq[right];
            if (is_good(base)) {
                add_base(match.state, base); // no need to trim off the end, the mask will handle that.
            } else {
                shift(match.state);
                match.bad.push_back(right);
            }

            ++match.position;
            ++right;
            if (match.bad.empty()) {
                full_match(match);
                return;
            }
        }

        match.finished = true;
        return;
    }

private:
    std::bitset<N> forward_ref, forward_mask;
    std::bitset<N> reverse_ref, reverse_mask;
    size_t length;
    int mismatches;
    bool forward, reverse;

    std::vector<std::pair<int, int> > forward_variables, reverse_variables;

    static int quick_match(const std::bitset<N>& current, const std::bitset<N>& ref, const std::bitset<N>& mask) {
        return ((current & mask) ^ ref).count();
    }

    void full_match(MatchDetails& match) {
        if (forward) {
            match.forward_mismatches = quick_match(match.state, forward_ref, forward_mask);
        }
        if (reverse) {
            match.reverse_mismatches = quick_match(match.state, reverse_ref, reverse_mask);
        }
    }
};

}

#endif
