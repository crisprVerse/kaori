#ifndef KAORI_EXACT_LIBRARY_HPP
#define KAORI_EXACT_LIBRARY_HPP

#include <unordered_map>
#include <bitset>
#include <vector>
#include <cmath>
#include <tuple>

namespace kaori {

template<size_t N>
struct SingleMismatchLibrary {
public:
    SingleMismatchLibrary(int seq_len, size_t expected_num_seq, bool forward, bool reverse, bool indels) :
        length(seq_len),
        allow_indels(indels),
        frag_len(seq_len / 2) // implicit floor.
    {
        if (seq_len * 2 > N) {
            throw std::runtime_error("sequence length must be no greater than the hard limit of " + std::to_string(N/2));
        }
        if (forward || reverse) {
            expected_num_seq *= 2; // we'll be creating two halves.
            if (forward && reverse) {
                expected_num_seq *= 2;
            }
            mappings.reserve(expected_num_seq);
        }
    }

    void add(const char* seq) {
        if (forward) {
            add(seq, false);
        }

        if (reverse) {
            // Storing the reverse complement. Note that the interval corresponding to the first half
            // becomes the second half here, because everything is reversed.
            std::vector<char> revcomp;
            revcomp.reserve(length);
            for (int i = 0; i < length; ++i) {
                revcomp.push_back(reverse_complement(seq[length - i - 1]));
            }
            add(revcomp.data(), true);
        }

        ++counter;
        return;
    }

    int hash_length() const {
        return frag_len;
    }

public:
    struct Match {
        Match() {}
        int id = -1;
        int edits = 0;
        bool reverse = false;
    };

    Match query(const std::bitset<N>& hash, const char* before, size_t before_len, const char* after, size_t after_len) const {
        Match output;
        output.edits = 2;

        auto it = mappings.find(hash);
        if (it != mappings.end()) {
            const auto& before_partials = (it->second).before;
            const auto& after_partials = (it->second).after;

            bool ambiguous = false;
            if (!allow_indels) {
                if (before_len && before_partials.size()) {
                    ambiguous = no_indel<true>(output, before, before_len, before_partials);
                }
                if (after_len && after_partials.size() && !ambiguous) {
                    no_indel<false>(output, after, after_len, after_partials);
                }
            } else {
                if (before_len && before_partials.size()) {
                    ambiguous = with_indel<true>(output, before, before_len, before_partials);
                }
                if (after_len && after_partials.size() && !ambiguous) {
                    with_indel<false>(output, after, after_len, after_partials);
                }
            }
        }

        return output;
    }

private:
    template<bool before>
    static char pick_base(const char* s, size_t len, size_t pos) {
        if constexpr(!before) {
            return s[pos];
        } else {
            return s[len - pos - 1];
        }
    }

    template<bool before>
    bool no_indel(Match& match, const char* seq, size_t len, const std::vector<std::pair<std::string, bool> >& partials) const {
        bool ambiguous = false;

        for (const auto& p : partials) {
            const auto& current = p.sequence;

            // Skipping if it's shorter, as that implies a deletion where it is not allowed.
            if (current.size() > len) {
                continue;
            }

            int mismatches = 0;
            size_t pos = 0;
            bool failed = false;

            for (; pos < current.size(); ++pos) {
                mismatches += (pick_base<before>(current.c_str(), current.size(), pos) != pick_base<before>(seq, len, pos));
                if (mismatches > match.edits) {
                    failed = true;
                    break;
                }
            }

            if (!failed) {
                if (mismatches < match.edits) {
                    match.id = p.id;
                    match.reverse = p.reverse;
                    match.edits = mismatches;
                } else if (mismatches == match.edits) {
                    // Instant failure. Because we're only allowing 1 mismatch,
                    // any further hits will be ambiguous at best, so we quit.
                    match.id = -1; 
                    ambiguous = true;
                    break;
                }
            }
        }

        return ambiguous;
    }

    template<bool before>
    bool with_indel(Match& match, const char* seq, size_t len, const std::vector<std::pair<std::string, bool> >& partials) const {
        bool ambiguous = false;
        for (const auto& p : partials) {
            const auto& current = p.sequence;

            // Skipping if it's much shorter, as that implies more deletions than allowed.
            if (current.size() + 1 > len) {
                continue;
            }

            // Doing a Smith-Waterman after restricting to only one indel or substitution.
            int mismatches = 0;
            size_t limit = std::min(current.size(), len); 
            size_t pos = 0;
            bool failed = false;
            int same = 0, ins = 1, del = 1;

            for (; pos < limit; ++pos) {
                auto match = (pick_base<before>(current.c_str(), current.size(), pos) != pick_base<before>(seq, len, pos));

                if (same == 0) {
                    ins = std::min(ins + match, same + 1);
                    same += match; // can't be ins + 1 or del + 1, because that would already be > 1.
                    del = std::min(same + 1, del + match);
                } else if (match) {
                    failed = true;
                    break;
                } else {
                    // 'same' is unchanged, 'ins' and 'del' can't afford to be
                    // 'same + 1' as this would push it over the limit, so they
                    // must be 'ins + match = ins + 0', i.e., unchanged.
                }
            }

            if (!failed) {
                int mismatches;
                if (limit == current.size()) {
                    mismatches = std::min({ ins, del, same });
                } else {
                    mismatches = std::min(ins + match, same + 1);
                }

                if (mismatches < match.edits) {
                    match.id = p.id;
                    match.reverse = p.reverse;
                    match.edits = mismatches;
                } else if (mismatches == match.edits) {
                    // Instant failure. Because we're only allowing 1 mismatch,
                    // any further hits will be ambiguous at best, so we quit.
                    match.id = -1; 
                    ambiguous = true;
                    break;
                }
            }
        }

        return ambiguous;
    }

private:
    struct Fragment {
        Fragment(std::string s, int i, bool r) : sequence(std::move(s)), id(i), reverse(r) {}
        std::string sequence;
        int id;
        bool reverse;
    };

    struct Partial {
        std::vector<Fragment> before, after;
    };

    std::unordered_map<std::bitset<N>, Partial> mappings;
    int counter = 0;

    int length;
    bool allow_indels;
    size_t frag_len;

private:
    void add(const char* seq, bool reverse) {
        for (int frag = 0; frag < 2; ++frag) {
            std::bitset<N> current;
            size_t start = (frag == 0 ? 0 : length - frag_len);
            for (int i = 0; i < frag_len; ++i) {
                add_base<N>(current, seq[start + i]);
            }

            auto& partials = mappings[current];
            auto& current = (frag == 0 ? partials.before : partials.after);
            current.emplace_back(std::string(seq + start, seq + frag_len), counter, reverse);
        }
    }
};

}

#endif
