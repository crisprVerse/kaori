#ifndef KAORI_MATCH_SEQUENCE_HPP
#define KAORI_MATCH_SEQUENCE_HPP

#include "ConstantTemplate.hpp"
#include "MismatchTrie.hpp"
#include "utils.hpp"

namespace kaori {

template<size_t N>
class MatchSequence {
public:
    MatchSequence(const char* s, size_t n, bool f, bool r, std::vector<int> cat, const std::vector<std::vector<const char*> >& var) : 
        forward(f), 
        reverse(r),
        constant(s, n, f, r), 
        forward_categories(std::move(cat))
    {
        const auto& fvar = constant.forward_variable_regions();
        size_t nvariable = fvar.size();
        forward_variable.resize(nvariable);
        forward_trie.resize(nvariable);
        reverse_variable.resize(nvariable);
        reverse_trie.resize(nvariable);

        // Check that all entries are within range. Note that the
        // reverse search needs to flip the categories.
        for (auto f : forward_categories) {
            if (f < 0 || static_cast<size_t>(f) >= nvariable) {
                throw std::runtime_error("categories for variable regions are out of range");
            }
            reverse_categories.push_back(nvariable - f - 1);
        }

        if (nvariable != var.size()) {
            throw std::runtime_error("number of vectors for variable sequences should be equal to number of variable regions"); 
        }

        for (size_t v = 0; v < nvariable; ++v) {
            const auto& curvar = var[v];
            const auto& curregion = fvar[v];
            size_t len = curregion.second - curregion.first;

            if (forward) {
                auto& fseq = forward_variable[v];
                auto& ftrie = forward_trie[v];
                ftrie = MismatchTrie(len);

                for (size_t i = 0; i < v.size(); ++i) {
                    auto ptr = v[i];
                    std::string current(ptr, ptr + len);
                    if (forward_variable.find(current) != forward_variable.end()) {
                        throw std::runtime_error("already present");
                    }
                    forward_variable[current] = i;
                    ftrie.add(current.c_str(), len);
                }
            }

            if (reverse) {
                auto& rseq = reverse_variable[v];
                auto& rtrie = reverse_trie[v];
                rtrie = MismatchTrie(len);

                for (size_t i = 0; i < v.size(); ++i) {
                    auto ptr = v[i];
                    std::string current;
                    for (int j = 0; j < len; ++j) {
                        current += reverse_complement(ptr[len - j - 1]);
                    }
                    if (reverse_variable.find(current) != reverse_variable.end()) {
                        throw std::runtime_error("already present");
                    } 
                    reverse_variable[current] = i;
                    rtrie.add(current.c_str(), len);
                }
            }
        }
    }

public:

private:
    bool forward, reverse;
    ConstantTemplate<N> constant;
    std::vector<int> forward_categories, reverse_categories;
    std::vector<std::unordered_map<std::string, int> > forward_variable, reverse_variable;
    std::vector<MismatchTrie> forward_trie, reverse_trie;
};

}

#endif
