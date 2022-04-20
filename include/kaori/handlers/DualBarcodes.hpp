#ifndef KAORI_DUAL_BARCODES_HPP
#define KAORI_DUAL_BARCODES_HPP

namespace kaori {

template<size_t N>
class DualBarcodes { 
public:
    DualBarcodes(
        const char* con1, size_t n1, bool rev1, const std::vector<const char*>& var1,
        const char* con2, size_t n2, bool rev2, const std::vector<const char*>& var2,
        int total_mismatches, bool random = false
    ) :
        reverse1(rev1),
        reverse2(rev2),
        constant1(con1, n1, !reverse1, reverse1),
        constant2(con2, n2, !reverse2, reverse2),
        constant_mismatches1(total_mismatches),
        constant_mismatches2(total_mismatches),
        max_mismatches(total_mismatches),
        randomized(random)
    {
        auto num_options = var1.size();
        if (num_options != var2.size()) {
            throw std::runtime_error("each read should contain the same number of choices for the variable region");
        }
        counts.resize(num_options);

        size_t len1;
        {
            const auto& regions = constant1.variable_regions();
            if (regions.size() != 1) { 
                throw std::runtime_error("expected one variable region in the first constant template");
            }
            len1 = regions[0].second - regions[0].first;
        }

        size_t len2;
        {
            const auto& regions = constant2.variable_regions();
            if (regions.size() != 1) { 
                throw std::runtime_error("expected one variable region in the second constant template");
            }
            len2 = regions[0].second - regions[0].first;
        }

        // Constructing the combined thing.
        std::vector<std::string> combined;
        combined.reserve(num_options);

        for (size_t i = 0; i < num_options; ++i) {
            std::string current;

            auto ptr1 = var1[i];
            if (reverse1) {
                for (int j = 0; j < len1; ++j) {
                    current += reverse_complement(ptr1[len1 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr1, ptr1 + len1);
            }

            auto ptr2 = var2[i];
            if (reverse2) {
                for (int j = 0; j < len2; ++j) {
                    current += reverse_complement(ptr2[len2 - j - 1]);
                }
            } else {
                current.insert(current.end(), ptr2, ptr2 + len2);
            }

            combined.push_back(std::move(current));
        }

        // Constructing the combined varlib.
        std::vector<const char*> ptrs;
        ptrs.reserve(num_options);

        for (size_t i =0 ; i <num_options; ++i) {
            ptrs.push_back(combined[i].c_str());
        }
        varlib = VariableLibrary(ptrs, len1 + len2, max_mismatches, false);
    }

public:
    struct State {
        State(size_t n = 0) : counts(n) {}
        std::vector<int> counts;
        int total = 0;

        /**
         * @cond
         */
        // Default constructors should be called in this case, so it should be fine.
        typename VariableLibrary::SearchState details;
        /**
         * @endcond
         */
    };

    State initialize() const {
        return State(num_options);
    }

    void reduce(State& s) {
        varlib.reduce(s.details);
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i] += s.counts[i];
        }
        total += s.total;
    }

private:
    void inner_process(
        bool reverse, 
        const ConstantTemplate<N>& constant, 
        int constant_mismatches,
        const std::pair<const char*, const char*>& against, 
        std::string& current) 
    {
        bool found = false;
        const deets = constant.initialize(against.first, against.second - against.first);

        while (!deets.finished) {
            constant.next(deets);
            if (reverse) {
                if (deets.reverse_mismatches <= constant_mismatches) {
                    const auto& reg = constant.variable_regions<true>();
                    auto start = against.first + deets.position;
                    search.insert(search.end(), start + reg.first, start + reg.second);
                    found = true;
                    break;
                }
            } else {
                if (deets.forward_mismatches <= constant_mismatches) {
                    const auto& reg = constant.variable_regions();
                    auto start = against.first + deets.position;
                    search.insert(search.end(), start + reg.first, start + reg.second);
                    found = true;
                    break;
                }
            }
        }

        return found;
    }

    bool process_(State& state, const std::pair<const char*, const char*>& against1, const std::pair<const char*, const char*>& against2) {
        std::string search;

        if (!inner_process(reverse1, constant1, constant_mismatches1, against1, search)) {
            return false;
        }
        if (!inner_process(reverse2, constant2, constant_mismatches2, against2, search)) {
            return false;
        }

        varlib.match(state.details, search.c_str());
        if (state.details.index >= 0 && state.details.mismatch <= max_mismatches) {
            ++(counts[state.details.index]);
        }

        // i.e., was any constant match found? We don't care so much about the 
        // successful recovery of the variable region when deciding whether to 
        // quit in the upstream function.
        return true;
    }

public:
    void process(State& state, const std::pair<const char*, const char*>& r1, const std::pair<const char*, const char*>& r2) {
        bool constant_found = process_(state, r1, r2);

        if (!constant_found && randomized) {
            process_(state, r2, r1);
        }

        ++total;
    }

private:
    bool reverse1, reverse2;

    ConstantTemplate<N> constant1, constant2;
    int constant_mismatches1, constant_mismatches2;

    VariableLibrary varlib;
    int max_mismatches;

    bool randomized;

    std::vector<int> counts;
    int total = 0;

public:
    const std::vector<int>& get_counts() const {
        return counts;
    }

    int get_total() const {
        return total;
    }
};

}

#endif
