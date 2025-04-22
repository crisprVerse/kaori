#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <algorithm>

inline std::string convert_to_fastq(const std::vector<std::string>& reads, std::string prefix = "READ") {
    std::string output;

    for (size_t i = 0; i < reads.size(); ++i) {
        output += "@" + prefix + std::to_string(i+1) + "\n";
        output += reads[i] + "\n";
        output += "+\n" + std::string(reads[i].size(), '!') + "\n";
    }

    return output;
}

inline std::pair<const char*, const char*> bounds(const std::string& s) {
    return std::make_pair(s.c_str(), s.c_str() + s.size());
}

template<int num_variable_, typename ResultMap_>
std::vector<std::pair<std::array<kaori::BarcodeIndex, num_variable_>, kaori::Count> flatten_results(const ResultMap_& res) {
    std::vector<std::pair<std::array<kaori::BarcodeIndex, num_variable_>, kaori::Count> output;
    output.reserve(res.size();
    for (const auto& r : res) {
        output.push_back(r);
    }
    std::sort(output.begin(), output.end());
    return output;
}

#endif
