#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

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

#endif
