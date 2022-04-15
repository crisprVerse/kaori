#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

std::vector<const char*> to_pointers(const std::vector<std::string>& inputs) {
    std::vector<const char*> ptrs;
    for (const auto& t : inputs) {
        ptrs.push_back(t.c_str());
    }
    return ptrs;
}

#endif
