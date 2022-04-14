#include <gtest/gtest.h>
#include "kaori/MismatchTrie.hpp"
#include <string>

TEST(MismatchTrie, Basic) {
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };

    std::vector<const char*> ptrs;
    for (const auto& t : things) {
        ptrs.push_back(t.c_str());
    }

    kaori::MismatchTrie stuff(ptrs, 4);
    {
        auto res = stuff.search_single("ACGT", 0);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search_single("AAAT", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search_single("CCAG", 2);
        EXPECT_EQ(res.first, 2);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search_single("AGTT", 0);
        EXPECT_EQ(res.first, 3);
        EXPECT_EQ(res.second, 0);
    }
}
