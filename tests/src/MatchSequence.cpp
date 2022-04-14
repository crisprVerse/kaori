#include <gtest/gtest.h>
#include "kaori/MatchSequence.hpp"
#include <string>
#include <vector>
#include "utils.h"

TEST(MatchSequence, BasicFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::MatchSequence<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

    // Perfect match.
    {
        std::string seq = "cagcatcgatcgtgaACGTAAAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 0);
    }

    // 1 mismatch in constant region.
    {
        std::string seq = "cagcatcgatcgtgaACGGAAAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 1);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 0);
    }

    // 1 mismatch in variable region.
    {
        std::string seq = "cagcatcgatcgtgaACGTATAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 1);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 1);
    }

    // Blocked for exceeding the number of mismatches.
    {
        std::string seq = "cagcatcgatcgtgaACGTATTATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 1, state));
    }

    // Blocked for ambiguity.
    {
        std::string seq = "cagcatcgatcgtgaACGTATTATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
    }
}
