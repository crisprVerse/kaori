#include <gtest/gtest.h>
#include "kaori/VariableLibrary.hpp"
#include <string>
#include <vector>
#include "utils.h"

TEST(SimpleVariableLibrary, Basic) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);

    kaori::SimpleVariableLibrary stuff(ptrs, 4);
    auto init = stuff.initialize();

    stuff.match("AAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.match("TTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.match("CCAC", init);
    EXPECT_EQ(init.index, -1);
}

TEST(SimpleVariableLibrary, ReverseComplement) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);

    kaori::SimpleVariableLibrary stuff(ptrs, 4, 0, true);
    auto init = stuff.initialize();

    stuff.match("AAAA", init);
    EXPECT_EQ(init.index, 3);
    stuff.match("TTTT", init);
    EXPECT_EQ(init.index, 0);
}

TEST(SimpleVariableLibrary, Mismatches) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);

    {
        kaori::SimpleVariableLibrary stuff(ptrs, 4, 1);
        auto init = stuff.initialize();

        stuff.match("AAAA", init);
        EXPECT_EQ(init.index, 0);
        EXPECT_EQ(init.mismatches, 0);
        
        stuff.match("TTTT", init);
        EXPECT_EQ(init.index, 3);
        EXPECT_EQ(init.mismatches, 0);

        stuff.match("CCAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 1);
        stuff.match("CGAC", init);
        EXPECT_EQ(init.index, -1);
    }

    {
        kaori::SimpleVariableLibrary stuff(ptrs, 4, 2);
        auto init = stuff.initialize();

        stuff.match("CCAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 1);

        stuff.match("CGAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 2);

        stuff.match("CGGC", init); // ambiguous.
        EXPECT_EQ(init.index, -1);
    }
}

TEST(SimpleVariableLibrary, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleVariableLibrary stuff(ptrs, 4, 1);

    auto state = stuff.initialize();

    // No cache when there is no mismatch.
    {
        stuff.match("AAAA", state);
        auto it = state.cache.find("AAAA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // Stored in cache for >1 mismatches.
    {
        stuff.match("AATA", state);
        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).first, 0);
        EXPECT_EQ((it->second).second, 1);
    }

    {
        stuff.match("ACTA", state);
        auto it = state.cache.find("ACTA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).first, -1);
    }
    
    // Checking that the reduction works correctly.
    state.cache["AATA"].first = 2;
    stuff.reduce(state);
    EXPECT_TRUE(state.cache.empty());

    {
        stuff.match("AATA", state);
        EXPECT_EQ(state.index, 2); // re-uses the cache value!
    }
}

TEST(SegmentedVariableLibrary, Basic) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    auto ptrs = to_pointers(variables);

    kaori::SegmentedVariableLibrary<2> stuff(ptrs, { 2, 4 }, { 0, 1 });
    auto init = stuff.initialize();

    stuff.match("AAAAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.match("AATTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.match("AACCAC", init); // 1 mismatch
    EXPECT_EQ(init.index, 1);

    stuff.match("AAccgg", init); // ambiguous.
    EXPECT_EQ(init.index, -1);
}

TEST(SegmentedVariableLibrary, ReverseComplement) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    auto ptrs = to_pointers(variables);

    kaori::SegmentedVariableLibrary<2> stuff(ptrs, { 2, 4 }, { 0, 1 }, true);
    auto init = stuff.initialize();

    stuff.match("AAAATT", init);
    EXPECT_EQ(init.index, 3);
    stuff.match("GGGGTT", init);
    EXPECT_EQ(init.index, 1);

    stuff.match("CCACTT", init);
    EXPECT_EQ(init.index, 2);

    stuff.match("GGCCTT", init); // ambiguous.
    EXPECT_EQ(init.index, -1);
}

TEST(SegmentedVariableLibrary, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SegmentedVariableLibrary<2> stuff(ptrs, {2, 2}, {1, 1});

    auto state = stuff.initialize();

    // No cache when there is no mismatch.
    {
        stuff.match("AAAA", state);
        auto it = state.cache.find("AAAA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // Stored in cache for >1 mismatches.
    {
        stuff.match("AATA", state);
        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, 0);
        EXPECT_EQ((it->second).total, 1);
    }

    {
        stuff.match("ACTA", state);
        auto it = state.cache.find("ACTA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, -1);
    }
    
    // Checking that the reduction works correctly.
    state.cache["AATA"].index = 2;
    stuff.reduce(state);
    EXPECT_TRUE(state.cache.empty());

    {
        stuff.match("AATA", state);
        EXPECT_EQ(state.index, 2); // re-uses the cache value!
    }
}
