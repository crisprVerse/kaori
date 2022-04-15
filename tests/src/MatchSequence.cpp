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

    // At the start or end.
    {
        std::string seq = "ACGTAAAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));
        EXPECT_EQ(state.position, 0);
    }
    {
        std::string seq = "acacacacacACGTAAAATGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));
        EXPECT_EQ(state.position, 10);
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

    // Only returns the first instance.
    {
        std::string seq = "cagACGTCCCCTGCAcacACGTAAAATGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));

        EXPECT_EQ(state.position, 3);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 1);
    }
}

TEST(MatchSequence, ReverseComplementFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::MatchSequence<64> forward_only(constant.c_str(), constant.size(), true, false, ptrs);
    kaori::MatchSequence<64> reverse_only(constant.c_str(), constant.size(), false, true, ptrs);
    kaori::MatchSequence<64> both(constant.c_str(), constant.size(), true, true, ptrs);

    // Forward only.
    {
        std::string seq = "tcgatcgtgaACGTGGGGTGCAcacggaggaga";

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), 0, stateB));
        EXPECT_EQ(stateB.position, 10);
        EXPECT_FALSE(stateB.reverse);
        ASSERT_EQ(stateB.identity.size(), 1);
        EXPECT_EQ(stateB.identity[0].first, 2);

        auto stateR = reverse_only.initialize();
        EXPECT_FALSE(reverse_only.search_first(seq.c_str(), seq.size(), 0, stateR));
    }

    // Reverse only.
    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacggaggaga";

        auto stateF = forward_only.initialize();
        EXPECT_FALSE(forward_only.search_first(seq.c_str(), seq.size(), 0, stateF));

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), 0, stateB));
        EXPECT_EQ(stateB.position, 10);
        EXPECT_TRUE(stateB.reverse);
        ASSERT_EQ(stateB.identity.size(), 1);
        EXPECT_EQ(stateB.identity[0].first, 2);

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_first(seq.c_str(), seq.size(), 0, stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        ASSERT_EQ(stateR.identity.size(), 1);
        EXPECT_EQ(stateR.identity[0].first, 2);
    }

    // Present in both.
    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacACGTTTTTTGCA";

        auto stateF = forward_only.initialize();
        EXPECT_TRUE(forward_only.search_first(seq.c_str(), seq.size(), 0, stateF));
        EXPECT_EQ(stateF.position, 25);
        EXPECT_FALSE(stateF.reverse);
        ASSERT_EQ(stateF.identity.size(), 1);
        EXPECT_EQ(stateF.identity[0].first, 3);

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), 0, stateB)); 
        EXPECT_EQ(stateB.position, 10);
        EXPECT_TRUE(stateB.reverse);

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_first(seq.c_str(), seq.size(), 0, stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        ASSERT_EQ(stateR.identity.size(), 1);
        EXPECT_EQ(stateR.identity[0].first, 2);
    }
}

TEST(MatchSequence, BasicBest) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::MatchSequence<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

    // Favors the perfect match.
    {
        std::string seq = "gatcgtgaACGTATAATGCAcacggagACGTGGGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), 1, state));

        EXPECT_EQ(state.position, 27);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 2);
        EXPECT_EQ(state.identity[0].second, 0);
    }

    {
        std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGCGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), 1, state));

        EXPECT_EQ(state.position, 8);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 0);
    }

    // ... unless it's ambiguous.
    {
        std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGGGGTGCA";
        auto state = stuff.initialize();
        EXPECT_FALSE(stuff.search_best(seq.c_str(), seq.size(), 0, state));
    }

    // Works at the start.
    {
        std::string seq = "ACGTAAAATGCAACGTGCGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), 1, state));

        EXPECT_EQ(state.position, 0);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
    }

    // Works on the other strands.
    kaori::MatchSequence<64> reverse_only(constant.c_str(), constant.size(), false, true, ptrs);
    kaori::MatchSequence<64> both(constant.c_str(), constant.size(), true, true, ptrs);

    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacACGTTTTTTGCA";

        auto stateF = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), 0, stateF));
        EXPECT_EQ(stateF.position, 25);
        EXPECT_FALSE(stateF.reverse);
        ASSERT_EQ(stateF.identity.size(), 1);
        EXPECT_EQ(stateF.identity[0].first, 3);

        auto stateB = both.initialize();
        EXPECT_FALSE(both.search_best(seq.c_str(), seq.size(), 0, stateB)); // ambiguous

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_best(seq.c_str(), seq.size(), 0, stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        ASSERT_EQ(stateR.identity.size(), 1);
        EXPECT_EQ(stateR.identity[0].first, 2);
    }
}

