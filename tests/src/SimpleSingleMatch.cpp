#include <gtest/gtest.h>
#include "kaori/SimpleSingleMatch.hpp"
#include <string>
#include <vector>
#include "utils.h"

TEST(SimpleSingleMatch, BasicFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

    // Perfect match.
    {
        std::string seq = "cagcatcgatcgtgaACGTAAAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.variable_mismatches, 0);
    }

    // At the start or end.
    {
        std::string seq = "ACGTAAAATGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), state));
        EXPECT_EQ(state.position, 0);
    }
    {
        std::string seq = "acacacacacACGTAAAATGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), state));
        EXPECT_EQ(state.position, 10);
    }
}

TEST(SimpleSingleMatch, MismatchFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> stuff1(constant.c_str(), constant.size(), true, false, ptrs, 1);
    kaori::SimpleSingleMatch<64> stuff2(constant.c_str(), constant.size(), true, false, ptrs, 2);

    // 1 mismatch in constant region.
    {
        std::string seq = "cagcatcgatcgtgaACGGAAAATGCAcacggaggaga";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 1);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.variable_mismatches, 0);
    }

    // Mismatches in variable region.
    {
        std::string seq = "cagcatcgatcgtgaACGTATAATGCAcacggaggaga";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 1);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.variable_mismatches, 1);
    }

    {
        std::string seq = "cagcatcgatcgtgaACGTTTACTGCAcacggaggaga";
        auto state = stuff2.initialize();
        EXPECT_TRUE(stuff2.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 3);
        EXPECT_EQ(state.variable_mismatches, 2);
    }

    // Mismatches in both.
    {
        std::string seq = "cagcatcgatcgtgaACCTTTATTGCAcacggaggaga";
        auto state = stuff2.initialize();
        EXPECT_TRUE(stuff2.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 3);
        EXPECT_EQ(state.variable_mismatches, 1);
    }

    // Blocked for exceeding the number of mismatches.
    {
        std::string seq = "cagcatcgatcgtgaACGTATTATGCAcacggaggaga";
        auto state = stuff1.initialize();
        EXPECT_FALSE(stuff1.search_first(seq.c_str(), seq.size(), state));
    }

    // Blocked for ambiguity.
    {
        std::string seq = "cagcatcgatcgtgaACGTATTATGCAcacggaggaga";
        auto state = stuff2.initialize();
        EXPECT_FALSE(stuff2.search_first(seq.c_str(), seq.size(), state));
    }

    // Only returns the first instance, even though the second is better.
    {
        std::string seq = "cagACGTCCCCTGCAcacACGTAAAATGCA";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_first(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.index, 1);
    }
}

TEST(SimpleSingleMatch, ReverseComplementFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> forward_only(constant.c_str(), constant.size(), true, false, ptrs);
    kaori::SimpleSingleMatch<64> reverse_only(constant.c_str(), constant.size(), false, true, ptrs);
    kaori::SimpleSingleMatch<64> both(constant.c_str(), constant.size(), true, true, ptrs);

    // Forward only.
    {
        std::string seq = "tcgatcgtgaACGTGGGGTGCAcacggaggaga";

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), stateB));
        EXPECT_EQ(stateB.position, 10);
        EXPECT_FALSE(stateB.reverse);
        EXPECT_EQ(stateB.index, 2);

        auto stateR = reverse_only.initialize();
        EXPECT_FALSE(reverse_only.search_first(seq.c_str(), seq.size(), stateR));
    }

    // Reverse only.
    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacggaggaga";

        auto stateF = forward_only.initialize();
        EXPECT_FALSE(forward_only.search_first(seq.c_str(), seq.size(), stateF));

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), stateB));
        EXPECT_EQ(stateB.position, 10);
        EXPECT_TRUE(stateB.reverse);
        EXPECT_EQ(stateB.index, 2);

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_first(seq.c_str(), seq.size(), stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        EXPECT_EQ(stateR.index, 2);
    }

    // Present in both.
    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacACGTTTTTTGCA";

        auto stateF = forward_only.initialize();
        EXPECT_TRUE(forward_only.search_first(seq.c_str(), seq.size(), stateF));
        EXPECT_EQ(stateF.position, 25);
        EXPECT_FALSE(stateF.reverse);
        EXPECT_EQ(stateF.index, 3);

        auto stateB = both.initialize();
        EXPECT_TRUE(both.search_first(seq.c_str(), seq.size(), stateB)); 
        EXPECT_EQ(stateB.position, 10);
        EXPECT_TRUE(stateB.reverse);

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_first(seq.c_str(), seq.size(), stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        EXPECT_EQ(stateR.index, 2);
    }
}

TEST(SimpleSingleMatch, BasicBest) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> stuff1(constant.c_str(), constant.size(), true, false, ptrs, 1);

    // Favors the perfect match.
    {
        std::string seq = "gatcgtgaACGTATAATGCAcacggagACGTGGGGTGCA";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 27);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 2);
        EXPECT_EQ(state.variable_mismatches, 0);
    }

    {
        std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGCGGTGCA";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 8);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.variable_mismatches, 0);
    }

    // ... unless it's ambiguous.
    {
        kaori::SimpleSingleMatch<64> stuff0(constant.c_str(), constant.size(), true, false, ptrs);
        std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGGGGTGCA";
        auto state = stuff0.initialize();
        EXPECT_FALSE(stuff0.search_best(seq.c_str(), seq.size(), state));
    }

    // Works at the start.
    {
        std::string seq = "ACGTAAAATGCAACGTGCGGTGCA";
        auto state = stuff1.initialize();
        EXPECT_TRUE(stuff1.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 0);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
    }

    // Handles ambiguity correctly within a single region.
    { 
        kaori::SimpleSingleMatch<64> stuff2(constant.c_str(), constant.size(), true, false, ptrs);
        std::string seq = "ACGTATTATGCA";
        auto state = stuff2.initialize();
        EXPECT_FALSE(stuff2.search_best(seq.c_str(), seq.size(), state));
    }
}

TEST(SimpleSingleMatch, ReverseComplementBest) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

    kaori::SimpleSingleMatch<64> reverse_only(constant.c_str(), constant.size(), false, true, ptrs);
    kaori::SimpleSingleMatch<64> both(constant.c_str(), constant.size(), true, true, ptrs);

    // Handles hits in both directions.
    {
        std::string seq = "tcgatcgtgaTGCACCCCACGTcacACGTTTTTTGCA";

        auto stateF = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), stateF));
        EXPECT_EQ(stateF.position, 25);
        EXPECT_FALSE(stateF.reverse);
        EXPECT_EQ(stateF.index, 3);

        auto stateB = both.initialize();
        EXPECT_FALSE(both.search_best(seq.c_str(), seq.size(), stateB)); // ambiguous

        auto stateR = reverse_only.initialize();
        EXPECT_TRUE(reverse_only.search_best(seq.c_str(), seq.size(), stateR));
        EXPECT_EQ(stateR.position, 10);
        EXPECT_TRUE(stateR.reverse);
        EXPECT_EQ(stateR.index, 2);
    }
}

TEST(SimpleSingleMatch, Caching) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::SimpleSingleMatch<64> stuff(constant.c_str(), constant.size(), true, true, ptrs, 1);

    auto state = stuff.initialize();

    // Fill it up with some cache hits.
    std::string seq = "tcgatcgtgaTGCACCTCACGTcacACGTTATTTGCA";
    EXPECT_FALSE(stuff.search_best(seq.c_str(), seq.size(), state));

    // This should reuse the hits.
    EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), state));

    // Get some coverage on the reduction method.
    stuff.reduce(state);
}

