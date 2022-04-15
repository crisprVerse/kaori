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
}

TEST(MatchSequence, MismatchFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::MatchSequence<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

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

    // Mismatches in variable region.
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

    {
        std::string seq = "cagcatcgatcgtgaACGTTTACTGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 2, state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 3);
        EXPECT_EQ(state.identity[0].second, 2);
    }

    // Mismatches in both.
    {
        std::string seq = "cagcatcgatcgtgaACCTTTATTGCAcacggaggaga";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 2, state));

        EXPECT_EQ(state.position, 15);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 3);
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

    // Only returns the first instance, even though the second is better.
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

TEST(MatchSequence, MultiFirst) {
    std::string constant = "AAAA----CGGC------TTTT";
    std::vector<std::string> variables1 { "AAAA", "CCCC", "GGGG", "TTTT" };
    std::vector<std::string> variables2 { "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" };
    kaori::MatchSequence<128> stuff(constant.c_str(), constant.size(), true, false, std::vector<int>{0, 1}, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) });

    // Perfect match.
    {
        std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacac";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));

        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 2);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 0);
        EXPECT_EQ(state.identity[1].first, 1);
        EXPECT_EQ(state.identity[1].second, 0);
    }

    // One mismatch.
    {
        std::string seq = "cagAAAAATAACGGCTGTGTGTTTTacac";
        {
            auto state1 = stuff.initialize();
            EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 0, state1));
        }

        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));
        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 1);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 2);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 1);
        EXPECT_EQ(state.identity[1].first, 1);
        EXPECT_EQ(state.identity[1].second, 0);
    }

    // Mismatches spread across two variable regions.
    {
        std::string seq = "cagAAAAATAACGGCTGTCTGTTTTacac";

        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 2);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 1);
        EXPECT_EQ(state.identity[1].first, 1);
        EXPECT_EQ(state.identity[1].second, 1);
    }

    {
        std::string seq = "cagAAAATCAACGGCTGTGTGTTTTacac";

        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 2);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 2);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 2);
        EXPECT_EQ(state.identity[1].first, 1);
        EXPECT_EQ(state.identity[1].second, 0);
    }

    {
        std::string seq = "cagAAAATTAACGGCTGTCTGTTTTacac";
        auto state = stuff.initialize();
        EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
    }

    // Mismatches in constant and variable regions.
    {
        std::string seq = "cagAATAATAACGGCTGTCTGTTTTacac";
        {
            auto state = stuff.initialize();
            EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 1, state));
        }
        {
            auto state = stuff.initialize();
            EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
        }

        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 3, state));
        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 3);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 2);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 1);
        EXPECT_EQ(state.identity[1].first, 1);
        EXPECT_EQ(state.identity[1].second, 1);
    }

    // Rejects on ambiguities for any variable region.
    {
        std::string seq = "cagAAAATTAACGGCTGTGTGTTTTacac";
        auto state = stuff.initialize();
        EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 2, state));
    }
}

TEST(MatchSequence, FragmentedFirst) {
    std::string constant = "AAAA----CGGC------TTTT";
    std::vector<std::string> variables { "AAAATTTTTT", "CCCCGGGGGG", "GGGGAAAAAA", "TTTTAAAAAA" };
    kaori::MatchSequence<128> stuff(constant.c_str(), constant.size(), true, false, to_pointers(variables));

    // Perfect match.
    {
        std::string seq = "cagAAAAAAAACGGCTTTTTTTTTTacac";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 0, state));

        EXPECT_EQ(state.position, 3);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 0);
        EXPECT_EQ(state.identity[0].second, 0);
    }
}

TEST(MatchSequence, Caching) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    auto ptrs = to_pointers(variables);
    kaori::MatchSequence<64> stuff(constant.c_str(), constant.size(), true, false, ptrs);

    auto state = stuff.initialize();

    // No cache when there is no mismatch.
    {
        std::string seq = "cagcatcgatcgtgaACGTAAAATGCAcacggaggaga";
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));

        const auto& cache = state.forward_cache[0];
        auto it = cache.find("AAAA");
        EXPECT_TRUE(it == cache.end());
    }

    // Stored in cache for >1 mismatches.
    {
        std::string seq = "cagcatcgatcgtgaACGTAATATGCAcacggaggaga";
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));

        const auto& cache = state.forward_cache[0];
        auto it = cache.find("AATA");
        EXPECT_TRUE(it != cache.end());
        EXPECT_EQ((it->second).first, 0);
        EXPECT_EQ((it->second).second, 1);
    }

    {
        std::string seq = "cagcatcgatcgtgaACGTACTATGCAcacggaggaga";
        EXPECT_FALSE(stuff.search_first(seq.c_str(), seq.size(), 1, state));

        const auto& cache = state.forward_cache[0];
        auto it = cache.find("ACTA");
        EXPECT_TRUE(it != cache.end());
        EXPECT_EQ((it->second).first, -1);
    }
    
    // Checking that the reduction works correctly.
    state.forward_cache[0]["AATA"].first = 2;
    stuff.reduce(state);
    EXPECT_TRUE(state.forward_cache[0].empty());

    {
        std::string seq = "cagcatcgatcgtgaACGTAATATGCAcacggaggaga";
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), 1, state));
        ASSERT_EQ(state.identity.size(), 1);
        EXPECT_EQ(state.identity[0].first, 2); // re-uses the cache value!
    }
}
