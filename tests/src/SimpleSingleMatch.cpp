#include <gtest/gtest.h>
#include "kaori/SimpleSingleMatch.hpp"
#include <string>
#include <vector>
#include "utils.h"

class SimpleSingleMatchTest : public ::testing::Test {
protected:
    template<size_t N>
    using Options = typename kaori::SimpleSingleMatch<N>::Options;
};

TEST_F(SimpleSingleMatchTest, BasicFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, Options<16>());

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

TEST_F(SimpleSingleMatchTest, MismatchFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SimpleSingleMatch<16> stuff1(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.max_mismatches = 1;
        return opt;
    }());
    kaori::SimpleSingleMatch<16> stuff2(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.max_mismatches = 2;
        return opt;
    }());

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

TEST_F(SimpleSingleMatchTest, ReverseComplementFirst) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SimpleSingleMatch<16> forward_only(constant.c_str(), constant.size(), ptrs, Options<16>());
    kaori::SimpleSingleMatch<16> reverse_only(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        return opt;
    }());
    kaori::SimpleSingleMatch<16> both(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.strand = kaori::SearchStrand::BOTH;
        return opt;
    }());

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

TEST_F(SimpleSingleMatchTest, BasicBest) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.max_mismatches = 1;
        return opt;
    }());

    // Favors the perfect match.
    {
        std::string seq = "gatcgtgaACGTATAATGCAcacggagACGTGGGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 27);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 2);
        EXPECT_EQ(state.variable_mismatches, 0);

        auto fstate = stuff.initialize();
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), fstate));

        EXPECT_EQ(fstate.position, 8);
        EXPECT_EQ(fstate.mismatches, 1);
        EXPECT_FALSE(fstate.reverse);
        EXPECT_EQ(fstate.index, 0);
        EXPECT_EQ(fstate.variable_mismatches, 1);
    }

    {
        std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGCGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 8);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.variable_mismatches, 0);

        auto fstate = stuff.initialize(); // just to check as a reference.
        EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), fstate));
        EXPECT_EQ(fstate.position, 8);
    }

    // ... unless it's ambiguous.
    {
        {
            kaori::SimpleSingleMatch<16> stuff2(constant.c_str(), constant.size(), ptrs, Options<16>());
            std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTGGGGTGCA";
            auto state = stuff2.initialize();
            EXPECT_FALSE(stuff2.search_best(seq.c_str(), seq.size(), state));
        }

        // ... unless the ambiguity refers to the same barcode!
        {
            kaori::SimpleSingleMatch<16> stuff2(constant.c_str(), constant.size(), ptrs, Options<16>());
            std::string seq = "gatcgtgaACGTAAAATGCAcacggagACGTAAAATGCA";
            auto state = stuff2.initialize();
            EXPECT_TRUE(stuff2.search_best(seq.c_str(), seq.size(), state));
            EXPECT_EQ(state.position, 8);
            EXPECT_EQ(state.mismatches, 0);
        }
    }

    // Works at the start.
    {
        std::string seq = "ACGTAAAATGCAACGTGCGGTGCA";
        auto state = stuff.initialize();
        EXPECT_TRUE(stuff.search_best(seq.c_str(), seq.size(), state));

        EXPECT_EQ(state.position, 0);
        EXPECT_EQ(state.mismatches, 0);
        EXPECT_FALSE(state.reverse);
        EXPECT_EQ(state.index, 0);
    }

    // Handles ambiguity correctly within a single region.
    { 
        kaori::SimpleSingleMatch<16> stuff2(constant.c_str(), constant.size(), ptrs, [&]{
            Options<16> opt;
            opt.max_mismatches = 2;
            return opt;
        }());

        std::string seq = "ACGTATTATGCA";
        auto state = stuff2.initialize();
        EXPECT_FALSE(stuff2.search_best(seq.c_str(), seq.size(), state));

        std::string control_seq = "ACGTACGCTGCA"; // Control to check that the max_mismatches is respected.
        EXPECT_TRUE(stuff2.search_best(control_seq.c_str(), control_seq.size(), state));
        EXPECT_EQ(state.position, 0);
        EXPECT_EQ(state.index, 1);
        EXPECT_EQ(state.variable_mismatches, 2);
    }
}

TEST_F(SimpleSingleMatchTest, ReverseComplementBest) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, Options<16>());
    kaori::SimpleSingleMatch<16> reverse_only(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        return opt;
    }());
    kaori::SimpleSingleMatch<16> both(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.strand = kaori::SearchStrand::BOTH;
        return opt;
    }());

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

TEST_F(SimpleSingleMatchTest, Caching) {
    std::string constant = "ACGT----TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, [&]{
        Options<16> opt;
        opt.strand = kaori::SearchStrand::BOTH;
        opt.max_mismatches = 1;
        return opt;
    }());

    auto state = stuff.initialize();

    // Fill it up with some cache hits.
    std::string seq = "tcgatcgtgaTGCACCTCACGTcacACGTTATTTGCA";
    EXPECT_FALSE(stuff.search_best(seq.c_str(), seq.size(), state));

    // This should reuse the hits.
    EXPECT_TRUE(stuff.search_first(seq.c_str(), seq.size(), state));

    // Get some coverage on the reduction method.
    stuff.reduce(state);
}

TEST_F(SimpleSingleMatchTest, Error) {
    std::string constant = "ACGT------TGCA";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    EXPECT_ANY_THROW({
        try {
            kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, Options<16>());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("should be the same") != std::string::npos);
            throw e;
        }
    });

    EXPECT_ANY_THROW({
        try {
            constant = "ACACACCAC";
            kaori::SimpleSingleMatch<16> stuff(constant.c_str(), constant.size(), ptrs, Options<16>());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("expected one variable region") != std::string::npos);
            throw e;
        }
    });
}
