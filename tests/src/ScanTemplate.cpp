#include <gtest/gtest.h>
#include "kaori/ScanTemplate.hpp"
#include <string>

TEST(ScanTemplate, Basic) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.variable_regions<true>();
    EXPECT_EQ(rvar.size(), 0); // not computed if reverse = false.

    {
        std::string seq = "ACGTAAAATTTT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_EQ(out.reverse_mismatches, -1);
        EXPECT_TRUE(out.finished);
    }

    {
        std::string seq = "GACGTAAAATTTTA";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_TRUE(out.forward_mismatches > 2);
        EXPECT_FALSE(out.finished);

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);

        stuff.next(out);
        EXPECT_TRUE(out.forward_mismatches > 2);
        EXPECT_TRUE(out.finished);
    }
}

TEST(ScanTemplate, TooShort) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    std::string seq = "ACGT";
    auto out = stuff.initialize(seq.c_str(), seq.size());

    EXPECT_EQ(out.forward_mismatches, -1);
    EXPECT_EQ(out.reverse_mismatches, -1);
    EXPECT_TRUE(out.finished);
}

TEST(ScanTemplate, ReverseComplement) {
    std::string thing = "ACGT----TTTT"; 

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::BOTH);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_TRUE(out.forward_mismatches > 0);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);
    }

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::REVERSE);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, -1);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);

        const auto& fvar = stuff.variable_regions();
        ASSERT_EQ(fvar.size(), 1); // still reported correctly.
        EXPECT_EQ(fvar.front().first, 4);
        EXPECT_EQ(fvar.front().second, 8);
    }

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_TRUE(out.forward_mismatches > 0);
        EXPECT_EQ(out.reverse_mismatches, -1);
        EXPECT_TRUE(out.finished);
    }
 }

TEST(ScanTemplate, Multiple) {
    std::string thing = "ACGT----TT-----GG"; 
    kaori::ScanTemplate<32> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::BOTH);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 2);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);
    EXPECT_EQ(fvar.back().first, 10);
    EXPECT_EQ(fvar.back().second, 15);

    auto rvar = stuff.variable_regions<true>();
    ASSERT_EQ(rvar.size(), 2);
    EXPECT_EQ(rvar.front().first, 2);
    EXPECT_EQ(rvar.front().second, 7);
    EXPECT_EQ(rvar.back().first, 9);
    EXPECT_EQ(rvar.back().second, 13);

    {
        std::string seq = "aaaaaACGTAAAATTTTGGGGGggggg";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        
        for (int i = 0; i < 5; ++i) {
            stuff.next(out);
            EXPECT_TRUE(out.forward_mismatches > 1);
            EXPECT_FALSE(out.finished);
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);

        for (int i = 0; i < 5; ++i) {
            stuff.next(out);
            EXPECT_TRUE(out.forward_mismatches > 1);
            EXPECT_EQ(out.finished, i == 4);
        }
    }
}

TEST(ScanTemplate, Mismatches) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.variable_regions<true>();
    EXPECT_EQ(rvar.size(), 0); // as we set reverse=false.

    {
        std::string seq = "aACGTAAAAGTTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 1);
        EXPECT_FALSE(out.finished);
    }
    
    {
        std::string seq = "aACGTAAAAGGTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);
    }

    {
        std::string seq = "aCCGTAAAAGGTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 3);
        EXPECT_FALSE(out.finished);
    }
}

TEST(ScanTemplate, BadBases) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    {
        std::string seq = "aACGTNNNNTTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
    }

    // Runs into an N later.
    {
        std::string seq = "aaaaaaaaaaACGNAAAATTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i < 10; ++i) {
            stuff.next(out);
            EXPECT_TRUE(out.forward_mismatches > 0);
            EXPECT_FALSE(out.finished);
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 1);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ(out.ambiguous.count(), 4);
        EXPECT_EQ(out.bad.size(), 1);
    }

    // Clears existing Ns.
    {
        std::string seq = "NNNNACGTAAAATTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i < 4; ++i) {
            stuff.next(out);
            EXPECT_TRUE(out.forward_mismatches > 0);
            EXPECT_FALSE(out.finished);
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
        EXPECT_TRUE(out.bad.empty());
    }
}
