#include <gtest/gtest.h>
#include "kaori/ConstantTemplate.hpp"
#include <string>

TEST(ConstantTemplate, Basic) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.variable_regions(false);
    ASSERT_EQ(rvar.size(), 1);
    EXPECT_EQ(rvar.front().first, 4);
    EXPECT_EQ(rvar.front().second, 8);

    {
        std::string seq = "ACGTAAAATTTT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_EQ(out.reverse_mismatches, -1);
        EXPECT_TRUE(out.finished);
    }

    {
        std::string seq = "GACGTAAAATTTTA";
        auto out = stuff.initialize(seq.c_str(), seq.size());
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

TEST(ConstantTemplate, TooShort) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);

    std::string seq = "ACGT";
    auto out = stuff.initialize(seq.c_str(), seq.size());
    EXPECT_EQ(out.forward_mismatches, -1);
    EXPECT_EQ(out.reverse_mismatches, -1);
    EXPECT_TRUE(out.finished);
}

TEST(ConstantTemplate, ReverseComplement) {
    std::string thing = "ACGT----TTTT"; 

    {
        kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, true);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        EXPECT_TRUE(out.forward_mismatches > 0);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);
    }

    {
        kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), false, true);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        EXPECT_EQ(out.forward_mismatches, -1);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);

        const auto& fvar = stuff.variable_regions();
        ASSERT_EQ(fvar.size(), 1); // still reported correctly.
        EXPECT_EQ(fvar.front().first, 4);
        EXPECT_EQ(fvar.front().second, 8);
    }

    {
        kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        EXPECT_TRUE(out.forward_mismatches > 0);
        EXPECT_EQ(out.reverse_mismatches, -1);
        EXPECT_TRUE(out.finished);
    }
 }

TEST(ConstantTemplate, Multiple) {
    std::string thing = "ACGT----TT-----GG"; 
    kaori::ConstantTemplate<128> stuff(thing.c_str(), thing.size(), true, true);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 2);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);
    EXPECT_EQ(fvar.back().first, 10);
    EXPECT_EQ(fvar.back().second, 15);

    auto rvar = stuff.variable_regions(true);
    ASSERT_EQ(rvar.size(), 2);
    EXPECT_EQ(rvar.front().first, 2);
    EXPECT_EQ(rvar.front().second, 7);
    EXPECT_EQ(rvar.back().first, 9);
    EXPECT_EQ(rvar.back().second, 13);

    {
        std::string seq = "aaaaaACGTAAAATTTTGGGGGggggg";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        
        for (int i = 0; i < 5; ++i) {
            EXPECT_TRUE(out.forward_mismatches > 1);
            EXPECT_FALSE(out.finished);
            stuff.next(out);
        }

        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
        stuff.next(out);

        for (int i = 0; i < 5; ++i) {
            EXPECT_TRUE(out.forward_mismatches > 1);
            EXPECT_EQ(out.finished, i == 4);
            if (i < 4) {
                stuff.next(out);
            }
        }
    }
}

TEST(ConstantTemplate, Mismatches) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);

    auto fvar = stuff.variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.variable_regions(true);
    EXPECT_EQ(rvar.size(), 0); // as we set reverse=false.

    {
        std::string seq = "aACGTAAAAGTTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 1);
        EXPECT_FALSE(out.finished);
    }
    
    {
        std::string seq = "aACGTAAAAGGTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);
    }

    {
        std::string seq = "aCCGTAAAAGGTTg";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 3);
        EXPECT_FALSE(out.finished);
    }
}

TEST(ConstantTemplate, BadBases) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);

    {
        std::string seq = "aACGTNNNNTTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
    }

    // Runs into an N later.
    {
        std::string seq = "aaaaaaaaaaACGNAAAATTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i < 10; ++i) {
            EXPECT_TRUE(out.forward_mismatches > 0);
            EXPECT_FALSE(out.finished);
            stuff.next(out);
        }

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
            EXPECT_TRUE(out.forward_mismatches > 0);
            EXPECT_FALSE(out.finished);
            stuff.next(out);
        }

        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
        EXPECT_TRUE(out.bad.empty());
    }
}
