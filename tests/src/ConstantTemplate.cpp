#include <gtest/gtest.h>
#include "kaori/ConstantTemplate.hpp"
#include <string>

TEST(ConstantTemplate, Basic) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ConstantTemplate<64> stuff(thing.c_str(), thing.size(), true, false);

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    ASSERT_EQ(fvar.front().first, 4);
    ASSERT_EQ(fvar.front().second, 8);

    auto rvar = stuff.reverse_variable_regions();
    ASSERT_EQ(rvar.size(), 1);
    ASSERT_EQ(rvar.front().first, 4);
    ASSERT_EQ(rvar.front().second, 8);

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
    kaori::ConstantTemplate<128> stuff(thing.c_str(), thing.size(), true, false);

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 2);
    ASSERT_EQ(fvar.front().first, 4);
    ASSERT_EQ(fvar.front().second, 8);
    ASSERT_EQ(fvar.back().first, 10);
    ASSERT_EQ(fvar.back().second, 15);

    auto rvar = stuff.reverse_variable_regions();
    ASSERT_EQ(rvar.size(), 2);
    ASSERT_EQ(rvar.front().first, 2);
    ASSERT_EQ(rvar.front().second, 7);
    ASSERT_EQ(rvar.back().first, 9);
    ASSERT_EQ(rvar.back().second, 13);

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

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    ASSERT_EQ(fvar.front().first, 4);
    ASSERT_EQ(fvar.front().second, 8);

    auto rvar = stuff.reverse_variable_regions();
    ASSERT_EQ(rvar.size(), 1);
    ASSERT_EQ(rvar.front().first, 4);
    ASSERT_EQ(rvar.front().second, 8);

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
