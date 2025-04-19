#include <gtest/gtest.h>
#include "kaori/ScanTemplate.hpp"
#include <string>

TEST(ScanTemplate, Basic) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.reverse_variable_regions();
    EXPECT_EQ(rvar.size(), 0); // not computed if reverse = false.

    {
        std::string seq = "ACGTAAAATTTT";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_TRUE(out.finished);
    }

    {
        std::string seq = "GACGTAAAATTTTA";
        auto out = stuff.initialize(seq.c_str(), seq.size());
        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);

        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 2);
        EXPECT_TRUE(out.finished);
    }
}

TEST(ScanTemplate, TooShort) {
    std::string thing = "ACGT----TTTT"; 
    kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);

    std::string seq = "ACGT";
    auto out = stuff.initialize(seq.c_str(), seq.size());
    EXPECT_TRUE(out.finished);
}

TEST(ScanTemplate, ReverseComplement) {
    std::string thing = "ACGT----TTTT"; 

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::BOTH);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 0);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);
    }

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::REVERSE);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_TRUE(out.finished);

        const auto& fvar = stuff.forward_variable_regions();
        ASSERT_EQ(fvar.size(), 1); // still reported correctly.
        EXPECT_EQ(fvar.front().first, 4);
        EXPECT_EQ(fvar.front().second, 8);
    }

    {
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::FORWARD);
        std::string seq = "AAAATTTTACGT";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 0);
        EXPECT_TRUE(out.finished);
    }
 }

TEST(ScanTemplate, Multiple) {
    std::string thing = "ACGT----TT-----GG"; 
    kaori::ScanTemplate<32> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::BOTH);

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 2);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);
    EXPECT_EQ(fvar.back().first, 10);
    EXPECT_EQ(fvar.back().second, 15);

    auto rvar = stuff.reverse_variable_regions();
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
            EXPECT_GT(out.forward_mismatches, 1);
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

    auto fvar = stuff.forward_variable_regions();
    ASSERT_EQ(fvar.size(), 1);
    EXPECT_EQ(fvar.front().first, 4);
    EXPECT_EQ(fvar.front().second, 8);

    auto rvar = stuff.reverse_variable_regions();
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

    std::bitset<16> amask;
    for (size_t i = 0; i < thing.size(); ++i) {
        amask.set(i);
    }

    // Runs into an N later.
    {
        std::string seq = "aaaaaaaaaaACGNAAAATTTTacgatcgatcagctag";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i < 10; ++i) {
            stuff.next(out);
            EXPECT_GT(out.forward_mismatches, 1);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ(out.any_ambiguous, i > 1); // at the third base (i.e., i = 2), the N comes into the window.
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 1);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 1);
        EXPECT_TRUE(out.any_ambiguous);

        for (int i = 0; i < 3; ++i) { // next three shifts still overlap the N.
            stuff.next(out);
            EXPECT_GT(out.forward_mismatches, 1);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ((out.ambiguous & amask).count(), 1);
            EXPECT_TRUE(out.any_ambiguous);
        }

        stuff.next(out); // past the N now, so ambiguity is now dropped.
        EXPECT_GT(out.forward_mismatches, 1);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 0);
        EXPECT_FALSE(out.any_ambiguous);
    }

    // Clears existing Ns.
    {
        std::string seq = "NNNNACGTAAAATTTTa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i < 4; ++i) {
            stuff.next(out);
            EXPECT_GT(out.forward_mismatches, 0);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ((out.ambiguous & amask).count(), 4 - i);
            EXPECT_TRUE(out.any_ambiguous);
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 0);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 0);
        EXPECT_FALSE(out.any_ambiguous);
    }

    // Works with separated N's.
    {
        std::string seq = "aaaaaaaANGTAAAATTTNaaaaaaaaaaaaaa";
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i <= 6; ++i) {
            stuff.next(out);
            EXPECT_GT(out.forward_mismatches, 2);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ((out.ambiguous & amask).count(), 1);
            EXPECT_TRUE(out.any_ambiguous);
        }

        stuff.next(out);
        EXPECT_EQ(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 2);
        EXPECT_TRUE(out.any_ambiguous);

        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 2);
        EXPECT_TRUE(out.any_ambiguous);

        for (int i = 0; i <= 9; ++i) {
            stuff.next(out);
            EXPECT_GT(out.forward_mismatches, 2);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ((out.ambiguous & amask).count(), 1);
            EXPECT_TRUE(out.any_ambiguous);
        }

        stuff.next(out);
        EXPECT_GT(out.forward_mismatches, 2);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 0);
        EXPECT_FALSE(out.any_ambiguous);
    }

    // Works in reverse.
    {
        std::string seq = "aaNaaAAAATTTTACGTaaNaaaaa";
        kaori::ScanTemplate<16> stuff(thing.c_str(), thing.size(), kaori::SearchStrand::REVERSE);
        auto out = stuff.initialize(seq.c_str(), seq.size());

        for (int i = 0; i <= 4; ++i) {
            stuff.next(out);
            EXPECT_GT(out.reverse_mismatches, 1);
            EXPECT_FALSE(out.finished);
            EXPECT_EQ((out.ambiguous & amask).count(), (i <= 2)); // i.e., before we pass the N.
            EXPECT_EQ(out.any_ambiguous, (i <= 2));
        }

        stuff.next(out);
        EXPECT_EQ(out.reverse_mismatches, 0);
        EXPECT_FALSE(out.finished);
        EXPECT_EQ((out.ambiguous & amask).count(), 0);
        EXPECT_FALSE(out.any_ambiguous);

        size_t remaining = seq.size() - thing.size() - out.position;
        for (size_t i = 0; i < remaining; ++i) {
            stuff.next(out);
            EXPECT_GT(out.reverse_mismatches, 1);
            EXPECT_EQ(out.finished, i + 1 == remaining);
            EXPECT_EQ((out.ambiguous & amask).count(), (i > 1)); // i.e., after we hit the next N.
            EXPECT_EQ(out.any_ambiguous, (i > 1));
        }
    }
}
