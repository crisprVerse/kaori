#include <gtest/gtest.h>
#include "kaori/MismatchTrie.hpp"
#include "kaori/BarcodePool.hpp"
#include <string>
#include "utils.h"

class AnyMismatchesTest : public ::testing::Test {
protected:
    static kaori::AnyMismatches populate(const kaori::BarcodePool& ptrs) {
        kaori::AnyMismatches output(ptrs.length, kaori::DuplicateAction::ERROR);
        for (auto p : ptrs.pool) {
            output.add(p);
        }
        return output;
    }
};

TEST_F(AnyMismatchesTest, Basic) {
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate(ptrs);

    {
        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search("AAAT", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("CCAG", 2);
        EXPECT_EQ(res.first, 2);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("AGTT", 0);
        EXPECT_EQ(res.first, 3);
        EXPECT_EQ(res.second, 0);
    }
}

TEST_F(AnyMismatchesTest, MoreMismatches) {
    std::vector<std::string> things { "ACGTACGTACGT", "TTTGGGCCCAAA" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate(ptrs);

    {
        auto res = stuff.search("ACGTACGTCCGT", 2);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("TCGTACGTCCGT", 2);
        EXPECT_EQ(res.first, 0);
        EXPECT_EQ(res.second, 2);
    }

    {
        auto res = stuff.search("TTTGGGGCCAAA", 2);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("TTGGGGGCCAAA", 2);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 2);
    }
}

TEST_F(AnyMismatchesTest, MismatchesWithNs) {
    std::vector<std::string> things { "ACGTACGTACGT", "TTTGGGCCCAAA" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate(ptrs);

    {
        auto res = stuff.search("ACGTACGTACGN", 0);
        EXPECT_EQ(res.first, -1);

        auto res2 = stuff.search("ACGTACGTACGN", 1);
        EXPECT_EQ(res2.first, 0);
        EXPECT_EQ(res2.second, 1);
    }

    {
        auto res = stuff.search("TTNGGGNCCAAA", 1);
        EXPECT_EQ(res.first, -1);

        auto res2 = stuff.search("TTNGGGNCCAAA", 2);
        EXPECT_EQ(res2.first, 1);
        EXPECT_EQ(res2.second, 2);
    }
}

TEST_F(AnyMismatchesTest, CappedMismatch) {
    // Force an early return.
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate(ptrs);

    {
        auto res = stuff.search("AAAT", 0);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }

    {
        auto res = stuff.search("CCAG", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 2);
    }
}

TEST_F(AnyMismatchesTest, Ambiguous) {
    std::vector<std::string> things { "AAAAGAAAA", "AAAACAAAA", "AAAAAAAAG", "AAAAAAAAC" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate(ptrs);

    {
        // Positive control first.
        auto res = stuff.search("AAAACAAAA", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search("AAAATAAAA", 1);
        EXPECT_EQ(res.first, -2);
        EXPECT_EQ(res.second, 1);
    }

    {
        // Handles ambiguity at the end of the sequence.
        auto res = stuff.search("AAAAAAAAT", 1);
        EXPECT_EQ(res.first, -2);
        EXPECT_EQ(res.second, 1);
    }
}

TEST_F(AnyMismatchesTest, Duplicates) {
    std::vector<std::string> things { "ACGT", "ACGT", "AGTT", "AGTT" };
    kaori::BarcodePool ptrs(things);

    EXPECT_ANY_THROW({
        try {
            populate(ptrs);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("duplicate") != std::string::npos);
            throw e;
        }
    });

    auto CHECK = [&](auto status, bool isdup, bool replaced, bool cleared) {
        EXPECT_FALSE(status.has_ambiguous);

        if (isdup) {
            EXPECT_TRUE(status.is_duplicate);
        } else {
            EXPECT_FALSE(status.is_duplicate);
        }

        if (replaced) {
            EXPECT_TRUE(status.duplicate_replaced);
        } else {
            EXPECT_FALSE(status.duplicate_replaced);
        }

        if (cleared) {
            EXPECT_TRUE(status.duplicate_cleared);
        } else {
            EXPECT_FALSE(status.duplicate_cleared);
        }
    };

    // Gets the first occurrence.
    {
        kaori::AnyMismatches stuff(ptrs.length, kaori::DuplicateAction::FIRST);
        CHECK(stuff.add(ptrs.pool[0]), false, false, false);
        CHECK(stuff.add(ptrs.pool[1]), true, false, false);
        CHECK(stuff.add(ptrs.pool[2]), false, false, false);
        CHECK(stuff.add(ptrs.pool[3]), true, false, false);

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 0);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, 2);
    }

    // Gets the last occurrence.
    {
        kaori::AnyMismatches stuff(ptrs.length, kaori::DuplicateAction::LAST);
        CHECK(stuff.add(ptrs.pool[0]), false, false, false);
        CHECK(stuff.add(ptrs.pool[1]), true, true, false);
        CHECK(stuff.add(ptrs.pool[2]), false, false, false);
        CHECK(stuff.add(ptrs.pool[3]), true, true, false);

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 1);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, 3);
    }

    // Gets nothing.
    {
        kaori::AnyMismatches stuff(ptrs.length, kaori::DuplicateAction::NONE);
        CHECK(stuff.add(ptrs.pool[0]), false, false, false);
        CHECK(stuff.add(ptrs.pool[1]), true, false, true);
        CHECK(stuff.add(ptrs.pool[2]), false, false, false);
        CHECK(stuff.add(ptrs.pool[3]), true, false, true);
        CHECK(stuff.add(ptrs.pool[3]), true, false, false); // next addition sets duplicate_cleared = false as it's already cleared.

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, -2);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, -2);
    }
}

TEST_F(AnyMismatchesTest, Iupac) {
    {
        std::vector<std::string> things { "rACsCGk", "YacWcgM" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate(ptrs);

        EXPECT_EQ(stuff.search("AacGcgT", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("GacCcgG", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("CacTcgA", 0), std::make_pair(1, 0));
        EXPECT_EQ(stuff.search("TacAcgC", 0), std::make_pair(1, 0));

        EXPECT_EQ(stuff.search("AacAcgA", 0).first, -1);
    }

    {
        std::vector<std::string> things { "Bcgt", "aDgt", "acHt", "acgV" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate(ptrs);

        EXPECT_EQ(stuff.search("ccgt", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("aagt", 0), std::make_pair(1, 0));
        EXPECT_EQ(stuff.search("actt", 0), std::make_pair(2, 0));
        EXPECT_EQ(stuff.search("acgg", 0), std::make_pair(3, 0));

        EXPECT_EQ(stuff.search("acgt", 0).first, -1);
    }

    {
        std::vector<std::string> things { "ANNA", "CNNC" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate(ptrs);

        EXPECT_EQ(stuff.search("acga", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("catc", 0), std::make_pair(1, 0));

        EXPECT_EQ(stuff.search("catg", 0).first, -1);
    }

    {
        std::vector<std::string> things { "__U__" };
        kaori::BarcodePool ptrs(things);

        EXPECT_ANY_THROW({
            try {
                populate(ptrs);
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("unknown base") != std::string::npos);
                throw e;
            }
        });
    }

    {
        kaori::AnyMismatches stuff(6, kaori::DuplicateAction::ERROR);
        EXPECT_FALSE(stuff.add("AAAAAA").has_ambiguous);
        EXPECT_TRUE(stuff.add("RYSWKM").has_ambiguous);
    }
}

TEST_F(AnyMismatchesTest, IupacAmbiguity) {
    // Avoids detecting ambiguity for IUPAC-induced mismatches to the same sequence.
    {
        std::vector<std::string> things { "AAAAAAB", "TTTVTTT" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate(ptrs);
        EXPECT_EQ(stuff.search("AAAAAAA", 1), std::make_pair(0, 1));
        EXPECT_EQ(stuff.search("AAAAAAC", 1), std::make_pair(0, 0)); // control
        EXPECT_EQ(stuff.search("TTTTTTT", 1), std::make_pair(1, 1));
        EXPECT_EQ(stuff.search("TTTGTTT", 1), std::make_pair(1, 0)); // control
    }

    // Respects ambiguity from mismatches.
    {
        std::vector<std::string> things { "AAAAAAB", "AAAAAAY", "TTTVTTT", "TTTRTTT" };
        kaori::BarcodePool ptrs(things);
        kaori::AnyMismatches stuff(ptrs.length, kaori::DuplicateAction::NONE);
        for (auto p : ptrs.pool) {
            stuff.add(p);
        }

        EXPECT_EQ(stuff.search("AAAAAAA", 1), std::make_pair(-2, 1)); // ambiguous
        EXPECT_EQ(stuff.search("AAAAAAC", 1), std::make_pair(-2, 0)); // still ambiguous, even with no mismatch
        EXPECT_EQ(stuff.search("AAAAAAG", 1), std::make_pair(0, 0)); // okay

        EXPECT_EQ(stuff.search("TTTTTTT", 1), std::make_pair(-2, 1)); // ambiguous
        EXPECT_EQ(stuff.search("TTTATTT", 1), std::make_pair(-2, 0)); // still ambiguous, even with no mismatch
        EXPECT_EQ(stuff.search("TTTCTTT", 1), std::make_pair(2, 0)); // okay

        // Unless we want to report duplicates.
        for (int i = 0; i < 2; ++i) {
            kaori::AnyMismatches stuff2(ptrs.length, (i == 0 ? kaori::DuplicateAction::FIRST : kaori::DuplicateAction::LAST));
            for (auto p : ptrs.pool) {
                stuff2.add(p);
            }

            if (i == 0) {
                EXPECT_EQ(stuff2.search("AAAAAAA", 1), std::make_pair(0, 1)); 
                EXPECT_EQ(stuff2.search("AAAAAAC", 1), std::make_pair(0, 0));
            } else {
                EXPECT_EQ(stuff2.search("AAAAAAA", 1), std::make_pair(1, 1)); 
                EXPECT_EQ(stuff2.search("AAAAAAC", 1), std::make_pair(1, 0)); 
            }
            EXPECT_EQ(stuff2.search("AAAAAAG", 1), std::make_pair(0, 0));

            if (i == 0) {
                EXPECT_EQ(stuff2.search("TTTTTTT", 1), std::make_pair(2, 1)); 
                EXPECT_EQ(stuff2.search("TTTATTT", 1), std::make_pair(2, 0));
            } else {
                EXPECT_EQ(stuff2.search("TTTTTTT", 1), std::make_pair(3, 1)); 
                EXPECT_EQ(stuff2.search("TTTATTT", 1), std::make_pair(3, 0));
            }
            EXPECT_EQ(stuff2.search("TTTCTTT", 1), std::make_pair(2, 0));
        }
    }
}

TEST_F(AnyMismatchesTest, Optimized) {
    // Deliberately not in any order, to check whether optimization behaves
    // correctly. We try it with and without optimization.
    std::vector<std::string> things { "ACCA", "TGCC", "CAAA", "AACT", "CGCG", "GGTG" };
    kaori::BarcodePool ptrs(things);

    for (size_t i = 0; i < 2; ++i) {
        auto stuff = populate(ptrs);
        if (i != 0) {
            stuff.optimize();
        }

        {
            auto res = stuff.search("AAAT", 1);
            EXPECT_EQ(res.first, 3);
            EXPECT_EQ(res.second, 1);
        }

        {
            auto res = stuff.search("CCCG", 1);
            EXPECT_EQ(res.first, 4);
            EXPECT_EQ(res.second, 1);
        }

        {
            auto res = stuff.search("GGTG", 1);
            EXPECT_EQ(res.first, 5);
            EXPECT_EQ(res.second, 0);
        }

        // Not found with the requested number of mismatches.
        {
            auto res = stuff.search("AGGA", 1);
            EXPECT_EQ(res.first, -1);
            EXPECT_EQ(res.second, 2);
        }

        {
            auto res = stuff.search("AGGA", 2);
            EXPECT_EQ(res.first, 0);
            EXPECT_EQ(res.second, 2);
        }

        // Ambiguous.
        {
            auto res = stuff.search("TTTT", 3);
            EXPECT_EQ(res.first, -2);
            EXPECT_EQ(res.second, 3);
        }
    }
}

class SegmentedMismatchesTest : public ::testing::Test {
protected:
    template<size_t num_segments>
    static kaori::SegmentedMismatches<num_segments> populate(const kaori::BarcodePool& ptrs, std::array<int, num_segments> segments) {
        kaori::SegmentedMismatches output(segments, kaori::DuplicateAction::ERROR);
        assert(ptrs.length == output.get_length());
        for (auto p : ptrs.pool) {
            output.add(p);
        }
        return output;
    }
};

TEST_F(SegmentedMismatchesTest, Segmented) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate<2>(ptrs, {4, 2});

    {
        auto res = stuff.search("AAAAAAA", { 0, 0 });
        EXPECT_EQ(res.index, 0);
        EXPECT_EQ(res.total, 0);
        EXPECT_EQ(res.per_segment[0], 0);
        EXPECT_EQ(res.per_segment[1], 0);
    }

    {
        auto res = stuff.search("TTTTTT", { 0, 0 });
        EXPECT_EQ(res.index, 3);
        EXPECT_EQ(res.total, 0);
        EXPECT_EQ(res.per_segment[0], 0);
        EXPECT_EQ(res.per_segment[1], 0);
    }

    // Fails on one mismatch.
    {
        auto res = stuff.search("CCCCCTC", { 0, 0 });
        EXPECT_EQ(res.index, -1);
    }
}

TEST_F(SegmentedMismatchesTest, Mismatches) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate<2>(ptrs, {4, 2});

    // Handles one mismatch.
    {
        auto res = stuff.search("CCCCTC", { 0, 1 });
        EXPECT_EQ(res.index, 1);
        EXPECT_EQ(res.total, 1);
        EXPECT_EQ(res.per_segment[0], 0);
        EXPECT_EQ(res.per_segment[1], 1);
    }

    // But not in the wrong place.
    {
        auto res = stuff.search("CCCCTC", { 1, 0 });
        EXPECT_EQ(res.index, -1);
    }

    // Testing the handling of mismatches at the end of each segment.
    {
        auto res = stuff.search("TTTTTA", { 0, 1 });
        EXPECT_EQ(res.index, 3);
        EXPECT_EQ(res.total, 1);
        EXPECT_EQ(res.per_segment[0], 0);
        EXPECT_EQ(res.per_segment[1], 1);
    }

    {
        auto res = stuff.search("AAACAA", { 1, 0 });
        EXPECT_EQ(res.index, 0);
        EXPECT_EQ(res.total, 1);
        EXPECT_EQ(res.per_segment[0], 1);
        EXPECT_EQ(res.per_segment[1], 0);
    }

    // Mismatches in both.
    {
        auto res = stuff.search("GGTGGC", { 1, 1 });
        EXPECT_EQ(res.index, 2);
        EXPECT_EQ(res.total, 2);
        EXPECT_EQ(res.per_segment[0], 1);
        EXPECT_EQ(res.per_segment[1], 1);
    }

    // More mismatches.
    {
        auto res = stuff.search("GGTGTC", { 1, 1 });
        EXPECT_EQ(res.index, -1);
    }

    {
        auto res = stuff.search("GGTGTC", { 2, 2 });
        EXPECT_EQ(res.index, 2);
        EXPECT_EQ(res.total, 3);
        EXPECT_EQ(res.per_segment[0], 1);
        EXPECT_EQ(res.per_segment[1], 2);
    }
}

TEST_F(SegmentedMismatchesTest, MismatchesWithNs) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    auto stuff = populate<2>(ptrs, {4, 2});

    {
        auto res = stuff.search("CCCCNC", { 0, 1 });
        EXPECT_EQ(res.index, 1);
        EXPECT_EQ(res.total, 1);
        EXPECT_EQ(res.per_segment[0], 0);
        EXPECT_EQ(res.per_segment[1], 1);
    }

    {
        auto res = stuff.search("GNGGGN", { 1, 2 });
        EXPECT_EQ(res.index, 2);
        EXPECT_EQ(res.total, 2);
        EXPECT_EQ(res.per_segment[0], 1);
        EXPECT_EQ(res.per_segment[1], 1);
    }

    // Not in the wrong place, though.
    {
        auto res = stuff.search("CCCCNC", { 1, 0 });
        EXPECT_EQ(res.index, -1);
    }
}

TEST_F(SegmentedMismatchesTest, Ambiguity) {
    {
        std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate<2>(ptrs, {4, 2});

        // Handles ambiguity properly.
        {
            auto res = stuff.search("TTGGTG", { 2, 1 });
            EXPECT_EQ(res.index, -2);
        }

        {
            auto res = stuff.search("TTTGGG", { 3, 2 });
            EXPECT_EQ(res.index, -2);
        }

        // Not ambiguous due to mismatch restrictions.
        {
            auto res = stuff.search("GGGTTT", { 2, 2 });
            EXPECT_EQ(res.index, 2);
        }
    }

    // Handles ambiguity properly at the end
    {
        std::vector<std::string> things { "AAAAAA", "AAAAAT" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate<2>(ptrs, {4, 2});

        {
            auto res = stuff.search("AAAAAC", { 0, 1 });
            EXPECT_EQ(res.index, -2);
        }

        {
            auto res = stuff.search("CAAAAG", { 1, 1 });
            EXPECT_EQ(res.index, -2);
        }
    }
}

TEST_F(SegmentedMismatchesTest, Optimized) {
    // Deliberately out of order.
    std::vector<std::string> things { "CCCCCC", "AAAAAA", "TTTTTT", "GGGGGG" };
    kaori::BarcodePool ptrs(things);

    for (size_t i = 0; i < 2; ++i) {
        auto stuff = populate<2>(ptrs, {3, 3});
        if (i != 0) {
            stuff.optimize();
        }

        {
            auto res = stuff.search("GGGGGG", { 1, 1 });
            EXPECT_EQ(res.index, 3);
            EXPECT_EQ(res.total, 0);
            EXPECT_EQ(res.per_segment[0], 0);
            EXPECT_EQ(res.per_segment[1], 0);
        }

        {
            auto res = stuff.search("TTTTTT", { 0, 0 });
            EXPECT_EQ(res.index, 2);
            EXPECT_EQ(res.total, 0);
            EXPECT_EQ(res.per_segment[0], 0);
            EXPECT_EQ(res.per_segment[1], 0);
        }

        // Handles a few mismatches.
        {
            auto res = stuff.search("GAAAAA", { 1, 1 });
            EXPECT_EQ(res.index, 1);
            EXPECT_EQ(res.total, 1);
            EXPECT_EQ(res.per_segment[0], 1);
            EXPECT_EQ(res.per_segment[1], 0);
        }

        {
            auto res = stuff.search("CACCTC", { 1, 1 });
            EXPECT_EQ(res.index, 0);
            EXPECT_EQ(res.total, 2);
            EXPECT_EQ(res.per_segment[0], 1);
            EXPECT_EQ(res.per_segment[1], 1);
        }

      {
          auto res = stuff.search("AACAAT", { 1, 1 });
          EXPECT_EQ(res.index, 1);
          EXPECT_EQ(res.total, 2);
          EXPECT_EQ(res.per_segment[0], 1);
          EXPECT_EQ(res.per_segment[1], 1);
      }

        {
            auto res = stuff.search("AAAAGG", { 2, 2 });
            EXPECT_EQ(res.index, 1);
            EXPECT_EQ(res.total, 2);
            EXPECT_EQ(res.per_segment[0], 0);
            EXPECT_EQ(res.per_segment[1], 2);
        }

        // But not in the wrong place.
        {
            auto res = stuff.search("CACCTC", { 1, 0 });
            EXPECT_EQ(res.index, -1);
        }

        // Ambiguous.
        {
            auto res = stuff.search("CCAAAC", { 2, 2 });
            EXPECT_EQ(res.index, -2);
            EXPECT_EQ(res.total, 3);
        }
    }
}

TEST_F(SegmentedMismatchesTest, Iupac) {
    // Avoids detecting ambiguity for IUPAC-induced mismatches to the same sequence.
    {
        std::vector<std::string> things { "AAAAAAB", "TTTVTTT" };
        kaori::BarcodePool ptrs(things);
        auto stuff = populate<2>(ptrs, { 3, 4 });

        auto res = stuff.search("AAAAAAA", {0, 1});
        EXPECT_EQ(res.index, 0);
        EXPECT_EQ(res.total, 1);

        res = stuff.search("AAAAAAC", {0, 1}); // control
        EXPECT_EQ(res.index, 0);
        EXPECT_EQ(res.total, 0);

        res = stuff.search("TTTTTTT", {0, 1});
        EXPECT_EQ(res.index, 1);
        EXPECT_EQ(res.total, 1);

        res = stuff.search("TTTGTTT", {0, 1}); // control
        EXPECT_EQ(res.index, 1);
        EXPECT_EQ(res.total, 0);
    }

    // Respects ambiguity from mismatches.
    {
        std::vector<std::string> things { "AAAAAAB", "AAAAAAY", "TTTVTTT", "TTTRTTT" };
        kaori::BarcodePool ptrs(things);
        kaori::SegmentedMismatches<2> stuff({ 4, 3 }, kaori::DuplicateAction::NONE);
        for (auto p : ptrs.pool) {
            stuff.add(p);
        }

        auto res = stuff.search("AAAAAAA", {0, 1}); // ambiguous.
        EXPECT_EQ(res.index, -2);
        EXPECT_EQ(res.total, 1);

        res = stuff.search("AAAAAAC", {0, 1}); // still ambiguous.
        EXPECT_EQ(res.index, -2);
        EXPECT_EQ(res.total, 0);

        res = stuff.search("AAAAAAG", {0, 1}); // okay
        EXPECT_EQ(res.index, 0);
        EXPECT_EQ(res.total, 0);

        res = stuff.search("TTTTTTT", {1, 0}); // ambiguous
        EXPECT_EQ(res.index, -2);
        EXPECT_EQ(res.total, 1);

        res = stuff.search("TTTATTT", {1, 0}); // still ambiguous
        EXPECT_EQ(res.index, -2);
        EXPECT_EQ(res.total, 0);

        res = stuff.search("TTTCTTT", {1, 0}); // okay
        EXPECT_EQ(res.index, 2);
        EXPECT_EQ(res.total, 0);

        // Unless we want to report duplicates.
        for (int i = 0; i < 2; ++i) {
            kaori::SegmentedMismatches<2> stuff2({ 4, 3 }, (i == 0 ? kaori::DuplicateAction::FIRST : kaori::DuplicateAction::LAST));
            for (auto p : ptrs.pool) {
                stuff2.add(p);
            }

            auto res = stuff2.search("AAAAAAA", {0, 1});
            if (i == 0) {
                EXPECT_EQ(res.index, 0);
            } else {
                EXPECT_EQ(res.index, 1);
            }
            EXPECT_EQ(res.total, 1);

            res = stuff2.search("AAAAAAC", {0, 1});
            if (i == 0) {
                EXPECT_EQ(res.index, 0);
            } else {
                EXPECT_EQ(res.index, 1);
            }
            EXPECT_EQ(res.total, 0);

            res = stuff2.search("AAAAAAG", {0, 1});
            EXPECT_EQ(res.index, 0);
            EXPECT_EQ(res.total, 0);

            res = stuff2.search("TTTTTTT", {1, 0});
            if (i == 0) {
                EXPECT_EQ(res.index, 2);
            } else {
                EXPECT_EQ(res.index, 3);
            }
            EXPECT_EQ(res.total, 1);

            res = stuff2.search("TTTATTT", {1, 0});
            if (i == 0) {
                EXPECT_EQ(res.index, 2);
            } else {
                EXPECT_EQ(res.index, 3);
            }
            EXPECT_EQ(res.total, 0);

            res = stuff2.search("TTTCTTT", {1, 0});
            EXPECT_EQ(res.index, 2);
            EXPECT_EQ(res.total, 0);
        }
    }
}
