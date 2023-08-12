#include <gtest/gtest.h>
#include "kaori/MismatchTrie.hpp"
#include <string>
#include "utils.h"

TEST(AnyMismatches, Basic) {
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };
    kaori::BarcodePool ptrs(things);
    kaori::AnyMismatches stuff(ptrs);

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

TEST(AnyMismatches, MoreMismatches) {
    std::vector<std::string> things { "ACGTACGTACGT", "TTTGGGCCCAAA" };
    kaori::BarcodePool ptrs(things);
    kaori::AnyMismatches stuff(ptrs);

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

TEST(AnyMismatches, MismatchesWithNs) {
    std::vector<std::string> things { "ACGTACGTACGT", "TTTGGGCCCAAA" };
    kaori::BarcodePool ptrs(things);
    kaori::AnyMismatches stuff(ptrs);

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

TEST(AnyMismatches, CappedMismatch) {
    // Force an early return.
    std::vector<std::string> things { "ACGT", "AAAA", "ACAA", "AGTT" };
    kaori::BarcodePool ptrs(things);
    kaori::AnyMismatches stuff(ptrs);

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

TEST(AnyMismatches, Ambiguous) {
    std::vector<std::string> things { "AAAAGAAAA", "AAAACAAAA", "AAAAAAAAG", "AAAAAAAAC" };
    kaori::BarcodePool ptrs(things);
    kaori::AnyMismatches stuff(ptrs);

    {
        // Positive control first.
        auto res = stuff.search("AAAACAAAA", 1);
        EXPECT_EQ(res.first, 1);
        EXPECT_EQ(res.second, 0);
    }

    {
        auto res = stuff.search("AAAATAAAA", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }

    {
        // Handles ambiguity at the end of the sequence.
        auto res = stuff.search("AAAAAAAAT", 1);
        EXPECT_EQ(res.first, -1);
        EXPECT_EQ(res.second, 1);
    }
}

TEST(AnyMismatches, Duplicates) {
    std::vector<std::string> things { "ACGT", "ACGT", "AGTT", "AGTT" };
    kaori::BarcodePool ptrs(things);

    EXPECT_ANY_THROW({
        try {
            kaori::AnyMismatches stuff(ptrs);
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
        kaori::AnyMismatches stuff(ptrs.length);
        CHECK(stuff.add(ptrs.pool[0], kaori::DuplicateAction::FIRST), false, false, false);
        CHECK(stuff.add(ptrs.pool[1], kaori::DuplicateAction::FIRST), true, false, false);
        CHECK(stuff.add(ptrs.pool[2], kaori::DuplicateAction::FIRST), false, false, false);
        CHECK(stuff.add(ptrs.pool[3], kaori::DuplicateAction::FIRST), true, false, false);

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 0);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, 2);
    }

    // Gets the last occurrence.
    {
        kaori::AnyMismatches stuff(ptrs.length);
        CHECK(stuff.add(ptrs.pool[0], kaori::DuplicateAction::LAST), false, false, false);
        CHECK(stuff.add(ptrs.pool[1], kaori::DuplicateAction::LAST), true, true, false);
        CHECK(stuff.add(ptrs.pool[2], kaori::DuplicateAction::LAST), false, false, false);
        CHECK(stuff.add(ptrs.pool[3], kaori::DuplicateAction::LAST), true, true, false);

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, 1);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, 3);
    }

    // Gets nothing.
    {
        kaori::AnyMismatches stuff(ptrs.length);
        CHECK(stuff.add(ptrs.pool[0], kaori::DuplicateAction::NONE), false, false, false);
        CHECK(stuff.add(ptrs.pool[1], kaori::DuplicateAction::NONE), true, false, true);
        CHECK(stuff.add(ptrs.pool[2], kaori::DuplicateAction::NONE), false, false, false);
        CHECK(stuff.add(ptrs.pool[3], kaori::DuplicateAction::NONE), true, false, true);
        CHECK(stuff.add(ptrs.pool[3], kaori::DuplicateAction::NONE), true, false, false); // next addition sets duplicate_cleared = false as it's already cleared.

        auto res = stuff.search("ACGT", 0);
        EXPECT_EQ(res.first, -1);

        auto res2 = stuff.search("AGTT", 0);
        EXPECT_EQ(res2.first, -1);
    }
}

TEST(AnyMismatches, Iupac) {
    {
        std::vector<std::string> things { "rACsCGk", "YacWcgM" };
        kaori::BarcodePool ptrs(things);
        kaori::AnyMismatches stuff(ptrs);

        EXPECT_EQ(stuff.search("AacGcgT", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("GacCcgG", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("CacTcgA", 0), std::make_pair(1, 0));
        EXPECT_EQ(stuff.search("TacAcgC", 0), std::make_pair(1, 0));

        EXPECT_EQ(stuff.search("AacAcgA", 0).first, -1);
    }

    {
        std::vector<std::string> things { "Bcgt", "aDgt", "acHt", "acgV" };
        kaori::BarcodePool ptrs(things);
        kaori::AnyMismatches stuff(ptrs);

        EXPECT_EQ(stuff.search("ccgt", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("aagt", 0), std::make_pair(1, 0));
        EXPECT_EQ(stuff.search("actt", 0), std::make_pair(2, 0));
        EXPECT_EQ(stuff.search("acgg", 0), std::make_pair(3, 0));

        EXPECT_EQ(stuff.search("acgt", 0).first, -1);
    }

    {
        std::vector<std::string> things { "ANNA", "CNNC" };
        kaori::BarcodePool ptrs(things);
        kaori::AnyMismatches stuff(ptrs);

        EXPECT_EQ(stuff.search("acga", 0), std::make_pair(0, 0));
        EXPECT_EQ(stuff.search("catc", 0), std::make_pair(1, 0));

        EXPECT_EQ(stuff.search("catg", 0).first, -1);
    }

    {
        std::vector<std::string> things { "__U__" };
        kaori::BarcodePool ptrs(things);

        EXPECT_ANY_THROW({
            try {
                kaori::AnyMismatches stuff(ptrs);
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("unknown base") != std::string::npos);
                throw e;
            }
        });
    }

    {
        kaori::AnyMismatches stuff(6);
        EXPECT_FALSE(stuff.add("AAAAAA", kaori::DuplicateAction::ERROR).has_ambiguous);
        EXPECT_TRUE(stuff.add("RYSWKM", kaori::DuplicateAction::ERROR).has_ambiguous);
    }
}

TEST(AnyMismatches, Optimized) {
    // Deliberately not in any order, to check whether optimization behaves
    // correctly. We try it with and without optimization.
    std::vector<std::string> things { "ACCA", "TGCC", "CAAA", "AACT", "CGCG", "GGTG" };
    kaori::BarcodePool ptrs(things);

    for (size_t i = 0; i < 2; ++i) {
        kaori::AnyMismatches stuff(ptrs);
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
            EXPECT_EQ(res.first, -1);
            EXPECT_EQ(res.second, 3);
        }
    }
}

TEST(SegmentedMismatches, Segmented) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    kaori::SegmentedMismatches<2> stuff(ptrs, {4, 2});

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

TEST(SegmentedMismatches, Mismatches) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    kaori::SegmentedMismatches<2> stuff(ptrs, {4, 2});

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

TEST(SegmentedMismatches, MismatchesWithNs) {
    std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
    kaori::BarcodePool ptrs(things);
    kaori::SegmentedMismatches<2> stuff(ptrs, {4, 2});

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

TEST(SegmentedMismatches, Ambiguity) {
    {
        std::vector<std::string> things { "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT" };
        kaori::BarcodePool ptrs(things);
        kaori::SegmentedMismatches<2> stuff(ptrs, {4, 2});

        // Handles ambiguity properly.
        {
            auto res = stuff.search("TTGGTG", { 2, 1 });
            EXPECT_EQ(res.index, -1);
        }

        {
            auto res = stuff.search("TTTGGG", { 3, 2 });
            EXPECT_EQ(res.index, -1);
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
        kaori::SegmentedMismatches<2> stuff(ptrs, {4, 2});

        {
            auto res = stuff.search("AAAAAC", { 0, 1 });
            EXPECT_EQ(res.index, -1);
        }

        {
            auto res = stuff.search("CAAAAG", { 1, 1 });
            EXPECT_EQ(res.index, -1);
        }
    }
}

TEST(SegmentedMismatches, Optimized) {
    // Deliberately out of order.
    std::vector<std::string> things { "CCCCCC", "AAAAAA", "TTTTTT", "GGGGGG" };
    kaori::BarcodePool ptrs(things);

    for (size_t i = 0; i < 2; ++i) {
      kaori::SegmentedMismatches<2> stuff(ptrs, {3, 3});
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
            EXPECT_EQ(res.index, -1);
            EXPECT_EQ(res.total, 3);
        }
    }
}
