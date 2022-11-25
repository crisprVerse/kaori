#include <gtest/gtest.h>
#include "kaori/BarcodeSearch.hpp"
#include <string>
#include <vector>
#include "utils.h"

TEST(SimpleBarcodeSearch, Basic) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs);
    auto init = stuff.initialize();

    stuff.search("AAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.search("TTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.search("CCAC", init);
    EXPECT_EQ(init.index, -1);
}

TEST(SimpleBarcodeSearch, ReverseComplement) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs, 0, true);
    auto init = stuff.initialize();

    stuff.search("AAAA", init);
    EXPECT_EQ(init.index, 3);
    stuff.search("TTTT", init);
    EXPECT_EQ(init.index, 0);
}

TEST(SimpleBarcodeSearch, Mismatches) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, 1);
        auto init = stuff.initialize();

        stuff.search("AAAA", init);
        EXPECT_EQ(init.index, 0);
        EXPECT_EQ(init.mismatches, 0);
        
        stuff.search("TTTT", init);
        EXPECT_EQ(init.index, 3);
        EXPECT_EQ(init.mismatches, 0);

        stuff.search("CCAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 1);
        stuff.search("CGAC", init);
        EXPECT_EQ(init.index, -1);
    }

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, 2);
        auto init = stuff.initialize();

        stuff.search("CCAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 1);

        stuff.search("CGAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 2);

        stuff.search("CGGC", init); // ambiguous.
        EXPECT_EQ(init.index, -1);
    }
}

TEST(SimpleBarcodeSearch, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs, 1);

    auto state = stuff.initialize();

    // No cache when there is no mismatch.
    {
        stuff.search("AAAA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 0);
        auto it = state.cache.find("AAAA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // No cache when the number of mismatches is lower than that in the constructor.
    {
        stuff.search("AATA", state, 0);
        EXPECT_EQ(state.index, -1);
        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // Stored in cache for >1 mismatches.
    {
        stuff.search("AATA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 1);

        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).first, 0);
        EXPECT_EQ((it->second).second, 1);
    }

    {
        stuff.search("ACTA", state);
        EXPECT_EQ(state.index, -1);

        auto it = state.cache.find("ACTA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).first, -1);
    }

    // Retrieval from cache respects a lower mismatch threshold.  This uses the
    // same barcode as above, which should be cached with a match that is
    // revoked at a lower mismatch threshold.
    {
        stuff.search("AATA", state, 0); 
        EXPECT_EQ(state.index, -1);
    }
 
    // Checking that the reduction works correctly.
    state.cache["AATA"].first = 2;
    stuff.reduce(state);
    EXPECT_TRUE(state.cache.empty());

    {
        stuff.search("AATA", state);
        EXPECT_EQ(state.index, 2); // re-uses the cache value!
    }
}

TEST(SimpleBarcodeSearch, Duplicates) {
    std::vector<std::string> things { "ACGT", "ACGT", "AGTT", "AGTT" };
    kaori::BarcodePool ptrs(things);

    EXPECT_ANY_THROW({
        try {
            kaori::SimpleBarcodeSearch stuff(ptrs, 0, false, false);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("duplicate") != std::string::npos);
            throw e;
        }
    });

    // Gets the first occurrence.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, 0, false, true);
        auto state = stuff.initialize();

        stuff.search("ACGT", state);
        EXPECT_EQ(state.index, 0);

        stuff.search("AGTT", state);
        EXPECT_EQ(state.index, 2);
    }

    // ... even with a mismatch.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, 1, false, true);
        auto state = stuff.initialize();

        stuff.search("ACGA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 1);

        stuff.search("AGTA", state);
        EXPECT_EQ(state.index, 2);
        EXPECT_EQ(state.mismatches, 1);
    }
}

TEST(SegmentedBarcodeSearch, Basic) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    kaori::BarcodePool ptrs(variables);
 
    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, { 2, 4 }, { 0, 1 });
    auto init = stuff.initialize();

    stuff.search("AAAAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.search("AATTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.search("AACCAC", init); // 1 mismatch
    EXPECT_EQ(init.index, 1);

    stuff.search("AAccgg", init); // ambiguous.
    EXPECT_EQ(init.index, -1);
}

TEST(SegmentedBarcodeSearch, ReverseComplement) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, { 2, 4 }, { 0, 1 }, true);
    auto init = stuff.initialize();

    stuff.search("AAAATT", init);
    EXPECT_EQ(init.index, 3);
    stuff.search("GGGGTT", init);
    EXPECT_EQ(init.index, 1);

    stuff.search("CCACTT", init);
    EXPECT_EQ(init.index, 2);

    stuff.search("GGCCTT", init); // ambiguous.
    EXPECT_EQ(init.index, -1);
}

TEST(SegmentedBarcodeSearch, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, {2, 2}, {1, 1});

    auto state = stuff.initialize();

    // No cache when there is no mismatch.
    {
        stuff.search("AAAA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 0);
        auto it = state.cache.find("AAAA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // No cache when the number of mismatches is less than that in the constructor.
    {
        stuff.search("AATA", state, { 0, 0 });
        EXPECT_EQ(state.index, -1);
        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it == state.cache.end());
    }

    // Stored in cache for >1 mismatches.
    {
        stuff.search("AATA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 1);
        auto it = state.cache.find("AATA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, 0);
        EXPECT_EQ((it->second).total, 1);
    }

    {
        stuff.search("ACCA", state);
        EXPECT_EQ(state.index, -1);
        auto it = state.cache.find("ACCA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, -1);
    }

    // Retrieval from cache respects a lower mismatch threshold.  This uses the
    // same barcode as above, which should be cached with a match that is
    // revoked at a lower mismatch threshold.
    {
        stuff.search("AATA", state, { 0, 1 }); //ok
        EXPECT_EQ(state.index, 0);

        stuff.search("AATA", state, { 0, 0 }); //fail
        EXPECT_EQ(state.index, -1);
    }

    // Checking that the reduction works correctly.
    state.cache["AATA"].index = 2;
    stuff.reduce(state);
    EXPECT_TRUE(state.cache.empty());

    {
        stuff.search("AATA", state);
        EXPECT_EQ(state.index, 2); // re-uses the cache value!
    }
}
