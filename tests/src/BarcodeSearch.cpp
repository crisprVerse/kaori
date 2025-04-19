#include <gtest/gtest.h>
#include "kaori/BarcodeSearch.hpp"
#include <string>
#include <vector>
#include "utils.h"

class SimpleBarcodeSearchTest : public ::testing::Test {
protected:
    typedef typename kaori::SimpleBarcodeSearch::Options Options;
};

TEST_F(SimpleBarcodeSearchTest, Basic) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs, Options());
    auto init = stuff.initialize();

    stuff.search("AAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.search("TTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.search("CCAC", init);
    EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);
}

TEST_F(SimpleBarcodeSearchTest, ReverseComplement) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs, [&]{ 
        Options opt;
        opt.reverse = true;
        return opt;
    }());
    auto init = stuff.initialize();

    stuff.search("AAAA", init);
    EXPECT_EQ(init.index, 3);
    stuff.search("TTTT", init);
    EXPECT_EQ(init.index, 0);
}

TEST_F(SimpleBarcodeSearchTest, Iupac) {
    std::vector<std::string> variables { "ARYA", "CSWC", "GKMG", "ABDC", "GHVT", "ANNT" };
    kaori::BarcodePool ptrs(variables);

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, Options());
        auto init = stuff.initialize();

        stuff.search("AATA", init);
        EXPECT_EQ(init.index, 0);
        stuff.search("CGAC", init);
        EXPECT_EQ(init.index, 1);
        stuff.search("GTCG", init);
        EXPECT_EQ(init.index, 2);
        stuff.search("AGGC", init);
        EXPECT_EQ(init.index, 3);
        stuff.search("GCCT", init);
        EXPECT_EQ(init.index, 4);
        stuff.search("ACGT", init);
        EXPECT_EQ(init.index, 5);

        stuff.search("AAAA", init);
        EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);
    }

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.reverse = true;
            return opt;
        }());
        auto init = stuff.initialize();

        stuff.search("TATT", init);
        EXPECT_EQ(init.index, 0);
        stuff.search("GTCG", init);
        EXPECT_EQ(init.index, 1);
        stuff.search("CGAC", init);
        EXPECT_EQ(init.index, 2);
        stuff.search("GCCT", init);
        EXPECT_EQ(init.index, 3);
        stuff.search("AGGC", init);
        EXPECT_EQ(init.index, 4);
        stuff.search("ACGT", init);
        EXPECT_EQ(init.index, 5);

        stuff.search("AATA", init);
        EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);
    }
}

TEST_F(SimpleBarcodeSearchTest, Mismatches) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.max_mismatches = 1;
            return opt;
        }());
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
        EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);
    }

    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.max_mismatches = 2;
            return opt;
        }());
        auto init = stuff.initialize();

        stuff.search("CCAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 1);

        stuff.search("CGAC", init);
        EXPECT_EQ(init.index, 1);
        EXPECT_EQ(init.mismatches, 2);

        stuff.search("CGGC", init); // ambiguous.
        EXPECT_EQ(init.index, kaori::STATUS_AMBIGUOUS);
    }
}

TEST_F(SimpleBarcodeSearchTest, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
        Options opt;
        opt.max_mismatches = 1;
        return opt;
    }());

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
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
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
        EXPECT_EQ((it->second).mismatches, 1);
    }

    {
        stuff.search("ACTA", state);
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);

        auto it = state.cache.find("ACTA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, kaori::STATUS_UNMATCHED);
    }

    // Retrieval from cache respects a lower mismatch threshold.  This uses the
    // same barcode as above, which should be cached with a match that is
    // revoked at a lower mismatch threshold.
    {
        stuff.search("AATA", state, 0); 
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
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

TEST_F(SimpleBarcodeSearchTest, Duplicates) {
    std::vector<std::string> things { "ACGT", "ACGT", "AGTT", "AGTT" };
    kaori::BarcodePool ptrs(things);

    EXPECT_ANY_THROW({
        try {
            kaori::SimpleBarcodeSearch stuff(ptrs, Options());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("duplicate") != std::string::npos);
            throw e;
        }
    });

    // Gets the first occurrence.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.duplicates = kaori::DuplicateAction::FIRST;
            return opt;
        }());
        auto state = stuff.initialize();

        stuff.search("ACGT", state);
        EXPECT_EQ(state.index, 0);

        stuff.search("AGTT", state);
        EXPECT_EQ(state.index, 2);
    }

    // ... even with a mismatch.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.max_mismatches = 1;
            opt.duplicates = kaori::DuplicateAction::FIRST;
            return opt;
        }());
        auto state = stuff.initialize();

        stuff.search("ACGA", state);
        EXPECT_EQ(state.index, 0);
        EXPECT_EQ(state.mismatches, 1);

        stuff.search("AGTA", state);
        EXPECT_EQ(state.index, 2);
        EXPECT_EQ(state.mismatches, 1);
    }

    // Gets the last occurrence.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.duplicates = kaori::DuplicateAction::LAST;
            return opt;
        }());
        auto state = stuff.initialize();

        stuff.search("ACGT", state);
        EXPECT_EQ(state.index, 1);

        stuff.search("AGTT", state);
        EXPECT_EQ(state.index, 3);
    }

    // Gets no occurrence.
    {
        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.duplicates = kaori::DuplicateAction::NONE;
            return opt;
        }());
        auto state = stuff.initialize();

        stuff.search("ACGT", state);
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);

        stuff.search("AGTT", state);
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
    }

    // Continues getting no occurrences, even with > 2 occurrences of the duplicates.
    {
        std::vector<std::string> things { "ACGT", "ACGT", "AGTT", "ACGT", "AGTT", "AGTT", "ACGT" };
        kaori::BarcodePool ptrs(things);

        kaori::SimpleBarcodeSearch stuff(ptrs, [&]{
            Options opt;
            opt.duplicates = kaori::DuplicateAction::NONE;
            return opt;
        }());
        auto state = stuff.initialize();

        stuff.search("ACGT", state);
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);

        stuff.search("AGTT", state);
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
    }
}

class SegmentedBarcodeSearchTest : public ::testing::Test {
protected:
    template<size_t num_segments>
    using Options = typename kaori::SegmentedBarcodeSearch<num_segments>::Options;
};

TEST_F(SegmentedBarcodeSearchTest, Basic) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    kaori::BarcodePool ptrs(variables);
 
    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, { 2, 4 }, [&]{
        Options<2> opt;
        opt.max_mismatches = { 0, 2 };
        return opt;
    }());
    auto init = stuff.initialize();

    stuff.search("AAAAAA", init);
    EXPECT_EQ(init.index, 0);
    stuff.search("AATTTT", init);
    EXPECT_EQ(init.index, 3);

    stuff.search("AACCAC", init); // 1 mismatch
    EXPECT_EQ(init.index, 1);

    stuff.search("ACCCCC", init); // 1 mismatch in the wrong place.
    EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);

    stuff.search("AAccgg", init); // ambiguous.
    EXPECT_EQ(init.index, kaori::STATUS_AMBIGUOUS);
}

TEST_F(SegmentedBarcodeSearchTest, ReverseComplement) {
    std::vector<std::string> variables { "AAAAAA", "AACCCC", "AAGGGG", "AATTTT" };
    kaori::BarcodePool ptrs(variables);

    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, { 2, 4 }, [&]{
        Options<2> opt;
        opt.max_mismatches = { 0, 2 };
        opt.reverse = true;
        return opt;
    }());
    auto init = stuff.initialize();

    stuff.search("AAAATT", init);
    EXPECT_EQ(init.index, 3);
    stuff.search("GGGGTT", init);
    EXPECT_EQ(init.index, 1);

    stuff.search("CCACTT", init);
    EXPECT_EQ(init.index, 2);

    stuff.search("AGGCTT", init); // two mismatches.
    EXPECT_EQ(init.index, 1);

    stuff.search("GGGGGG", init); // two mismatches in the wrong place.
    EXPECT_EQ(init.index, kaori::STATUS_UNMATCHED);

    stuff.search("GGCCTT", init); // ambiguous.
    EXPECT_EQ(init.index, kaori::STATUS_AMBIGUOUS);
}

TEST_F(SegmentedBarcodeSearchTest, Caching) {
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };
    kaori::BarcodePool ptrs(variables);
    kaori::SegmentedBarcodeSearch<2> stuff(ptrs, {2, 2}, [&]{
        Options<2> opt;
        opt.max_mismatches = {1, 1};
        return opt;
    }());
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
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
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
        EXPECT_EQ((it->second).mismatches, 1);
    }

    // Stored in cache if it's ambiguous.
    {
        stuff.search("ACCA", state);
        EXPECT_EQ(state.index, kaori::STATUS_AMBIGUOUS);
        auto it = state.cache.find("ACCA");
        EXPECT_TRUE(it != state.cache.end());
        EXPECT_EQ((it->second).index, kaori::STATUS_AMBIGUOUS);
    }

    // Retrieval from cache respects a lower mismatch threshold.  This uses the
    // same barcode as above, which should be cached with a match that is
    // revoked at a lower mismatch threshold.
    {
        stuff.search("AATA", state, { 0, 1 }); //ok
        EXPECT_EQ(state.index, 0);

        stuff.search("AATA", state, { 0, 0 }); //fail
        EXPECT_EQ(state.index, kaori::STATUS_UNMATCHED);
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
