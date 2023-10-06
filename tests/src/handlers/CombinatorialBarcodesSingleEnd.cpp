#include <gtest/gtest.h>
#include "kaori/handlers/CombinatorialBarcodesSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class CombinatorialBarcodesSingleEndTest : public testing::Test {
protected:
    template<size_t max_size, size_t num_variable>
    using Options = typename kaori::CombinatorialBarcodesSingleEnd<max_size, num_variable>::Options;

    CombinatorialBarcodesSingleEndTest() : 
        constant("AAAA----CGGC------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::array<kaori::BarcodePool, 2> make_pointers() const {
        return std::array<kaori::BarcodePool, 2>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) };
    }

    std::string constant;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(CombinatorialBarcodesSingleEndTest, BasicFirst) {
    kaori::CombinatorialBarcodesSingleEnd<128, 2> stuff(constant.c_str(), constant.size(), make_pointers(), Options<128, 2>());

    // Perfect match.
    {
        std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacac";
        auto state = stuff.initialize();
        stuff.process(state, bounds(seq));
        stuff.reduce(state);

        const auto& combos = stuff.get_combinations();
        ASSERT_EQ(combos.size(), 1);
        EXPECT_EQ(combos.front()[0], 0);
        EXPECT_EQ(combos.front()[1], 1);
    }

    // Works when wrapped inside a process_single_end_data.
    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaAAAACCCCCGGCAGAGAGTTTTacggaggaga", 
        "AAAAGGGGCGGCCTCTCTTTTTacg", 
        "aaaaaAAAATTTTCGGCACACACTTTT" 
    };
    std::string fq = convert_to_fastq(seq);

    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, stuff);

        const auto& combos = stuff.get_combinations();
        ASSERT_EQ(combos.size(), 4); // from the previous invokation of process().
        EXPECT_EQ(stuff.get_total(), 4);

        EXPECT_EQ(combos[1][0], 1);
        EXPECT_EQ(combos[1][1], 2);
        EXPECT_EQ(combos[2][0], 2);
        EXPECT_EQ(combos[2][1], 3);
        EXPECT_EQ(combos[3][0], 3);
        EXPECT_EQ(combos[3][1], 0);
    }
}

TEST_F(CombinatorialBarcodesSingleEndTest, ReverseComplementFirst) {
    auto ptrs = make_pointers();
    kaori::CombinatorialBarcodesSingleEnd<128, 2> forward(constant.c_str(), constant.size(), ptrs, Options<128, 2>());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> reverse(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        return opt;
    }());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> both(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.strand = kaori::SearchStrand::BOTH;
        return opt;
    }());

    std::vector<std::string> seq{ 
        "AAAACACACAGCCGCCCCTTTTccccc", // (GGGG = 2, TGTGTG = 1), and then reverse complemented.
        "aaaaaAAAATTTTCGGCACACACTTTT"  // (3, 0)
    };
    std::string fq = convert_to_fastq(seq);

    // Forward only
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, forward);
        EXPECT_EQ(forward.get_total(), 2);

        const auto& combos = forward.get_combinations();
        ASSERT_EQ(combos.size(), 1);
        EXPECT_EQ(combos.front()[0], 3);
        EXPECT_EQ(combos.front()[1], 0);
    }

    // Reverse only
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, reverse);
        EXPECT_EQ(reverse.get_total(), 2);

        const auto& combos = reverse.get_combinations();
        ASSERT_EQ(combos.size(), 1);
        EXPECT_EQ(combos.front()[0], 2);
        EXPECT_EQ(combos.front()[1], 1);
    }

    // Both.
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, both);
        const auto& combos = both.get_combinations();
        ASSERT_EQ(combos.size(), 2);
    }
}

TEST_F(CombinatorialBarcodesSingleEndTest, MismatchFirst) {
    auto ptrs = make_pointers();
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm0(constant.c_str(), constant.size(), ptrs, Options<128, 2>());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm1(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.max_mismatches = 1;
        return opt;
    }());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm2(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.max_mismatches = 2;
        return opt;
    }());

    // One mismatch.
    {
        std::string seq = "cagAAAAATAACGGCTGTGTGTTTTacac";

        auto state0 = mm0.initialize();
        mm0.process(state0, bounds(seq));
        EXPECT_EQ(state0.collected.size(), 0);

        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        ASSERT_EQ(state1.collected.size(), 1);
        EXPECT_EQ(state1.collected.front()[0], 0);
        EXPECT_EQ(state1.collected.front()[1], 1);
    }

    // Mismatches spread across two variable regions.
    {
        std::string seq = "cagAAAAATAACGGCTGTCTGTTTTacac";

        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        EXPECT_EQ(state1.collected.size(), 0);

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        ASSERT_EQ(state2.collected.size(), 1);
        EXPECT_EQ(state2.collected.front()[0], 0);
        EXPECT_EQ(state2.collected.front()[1], 1);
    }

    // Two mismatches in one variable region.
    {
        std::string seq = "cagAAAATCAACGGCTGTGTGTTTTacac";

        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        EXPECT_EQ(state1.collected.size(), 0);

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        ASSERT_EQ(state2.collected.size(), 1);
        EXPECT_EQ(state2.collected.front()[0], 0);
        EXPECT_EQ(state2.collected.front()[1], 1);
    }

    // More than two msimatches in the variable regions.
    {
        std::string seq = "cagAAAATTAACGGCTGTCTGTTTTacac";

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        EXPECT_EQ(state2.collected.size(), 0);
    }

    // Mismatches in constant and variable regions.
    {
        std::string seq = "cagAATAATAACGGCTGTGTGTTTTacac";

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        ASSERT_EQ(state2.collected.size(), 1);
        EXPECT_EQ(state2.collected.front()[0], 0);
        EXPECT_EQ(state2.collected.front()[1], 1);
    }

    {
        std::string seq = "cagAATAATAACGGCTGTCTGTTTTacac";

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        EXPECT_EQ(state2.collected.size(), 0);
    }

    // Rejects on ambiguities for any variable region.
    {
        std::string seq = "cagAAAATTAACGGCTGTGTGTTTTacac";

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        EXPECT_EQ(state2.collected.size(), 0);
    }
}

TEST_F(CombinatorialBarcodesSingleEndTest, ReverseComplementMismatchFirst) {
    auto ptrs = make_pointers();
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm0(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        return opt;
    }());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm1(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        opt.max_mismatches = 1;
        return opt;
    }());
    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm2(constant.c_str(), constant.size(), ptrs, [&]{
        Options<128, 2> opt;
        opt.strand = kaori::SearchStrand::REVERSE;
        opt.max_mismatches = 2;
        return opt;
    }());

    // One mismatch.
    {
        std::string seq = "cagAAAAAGAGAGGCCGTCTTTTTT"; // (AAAA = 0, CTCTCT = 3)

        auto state0 = mm0.initialize();
        mm0.process(state0, bounds(seq));
        EXPECT_EQ(state0.collected.size(), 0);

        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        ASSERT_EQ(state1.collected.size(), 1);
        EXPECT_EQ(state1.collected.front()[0], 0);
        EXPECT_EQ(state1.collected.front()[1], 3);
    }

    // Mismatches spread across two variable regions.
    {
        std::string seq = "AAAAATAGAGGCCGTCTTTTTTccac"; 

        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        EXPECT_EQ(state1.collected.size(), 0);

        auto state2 = mm2.initialize();
        mm2.process(state2, bounds(seq));
        ASSERT_EQ(state2.collected.size(), 1);
        EXPECT_EQ(state2.collected.front()[0], 0);
        EXPECT_EQ(state2.collected.front()[1], 3);
    }
}

TEST_F(CombinatorialBarcodesSingleEndTest, Best) {
    kaori::CombinatorialBarcodesSingleEnd<128, 2> bst(constant.c_str(), constant.size(), make_pointers(), [&]{
        Options<128, 2> opt;
        opt.use_first = false;
        opt.max_mismatches = 1;
        return opt;
    }());

    kaori::CombinatorialBarcodesSingleEnd<128, 2> mm1(constant.c_str(), constant.size(), make_pointers(), [&]{
        Options<128, 2> opt;
        opt.max_mismatches = 1;
        return opt;
    }());

    // Overrides the first mismatch.
    {
        std::string seq = "cagAAAATAAACGGCTGTGTGTTTTacacAAAATTTTCGGCCTCTCTTTTT";

        auto state = bst.initialize();
        bst.process(state, bounds(seq));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 3);
        EXPECT_EQ(state.collected.front()[1], 3);
        
        auto state1 = mm1.initialize();
        mm1.process(state1, bounds(seq));
        ASSERT_EQ(state1.collected.size(), 1);
        EXPECT_EQ(state1.collected.front()[0], 0);
        EXPECT_EQ(state1.collected.front()[1], 1);
    }

    // Two perfect matches => ambiguous.
    {
        {
            std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacacAAAATTTTCGGCCTCTCTTTTT";

            auto state = bst.initialize();
            bst.process(state, bounds(seq));
            EXPECT_EQ(state.collected.size(), 0);
        }

        // Unless they're the same match!
        {
            std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacacAAAAAAAACGGCTGTGTGTTTT";

            auto state = bst.initialize();
            bst.process(state, bounds(seq));
            ASSERT_EQ(state.collected.size(), 1);
            EXPECT_EQ(state.collected.front()[0], 0);
            EXPECT_EQ(state.collected.front()[1], 1);
        }
    }

    // Works across strands.
    {
        kaori::CombinatorialBarcodesSingleEnd<128, 2> reverse(constant.c_str(), constant.size(), make_pointers(), [&]{
            Options<128, 2> opt;
            opt.use_first = false;
            opt.strand = kaori::SearchStrand::REVERSE;
            return opt;
        }());

        kaori::CombinatorialBarcodesSingleEnd<128, 2> both(constant.c_str(), constant.size(), make_pointers(), [&]{
            Options<128, 2> opt;
            opt.use_first = false;
            opt.strand = kaori::SearchStrand::BOTH;
            return opt;
        }());

        std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacacAAAAGTGTGTGCCGGGGGTTTT"; // second one is (CCCC = 1, ACACAC = 0)

        auto state = bst.initialize(); // only considers the forward strand.
        bst.process(state, bounds(seq));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 0);
        EXPECT_EQ(state.collected.front()[1], 1);

        auto rstate = reverse.initialize(); // only considers the reverse strand.
        reverse.process(rstate, bounds(seq));
        ASSERT_EQ(rstate.collected.size(), 1);
        EXPECT_EQ(rstate.collected.front()[0], 1);
        EXPECT_EQ(rstate.collected.front()[1], 0);

        auto bstate = both.initialize(); // ambiguous on both strands.
        both.process(bstate, bounds(seq));
        EXPECT_EQ(bstate.collected.size(), 0);
    }
}

TEST_F(CombinatorialBarcodesSingleEndTest, Sorting) {
    kaori::CombinatorialBarcodesSingleEnd<128, 2> x(constant.c_str(), constant.size(), make_pointers(), Options<128, 2>());

    auto state = x.initialize(); 
    state.collected.push_back(std::array<int, 2>{ 3, 1 });
    state.collected.push_back(std::array<int, 2>{ 1, 3 });
    state.collected.push_back(std::array<int, 2>{ 2, 3 });
    state.collected.push_back(std::array<int, 2>{ 3, 2 });
    state.collected.push_back(std::array<int, 2>{ 3, 1 });
    state.collected.push_back(std::array<int, 2>{ 0, 2 });
    state.collected.push_back(std::array<int, 2>{ 1, 3 });

    auto copy = state.collected;
    std::sort(copy.begin(), copy.end());
    x.reduce(state);
    x.sort();

    EXPECT_EQ(x.get_combinations().front()[0], 0);
    EXPECT_EQ(x.get_combinations().front()[1], 2);
    EXPECT_EQ(x.get_combinations().back()[0], 3);
    EXPECT_EQ(x.get_combinations().back()[1], 2);

    EXPECT_EQ(x.get_combinations(), copy);
}

TEST_F(CombinatorialBarcodesSingleEndTest, Error) {
    typedef kaori::CombinatorialBarcodesSingleEnd<128, 2> Thing; // Type causes problems inside the macro... who knows what's going on.

    EXPECT_ANY_THROW({
        try {
            std::string constant2 = "AAAA----CGGC----TTTT";
            Thing stuff(constant2.c_str(), constant2.size(), make_pointers(), Options<128, 2>());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("should be the same") != std::string::npos);
            throw e;
        }
    });

    EXPECT_ANY_THROW({
        try {
            std::string constant2 = "AAAA----CGGC----TTTT";
            Thing stuff(constant2.c_str(), constant2.size(), std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables2) }, Options<128, 2>());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("length of 'barcode_pools'") != std::string::npos);
            throw e;
        }
    });

    EXPECT_ANY_THROW({
        try {
            std::string constant2 = "AAAA----CGGCTTTT";
            Thing stuff(constant2.c_str(), constant2.size(), make_pointers(), Options<128, 2>());
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("expected 2 variable regions") != std::string::npos);
            throw e;
        }
    });
}
