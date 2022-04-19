#include <gtest/gtest.h>
#include "kaori/handlers/CombinatorialBarcodes.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class CombinatorialBarcodesTest : public testing::Test {
protected:
    CombinatorialBarcodesTest() : 
        constant("AAAA----CGGC------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::vector<std::vector<const char*> > make_pointers() const {
        return std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) };
    }

    std::pair<const char*, const char*> bounds(const std::string& s) {
        return std::make_pair(s.c_str(), s.c_str() + s.size());
    }

    std::string constant;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(CombinatorialBarcodesTest, BasicFirst) {
    kaori::CombinatorialBarcodes<128, 2> stuff(constant.c_str(), constant.size(), 0, make_pointers());

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

TEST_F(CombinatorialBarcodesTest, ReverseComplementFirst) {
    auto ptrs = make_pointers();
    kaori::CombinatorialBarcodes<128, 2> forward(constant.c_str(), constant.size(), 0, ptrs);
    kaori::CombinatorialBarcodes<128, 2> reverse(constant.c_str(), constant.size(), 1, ptrs);
    kaori::CombinatorialBarcodes<128, 2> both(constant.c_str(), constant.size(), 2, ptrs);

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

TEST_F(CombinatorialBarcodesTest, MismatchFirst) {
    kaori::CombinatorialBarcodes<128, 2> mm0(constant.c_str(), constant.size(), 0, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 0);

    kaori::CombinatorialBarcodes<128, 2> mm1(constant.c_str(), constant.size(), 0, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 1);

    kaori::CombinatorialBarcodes<128, 2> mm2(constant.c_str(), constant.size(), 0, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 2);

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

TEST_F(CombinatorialBarcodesTest, ReverseComplementMismatchFirst) {
    kaori::CombinatorialBarcodes<128, 2> mm0(constant.c_str(), constant.size(), 1, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 0);

    kaori::CombinatorialBarcodes<128, 2> mm1(constant.c_str(), constant.size(), 1, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 1);

    kaori::CombinatorialBarcodes<128, 2> mm2(constant.c_str(), constant.size(), 1, 
        std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) }, 2);

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

TEST_F(CombinatorialBarcodesTest, Best) {
    kaori::CombinatorialBarcodes<128, 2> bst(constant.c_str(), constant.size(), 0, make_pointers());
    bst.set_first(false);
    
    kaori::CombinatorialBarcodes<128, 2> mm1(constant.c_str(), constant.size(), 0, make_pointers(), 1);

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
        std::string seq = "cagAAAAAAAACGGCTGTGTGTTTTacacAAAATTTTCGGCCTCTCTTTTT";

        auto state = bst.initialize();
        bst.process(state, bounds(seq));
        EXPECT_EQ(state.collected.size(), 0);
    }

    // Works across strands.
    {
        kaori::CombinatorialBarcodes<128, 2> reverse(constant.c_str(), constant.size(), 1, make_pointers());
        reverse.set_first(false);

        kaori::CombinatorialBarcodes<128, 2> both(constant.c_str(), constant.size(), 2, make_pointers());
        both.set_first(false);

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

TEST_F(CombinatorialBarcodesTest, Sorting) {
    kaori::CombinatorialBarcodes<128, 2> x(constant.c_str(), constant.size(), 0, make_pointers());

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
