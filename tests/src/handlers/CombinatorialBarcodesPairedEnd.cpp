#include <gtest/gtest.h>
#include "kaori/handlers/CombinatorialBarcodesPairedEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class CombinatorialBarcodesPairedEndTest : public testing::Test {
protected:
    CombinatorialBarcodesPairedEndTest() : 
        constant1("AAAA----CGGC"),
        constant2("AGCT------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::vector<std::vector<const char*> > make_pointers() const {
        return std::vector<std::vector<const char*> >{ to_pointers(variables1), to_pointers(variables2) };
    }

    std::string constant1, constant2;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(CombinatorialBarcodesPairedEndTest, BasicFirst) {
    kaori::CombinatorialBarcodesPairedEnd<128> stuff(
        constant1.c_str(), constant1.size(), false, to_pointers(variables1), 0,
        constant2.c_str(), constant2.size(), false, to_pointers(variables2), 0
    );

    // Works in the simple case.
    {
        auto state = stuff.initialize();
        std::string seq1 = "AAAATTTTCGGC", seq2 = "AGCTACACACTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 3);
        EXPECT_EQ(state.collected.front()[1], 0);
    }

    // Making it work for it.
    {
        auto state = stuff.initialize();
        std::string seq1 = "cacacacAAAAAAAACGGC", seq2 = "ggggAGCTAGAGAGTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 0);
        EXPECT_EQ(state.collected.front()[1], 2);
    }

    // Integrated.
    {
        std::vector<std::string> seq1{ 
            "cagcatcgatcgtgaAAAACCCCCGGCacggaggaga", 
            "AAAAGGGGCGGCaaaaccccggg"
        };
        std::string fq1 = convert_to_fastq(seq1);

        std::vector<std::string> seq2{ 
            "cagcatcgatcgtgaAGCTTGTGTGTTTT", 
            "AGCTAGAGAGTTTTaaaaccccggg"
        };
        std::string fq2 = convert_to_fastq(seq2);

        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, stuff);

        const auto& out = stuff.get_combinations();
        ASSERT_EQ(out.size(), 2);
        EXPECT_EQ(out[0][0], 1);
        EXPECT_EQ(out[0][1], 1);
        EXPECT_EQ(out[1][0], 2);
        EXPECT_EQ(out[1][1], 2);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, ReverseComplementFirst) {
    {
        kaori::CombinatorialBarcodesPairedEnd<128> stuff(
            constant1.c_str(), constant1.size(), true, to_pointers(variables1), 0,
            constant2.c_str(), constant2.size(), false, to_pointers(variables2), 0
        );

        auto state = stuff.initialize();
        std::string seq1 = "atacacacaGCCGAAAATTTT", seq2 = "AGCTACACACTTTTtgaaca";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 3);
        EXPECT_EQ(state.collected.front()[1], 0);
    }

    {
        kaori::CombinatorialBarcodesPairedEnd<128> stuff(
            constant1.c_str(), constant1.size(), false, to_pointers(variables1), 0,
            constant2.c_str(), constant2.size(), true, to_pointers(variables2), 0
        );

        auto state = stuff.initialize();
        std::string seq1 = "AAAACCCCCGGCacacacac", seq2 = "gacgacgaAAAAAGAGAGAGCT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 1);
        EXPECT_EQ(state.collected.front()[1], 3);
    }
}
