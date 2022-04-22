#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodesWithDiagnostics.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesWithDiagnosticsTest : public testing::Test {
protected:
    DualBarcodesWithDiagnosticsTest() : 
        constant1("AAAA----CGGC"),
        constant2("AGCT------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::string constant1, constant2;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(DualBarcodesWithDiagnosticsTest, BasicFirst) {
    std::vector<std::string> seq1{ 
        "cagcatcgatcgtgaAAAACCCCCGGCacggaggaga",  // index 1
        "AAAAGGGGCGGCaaaaccccggg", // index 2
        "AAAAGGGGCGGCaaaaccccggg", // invalid 
        "AAAAGGGGCGGCaaaaccccggg", // read 1 only
        "acacacacacac" // read 2 only 
    };
    std::string fq1 = convert_to_fastq(seq1);

    std::vector<std::string> seq2{ 
        "cagcatcgatcgtgaAGCTTGTGTGTTTT", 
        "AGCTAGAGAGTTTTaaaaccccggg",
        "cagcatcgatcgtgaAGCTTGTGTGTTTT",
        "acacacacacaca",
        "cagcatcgatcgtgaAGCTTGTGTGTTTT"
    };
    std::string fq2 = convert_to_fastq(seq2);

    byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
    byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());

    kaori::DualBarcodesWithDiagnostics<128> stuff(
        constant1.c_str(), constant1.size(), false, to_pointers(variables1), 0,
        constant2.c_str(), constant2.size(), false, to_pointers(variables2), 0
    );
    kaori::process_paired_end_data(&reader1, &reader2, stuff);

    EXPECT_EQ(stuff.get_total(), 5);
    EXPECT_EQ(stuff.get_counts()[0], 0);
    EXPECT_EQ(stuff.get_counts()[1], 1);
    EXPECT_EQ(stuff.get_counts()[2], 1);
    EXPECT_EQ(stuff.get_counts()[3], 0);

    stuff.sort();
    const auto& combos = stuff.get_combinations();
    ASSERT_EQ(combos.size(), 1);
    EXPECT_EQ(combos.front()[0], 2);
    EXPECT_EQ(combos.front()[1], 1);

    EXPECT_EQ(stuff.get_read1_only(), 1);
    EXPECT_EQ(stuff.get_read2_only(), 1);
}
