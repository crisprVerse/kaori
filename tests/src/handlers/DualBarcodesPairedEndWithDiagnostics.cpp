#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodesPairedEndWithDiagnostics.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesPairedEndWithDiagnosticsTest : public testing::Test {
protected:
    DualBarcodesPairedEndWithDiagnosticsTest() : 
        constant1("AAAA----CGGC"),
        constant2("AGCT------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::string constant1, constant2;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;

    template<size_t max_size>
    using Options = typename kaori::DualBarcodesPairedEnd<max_size>::Options;
};

TEST_F(DualBarcodesPairedEndWithDiagnosticsTest, BasicFirst) {
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

    kaori::DualBarcodesPairedEndWithDiagnostics<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );
    kaori::process_paired_end_data(&reader1, &reader2, stuff, {});

    EXPECT_EQ(stuff.get_total(), 5);
    EXPECT_EQ(stuff.get_counts()[0], 0);
    EXPECT_EQ(stuff.get_counts()[1], 1);
    EXPECT_EQ(stuff.get_counts()[2], 1);
    EXPECT_EQ(stuff.get_counts()[3], 0);

    auto combos = flatten_results<2>(stuff.get_combinations());
    ASSERT_EQ(combos.size(), 1);
    EXPECT_EQ(combos[0].first[0], 2);
    EXPECT_EQ(combos[0].first[1], 1);
    EXPECT_EQ(combos[0].second, 1);

    EXPECT_EQ(stuff.get_barcode1_only(), 1);
    EXPECT_EQ(stuff.get_barcode2_only(), 1);
}

TEST_F(DualBarcodesPairedEndWithDiagnosticsTest, WithDuplicates) {
    // Inserting duplicate entries, even though the combinations are unique.
    variables1.push_back("AAAA");
    EXPECT_EQ(variables1.front(), variables1.back());
    variables2.push_back("CTCTCT");
    EXPECT_EQ(variables2[3], variables2[4]);

    std::vector<std::string> seq1{ 
        "cagcatcgatcgtgaAAAAAAAACGGCacggaggaga",  
        "cagcatcgatcgtgaAAAAAAAACGGCacggaggaga",  
        "cagcatcgatcgtgaAAAAAAAACGGCacggaggaga",  
    };
    std::string fq1 = convert_to_fastq(seq1);

    std::vector<std::string> seq2{ 
        "cagcatcgatcgtgaAGCTACACACTTTT", 
        "cagcatcgatcgtgaAGCTCTCTCTTTTT", 
        "cagcatcgatcgtgaAGCTTGTGTGTTTT", 
    };
    std::string fq2 = convert_to_fastq(seq2);

    byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
    byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());

    kaori::DualBarcodesPairedEndWithDiagnostics<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );
    kaori::process_paired_end_data(&reader1, &reader2, stuff, {});

    EXPECT_EQ(stuff.get_total(), 3);
    EXPECT_EQ(stuff.get_counts()[0], 1);
    EXPECT_EQ(stuff.get_counts()[1], 0);
    EXPECT_EQ(stuff.get_counts()[2], 0);
    EXPECT_EQ(stuff.get_counts()[3], 0);
    EXPECT_EQ(stuff.get_counts()[4], 1);

    auto combos = flatten_results<2>(stuff.get_combinations());
    ASSERT_EQ(combos.size(), 1);
    EXPECT_EQ(combos[0].first[0], 0);
    EXPECT_EQ(combos[0].first[1], 1);
    EXPECT_EQ(combos[0].second, 1);
}
