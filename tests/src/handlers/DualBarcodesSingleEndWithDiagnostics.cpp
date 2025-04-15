#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodesSingleEndWithDiagnostics.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesSingleEndWithDiagnosticsTest : public testing::Test {
protected:
    DualBarcodesSingleEndWithDiagnosticsTest() : 
        constant("AAAA----CGGCAGCT------TTTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::string constant1, constant;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;

    template<size_t max_size>
    using Options = typename kaori::DualBarcodesSingleEnd<max_size>::Options;
};

TEST_F(DualBarcodesSingleEndWithDiagnosticsTest, BasicFirst) {
    std::vector<std::string> seq { 
        "cagcatcgatcgtgaAAAACCCCCGGCAGCTTGTGTGTTTTacggaggaga",  // index 1
        "AAAAGGGGCGGCAGCTAGAGAGTTTTaaaaccccggg", // index 2
        "AAAAGGGGCGGCAGCTTGTGTGTTTTaaaaccccggg", // invalid: (2, 1)
        "aaccaccaAAAATTTTCGGCAGCTACACACTTTTaaaaccccggg", // invalid (3, 0)
        "cagacgagcagcgagcagcatcagca" // matches nothing.
    };
    std::string fq = convert_to_fastq(seq);

    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());

    kaori::DualBarcodesSingleEndWithDiagnostics<32, 2> stuff(
        constant.c_str(), constant.size(), 
        { kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        Options<32>()
    );
    kaori::process_single_end_data(&reader, stuff, {});

    EXPECT_EQ(stuff.get_total(), 5);
    EXPECT_EQ(stuff.get_counts()[0], 0);
    EXPECT_EQ(stuff.get_counts()[1], 1);
    EXPECT_EQ(stuff.get_counts()[2], 1);
    EXPECT_EQ(stuff.get_counts()[3], 0);

    stuff.sort();
    const auto& combos = stuff.get_combinations();
    ASSERT_EQ(combos.size(), 2);
    EXPECT_EQ(combos[0][0], 2);
    EXPECT_EQ(combos[0][1], 1);

    EXPECT_EQ(combos[1][0], 3);
    EXPECT_EQ(combos[1][1], 0);
}

TEST_F(DualBarcodesSingleEndWithDiagnosticsTest, WithDuplicates) {
    // Inserting duplicate entries, even though the combinations are unique.
    variables1.push_back("AAAA");
    EXPECT_EQ(variables1.front(), variables1.back());
    variables2.push_back("CTCTCT");
    EXPECT_EQ(variables2[3], variables2[4]);

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaAAAAAAAACGGCAGCTACACACTTTTcagcatcgatcgtga", // ok, index 1
        "cagcatcgatcgtgaAAAAAAAACGGCAGCTCTCTCTTTTTcagcatcgatcgtga", // ok, index 4 
        "cagcatcgatcgtgaAAAAAAAACGGCAGCTTGTGTGTTTTcagcatcgatcgtga"  // invalid (0, 1), as the first hit is reported.
    };
    std::string fq = convert_to_fastq(seq);

    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());

    kaori::DualBarcodesSingleEndWithDiagnostics<32, 2> stuff(
        constant.c_str(), constant.size(), 
        { kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        Options<32>()
    );
    kaori::process_single_end_data(&reader, stuff, {});

    EXPECT_EQ(stuff.get_total(), 3);
    EXPECT_EQ(stuff.get_counts()[0], 1);
    EXPECT_EQ(stuff.get_counts()[1], 0);
    EXPECT_EQ(stuff.get_counts()[2], 0);
    EXPECT_EQ(stuff.get_counts()[3], 0);
    EXPECT_EQ(stuff.get_counts()[4], 1);

    stuff.sort();
    const auto& combos = stuff.get_combinations();
    ASSERT_EQ(combos.size(), 1);
    EXPECT_EQ(combos.front()[0], 0);
    EXPECT_EQ(combos.front()[1], 1);
}
