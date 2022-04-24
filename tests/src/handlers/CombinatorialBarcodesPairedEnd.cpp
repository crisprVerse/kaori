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

    std::string constant1, constant2;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;
};

TEST_F(CombinatorialBarcodesPairedEndTest, BasicFirst) {
    kaori::CombinatorialBarcodesPairedEnd<128> stuff(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
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

        stuff.sort(); // for some coverage.

        EXPECT_EQ(stuff.get_total(), 2);
        EXPECT_EQ(stuff.get_barcode1_only(), 0);
        EXPECT_EQ(stuff.get_barcode2_only(), 0);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, ReverseComplementFirst) {
    {
        kaori::CombinatorialBarcodesPairedEnd<128> stuff(
            constant1.c_str(), constant1.size(), true, kaori::BarcodePool(variables1), 0,
            constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
        );

        auto state = stuff.initialize();
        std::string seq1 = "atacacacaGCCGAAAATTTT", seq2 = "AGCTACACACTTTTtgaaca";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 3);
        EXPECT_EQ(state.collected.front()[1], 0);
    }

    // Just some due diligence here...
    {
        kaori::CombinatorialBarcodesPairedEnd<128> stuff(
            constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
            constant2.c_str(), constant2.size(), true, kaori::BarcodePool(variables2), 0
        );

        auto state = stuff.initialize();
        std::string seq1 = "AAAACCCCCGGCacacacac", seq2 = "gacgacgaAAAAAGAGAGAGCT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 1);
        EXPECT_EQ(state.collected.front()[0], 1);
        EXPECT_EQ(state.collected.front()[1], 3);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, MismatchesFirst) {
    kaori::CombinatorialBarcodesPairedEnd<128> stuff(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );

    kaori::CombinatorialBarcodesPairedEnd<128> stuff10(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );

    kaori::CombinatorialBarcodesPairedEnd<128> stuff20(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 2,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );

    // Works in the simple case.
    {
        std::string seq1 = "AAAATTATCGGC", seq2 = "AGCTACACACTTTT";

        auto state = stuff.initialize();
        stuff.process(state, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state.collected.size(), 0);

        auto state10 = stuff10.initialize();
        stuff10.process(state10, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state10.collected.size(), 1);
        EXPECT_EQ(state10.collected.front()[0], 3);
        EXPECT_EQ(state10.collected.front()[1], 0);
    }

    // Handles ambiguity.
    {
        std::string seq1 = "AAAATAATCGGC", seq2 = "AGCTACACACTTTT";

        auto state20 = stuff20.initialize();
        stuff.process(state20, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state20.collected.size(), 0);

        seq1[5] = 'C';
        stuff20.process(state20, bounds(seq1), bounds(seq2));
        ASSERT_EQ(state20.collected.size(), 1);
        EXPECT_EQ(state20.collected.front()[0], 3);
        EXPECT_EQ(state20.collected.front()[1], 0);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, RandomizedFirst) {
    kaori::CombinatorialBarcodesPairedEnd<128> nonrandom(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );

    kaori::CombinatorialBarcodesPairedEnd<128> randomized(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0,
        true
    );

    std::string seq1 = "ctagcgaAGCTTGTGTGTTTTagaga", seq2 = "acacAAAAGGGGCGGCacac";

    {
        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nrstate.collected.size(), 0);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(rstate.collected.size(), 1);
        EXPECT_EQ(rstate.collected.front()[0], 2);
        EXPECT_EQ(rstate.collected.front()[1], 1);
    }

    // Positive control, check that the sequences were right.
    {
        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq2), bounds(seq1));
        ASSERT_EQ(nrstate.collected.size(), 1);
        EXPECT_EQ(nrstate.collected.front()[0], 2);
        EXPECT_EQ(nrstate.collected.front()[1], 1);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq2), bounds(seq1));
        ASSERT_EQ(rstate.collected.size(), 1);
        EXPECT_EQ(rstate.collected.front()[0], 2);
        EXPECT_EQ(rstate.collected.front()[1], 1);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, DiagnosticsFirst) {
    kaori::CombinatorialBarcodesPairedEnd<128> nonrandom(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );

    kaori::CombinatorialBarcodesPairedEnd<128> randomized(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0,
        true
    );

    // Only read 1.
    {
        std::string seq1 = "AAAATTTTCGGC", seq2 = "acacaaacca";

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nrstate.barcode1_only, 1);
        EXPECT_EQ(nrstate.barcode2_only, 0);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.barcode1_only, 1);
        EXPECT_EQ(rstate.barcode2_only, 0);
    }

    // Only read 2.
    {
        std::string seq1 = "acacacacacacca", seq2 = "ctagcgaAGCTTGTGTGTTTTagaga";

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nrstate.barcode1_only, 0);
        EXPECT_EQ(nrstate.barcode2_only, 1);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.barcode1_only, 0);
        EXPECT_EQ(rstate.barcode2_only, 1);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, BasicBest) {
    kaori::CombinatorialBarcodesPairedEnd<128> best(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 1
    );
    best.set_first(false);

    kaori::CombinatorialBarcodesPairedEnd<128> first(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 1
    );

    {
        std::string seq1 = "AAAATTATCGGCacacatAAAACCCCCGGC", seq2 = "AGCTCCCTCTTTTTtttttttAGCTAGAGAGTTTT";

        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 1);
        EXPECT_EQ(bstate.collected.front()[0], 1);
        EXPECT_EQ(bstate.collected.front()[1], 2);

        // Positive control.
        auto fstate = first.initialize();
        first.process(fstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(fstate.collected.size(), 1);
        EXPECT_EQ(fstate.collected.front()[0], 3);
        EXPECT_EQ(fstate.collected.front()[1], 3);
    }

    // Recognizes ambiguity.
    {
        std::string seq1 = "AAAATTTTCGGCacacatAAAACCCCCGGC", seq2 = "AGCTCTCTCTTTTTtttttttAGCTAGAGAGTTTT";
        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 0);

        // Positive control.
        auto fstate = first.initialize();
        first.process(fstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(fstate.collected.size(), 1);
        EXPECT_EQ(fstate.collected.front()[0], 3);
        EXPECT_EQ(fstate.collected.front()[1], 3);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, RandomizedBest) {
    kaori::CombinatorialBarcodesPairedEnd<128> best(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 1,
        true
    );
    best.set_first(false);

    kaori::CombinatorialBarcodesPairedEnd<128> nonrandom(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 1
    );
    nonrandom.set_first(false);

    kaori::CombinatorialBarcodesPairedEnd<128> first(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 1,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 1,
        true
    );

    // Favors the other orientation.
    {
        std::string seq1 = "AAAATTATCGGCacacaAGCTCTCTCTTTTT", seq2 = "AGCTACACTCTTTTtttttttAAAAAAAACGGC";

        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 1);
        EXPECT_EQ(bstate.collected.front()[0], 0);
        EXPECT_EQ(bstate.collected.front()[1], 3);

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(nrstate.collected.size(), 1);
        EXPECT_EQ(nrstate.collected.front()[0], 3);
        EXPECT_EQ(nrstate.collected.front()[1], 0);

        // Positive control.
        auto fstate = first.initialize();
        first.process(fstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(fstate.collected.size(), 1);
        EXPECT_EQ(fstate.collected.front()[0], 3);
        EXPECT_EQ(fstate.collected.front()[1], 0);
    }

    // Testing when the original orientation is preferred.
    {
        std::string seq1 = "AAAATTTTCGGCacacaAGCTCGCTCTTTTT", seq2 = "AGCTACACACTTTTtttttttAAAAAAAACGGC";
        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 1);
        EXPECT_EQ(bstate.collected.front()[0], 3);
        EXPECT_EQ(bstate.collected.front()[1], 0);
    }

    // Recognizes ambiguity.
    {
        std::string seq1 = "AAAATTTTCGGCacacaAGCTCTCTCTTTTT", seq2 = "AGCTACACACTTTTtttttttAAAAAAAACGGC";
        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 0);

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(nrstate.collected.size(), 1);
        EXPECT_EQ(nrstate.collected.front()[0], 3);
        EXPECT_EQ(nrstate.collected.front()[1], 0);

        // Positive control.
        auto fstate = first.initialize();
        first.process(fstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(fstate.collected.size(), 1);
        EXPECT_EQ(fstate.collected.front()[0], 3);
        EXPECT_EQ(fstate.collected.front()[1], 0);
    }

    // ...unless it's just a duplicate of the same pair.
    {
        std::string seq1 = "AAAATTTTCGGCacacaAGCTACACACTTTT", seq2 = "AGCTACACACTTTTtttttttAAAATTTTCGGC";
        auto bstate = best.initialize();
        best.process(bstate, bounds(seq1), bounds(seq2));
        ASSERT_EQ(bstate.collected.size(), 1);
        EXPECT_EQ(bstate.collected.front()[0], 3);
        EXPECT_EQ(bstate.collected.front()[1], 0);
    }
}

TEST_F(CombinatorialBarcodesPairedEndTest, DiagnosticsBest) {
    kaori::CombinatorialBarcodesPairedEnd<128> nonrandom(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0
    );
    nonrandom.set_first(false);

    kaori::CombinatorialBarcodesPairedEnd<128> randomized(
        constant1.c_str(), constant1.size(), false, kaori::BarcodePool(variables1), 0,
        constant2.c_str(), constant2.size(), false, kaori::BarcodePool(variables2), 0,
        true
    );
    randomized.set_first(false);

    // Only read 1.
    {
        std::string seq1 = "AAAATTTTCGGC", seq2 = "acacaaacca";

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nrstate.barcode1_only, 1);
        EXPECT_EQ(nrstate.barcode2_only, 0);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.barcode1_only, 1);
        EXPECT_EQ(rstate.barcode2_only, 0);
    }

    // Only read 2.
    {
        std::string seq1 = "acacacacacacca", seq2 = "ctagcgaAGCTTGTGTGTTTTagaga";

        auto nrstate = nonrandom.initialize();
        nonrandom.process(nrstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nrstate.barcode1_only, 0);
        EXPECT_EQ(nrstate.barcode2_only, 1);

        auto rstate = randomized.initialize();
        randomized.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.barcode1_only, 0);
        EXPECT_EQ(rstate.barcode2_only, 1);
    }
}

