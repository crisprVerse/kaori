#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodesPairedEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesPairedEndTest : public testing::Test {
protected:
    DualBarcodesPairedEndTest() : 
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

TEST_F(DualBarcodesPairedEndTest, BasicFirst) {
    kaori::DualBarcodesPairedEnd<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );

    // Works in the simple case.
    {
        auto state = stuff.initialize();
        std::string seq1 = "AAAATTTTCGGC", seq2 = "AGCTCTCTCTTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 1);
        EXPECT_EQ(state.total, 1);
    }

    // Making it work for it.
    {
        auto state = stuff.initialize();
        std::string seq1 = "cacacacAAAAAAAACGGC", seq2 = "ggggAGCTACACACTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 1);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 0);
        EXPECT_EQ(state.total, 1);
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

        EXPECT_EQ(stuff.get_counts()[0], 0);
        EXPECT_EQ(stuff.get_counts()[1], 1);
        EXPECT_EQ(stuff.get_counts()[2], 1);
        EXPECT_EQ(stuff.get_counts()[3], 0);
        EXPECT_EQ(stuff.get_total(), 2);
    }
}

TEST_F(DualBarcodesPairedEndTest, ReverseComplementFirst) {
    // Works in the simple case.
    {
        kaori::DualBarcodesPairedEnd<32> stuff(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
            [&]{
                Options<32> opt;
                opt.strand1 = kaori::SearchStrand::REVERSE;
                return opt;
            }()
        );

        auto state = stuff.initialize();
        std::string seq1 = "GCCGAAAATTTTaacacac", seq2 = "ccacacAGCTCTCTCTTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 1);
        EXPECT_EQ(state.total, 1);
    }

    // And now the other one.
    {
        kaori::DualBarcodesPairedEnd<32> stuff(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            [&]{
                Options<32> opt;
                opt.strand2 = kaori::SearchStrand::REVERSE;
                return opt;
            }()
        );

        auto state = stuff.initialize();
        std::string seq1 = "cgataAAAAGGGGCGGC", seq2 = "ccacacAAAACTCTCTAGCT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 1);
        EXPECT_EQ(state.counts[3], 0);
        EXPECT_EQ(state.total, 1);
    }
}

TEST_F(DualBarcodesPairedEndTest, Iupac) {
    variables1 = std::vector<std::string>{ "ARRA", "CYYC", "GSSG", "TWWT" };
    variables2 = std::vector<std::string>{ "ABACAC", "TVTGTG", "AGHGAG", "CDCTCT" };

    // Forward.
    {
        kaori::DualBarcodesPairedEnd<32> stuff(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            Options<32>()
        );

        auto state = stuff.initialize();
        std::string seq1 = "AAAATTTTCGGC", seq2 = "AGCTCTCTCTTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 1);
        EXPECT_EQ(state.total, 1);
    }

    // Reverse-complement.
    {
        kaori::DualBarcodesPairedEnd<32> stuff(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            [&]{
                Options<32> opt;
                opt.strand1 = kaori::SearchStrand::REVERSE;
                opt.strand2 = kaori::SearchStrand::REVERSE;
                return opt;
            }()
        );

        auto state = stuff.initialize();
        std::string seq1 = "cgataGCCGCCCCTTTTacttc", seq2 = "ccacacAAAACTCTCTAGCT";
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 1);
        EXPECT_EQ(state.counts[3], 0);
        EXPECT_EQ(state.total, 1);
    }
}


TEST_F(DualBarcodesPairedEndTest, MismatchFirst) {
    kaori::DualBarcodesPairedEnd<32> stuff00(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );

    kaori::DualBarcodesPairedEnd<32> stuff10(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.max_mismatches1 = 1;
            return opt;
        }()
    );

    kaori::DualBarcodesPairedEnd<32> stuff11(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
        [&]{
            Options<32> opt;
            opt.max_mismatches1 = 1;
            opt.max_mismatches2 = 1;
            return opt;
        }()
    );

    kaori::DualBarcodesPairedEnd<32> stuff20(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.max_mismatches1 = 2;
            return opt;
        }()
    );

    // One mismatch.
    {
        std::string seq1 = "AAAATTATCGGC", seq2 = "AGCTCTCTCTTTTT";

        auto state00 = stuff00.initialize();
        stuff00.process(state00, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state00.counts[3], 0);

        auto state10 = stuff10.initialize();
        stuff10.process(state10, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state10.counts[3], 1);
    }

    // One mismatch, each.
    {
        std::string seq1 = "AAAATTATCGGC", seq2 = "AGCTCTCGCTTTTT";

        auto state10 = stuff10.initialize();
        stuff10.process(state10, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state10.counts[3], 0);

        auto state11 = stuff11.initialize();
        stuff11.process(state11, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state11.counts[3], 1);
    }

    // Two mismatches, in constant and variable regions.
    {
        std::string seq1 = "AAAATTATCGGG", seq2 = "AGCTCTCTCTTTTT";

        auto state10 = stuff10.initialize();
        stuff10.process(state10, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state10.counts[3], 0);

        auto state20 = stuff20.initialize();
        stuff20.process(state20, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state20.counts[3], 1);
    }
}

TEST_F(DualBarcodesPairedEndTest, AmbiguityFirst) {
    variables1.push_back("AAAA");
    variables2.push_back("AAAAAG");
    variables1.push_back("AAAA");
    variables2.push_back("AAAAAT");

    kaori::DualBarcodesPairedEnd<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.max_mismatches2 = 1;
            return opt;
        }()
    );

    {
        std::string seq1 = "AAAAAAAACGGC", seq2 = "AGCTAAAAACTTTT";

        auto state = stuff.initialize();
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts, std::vector<int>(variables1.size())); // nothing detected.

        seq2 = "AGCTAAAAATTTTT"; // as a control.
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[5], 1);
    }
}

TEST_F(DualBarcodesPairedEndTest, RandomizedFirst) {
    kaori::DualBarcodesPairedEnd<32> nonrandom(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );

    kaori::DualBarcodesPairedEnd<32> random(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
        [&]{
            Options<32> opt;
            opt.random = true;
            return opt;
        }()
    );

    std::string seq1 = "AAAATTTTCGGCcacacacaAGCTTGTGTGTTTT";
    std::string seq2 = "AGCTCTCTCTTTTTcgtacgacAAAACCCCCGGC";

    auto nstate = nonrandom.initialize();
    nonrandom.process(nstate, bounds(seq1), bounds(seq2));
    EXPECT_EQ(nstate.counts[3], 1);

    // Compromise the first hit, force it to look elsewhere.    
    seq2[0] = 'T';
    auto rstate = random.initialize();
    random.process(rstate, bounds(seq1), bounds(seq2));
    EXPECT_EQ(rstate.counts[1], 1);
    EXPECT_EQ(rstate.counts[3], 0);
}

TEST_F(DualBarcodesPairedEndTest, BasicBest) {
    kaori::DualBarcodesPairedEnd<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
        [&]{
            Options<32> opt;
            opt.max_mismatches1 = 1;
            opt.max_mismatches2 = 1;
            opt.use_first = false;
            return opt;
        }()
    );

    kaori::DualBarcodesPairedEnd<32> fstuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.max_mismatches1 = 1;
            opt.max_mismatches2 = 1;
            return opt;
        }()
    );

    kaori::DualBarcodesPairedEnd<32> fstuff0(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        Options<32>()
    );

    // Keeps searching for the best.
    {
        std::string seq1 = "AAAATTATCGGCcacacacaAAAACCCCCGGC", seq2 = "AGCTCTCTCTTTTTcgtacgactAGCTTGTGTGTTTT";

        auto state = stuff.initialize();
        stuff.process(state, bounds(seq1), bounds(seq2));
        EXPECT_EQ(state.counts[1], 1);

        auto fstate = fstuff.initialize();
        fstuff.process(fstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(fstate.counts[3], 1);

        auto fstate0 = fstuff0.initialize();
        fstuff0.process(fstate0, bounds(seq1), bounds(seq2));
        EXPECT_EQ(fstate0.counts[1], 1);
    }

    // Handles ambiguity properly.
    {
        auto state = stuff.initialize();
        std::string seq1 = "AAAATTATCGGCcacacacaAAAACCCCCGGC", seq2 = "AGCTCTCTCTTTTTcgtacgactAGCTTGTCTGTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2)); // ambiguous
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[3], 0);

        auto fstate = fstuff.initialize(); // takes the first
        fstuff.process(fstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(fstate.counts[3], 1);

        auto fstate0 = fstuff0.initialize(); // no match anyway
        fstuff0.process(fstate0, bounds(seq1), bounds(seq2));
        EXPECT_EQ(fstate0.counts[1], 0);
        EXPECT_EQ(fstate0.counts[3], 0);
    }

    // ... unless the ambiguous regions are the same.
    {
        auto state = stuff.initialize();
        std::string seq1 = "AAAATTTTCGGCcacacacaAAAATTTTCGGC", seq2 = "AGCTCTCTCTTTTTcgtacgactAGCTCTCTCTTTTT";
        stuff.process(state, bounds(seq1), bounds(seq2)); // ambiguous
        EXPECT_EQ(state.counts[3], 1);
    }
}

TEST_F(DualBarcodesPairedEndTest, ReverseComplementBest) {
    kaori::DualBarcodesPairedEnd<32> stuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.strand1 = kaori::SearchStrand::REVERSE;
            opt.strand2 = kaori::SearchStrand::REVERSE;
            opt.max_mismatches1 = 1;
            opt.max_mismatches2 = 1;
            opt.use_first = false;
            return opt;
        }()
    );

    kaori::DualBarcodesPairedEnd<32> fstuff(
        constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
        constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
        [&]{
            Options<32> opt;
            opt.strand1 = kaori::SearchStrand::REVERSE;
            opt.strand2 = kaori::SearchStrand::REVERSE;
            opt.max_mismatches1 = 1;
            opt.max_mismatches2 = 1;
            return opt;
        }()
    );

    std::string seq1 = "GCCGCCACTTTTcacacacacGCCGAAAATTTT"; // (GGGG, TTTT)
    std::string seq2 = "AAAACTCTCTAGCTcgtacgactAAAAAGAGAGAGCT"; // (AGAGAG, CTCTCT)

    auto state = stuff.initialize();
    stuff.process(state, bounds(seq1), bounds(seq2));
    EXPECT_EQ(state.counts[3], 1);

    auto fstate = fstuff.initialize();
    fstuff.process(fstate, bounds(seq1), bounds(seq2));
    EXPECT_EQ(fstate.counts[2], 1);
}

TEST_F(DualBarcodesPairedEndTest, RandomizedBest) {
    // No mismatches.
    {
        kaori::DualBarcodesPairedEnd<32> nonrandom(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            [&]{
                Options<32> opt;
                opt.use_first = false;
                return opt;
            }()
        );

        kaori::DualBarcodesPairedEnd<32> random(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            [&]{
                Options<32> opt;
                opt.use_first = false;
                opt.random = true;
                return opt;
            }()
        );

        std::string seq1 = "AAAATTTTCGGCcacacacaAGCTTGTGTGTTTT";
        std::string seq2 = "AGCTCTCTCTTTTTcgtacgacAAAACCCCCGGC";

        auto nstate = nonrandom.initialize();
        nonrandom.process(nstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nstate.counts[3], 1);
                
        auto rstate = random.initialize();
        random.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.counts, std::vector<int>(variables1.size())); // ambiguous.
    }

    // One mismatch.
    {
        kaori::DualBarcodesPairedEnd<32> nonrandom(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1), 
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
            [&]{
                Options<32> opt;
                opt.max_mismatches1 = 1;
                opt.max_mismatches2 = 1;
                opt.use_first = false;
                return opt;
            }()
        );

        kaori::DualBarcodesPairedEnd<32> random(
            constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
            constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
            [&]{
                Options<32> opt;
                opt.max_mismatches1 = 1;
                opt.max_mismatches2 = 1;
                opt.use_first = false;
                opt.random = true;
                return opt;
            }()
        );

        std::string seq1 = "AAAATTATCGGCcacacacaAGCTTGTGTGTTTT";
        std::string seq2 = "AGCTCTCTCTTTTTcgtacgacAAAACCCCCGGC";

        auto nstate = nonrandom.initialize();
        nonrandom.process(nstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(nstate.counts[3], 1);
                
        auto rstate = random.initialize();
        random.process(rstate, bounds(seq1), bounds(seq2));
        EXPECT_EQ(rstate.counts[1], 1);
        EXPECT_EQ(rstate.counts[3], 0);
    }
}

TEST_F(DualBarcodesPairedEndTest, Error) {
    constant2 = "ACGT---TGCA";
    EXPECT_ANY_THROW({
        try {
            kaori::DualBarcodesPairedEnd<32> stuff(
                constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
                constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2),
                Options<32>()
            );
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("should be the same") != std::string::npos);
            throw e;
        }
    });

    constant2 = "ACGTTGCA";
    EXPECT_ANY_THROW({
        try {
            kaori::DualBarcodesPairedEnd<32> stuff(
                constant1.c_str(), constant1.size(), kaori::BarcodePool(variables1),
                constant2.c_str(), constant2.size(), kaori::BarcodePool(variables2), 
                Options<32>()
            );
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("expected one variable region") != std::string::npos);
            throw e;
        }
    });
}
