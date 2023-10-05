#include <gtest/gtest.h>
#include "kaori/handlers/DualBarcodesSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class DualBarcodesSingleEndTest : public testing::Test {
protected:
    DualBarcodesSingleEndTest() : 
        constant("AAAA----CGGC------TTT"),
        variables1(std::vector<std::string>{ "AAAA", "CCCC", "GGGG", "TTTT" }),
        variables2(std::vector<std::string>{ "ACACAC", "TGTGTG", "AGAGAG", "CTCTCT" })
    {}

    std::string constant;
    std::vector<std::string> variables1;
    std::vector<std::string> variables2;

    template<size_t max_size>
    using Options = typename kaori::DualBarcodesSingleEnd<max_size>::Options;
};

TEST_F(DualBarcodesSingleEndTest, BasicFirst) {
    kaori::DualBarcodesSingleEnd<32> stuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        Options<32>()
    );

    // Works in the simple case.
    {
        auto state = stuff.initialize();
        std::string seq = "AAAATTTTCGGCCTCTCTTTT";
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 1);
        EXPECT_EQ(state.total, 1);
    }

    // Making it work for it.
    {
        auto state = stuff.initialize();
        std::string seq = "cacacacAAAAAAAACGGCACACACTTTacagcat";
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[0], 1);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 0);
        EXPECT_EQ(state.total, 1);
    }

    // Integrated.
    {
        std::vector<std::string> seq{ 
            "cagcatcgatcgtgaAAAACCCCCGGCTGTGTGTTTacggaggaga", 
            "AAAAGGGGCGGCAGAGAGTTTaaaaccccggg"
        };
        std::string fq = convert_to_fastq(seq);
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, stuff);

        EXPECT_EQ(stuff.get_counts()[0], 0);
        EXPECT_EQ(stuff.get_counts()[1], 1);
        EXPECT_EQ(stuff.get_counts()[2], 1);
        EXPECT_EQ(stuff.get_counts()[3], 0);
        EXPECT_EQ(stuff.get_total(), 2);
    }
}

TEST_F(DualBarcodesSingleEndTest, ReverseComplementFirst) {
    kaori::DualBarcodesSingleEnd<32> stuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.strand = kaori::SearchStrand::REVERSE;
            return opt;
        }()
    );

    auto state = stuff.initialize();
    std::string seq = "AAAAGAGAGGCCGAAAATTTTaacacac"; 
    stuff.process(state, bounds(seq));
    EXPECT_EQ(state.counts[0], 0);
    EXPECT_EQ(state.counts[1], 0);
    EXPECT_EQ(state.counts[2], 0);
    EXPECT_EQ(state.counts[3], 1);
    EXPECT_EQ(state.total, 1);

    seq = "ccacacAAACTCTCTGCCGCCCCTTTTcgata";
    stuff.process(state, bounds(seq));
    EXPECT_EQ(state.counts[2], 1);
    EXPECT_EQ(state.total, 2);
}

TEST_F(DualBarcodesSingleEndTest, Iupac) {
    variables1 = std::vector<std::string>{ "ARRA", "CYYC", "GSSG", "TWWT" };
    variables2 = std::vector<std::string>{ "ABACAC", "TVTGTG", "AGHGAG", "CDCTCT" };

    // Forward.
    {
        kaori::DualBarcodesSingleEnd<32> stuff(
            constant.c_str(), constant.size(),
            std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
            Options<32>()
        );

        auto state = stuff.initialize();
        std::string seq = "AAAATTTTCGGCCTCTCTTTTT";
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 0);
        EXPECT_EQ(state.counts[3], 1);
        EXPECT_EQ(state.total, 1);
    }

    // Reverse-complement.
    {
        kaori::DualBarcodesSingleEnd<32> stuff(
            constant.c_str(), constant.size(),
            std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
            [&]{
                Options<32> opt;
                opt.strand = kaori::SearchStrand::REVERSE;
                return opt;
            }()
        );

        auto state = stuff.initialize();
        std::string seq = "cgatAAACTCTCTGCCGCCCCTTTTacttc";
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[0], 0);
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[2], 1);
        EXPECT_EQ(state.counts[3], 0);
        EXPECT_EQ(state.total, 1);
    }
}

TEST_F(DualBarcodesSingleEndTest, MismatchFirst) {
    kaori::DualBarcodesSingleEnd<32> stuff0(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        Options<32>()
    );

    kaori::DualBarcodesSingleEnd<32> stuff1(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.max_mismatches = 1;
            return opt;
        }()
    );

    kaori::DualBarcodesSingleEnd<32> stuff2(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.max_mismatches = 2;
            return opt;
        }()
    );

    // One mismatch.
    {
        std::string seq = "AAAATTATCGGCCTCTCTTTT";

        auto state0 = stuff0.initialize();
        stuff0.process(state0, bounds(seq));
        EXPECT_EQ(state0.counts[3], 0);

        auto state1 = stuff1.initialize();
        stuff1.process(state1, bounds(seq));
        EXPECT_EQ(state1.counts[3], 1);
    }

    // Two mismatches, in constant and variable regions.
    {
        std::string seq = "AAAATTATCGGGCTCTCTTTT";

        auto state1 = stuff1.initialize();
        stuff1.process(state1, bounds(seq));
        EXPECT_EQ(state1.counts[3], 0);

        auto state2 = stuff2.initialize();
        stuff2.process(state2, bounds(seq));
        EXPECT_EQ(state2.counts[3], 1);
    }

    // Two mismatches in variable regions.
    {
        std::string seq = "AAAATTATCGGCCTCTCATTT";

        auto state1 = stuff1.initialize();
        stuff1.process(state1, bounds(seq));
        EXPECT_EQ(state1.counts[3], 0);

        auto state2 = stuff2.initialize();
        stuff2.process(state2, bounds(seq));
        EXPECT_EQ(state2.counts[3], 1);
    }
}

TEST_F(DualBarcodesSingleEndTest, AmbiguityFirst) {
    variables1.push_back("AAAA");
    variables2.push_back("AAAAAG");
    variables1.push_back("AAAA");
    variables2.push_back("AAAAAT");

    kaori::DualBarcodesSingleEnd<32> stuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.max_mismatches = 1;
            return opt;
        }()
    );

    {
        std::string seq = "AAAAAAAACGGCAAAAACTTTT";

        auto state = stuff.initialize();
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts, std::vector<int>(variables1.size())); // nothing detected.

        seq = "AAAAAAAAACGGCAAAAATTTTT"; // as a control.
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[5], 1);
    }
}

TEST_F(DualBarcodesSingleEndTest, BasicBest) {
    kaori::DualBarcodesSingleEnd<32> stuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.max_mismatches = 1;
            opt.use_first = false;
            return opt;
        }()
    );

    kaori::DualBarcodesSingleEnd<32> fstuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.max_mismatches = 1;
            return opt;
        }()
    );

    kaori::DualBarcodesSingleEnd<32> fstuff0(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        Options<32>()
    );

    // Keeps searching for the best.
    {
        std::string seq = "atatataAAAATTATCGGCCTCTCTTTTcacacacaAAAACCCCCGGCTGTGTGTTTacaca";

        auto state = stuff.initialize();
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[1], 1);

        auto fstate = fstuff.initialize();
        fstuff.process(fstate, bounds(seq));
        EXPECT_EQ(fstate.counts[3], 1);

        auto fstate0 = fstuff0.initialize();
        fstuff0.process(fstate0, bounds(seq));
        EXPECT_EQ(fstate0.counts[1], 1);
    }

    // Handles ambiguity properly.
    {
        auto state = stuff.initialize();
        std::string seq = "AAAATTATCGGCCTCTCTTTTcacacacaAAAACCCCCGGCTGTCTGTTT";
        stuff.process(state, bounds(seq)); // ambiguous
        EXPECT_EQ(state.counts[1], 0);
        EXPECT_EQ(state.counts[3], 0);

        auto fstate = fstuff.initialize(); // takes the first
        fstuff.process(fstate, bounds(seq));
        EXPECT_EQ(fstate.counts[3], 1);

        auto fstate0 = fstuff0.initialize(); // no match anyway
        fstuff0.process(fstate0, bounds(seq));
        EXPECT_EQ(fstate0.counts[1], 0);
        EXPECT_EQ(fstate0.counts[3], 0);
    }

    // ... unless the ambiguous regions are the same.
    {
        auto state = stuff.initialize();
        std::string seq = "AAAATTTTCGGCCTCTCTTTTcacacacaAAAATTTTCGGCCTCTCTTTT";
        stuff.process(state, bounds(seq));
        EXPECT_EQ(state.counts[3], 1);
    }
}

TEST_F(DualBarcodesSingleEndTest, ReverseComplementBest) {
    kaori::DualBarcodesSingleEnd<32> stuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.strand = kaori::SearchStrand::REVERSE;
            opt.max_mismatches = 1;
            opt.use_first = false;
            return opt;
        }()
    );

    kaori::DualBarcodesSingleEnd<32> fstuff(
        constant.c_str(), constant.size(),
        std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
        [&]{
            Options<32> opt;
            opt.strand = kaori::SearchStrand::REVERSE;
            opt.max_mismatches = 1;
            return opt;
        }()
    );

    //                  (G^TGG, AGAGAG)               (TTTT, CTCTCT)
    std::string seq = "AAACTCTCTGCCGCCACTTTTcacacacacAAAAGAGAGGCCGAAAATTTT"; 

    auto state = stuff.initialize();
    stuff.process(state, bounds(seq));
    EXPECT_EQ(state.counts[3], 1);

    auto fstate = fstuff.initialize();
    fstuff.process(fstate, bounds(seq));
    EXPECT_EQ(fstate.counts[2], 1);
}

TEST_F(DualBarcodesSingleEndTest, Error) {
    constant = "ACGT---TGCA";
    EXPECT_ANY_THROW({
        try {
            kaori::DualBarcodesSingleEnd<32> stuff(
                constant.c_str(), constant.size(), 
                std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
                Options<32>()
            );
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("number of variable regions") != std::string::npos);
            throw e;
        }
    });

    constant = "ACGT---TGCA---ACT";
    EXPECT_ANY_THROW({
        try {
            kaori::DualBarcodesSingleEnd<32> stuff(
                constant.c_str(), constant.size(), 
                std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
                Options<32>()
            );
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("length of variable region") != std::string::npos);
            throw e;
        }
    });


    constant = "ACGT----TGCA------ACT";
    variables1.pop_back();
    EXPECT_ANY_THROW({
        try {
            kaori::DualBarcodesSingleEnd<32> stuff(
                constant.c_str(), constant.size(), 
                std::vector<kaori::BarcodePool>{ kaori::BarcodePool(variables1), kaori::BarcodePool(variables2) },
                Options<32>()
            );
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("same length") != std::string::npos);
            throw e;
        }
    });
}
