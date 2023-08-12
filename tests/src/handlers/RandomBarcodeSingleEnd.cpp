#include <gtest/gtest.h>
#include "kaori/handlers/RandomBarcodeSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class RandomBarcodeSingleEndTest : public ::testing::Test {
protected:
    template<size_t max_size>
    using Options = typename kaori::RandomBarcodeSingleEnd<max_size>::Options;
};

TEST_F(RandomBarcodeSingleEndTest, ForwardOnly) {
    std::string thing = "ACGT----TTTT";

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", 
        "ACGTCCCCTTTTaaaaccccggg",
        "ccacacacaaaaaACGTAATATTGT", // 1 mismatch
        "cAGGTAATATCTTtttttt" // 2 mismatches
    };
    std::string fq = convert_to_fastq(seq);

    // No mismatches allowed.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), Options<16>());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 2);
        EXPECT_EQ(counts["AAAA"], 1);
        EXPECT_EQ(counts["CCCC"], 1);

        EXPECT_EQ(handler.get_total(), 4);
    }

    // Okay, 1 mismatch.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.max_mismatches = 1;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 3);
        EXPECT_EQ(counts["AAAA"], 1);
        EXPECT_EQ(counts["CCCC"], 1);
        EXPECT_EQ(counts["AATA"], 1);

        EXPECT_EQ(handler.get_total(), 4);
    }

    // Okay, 2 mismatches.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.max_mismatches = 2;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 3);
        EXPECT_EQ(counts["AAAA"], 1);
        EXPECT_EQ(counts["CCCC"], 1);
        EXPECT_EQ(counts["AATA"], 2);

        EXPECT_EQ(handler.get_total(), 4);
    }
}

TEST_F(RandomBarcodeSingleEndTest, Stranded) {
    std::string thing = "ACGT----TTTT";

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaACGTCCCCTTTTacggaggaga", 
        "AAAAAAAAACGTaaaaccccggg",
        "accgggAATATTCTACGTacaca", // 1 mismatch.
        "accgggAGGTCCGCTTTT" // 1 mismatch.
    };
    std::string fq = convert_to_fastq(seq);

    // Forward only.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), Options<16>());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 1);
        EXPECT_EQ(counts["CCCC"], 1);
    }

    // Reverse only.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.search_forward = false;
            opt.search_reverse = true;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 1);
        EXPECT_EQ(counts["TTTT"], 1);
    }

    // Both.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 2);
        EXPECT_EQ(counts["CCCC"], 1);
        EXPECT_EQ(counts["TTTT"], 1);
    }

    // Both plus mismatches.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{ 
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 2;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 4);
        EXPECT_EQ(counts["CCCC"], 1);
        EXPECT_EQ(counts["TTTT"], 1);
        EXPECT_EQ(counts["AGAA"], 1);
        EXPECT_EQ(counts["CCGC"], 1);
    }
}

TEST_F(RandomBarcodeSingleEndTest, Best) {
    std::string thing = "ACGT----TTTT";

    std::vector<std::string> seq{ 
        "accgggAAGATTCTACGTacacaACGTTTTTTTTT", // mismatch reverse, followed by perfect forward
        "accgggACGTCCGCTTTTcacacaAAAACCCCACGAttt", // perfect forward, followed by mismatch reverse
        "accgggAAAACCCCACGTacacaACGTGCCGTTTT", // two perfects: ambiguous and ignored in the best strat.
    };
    std::string fq = convert_to_fastq(seq);

    // First only.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 3);
        EXPECT_EQ(counts["AGAA"], 1);
        EXPECT_EQ(counts["CCGC"], 1);
        EXPECT_EQ(counts["GGGG"], 1);
    }

    // Best.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{;
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            opt.use_first = false;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 2);
        EXPECT_EQ(counts["TTTT"], 1);
        EXPECT_EQ(counts["CCGC"], 1);
    }
}

TEST_F(RandomBarcodeSingleEndTest, BadBases) {
    std::string thing = "ACGT----TTTT";

    // Passes N's through correctly.
    {
        std::vector<std::string> seq{ 
            "accgggAAAATNCTACGTacaca",
            "accACGTAANATTTTactgcgact",
            "accACGTNNNNTTTTactgcgact"
        };
        std::string fq = convert_to_fastq(seq);

        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 3);
        EXPECT_EQ(counts["AGNA"], 1);
        EXPECT_EQ(counts["AANA"], 1);
        EXPECT_EQ(counts["NNNN"], 1);
    }

    // Fails on seeing bad bases.
    {
        std::vector<std::string> seq{ 
            "accgggAAAATXCACGTacaca"
        };

        std::string fq = convert_to_fastq(seq);

        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());

        std::string failed;
        try {
            kaori::process_single_end_data(&reader, handler);
        } catch(std::exception& e) {
            failed = e.what();
        }
        EXPECT_TRUE(failed.find("X") != std::string::npos);
    }
}
