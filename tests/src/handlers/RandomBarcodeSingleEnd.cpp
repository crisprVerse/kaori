#include <gtest/gtest.h>
#include "kaori/handlers/RandomBarcodeSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

TEST(RandomBarcodeSingleEnd, ForwardOnly) {
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
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 0);
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
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 0, 1);
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
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 0, 2);
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

TEST(RandomBarcodeSingleEnd, Stranded) {
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
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 0);
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 1);
        EXPECT_EQ(counts["CCCC"], 1);
    }

    // Reverse only.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 1);
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 1);
        EXPECT_EQ(counts["TTTT"], 1);
    }

    // Both.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 2);
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 2);
        EXPECT_EQ(counts["CCCC"], 1);
        EXPECT_EQ(counts["TTTT"], 1);
    }

    // Both plus mismatches.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 2, 2);
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

TEST(RandomBarcodeSingleEnd, Best) {
    std::string thing = "ACGT----TTTT";

    std::vector<std::string> seq{ 
        "accgggAAGATTCTACGTacacaACGTTTTTTTTT", // mismatch reverse, followed by perfect forward
        "accgggACGTCCGCTTTTcacacaAAAACCCCACGAttt", // perfect forward, followed by mismatch reverse
        "accgggAAAACCCCACGTacacaACGTGCCGTTTT", // two perfects: ambiguous and ignored in the best strat.
    };
    std::string fq = convert_to_fastq(seq);

    // First only.
    {
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 2, 1);
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
        kaori::RandomBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), 2, 1);
        handler.set_first(false);
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        auto counts = handler.get_counts();
        EXPECT_EQ(counts.size(), 2);
        EXPECT_EQ(counts["TTTT"], 1);
        EXPECT_EQ(counts["CCGC"], 1);
    }
}
