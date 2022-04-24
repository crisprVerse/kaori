#include <gtest/gtest.h>
#include "kaori/handlers/SingleBarcodePairedEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

TEST(SingleBarcodePairedEnd, Forward) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq1{ 
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", 
        "ccacacacaaaaaACGTAATATTTT", // 1 mismatch
        "cAGGTAATATTTTtttttt" // 2 mismatches
    };
    std::string fq1 = convert_to_fastq(seq1);

    std::vector<std::string> seq2{ 
        "cagcatcgatcgtgaACGTCCCCTTTTacggaggaga", 
        "ccacacacaaaaaACGTCCCCTTTT", 
        "cACGTCCCCTTTTtttttt" 
    };
    std::string fq2 = convert_to_fastq(seq2);

    // No mismatches allowed.
    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 0, kaori::BarcodePool(variables));
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 1);
        EXPECT_EQ(counts[1], 2);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);

        EXPECT_EQ(handler.get_total(), 3);
    }

    // Okay, 2 mismatches, in which case the search on the first read is always favored.
    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 0, kaori::BarcodePool(variables), 2);
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 3);
        EXPECT_EQ(counts[1], 0);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);

        EXPECT_EQ(handler.get_total(), 3);
    }
}

TEST(SingleBarcodePairedEnd, Reverse) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq1{ 
        "cagcatcgatcgtgaAAAAAAAAACGTacggaggaga", 
        "ccacacacaaaaaAAAAAAATACGT", // 1 mismatch
        "cAAAAACAAAGGTtttttt" // 2 mismatches
    };
    std::string fq1 = convert_to_fastq(seq1);

    std::vector<std::string> seq2{ 
        "cAAAACCCCACGTacggaggaga", 
        "agagagaAAAACCCCACGTaaa", 
        "agagagagaAAAACCCCACGT" 
    };
    std::string fq2 = convert_to_fastq(seq2);

    // No mismatches allowed.
    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 1, kaori::BarcodePool(variables));
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 0);
        EXPECT_EQ(counts[2], 2);
        EXPECT_EQ(counts[3], 1);
    }

    // Okay, 2 mismatches, in which case the search on the first read is always favored.
    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 1, kaori::BarcodePool(variables), 2);
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 0);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 3);
    }
}

TEST(SingleBarcodePairedEnd, Best) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq1{ 
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", // 1a
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", // 1b
        "ccacacacaaaaaACGTAATATTTT", // 2
        "cACGTGGGGTTTTtttttt" // 3
    };
    std::string fq1 = convert_to_fastq(seq1);

    std::vector<std::string> seq2{ 
        "cagcatcgatcgtgaACGTCCCCTTTTacggaggaga", // 1a
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", // 1b
        "ccacacacaaaaaACGTCCCCTTTT", // 2
        "cAGGTAATATTTTtttttt" // 3
    };
    std::string fq2 = convert_to_fastq(seq2);

    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 0, kaori::BarcodePool(variables));
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        /*
         * 1a, 1b = first read, AAAA
         * 2 = second read wins, CCCC
         * 3 = first read wins, GGGG.
         */
        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 2);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 1);
        EXPECT_EQ(counts[3], 0);
    }

    {
        kaori::SingleBarcodePairedEnd<64> handler(thing.c_str(), thing.size(), 0, kaori::BarcodePool(variables));
        handler.set_first(false);
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fq1.c_str()), fq1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fq2.c_str()), fq2.size());
        kaori::process_paired_end_data(&reader1, &reader2, handler);

        /*
         * 1a = ambiguous and ignored.
         * 1b = match in both reads, but it's to the same ID, so whatever, we'll give it AAAA.
         * 2 = second read wins, CCCC
         * 3 = first read wins, GGGG.
         */
        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 1);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 1);
        EXPECT_EQ(counts[3], 0);
    }
}
