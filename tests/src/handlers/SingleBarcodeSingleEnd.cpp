#include <gtest/gtest.h>
#include "kaori/handlers/SingleBarcodeSingleEnd.hpp"
#include "kaori/process_data.hpp"
#include "byteme/RawBufferReader.hpp"
#include "../utils.h"
#include <string>

class SingleBarcodeSingleEndTest : public ::testing::Test {
protected:
    template<size_t max_size>
    using Options = typename kaori::SingleBarcodeSingleEnd<max_size>::Options;
};

TEST_F(SingleBarcodeSingleEndTest, ForwardOnly) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", 
        "ACGTCCCCTTTTaaaaccccggg",
        "ccacacacaaaaaACGTAATATTTT", // 1 mismatch
        "cAGGTAATATTTTtttttt" // 2 mismatches
    };
    std::string fq = convert_to_fastq(seq);

    // No mismatches allowed.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), Options<16>());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 1);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);

        EXPECT_EQ(handler.get_total(), 4);
    }

    // Okay, 2 mismatches.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.max_mismatches = 2;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 3);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);

        EXPECT_EQ(handler.get_total(), 4);
    }
}

TEST_F(SingleBarcodeSingleEndTest, Stranded) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaACGTCCCCTTTTacggaggaga", 
        "AAAAAAAAACGTaaaaccccggg",
        "accgggAAAATTCTACGTacaca", // 1 mismatch.
        "accgggACGTCCGCTTTT" // 1 mismatch.
    };
    std::string fq = convert_to_fastq(seq);

    // Forward only.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), Options<16>());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);
    }

    // Reverse only.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.search_forward = false;
            opt.search_reverse = true;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 0);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 1);
    }

    // Both.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 1);
    }

    // Both plus mismatches.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 2;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 1);
        EXPECT_EQ(counts[1], 2);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 1);
    }
}

TEST_F(SingleBarcodeSingleEndTest, Best) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "AAAA", "CCCC", "GGGG", "TTTT" };

    std::vector<std::string> seq{ 
        "accgggAAAATTCTACGTacacaACGTTTTTTTTT", 
        "accgggACGTCCGCTTTTcacacaAAAACCCCACGT" 
    };
    std::string fq = convert_to_fastq(seq);

    // First only.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 1);
        EXPECT_EQ(counts[1], 1);
        EXPECT_EQ(counts[2], 0);
        EXPECT_EQ(counts[3], 0);
    }

    // Best.
    {
        kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
            Options<16> opt;
            opt.search_forward = true;
            opt.search_reverse = true;
            opt.max_mismatches = 1;
            opt.use_first = false;
            return opt;
        }());
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
        kaori::process_single_end_data(&reader, handler);

        const auto& counts = handler.get_counts();
        EXPECT_EQ(counts[0], 0);
        EXPECT_EQ(counts[1], 0);
        EXPECT_EQ(counts[2], 1);
        EXPECT_EQ(counts[3], 1);
    }
}

TEST_F(SingleBarcodeSingleEndTest, Iupac) {
    std::string thing = "ACGT----TTTT";
    std::vector<std::string> variables { "ARRA", "CSSC", "GKKG", "TNNT" };

    std::vector<std::string> seq{ 
        "cagcatcgatcgtgaACGTAAAATTTTacggaggaga", 
        "acaccacaAAAACCCCACGTcacacacaca",
        "acaccacaAAAAACCAACGTcacacacaca",
        "cagcatcgatcgtgaACGTTCGTTTTTacggaggaga", 
        "acaccacaACGTCAACTTTTcacacacaca" // doesn't match anything
    };
    std::string fq = convert_to_fastq(seq);

    kaori::SingleBarcodeSingleEnd<16> handler(thing.c_str(), thing.size(), kaori::BarcodePool(variables), [&]{
        Options<16> opt;
        opt.search_forward = true;
        opt.search_reverse = true;
        return opt;
    }());
    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fq.c_str()), fq.size());
    kaori::process_single_end_data(&reader, handler);

    const auto& counts = handler.get_counts();
    EXPECT_EQ(counts[0], 1);
    EXPECT_EQ(counts[1], 0);
    EXPECT_EQ(counts[2], 1);
    EXPECT_EQ(counts[3], 2);
}
