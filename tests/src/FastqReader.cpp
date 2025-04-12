#include <gtest/gtest.h>
#include "kaori/FastqReader.hpp"
#include "byteme/RawBufferReader.hpp"
#include "byteme/RawFileReader.hpp"
#include <fstream>

TEST(BasicTests, Single) {
    std::string buffer = "@FOO\nACGT\n+\n!!!!"; // check it works without a terminating newline.
    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
    kaori::FastqReader fq(&reader);

    EXPECT_TRUE(fq());
    const auto& name = fq.get_name();
    EXPECT_EQ(std::string(name.begin(), name.end()), "FOO");
    const auto& seq = fq.get_sequence();
    EXPECT_EQ(std::string(seq.begin(), seq.end()), "ACGT");

    EXPECT_FALSE(fq());
}

TEST(BasicTests, Empty) {
    std::string buffer = "@FOO\n\n+\n";
    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
    kaori::FastqReader fq(&reader);

    EXPECT_TRUE(fq());
    const auto& name = fq.get_name();
    EXPECT_EQ(std::string(name.begin(), name.end()), "FOO");
    const auto& seq = fq.get_sequence();
    EXPECT_EQ(std::string(seq.begin(), seq.end()), "");

    EXPECT_FALSE(fq());
}

TEST(BasicTests, MultipleEntries) {
    std::string buffer = "@FOO and more info\nACGT\n+\n!!!!\n@WHEE\nTGCA\n+asdasd\naaaa\n";
    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
    kaori::FastqReader fq(&reader);

    EXPECT_TRUE(fq());
    const auto& name = fq.get_name();
    EXPECT_EQ(std::string(name.begin(), name.end()), "FOO");
    const auto& seq = fq.get_sequence();
    EXPECT_EQ(std::string(seq.begin(), seq.end()), "ACGT");

    EXPECT_TRUE(fq());
    EXPECT_EQ(std::string(name.begin(), name.end()), "WHEE");
    EXPECT_EQ(std::string(seq.begin(), seq.end()), "TGCA");

    EXPECT_FALSE(fq());
}

TEST(BasicTests, MultiLineEntries) {
    for (size_t i = 0; i < 2; ++i) {
        std::string buffer = "@FOO\nA\nCG\nTGCA\n+\n!!\n!!\n!!!\n@ARG\nACACGGT\nC\n+\n@@@\n@\n@@@@";
        if (i == 2) {
            buffer += '\n';
        }

        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);

        EXPECT_TRUE(fq());
        const auto& name = fq.get_name();
        EXPECT_EQ(std::string(name.begin(), name.end()), "FOO");
        const auto& seq = fq.get_sequence();
        EXPECT_EQ(std::string(seq.begin(), seq.end()), "ACGTGCA");

        EXPECT_TRUE(fq());
        const auto& name2 = fq.get_name();
        EXPECT_EQ(std::string(name2.begin(), name2.end()), "ARG");
        const auto& seq2 = fq.get_sequence();
        EXPECT_EQ(std::string(seq2.begin(), seq2.end()), "ACACGGTC");

        EXPECT_FALSE(fq());
    }
}

TEST(BasicTests, Errors) {
    {
        std::string buffer = "FOO";
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);

        EXPECT_ANY_THROW({
            try {
                fq();
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("read name should start") != std::string::npos);
                throw e;
            }
        });
    }

    {
        std::string buffer = "@FOO\nAC\n+";
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);

        EXPECT_ANY_THROW({
            try {
                fq();
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("premature end") != std::string::npos);
                throw e;
            }
        });
    }

    {
        std::string buffer = "@FOO\nAC\n+\n!!\n@WHEE\nACGT\n+\n!!"; // too short.
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);
        fq();

        EXPECT_ANY_THROW({
            try {
                fq();
            } catch (std::exception& e) {
                std::string msg(e.what());
                EXPECT_TRUE(msg.find("non-equal lengths") != std::string::npos);
                EXPECT_TRUE(msg.find("line 5") != std::string::npos);
                throw e;
            }
        });
    }

    {
        std::string buffer = "@FOO\nAC\n+\n!!\n@WHEE\nACGT\n+\n!!!@@!!@\n"; // too long.
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);
        fq();

        EXPECT_ANY_THROW({
            try {
                fq();
            } catch (std::exception& e) {
                std::string msg(e.what());
                EXPECT_TRUE(msg.find("non-equal lengths") != std::string::npos);
                EXPECT_TRUE(msg.find("line 5") != std::string::npos);
                throw e;
            }
        });
    }

    {
        std::string buffer = "@FOO\nAC\n+\n!!\nWHEE";
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.size());
        kaori::FastqReader fq(&reader);
        fq();

        EXPECT_ANY_THROW({
            try {
                fq();
            } catch (std::exception& e) {
                std::string msg(e.what());
                EXPECT_TRUE(msg.find("should start") != std::string::npos);
                EXPECT_TRUE(msg.find("line 5") != std::string::npos);
                throw e;
            }
        });
    }
}

class FastqReaderFileTest : public testing::TestWithParam<int> {};

TEST_P(FastqReaderFileTest, LongStrings) {
    std::string path = "TEST_reader.fastq";
    {
        std::ofstream out(path);
        out << "@FOOcagctacgt cagtcgact acgacactgactacg\nACGATCGATACGCATCGATCGATCGATCGATTTTAT\n+actgctacgatcgatcagcta\n!#!#!###!!#!#!#!!!!#!#!#!#!#!#!#!#!!";
    }
    byteme::RawFileReader reader(path.c_str(), [&]{
        byteme::RawFileReaderOptions ropt;
        ropt.buffer_size = GetParam();
        return ropt;
    }());
    kaori::FastqReader fq(&reader);

    EXPECT_TRUE(fq());
    const auto& name = fq.get_name();
    EXPECT_EQ(std::string(name.begin(), name.end()), "FOOcagctacgt");
    const auto& seq = fq.get_sequence();
    EXPECT_EQ(std::string(seq.begin(), seq.end()), "ACGATCGATACGCATCGATCGATCGATCGATTTTAT");

    EXPECT_FALSE(fq());
}

TEST_P(FastqReaderFileTest, StressTest) {
    std::string path = "TEST_reader.fastq";
    {
        std::ofstream out(path);
        for (size_t i = 0; i < 1000; ++i) {
            out << "@" << "READ_" << i << "\n";
            out << "AAAAAAAAAAAAAAAaaaaaaaaaaaaaa\n";
            out << "CCCCCCCCCCCCCCCcccccccccccccc\n";
            out << "GGGGGGGGGGGGGGGgggggggggggggg\n";
            out << "TTTTTTTTTTTTTTTtttttttttttttt\n";
            out << "+" << "\n";
            for (int i = 0; i < 4; ++i) {
                out << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            }
        }
    }
    byteme::RawFileReader reader(path.c_str(), [&]{
        byteme::RawFileReaderOptions ropt;
        ropt.buffer_size = GetParam();
        return ropt;
    }());
    kaori::FastqReader fq(&reader);

    std::string ref = "AAAAAAAAAAAAAAAaaaaaaaaaaaaaa";
    ref += "CCCCCCCCCCCCCCCcccccccccccccc";
    ref += "GGGGGGGGGGGGGGGgggggggggggggg";
    ref += "TTTTTTTTTTTTTTTtttttttttttttt";

    for (size_t i = 0; i < 1000; ++i) {
        EXPECT_TRUE(fq());
        const auto& name = fq.get_name();
        EXPECT_EQ(std::string(name.begin(), name.end()), "READ_" + std::to_string(i));
        const auto& seq = fq.get_sequence();
        EXPECT_EQ(std::string(seq.begin(), seq.end()), ref);
    }

    EXPECT_FALSE(fq());
}

INSTANTIATE_TEST_SUITE_P(
    FastqReader,
    FastqReaderFileTest, 
    ::testing::Values(5, 10, 50, 1000)
);
