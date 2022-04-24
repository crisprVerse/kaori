#include <gtest/gtest.h>
#include "kaori/process_data.hpp"
#include <random>
#include "byteme/RawBufferReader.hpp"
#include "utils.h"

template<bool unames, bool failtest = false>
class SingleEndCollector {
public:
    struct State {
        std::vector<std::string> reads, names;
    };

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        if constexpr(failtest) {
            throw std::runtime_error("I want a burger");
        }
        state.reads.emplace_back(x.first, x.second);
    }

    void process(State& state, const std::pair<const char*, const char*>& x, const std::pair<const char*, const char*>& y) const {
        state.names.emplace_back(x.first, x.second);
        state.reads.emplace_back(y.first, y.second);
    };

    State initialize() {
        return State();
    }

    void reduce(State& x) {
        collected_reads.insert(collected_reads.end(), x.reads.begin(), x.reads.end());
        if constexpr(use_names) {
            collected_names.insert(collected_names.end(), x.names.begin(), x.names.end());
        }
    }

    static constexpr bool use_names = unames;

    std::vector<std::string> collected_reads, collected_names;
};

class ProcessDataTester : public testing::TestWithParam<std::tuple<int, int> > {
protected:
    std::vector<std::string> simulate_reads(int seed) {
        const size_t n = 1000;
        std::mt19937_64 rng(seed);
        std::vector<std::string> output;
        output.reserve(n);

        for (size_t i = 1; i <= n; ++i) {
            std::string current;
            size_t n = rng() % 20 + 10;
            for (size_t j = 0; j < n; ++j) {
                switch (rng() % 4) {
                    case 0:
                        current += "A";
                        break;
                    case 1:
                        current += "C";
                        break;
                    case 2:
                        current += "G";
                        break;
                    case 3:
                        current += "T";
                        break;
                }
            }
            output.push_back(current);
        }

        return output;
    }
};

TEST_P(ProcessDataTester, SingleEnd) {
    auto param = GetParam();
    auto nthreads = std::get<0>(param);
    auto blocksize = std::get<1>(param);

    auto reads = simulate_reads(nthreads + blocksize);
    auto fastq_str = convert_to_fastq(reads);

    // Without names.
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
        SingleEndCollector<false> task;
        kaori::process_single_end_data(&reader, task, nthreads, blocksize);
        EXPECT_EQ(task.collected_reads, reads);
        EXPECT_TRUE(task.collected_names.empty());
    }

    // Plus names.
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
        SingleEndCollector<true> task;
        kaori::process_single_end_data(&reader, task, nthreads, blocksize);
        EXPECT_EQ(task.collected_reads, reads);
        EXPECT_EQ(task.collected_names.size(), reads.size());

        bool all_okay = true;
        for (size_t i = 0; i < task.collected_names.size(); ++i) {
            auto current = task.collected_names[i];
            std::string expected = "READ" + std::to_string(i + 1);
            if (current != expected) {
                all_okay = false;
            }
        }
        EXPECT_TRUE(all_okay);
    }
}

TEST_P(ProcessDataTester, SingleEndErrors) {
    // Errors in the processing are caught and handled correctly,
    // especially with respect to closing down all the threads.
    auto param = GetParam();
    auto nthreads = std::get<0>(param);
    auto blocksize = std::get<1>(param);

    auto reads = simulate_reads(nthreads + blocksize);
    auto fastq_str = convert_to_fastq(reads, "FOO");

    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
    SingleEndCollector<false, true> task;

    EXPECT_ANY_THROW({
        try {
            kaori::process_single_end_data(&reader, task, nthreads, blocksize);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()) == "I want a burger");
            throw e;
        }
    });
}

template<bool unames, bool failtest = false>
class PairedEndCollector {
public:
    SingleEndCollector<unames> read1, read2;

    struct State {
        typename SingleEndCollector<unames>::State read1, read2;
    };

    void process(State& state, const std::pair<const char*, const char*>& x1, const std::pair<const char*, const char*>& x2) const {
        if constexpr(failtest) {
            throw std::runtime_error("I want some fries");
        }
        read1.process(state.read1, x1);
        read2.process(state.read2, x2);
    }

    void process(State& state, 
        const std::pair<const char*, const char*>& x1,
        const std::pair<const char*, const char*>& y1,
        const std::pair<const char*, const char*>& x2,
        const std::pair<const char*, const char*>& y2)
    const {
        read1.process(state.read1, x1, y1);
        read2.process(state.read2, x2, y2);
    }

    State initialize() {
        return State();
    }

    void reduce(State& x) {
        read1.reduce(x.read1);
        read2.reduce(x.read2);
    }

    static constexpr bool use_names = unames;
};

TEST_P(ProcessDataTester, PairedEnd) {
    auto param = GetParam();
    auto nthreads = std::get<0>(param);
    auto blocksize = std::get<1>(param);

    auto reads1 = simulate_reads(nthreads + blocksize);
    auto reads2 = simulate_reads((nthreads + blocksize) * 2);
    auto fastq_str1 = convert_to_fastq(reads1, "FOO");
    auto fastq_str2 = convert_to_fastq(reads2, "BAR");

    // Without names.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<false> task;
        kaori::process_paired_end_data(&reader1, &reader2, task, nthreads, blocksize);

        EXPECT_EQ(task.read1.collected_reads, reads1);
        EXPECT_EQ(task.read2.collected_reads, reads2);
        EXPECT_TRUE(task.read1.collected_names.empty());
        EXPECT_TRUE(task.read2.collected_names.empty());
    }

    // Plus names.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<true> task;
        kaori::process_paired_end_data(&reader1, &reader2, task, nthreads, blocksize);

        EXPECT_EQ(task.read1.collected_reads, reads1);
        EXPECT_EQ(task.read2.collected_reads, reads2);
        EXPECT_EQ(task.read1.collected_names.size(), reads1.size());
        EXPECT_EQ(task.read2.collected_names.size(), reads2.size());

        bool all_okay = true;
        for (size_t i = 0; i < reads1.size(); ++i) {
            if (task.read1.collected_names[i] != "FOO" + std::to_string(i + 1)) {
                all_okay = false;
            }
            if (task.read2.collected_names[i] != "BAR" + std::to_string(i + 1)) { 
                all_okay = false;
            }
        }
        EXPECT_TRUE(all_okay);
    }
}

TEST_P(ProcessDataTester, PairedEndErrors) {
    auto param = GetParam();
    auto nthreads = std::get<0>(param);
    auto blocksize = std::get<1>(param);

    auto reads1 = simulate_reads(nthreads + blocksize);
    auto reads2 = simulate_reads((nthreads + blocksize) * 2);
    auto fastq_str1 = convert_to_fastq(reads1, "FOO");
    auto fastq_str2 = convert_to_fastq(reads2, "BAR");

    // Errors out correctly due to the process.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<false, true> task;
        EXPECT_ANY_THROW({
            try {
                kaori::process_paired_end_data(&reader1, &reader2, task, nthreads, blocksize);
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("I want some fries") != std::string::npos);
                throw e;
            }
        });
    }

    // Errors out correctly due to the read number.
    std::vector<size_t> size_options { 10, static_cast<size_t>(blocksize), reads2.size() - 1 };
    for (auto resized : size_options) {
        auto reads3 = reads2;
        reads3.resize(resized);
        auto fastq_str3 = convert_to_fastq(reads3, "FOO");

        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader3(reinterpret_cast<const unsigned char*>(fastq_str3.c_str()), fastq_str3.size());

        PairedEndCollector<false> task;
        EXPECT_ANY_THROW({
            try {
                kaori::process_paired_end_data(&reader1, &reader3, task, nthreads, blocksize);
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("different number of reads") != std::string::npos);
                throw e;
            }
        });
    }
}

INSTANTIATE_TEST_SUITE_P(
    ProcessData,
    ProcessDataTester, 
    ::testing::Combine(
        ::testing::Values(1, 2, 3), // number of threads
        ::testing::Values(10, 33, 71) // block size
    )
);
