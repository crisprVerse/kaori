#include <gtest/gtest.h>
#include "kaori/process_data.hpp"
#include <random>
#include "byteme/RawBufferReader.hpp"
#include "utils.h"

template<bool unames_, bool failtest_ = false>
class SingleEndCollector {
public:
    struct State {
        std::vector<std::string> reads, names;
    };

    void process(State& state, const std::pair<const char*, const char*>& x) const {
        if constexpr(failtest_) {
            throw std::runtime_error("I want a burger");
        }
        state.reads.emplace_back(x.first, x.second);
    }

    void process(State& state, const std::pair<const char*, const char*>& x, const std::pair<const char*, const char*>& y) const {
        state.names.emplace_back(x.first, x.second);
        state.reads.emplace_back(y.first, y.second);
    };

    State initialize() const {
        return State();
    }

    void reduce(State& x) {
        my_collected_reads.insert(my_collected_reads.end(), x.reads.begin(), x.reads.end());
        if constexpr(use_names) {
            my_collected_names.insert(my_collected_names.end(), x.names.begin(), x.names.end());
        }
    }

    static constexpr bool use_names = unames_;

private:
    std::vector<std::string> my_collected_reads, my_collected_names;

public:
    const auto& reads() const {
        return my_collected_reads;
    }

    const auto& names() const {
        return my_collected_names;
    }
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
    kaori::ProcessSingleEndDataOptions popt;
    popt.num_threads = std::get<0>(param);
    popt.block_size = std::get<1>(param);

    auto reads = simulate_reads(popt.num_threads + popt.block_size);
    auto fastq_str = convert_to_fastq(reads);

    // Without names.
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
        SingleEndCollector<false> task;
        kaori::process_single_end_data(&reader, task, popt);
        EXPECT_EQ(task.reads(), reads);
        EXPECT_TRUE(task.names().empty());
    }

    // Plus names.
    {
        byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
        SingleEndCollector<true> task;
        kaori::process_single_end_data(&reader, task, popt);
        EXPECT_EQ(task.reads(), reads);
        const auto& names = task.names();
        EXPECT_EQ(names.size(), reads.size());

        bool all_okay = true;
        for (size_t i = 0, end = names.size(); i < end; ++i) {
            auto current = names[i];
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
    kaori::ProcessSingleEndDataOptions popt;
    popt.num_threads = std::get<0>(param);
    popt.block_size = std::get<1>(param);

    auto reads = simulate_reads(popt.num_threads + popt.block_size);
    auto fastq_str = convert_to_fastq(reads, "FOO");

    byteme::RawBufferReader reader(reinterpret_cast<const unsigned char*>(fastq_str.c_str()), fastq_str.size());
    SingleEndCollector<false, true> task;

    EXPECT_ANY_THROW({
        try {
            kaori::process_single_end_data(&reader, task, popt);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()) == "I want a burger");
            throw e;
        }
    });
}

template<bool unames_, bool failtest_ = false>
class PairedEndCollector {
public:
    struct State {
        typename SingleEndCollector<unames_>::State read1, read2;
    };

    void process(State& state, const std::pair<const char*, const char*>& x1, const std::pair<const char*, const char*>& x2) const {
        if constexpr(failtest_) {
            throw std::runtime_error("I want some fries");
        }
        my_read1.process(state.read1, x1);
        my_read2.process(state.read2, x2);
    }

    void process(State& state, 
        const std::pair<const char*, const char*>& x1,
        const std::pair<const char*, const char*>& y1,
        const std::pair<const char*, const char*>& x2,
        const std::pair<const char*, const char*>& y2)
    const {
        my_read1.process(state.read1, x1, y1);
        my_read2.process(state.read2, x2, y2);
    }

    State initialize() const {
        return State();
    }

    void reduce(State& x) {
        my_read1.reduce(x.read1);
        my_read2.reduce(x.read2);
    }

    static constexpr bool use_names = unames_;

private:
    SingleEndCollector<unames_> my_read1, my_read2;

public:
    const auto& first_reads() const {
        return my_read1.reads();
    }

    const auto& second_reads() const {
        return my_read2.reads();
    }

    const auto& first_names() const {
        return my_read1.names();
    }

    const auto& second_names() const {
        return my_read2.names();
    }
};

TEST_P(ProcessDataTester, PairedEnd) {
    auto param = GetParam();
    kaori::ProcessPairedEndDataOptions popt;
    popt.num_threads = std::get<0>(param);
    popt.block_size = std::get<1>(param);

    auto reads1 = simulate_reads(popt.num_threads + popt.block_size);
    auto reads2 = simulate_reads((popt.num_threads + popt.block_size) * 2);
    auto fastq_str1 = convert_to_fastq(reads1, "FOO");
    auto fastq_str2 = convert_to_fastq(reads2, "BAR");

    // Without names.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<false> task;
        kaori::process_paired_end_data(&reader1, &reader2, task, popt);
        EXPECT_EQ(task.first_reads(), reads1);
        EXPECT_EQ(task.second_reads(), reads2);
        EXPECT_TRUE(task.first_names().empty());
        EXPECT_TRUE(task.second_names().empty());
    }

    // Plus names.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<true> task;
        kaori::process_paired_end_data(&reader1, &reader2, task, popt);
        EXPECT_EQ(task.first_reads(), reads1);
        EXPECT_EQ(task.second_reads(), reads2);

        const auto& first_names = task.first_names();
        const auto& second_names = task.second_names();
        EXPECT_EQ(first_names.size(), reads1.size());
        EXPECT_EQ(second_names.size(), reads2.size());

        bool all_okay = true;
        for (size_t i = 0, end = first_names.size(); i < end; ++i) {
            if (first_names[i] != "FOO" + std::to_string(i + 1)) {
                all_okay = false;
            }
            if (second_names[i] != "BAR" + std::to_string(i + 1)) { 
                all_okay = false;
            }
        }
        EXPECT_TRUE(all_okay);
    }
}

TEST_P(ProcessDataTester, PairedEndErrors) {
    auto param = GetParam();
    kaori::ProcessPairedEndDataOptions popt;
    popt.num_threads = std::get<0>(param);
    popt.block_size = std::get<1>(param);

    auto reads1 = simulate_reads(popt.num_threads + popt.block_size);
    auto reads2 = simulate_reads((popt.num_threads + popt.block_size) * 2);
    auto fastq_str1 = convert_to_fastq(reads1, "FOO");
    auto fastq_str2 = convert_to_fastq(reads2, "BAR");

    // Errors out correctly due to the process.
    {
        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader2(reinterpret_cast<const unsigned char*>(fastq_str2.c_str()), fastq_str2.size());

        PairedEndCollector<false, true> task;
        EXPECT_ANY_THROW({
            try {
                kaori::process_paired_end_data(&reader1, &reader2, task, popt);
            } catch (std::exception& e) {
                EXPECT_TRUE(std::string(e.what()).find("I want some fries") != std::string::npos);
                throw e;
            }
        });
    }

    // Errors out correctly due to the read number.
    std::vector<size_t> size_options { 10, static_cast<size_t>(popt.block_size), reads2.size() - 1 };
    for (auto resized : size_options) {
        auto reads3 = reads2;
        reads3.resize(resized);
        auto fastq_str3 = convert_to_fastq(reads3, "FOO");

        byteme::RawBufferReader reader1(reinterpret_cast<const unsigned char*>(fastq_str1.c_str()), fastq_str1.size());
        byteme::RawBufferReader reader3(reinterpret_cast<const unsigned char*>(fastq_str3.c_str()), fastq_str3.size());

        PairedEndCollector<false> task;
        EXPECT_ANY_THROW({
            try {
                kaori::process_paired_end_data(&reader1, &reader3, task, popt);
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
