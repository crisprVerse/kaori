#include <gtest/gtest.h>
#include "kaori/process_data.hpp"
#include <random>
#include "byteme/RawBufferReader.hpp"

template<bool unames>
class SingleEndCollector {
public:
    struct State {
        std::vector<std::string> reads, names;

        void process(const std::pair<const char*, const char*>& x) {
            reads.emplace_back(x.first, x.second);
        }

        void process(const std::pair<const char*, const char*>& x, const std::pair<const char*, const char*>& y) {
            names.emplace_back(x.first, x.second);
            reads.emplace_back(y.first, y.second);
        }
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

    std::string convert_to_fastq(const std::vector<std::string>& reads, std::string prefix = "READ") {
        std::string output;

        for (size_t i = 0; i < reads.size(); ++i) {
            output += "@" + prefix + std::to_string(i+1) + "\n";
            output += reads[i] + "\n";
            output += "+\n" + std::string(reads[i].size(), '!') + "\n";
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

template<bool unames>
class PairedEndCollector {
public:
    SingleEndCollector<unames> read1, read2;

    struct State {
        typename SingleEndCollector<unames>::State read1, read2;

        void process(const std::pair<const char*, const char*>& x1, const std::pair<const char*, const char*>& x2) {
            read1.process(x1);
            read2.process(x2);
        }

        void process(const std::pair<const char*, const char*>& x1, const std::pair<const char*, const char*>& y1,
                     const std::pair<const char*, const char*>& x2, const std::pair<const char*, const char*>& y2) {
            read1.process(x1, y1);
            read2.process(x2, y2);
        }
    };

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

    // Errors out correctly.
    {
        auto reads3 = reads2;
        reads3.resize(10);
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
