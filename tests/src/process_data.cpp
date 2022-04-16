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
            reads.emplace_back(x.first, x.second);
            names.emplace_back(x.first, x.second);
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

    std::string convert_to_fastq(const std::vector<std::string>& reads) {
        std::string output;

        for (size_t i = 0; i < reads.size(); ++i) {
            output += "READ" + std::to_string(i+1) + "\n";
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
        EXPECT_TRUE(task.collected_reads.empty());
    }
}

INSTANTIATE_TEST_SUITE_P(
    ProcessData,
    ProcessDataTester, 
    ::testing::Combine(
        ::testing::Values(1, 2, 3),
        ::testing::Values(10, 50, 100)
    )
);
