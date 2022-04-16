#ifndef KAORI_PROCESS_DATA_HPP
#define KAORI_PROCESS_DATA_HPP

#include <thread>
#include "FastqReader.hpp"
#include "byteme/Reader.hpp"

namespace kaori {

template<bool use_names>
struct ChunkOfReads {
    ChunkOfReads() : sequence_offset(1), name_offset(1) {} // zero is always the first element.

    void clear() {
        sequence_buffer.clear();
        sequence_offset.resize(1);
        if constexpr(use_names) {
            name_buffer.clear();
            name_offset.resize(1);
        }
    }

    void add_read_sequence(const std::vector<char>& sequence) {
        add_read_details(sequence, sequence_buffer, sequence_offset);
    }

    void add_read_name(const std::vector<char>& name) {
        add_read_details(name, name_buffer, name_offset);
    }

    size_t size() const {
        return sequence_offset.size() - 1;
    }

    std::pair<const char*, const char*> get_sequence(size_t i) const {
        return get_details(i, sequence_buffer, sequence_offset);
    }

    std::pair<const char*, const char*> get_name(size_t i) const {
        return get_details(i, name_buffer, name_offset);
    }

private:
    std::vector<char> sequence_buffer;
    std::vector<size_t> sequence_offset;
    std::vector<char> name_buffer;
    std::vector<size_t> name_offset;

    static void add_read_details(const std::vector<char>& src, std::vector<char>& dst, std::vector<size_t>& offset) {
        dst.insert(dst.end(), src.begin(), src.end());
        auto last = offset.back();
        offset.push_back(last + src.size());
    }

    static std::pair<const char*, const char*> get_details(size_t i, const std::vector<char>& dest, const std::vector<size_t>& offset) {
        const char * base = dest.data();
        return std::make_pair(base + offset[i], base + offset[i + 1]);
    }
};

template<class Task>
void process_single_end_data(byteme::Reader* input, Task& task, int num_threads = 1, int block_size = 100000) {
    FastqReader fastq(input);
    bool finished = false;

    std::vector<ChunkOfReads<Task::use_names> > reads(num_threads);
    std::vector<std::thread> jobs(num_threads);
    std::vector<decltype(task.initialize())> states(num_threads);

    while (!finished) {
        for (int t = 0; t < num_threads; ++t) {
            if (jobs[t].joinable()) {
                jobs[t].join();
                task.reduce(states[t]);
                reads[t].clear();
            }

            auto& curreads = reads[t];
            for (int b = 0; b < block_size; ++b) {
                if (!fastq()) {
                    finished = true;
                    break;
                }

                curreads.add_read_sequence(fastq.get_sequence());
                if constexpr(Task::use_names) {
                    curreads.add_read_name(fastq.get_name());
                }
            }

            states[t] = task.initialize();
            jobs[t] = std::thread([&](int i) -> void {
                auto& state = states[i];
                const auto& curreads = reads[i];
                size_t nreads = curreads.size();

                if constexpr(!Task::use_names) {
                    for (size_t b = 0; b < nreads; ++b) {
                        state.process(curreads.get_sequence(b));
                    }
                } else {
                    for (size_t b = 0; b < nreads; ++b) {
                        state.process(curreads.get_name(b), curreads.get_sequence(b));
                    }
                }
            }, t);

            if (finished) {
                // We won't get a future iteration to join the previous threads
                // that we kicked off, so we do it now.
                for (int u = 0; u < num_threads; ++u) {
                    auto pos = (u + t + 1) % num_threads;
                    if (jobs[pos].joinable()) {
                        jobs[pos].join();
                        task.reduce(states[pos]);
                    }
                }
                break;
            }
        }
    }

    return;
}

template<class Task>
void process_paired_end_data(byteme::Reader* input1, byteme::Reader* input2, Task& task, int num_threads = 1, int block_size = 100000) {
    FastqReader fastq1(input1);
    FastqReader fastq2(input2);
    bool finished = false;

    std::vector<ChunkOfReads<Task::use_names> > reads1(num_threads), reads2(num_threads);
    std::vector<std::thread> jobs(num_threads);
    std::vector<decltype(task.initialize())> states(num_threads);

    while (!finished) {
        for (int t = 0; t < num_threads; ++t) {
            if (jobs[t].joinable()) {
                jobs[t].join();
                task.reduce(states[t]);
                reads1[t].clear();
                reads2[t].clear();
            }

            bool finished1 = false;
            {
                auto& curreads = reads1[t];
                for (int b = 0; b < block_size; ++b) {
                    if (!fastq1()) {
                        finished1 = true;
                        break;
                    }

                    curreads.add_read_sequence(fastq1.get_sequence());
                    if constexpr(Task::use_names) {
                        curreads.add_read_name(fastq1.get_name());
                    }
                }
            }

            bool finished2 = false;
            {
                auto& curreads = reads2[t];
                for (int b = 0; b < block_size; ++b) {
                    if (!fastq2()) {
                        finished2 = true;
                        break;
                    }

                    curreads.add_read_sequence(fastq2.get_sequence());
                    if constexpr(Task::use_names) {
                        curreads.add_read_name(fastq2.get_name());
                    }
                }
            }

            if (finished1 != finished2) {
                // Wait for all remaining threads to join before throwing,
                // otherwise we get weird errors.
                for (int u = 0; u < num_threads; ++u) {
                    if (jobs[u].joinable()) {
                        jobs[u].join();
                    }
                }
                throw std::runtime_error("different number of reads in paired data");
            } else if (finished1) {
                finished = true;
            }

            states[t] = task.initialize();
            jobs[t] = std::thread([&](int i) -> void {
                auto& state = states[i];
                const auto& curreads1 = reads1[i];
                const auto& curreads2 = reads2[i];
                size_t nreads = curreads1.size();

                if constexpr(!Task::use_names) {
                    for (size_t b = 0; b < nreads; ++b) {
                        state.process(curreads1.get_sequence(b), curreads2.get_sequence(b));
                    }
                } else {
                    for (size_t b = 0; b < nreads; ++b) {
                        state.process(
                            curreads1.get_name(b), 
                            curreads1.get_sequence(b),
                            curreads2.get_name(b), 
                            curreads2.get_sequence(b)
                        );
                    }
                }
            }, t);

            if (finished) {
                // We won't get a future iteration to join the previous threads
                // that we kicked off, so we do it now.
                for (int u = 0; u < num_threads; ++u) {
                    auto pos = (u + t + 1) % num_threads;
                    if (jobs[pos].joinable()) {
                        jobs[pos].join();
                        task.reduce(states[pos]);
                    }
                }
                break;
            }
        }
    }

    return;
}


}

#endif
