#ifndef KAORI_THREAD_POOL_HPP
#define KAORI_THREAD_POOL_HPP

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <stdexcept>

namespace kaori {

template<typename Run_, typename Finish_>
class ThreadPool {
public:
    ThreadPool(Run_ run_fun, Finish_ finish_fun, int num_threads) : 
        my_run_fun(std::move(run_fun)),
        my_finish_fun(std::move(finish_fun)),
        my_threads_available(num_threads)
    {
        my_threads.reserve(num_threads);
        for (int t = 0; t < num_threads; ++t) {
            my_threads.emplace_back([&]() -> void { thread_loop(); });
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard lck(my_lock);
            my_terminated = true;
            my_input_ready = true;
        }
        my_input_cv.notify_all();
        for (auto& thread : my_threads) {
            thread.join();
        }
    }

public:
    void submit(int job) {
        std::unique_lock lck(my_lock);
        my_output_cv.wait(lck, [&]() -> bool { return my_threads_available > 0; });

        check_rethrow();
        my_input_ready = true;
        my_input_job = job;
        --my_threads_available; // need to do this here to ensure join_all() doesn't immediately return.

        lck.unlock();
        my_input_cv.notify_one();
    }

    void join_all() {
        std::unique_lock lck(my_lock);
        my_output_cv.wait(lck, [&]() -> bool { return my_threads_available == my_threads.size(); });
        check_rethrow();
    }

private:
    Run_ my_run_fun;
    Finish_ my_finish_fun;
    std::vector<std::thread> my_threads;

    std::condition_variable my_input_cv, my_output_cv;
    std::mutex my_lock;
    bool my_input_ready = false;
    bool my_terminated = false;
    int my_input_job = 0;
    typename decltype(my_threads)::size_type my_threads_available;
    std::exception_ptr my_error;

    void thread_loop() {
        while (1) {
            int cur_job;
            {
                std::unique_lock lck(my_lock);
                my_input_cv.wait(lck, [&]() -> bool { return my_input_ready; });
                if (my_terminated) {
                    return;
                }

                my_input_ready = false;
                cur_job = my_input_job;
            }

            bool okay = true;
            try {
                my_run_fun(cur_job);
            } catch (std::exception&) {
                my_error = std::current_exception();
                okay = false;
            }

            {
                std::lock_guard lck(my_lock);
                if (okay) {
                    my_finish_fun(cur_job);
                }
                ++my_threads_available;
            }
            my_output_cv.notify_one();
        }
    }

    void check_rethrow() {
        if (my_error) {
            auto copy = my_error;
            my_error = nullptr; // in case we want to continue submitting jobs.
            std::rethrow_exception(std::move(copy));
        }
    }
};

}

#endif
