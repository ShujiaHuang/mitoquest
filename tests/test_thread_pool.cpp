#include <iostream>
#include <cstdlib>

#include <htslib/thread_pool.h>

// Task function
void* task_function(void* arg) {
    int task_id = *(int*)arg;
    std::cout << "Executing task " << task_id << std::endl;
    return nullptr;
}

int main() {
    // Create a thread pool with 2 threads
    hts_tpool* pool = hts_tpool_init(2);
    if (!pool) {
        std::cerr << "Failed to create thread pool" << std::endl;
        return 1;
    }

    // Create a process queue with a size of 3, no output queue needed
    hts_tpool_process* process = hts_tpool_process_init(pool, 3, 1);
    if (!process) {
        std::cerr << "Failed to create process queue" << std::endl;
        hts_tpool_destroy(pool);
        return 1;
    }

    // Submit tasks
    for (int i = 0; i < 3; ++i) {
        int* task_id = new int;
        *task_id = i;
        if (hts_tpool_dispatch(pool, process, task_function, task_id) != 0) {
            std::cerr << "Failed to submit task " << i << std::endl;
            delete task_id;
        }
    }

    // Wait for all tasks to complete (no results to fetch since in_only is 1)
    hts_tpool_process_flush(process);

    // Cleanup
    hts_tpool_process_destroy(process);
    hts_tpool_destroy(pool);

    std::cout << "All tasks completed" << std::endl;
    return 0;
}
