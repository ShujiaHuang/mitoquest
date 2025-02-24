#include <iostream>
#include <vector>
#include <chrono>

#include "external/thread_pool.h"

int main(int argc, char *argv[]) {
    
    ThreadPool pool(4);
    std::vector< std::future<int> > results;

    for(int i = 0; i < 8; ++i) {
        results.emplace_back(
            pool.submit([i] {
                std::cout << "hello " << i << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(1));
                std::cout << "world " << i << std::endl;
                return i*i;
            })
        );
    }

    for(auto && result: results) {
        if (result.valid()) std::cout << result.get() << ' ';
    }
    std::cout << std::endl;
    
    return 0;
}
