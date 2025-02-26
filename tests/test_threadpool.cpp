#include <iostream>
#include <vector>
#include <chrono>

#include "external/thread_pool.h"

int main(int argc, char *argv[]) {
    
    ThreadPool pool(4);
    std::vector< std::future<int> > results;

    for(int i = 0; i < 20; ++i) {
        results.emplace_back(
            pool.submit([i] {
                //std::cout << "- hello " << i << std::endl;
                //std::this_thread::sleep_for(std::chrono::seconds(1));
                return i*i;
            })
        );
    }
/*
    // 收集所有任务的结果到vector中
    std::vector<int> output;
    for (auto& result : results) {
        output.push_back(result.get());  // 获取每个任务的返回值
    }

    // 输出结果
    std::cout << "Results: ";
    for (int value : output) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
*/

    std::cout << "Results: ";
    for(auto && result: results) {
        if (result.valid()) std::cout << result.get() << " ";
    }
    std::cout << std::endl;
    
    return 0;
}
