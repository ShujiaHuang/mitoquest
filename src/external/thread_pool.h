/**
 * A simple C++11 Thread Pool implementation.
 * 
 * This implementation is based on the example provided in the following link:
 * My Poe: https://poe.com/chat/35dqbwaajb24b8wotri 
 * 
 * Author: Shujia Huang (hshujia@qq.com)
 * Date: 2025-02-11
 * Version: 1.0.0 
 * 
 */
#ifndef __INCLUDE_THREAD_POOL_H__
#define __INCLUDE_THREAD_POOL_H__

#include <thread>
#include <queue>
#include <functional>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <future>
#include <memory>

class ThreadPool {
public:
    // 构造函数，初始化线程池
    explicit ThreadPool(size_t threadCount = std::thread::hardware_concurrency()): stop(false) {
        for (size_t i = 0; i < threadCount; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        // 等待直到有任务或需要停止
                        condition.wait(lock, [this] {
                            return this->stop || !this->tasks.empty();
                        });
                        
                        // 如果线程池停止且任务队列为空，退出线程
                        if (this->stop && this->tasks.empty()) {
                            return;
                        }
                        
                        // 获取任务队列中的第一个任务
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    } // ← 锁在此释放，{} 是必须的。

                    // 执行任务。submit() 用 packaged_task 包装，任务抛出的异常
                    // 会被捕获并存入 future，由调用方 f.get() 重抛 —— 这里的
                    // try/catch 对 packaged 任务永远不会触发（已移除旧的死分支）。
                    task();
                }
            });
        }
    }
    
    // 析构函数，确保所有线程正确退出
    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }  // 必须先释放锁，否则 worker 无法获取锁来检查 stop 标志 → 死锁，这里的 {} 是必须的。

        condition.notify_all();
        for (std::thread& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }
    
    // 提交任务到线程池，返回future对象
    template<class F, class... Args>
    auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        // 获取函数返回值类型
        using return_type = decltype(f(args...));
        
        // 创建任务包装器
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
        // 获取future对象
        std::future<return_type> result = task->get_future();
        
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            
            // 如果线程池已停止，抛出异常
            if (stop) {
                throw std::runtime_error("Cannot submit task to stopped ThreadPool");
            }
            
            // 将任务添加到队列
            tasks.emplace([task]() { (*task)(); });
        } // ← 锁在此释放
        
        // 通知一个等待中的线程
        condition.notify_one();
        return result;
    }
    
    // 获取当前等待执行的任务数量
    size_t get_task_count() const {
        std::unique_lock<std::mutex> lock(queueMutex);
        return tasks.size();
    }
    
    // 获取线程池中的线程数量
    size_t get_thread_count() const {
        return workers.size();
    }

private:
    std::vector<std::thread> workers;                // 工作线程容器
    std::queue<std::function<void()>> tasks;         // 任务队列
    mutable std::mutex queueMutex;                   // 互斥锁
    std::condition_variable condition;               // 条件变量
    bool stop;                                       // 停止标志
};

#endif
