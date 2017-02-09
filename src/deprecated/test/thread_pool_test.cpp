//
// Created by cheyulin on 8/8/16.
//
#include <thread>
#include "deprecated/parallel_utils/detail/thread_pool_base.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
    using namespace yche;
    ThreadPoolBase<int> pool(20);
    for (auto j = 0; j < 3000; j++) {
        cout << "Round:" << j << endl;
        for (int i = 0; i < 5000; ++i) {
            pool.AddTask([i]() {
                return i * i;
            });
        }
        cout << "Finish Add" << endl;
        pool.WaitAll();
    }
}