//
// Created by cheyulin on 8/9/16.
//

#include <iostream>
#include "deprecated/parallel_utils/detail/thread_pool_breakable.h"

using namespace yche;
using namespace std;

int main() {
    ThreadPoolBreakable breakable_pool(20);
    bool is_break = false;
    for (auto j = 0; j < 3000; j++) {
        cout << "Round:" << j << endl;
        auto integer = 0;
        for (auto i = 0; i < 5000; i++) {
            integer++;
            std::function<BreakWithCallBackRetType(void)> task_function = [integer]() -> BreakWithCallBackRetType {
                if (integer == 10) {
//                    return BreakWithCallBackRetType(true, []() { cout << "Cur Break" << endl; });
                    return BreakWithCallBackRetType();
                } else
                    return BreakWithCallBackRetType();
            };
            breakable_pool.AddTask(task_function);
        }
        cout << "Finish Add" << endl;
        breakable_pool.WaitForBreakOrTerminate(is_break);
    }
}