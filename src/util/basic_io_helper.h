//
// Created by cheyulin on 12/25/16.
//

#ifndef CODES_YCHE_BASIC_IO_HELPER_H
#define CODES_YCHE_BASIC_IO_HELPER_H

#include <vector>

namespace yche {
    using namespace std;

    auto Map2DArrWithDict(auto &arr_2d, auto &name_dict) {
        auto name_arr_2d = vector<vector<int>>();
        transform(arr_2d.begin(), arr_2d.end(), back_inserter(name_arr_2d), [&name_dict](auto &range) {
            auto tmp_vec = vector<int>();
            transform(std::begin(range), std::end(range), back_inserter(tmp_vec), [&name_dict](auto ele) {
                return name_dict[ele];
            });
            return tmp_vec;
        });
        return name_arr_2d;
    }
}

#endif //CODES_YCHE_BASIC_IO_HELPER_H
