//
// Created by cheyulin on 1/5/17.
//

#ifndef CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
#define CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H

#include <cmath>
#include <cassert>

#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace yche {
    using namespace std;

    struct dict_wrapper {
        unordered_map<size_t, double> weight_map_;

        double get(size_t index, double default_value = 0.0) {
            auto it = weight_map_.find(index);
            return it == weight_map_.end() ? default_value : it->second;
        }

        double sum() {
            return accumulate(weight_map_.begin(), weight_map_.end(), 0,
                              [](auto &&left, auto &&right) { return left + right.second; });
        }

        size_t max_index() {
            return weight_map_.size() == 0 ? 0 :
                   max_element(weight_map_.begin(), weight_map_.end(),
                               [](auto &&left, auto &&right) { return left.second < right.second; })->first;
        }
    };

    struct sparse_row {
        size_t n_, m_;
        size_t *vertices_;
        size_t *edges_;
        double *weight_;

        size_t sr_degree(size_t u) {
            return vertices_[u + 1] - vertices_[u];
        }
    };

    struct local_hkpr_stats {
        double conductance;
        double volume;
        double support;
        double steps;
        double cut;
    };

    class HKGrow {
    public:
        size_t ExpandSeed(sparse_row &graph, dict_wrapper &residual_dict, dict_wrapper &x_dict,
                          double t, double eps, size_t max_push_count);

        void SweepCut(sparse_row &G, dict_wrapper &x_dict, vector<size_t> &cluster,
                      double *out_cond, double *out_volume, double *out_cut);

        int HyperCluster(sparse_row &G, const vector<size_t> &seed_set, double t, double eps,
                         dict_wrapper &x_dict, dict_wrapper &residual_dict, vector<size_t> &cluster,
                         local_hkpr_stats *stats);

        void ExecuteHRGRow(sparse_row &G, vector<size_t> &seeds, double t, double eps, double &f_cond, double &f_cut,
                           double &f_vol, dict_wrapper &p, double &num_push);

    private:
        size_t GetTaylorDegree(double t, double eps);

        vector<double> ComputePsiVec(size_t taylor_deg, double t);

        vector<double> ComputePushCOefficientVec(size_t taylor_deg, double eps, double t, vector<double> &psi_vec);
    };


}
#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
