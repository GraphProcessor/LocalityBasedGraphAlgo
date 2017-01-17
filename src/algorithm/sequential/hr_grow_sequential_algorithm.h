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

    struct sparse_vec {
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
        double eps;
        double cut;
    };

    class HKGrow {
    public:
        size_t get_taylor_degree(double t, double eps);

        vector<double> compute_psi_vec(size_t taylor_deg, double t);

        vector<double> compute_threshold_vec(size_t taylor_deg, double eps, double t, vector<double> psi_vec);

        size_t gs_qexpm_seed(sparse_row &graph, sparse_vec &set, sparse_vec &y, double t, double eps,
                             size_t max_push_count);

        void cluster_from_sweep(sparse_row &G, sparse_vec &p, vector<size_t> &cluster, double *out_cond,
                                double *out_volume, double *out_cut);

        int hyper_cluster_heat_kernel_multiple(sparse_row &G, const vector<size_t> &seed_set, double t, double eps,
                                               sparse_vec &p, sparse_vec &r, vector<size_t> &cluster,
                                               local_hkpr_stats *stats);

        void hk_grow(sparse_row &G, vector<size_t> &seeds, double t, double eps, double &fcond, double &fcut,
                     double &fvol, sparse_vec &p, double &pushes);
    };


}
#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
