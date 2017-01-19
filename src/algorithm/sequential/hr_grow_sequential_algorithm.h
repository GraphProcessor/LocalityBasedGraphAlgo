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
#include <memory>
#include <unordered_set>
#include <unordered_map>

namespace yche {
    using namespace std;

    using SpareseVec=unordered_map<size_t, double>;

    struct SparseRow {
        size_t n_, m_;
        size_t *vertices_;
        size_t *edges_;
        double *weight_;

        size_t sr_degree(size_t u) {
            return vertices_[u + 1] - vertices_[u];
        }
    };

    struct LocalSweepCutStatus {
        double conductance_;
        double volume_;
        double support_;
        double steps_;
        double cut_;
    };

    class HKGrow {
    public:
        HKGrow(unique_ptr<SparseRow> graph_ptr, double t, double eps);

        void ExecuteHRGRow(vector<size_t> &seeds, double &f_cond, double &f_cut,
                           double &f_vol, SpareseVec &x_dict, double &num_push);

    private:
        unique_ptr<SparseRow> graph_ptr_;
        double t_;
        size_t taylor_deg_;
        vector<double> psi_vec_;
        vector<double> push_coefficient_vec_;

        static size_t GetTaylorDegree(double t, double eps);

        static vector<double> ComputePsiVec(size_t taylor_deg, double t);

        static vector<double> ComputePushCoefficientVec(size_t taylor_deg, double eps, double t,
                                                        const vector<double> &psi_vec);

        size_t DiffuseWeight(const SpareseVec &seed_dict, SpareseVec &x_dict, size_t max_push_count) const;

        LocalSweepCutStatus SweepCut(SpareseVec &x_dict, vector<size_t> &cluster) const;

        LocalSweepCutStatus HyperCluster(const vector<size_t> &seed_set, SpareseVec &x_dict,
                                         vector<size_t> &cluster) const;

    };


}
#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
