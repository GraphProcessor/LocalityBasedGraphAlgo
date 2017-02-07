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

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>

namespace yche {
    using namespace std;
    using namespace boost;
    using SpareseVec=unordered_map<size_t, double>;
    constexpr double DOUBLE_ACCURACY = 0.00001;

    struct LocalSweepCutStatus {
        double conductance_;
        double volume_;
        double support_;
        double steps_;
        double cut_;
    };

    class HKGrow {
    private:
        using CommunityVec=vector<vector<size_t>>;

    public:
        using Graph=compressed_sparse_row_graph<>;

        HKGrow(unique_ptr<Graph> graph_ptr, double t, double eps);

        CommunityVec ExecuteHRGRow();

    private:
        CommunityVec overlap_community_vec_;
        unique_ptr<Graph> graph_ptr_;
        double t_;
        size_t taylor_deg_;
        vector<double> psi_vec_;
        vector<double> push_coefficient_vec_;

        static double GetIntersectRatio(const vector<size_t> &left_community, const vector<size_t> &right_community);

        static vector<size_t> GetUnion(const vector<size_t> &left_community, const vector<size_t> &right_community);

        static size_t GetTaylorDegree(double t, double eps);

        static vector<double> ComputePsiVec(size_t taylor_deg, double t);

        static vector<double> ComputePushCoefficientVec(size_t taylor_deg, double eps, double t,
                                                        const vector<double> &psi_vec);

        size_t DiffuseWeight(const SpareseVec &seed_dict, SpareseVec &x_dict, size_t max_push_count) const;

        auto SweepCut(SpareseVec &x_dict) const;

        auto HyperCluster(const vector<size_t> &seed_set) const;

        void MergeCommToGlobal(vector<size_t> &result_community);
    };
}
#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
