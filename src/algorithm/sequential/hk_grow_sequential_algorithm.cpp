//
// Created by cheyulin on 1/5/17.
//

#include "hk_grow_sequential_algorithm.h"

namespace yche {
    HKGrow::HKGrow(unique_ptr<Graph> graph_ptr, double t, double eps) : t_(t), graph_ptr_(std::move(graph_ptr)) {
        taylor_deg_ = GetTaylorDegree(t, eps);
        psi_vec_ = ComputePsiVec(taylor_deg_, t);
        push_coefficient_vec_ = ComputePushCoefficientVec(taylor_deg_, eps, t, psi_vec_);
    }

    size_t HKGrow::GetTaylorDegree(double t, double eps) {
        auto eps_exp_t = eps * exp(t);
        auto error = exp(t) - 1;
        auto last = 1.0;
        auto k = 0u;
        while (error > eps_exp_t) {
            k++;
            last *= t / k;
            error -= last;
        }
        return max(k, 1u);
    }

    vector<double> HKGrow::ComputePsiVec(size_t taylor_deg, double t) {
        auto psi_vec = vector<double>(taylor_deg + 1, 0);
        psi_vec[taylor_deg] = 1;
        for (int k = 1; k <= taylor_deg; k++) {
            psi_vec[taylor_deg - k] = 1 + psi_vec[taylor_deg - k + 1] * t / (double) (taylor_deg - k + 1);
        }
        return psi_vec;
    }

    vector<double> HKGrow::ComputePushCoefficientVec(size_t taylor_deg, double eps, double t,
                                                     const vector<double> &psi_vec) {
        auto push_coefficient_vec = vector<double>(taylor_deg + 1, 0);
        push_coefficient_vec[0] = ((exp(t) * eps) / (double) taylor_deg) / psi_vec[0];
        for (int k = 1; k <= taylor_deg; k++) {
            push_coefficient_vec[k] = push_coefficient_vec[k - 1] * (psi_vec[k - 1] / psi_vec[k]);
        }
        return push_coefficient_vec;
    }

    size_t HKGrow::DiffuseWeight(const SpareseVec &seed_dict, SpareseVec &x_dict, size_t max_push_count) const {
        auto task_queue = queue<pair<size_t, size_t>>();
        auto ri = 0ul;
        auto rij = 0.0;
        auto push_num = 0ul;
        auto r_dict = SpareseVec();

#define rentry(i, j, n) ((i)+(j)*(n))
        for (auto &ele: seed_dict) {
            tie(ri, rij) = ele;
            r_dict[rentry(ri, 0, num_vertices(*graph_ptr_))] += rij;
            task_queue.emplace(ri, 0);
        }

        while (push_num < max_push_count) {
            size_t i, j;
            tie(i, j) = task_queue.front();
            task_queue.pop();

            auto deg_of_i = out_degree(i, *graph_ptr_);
            rij = r_dict[ri];
            x_dict[i] += rij;
            r_dict[ri] = 0;

            auto rijs = t_ * rij / (j + 1);
            auto ajv = 1.0 / deg_of_i;
            auto update = rijs * ajv;

            if (j == taylor_deg_ - 1) {
                for (auto vp = adjacent_vertices(i, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                    auto dst_v = *vp.first;
                    x_dict[dst_v] += update;
                }
            } else {
                for (auto vp = adjacent_vertices(i, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                    auto dst_v = *vp.first;
                    auto re = rentry(dst_v, j + 1, num_vertices(*graph_ptr_));
                    auto re_old = r_dict[re];
                    auto re_new = re_old + update;
                    r_dict[re] = re_new;
                    if (re_new >= out_degree(dst_v, *graph_ptr_) * push_coefficient_vec_[j + 1] &&
                        re_old < out_degree(dst_v, *graph_ptr_) * push_coefficient_vec_[j + 1]) {
                        task_queue.emplace(dst_v, j + 1);
                    }
                }
            }
            push_num += deg_of_i;
            if (task_queue.size() == 0) { return push_num; }
        }
#undef rentry
        return push_num;
    }

    auto HKGrow::SweepCut(SpareseVec &x_dict) const {
        auto cluster = vector<size_t>();
        for (auto &ele:x_dict) { ele.second *= (1.0 / max(out_degree(ele.first, *graph_ptr_), 1ul)); }
        auto pr_pairs = vector<pair<size_t, double >>(x_dict.begin(), x_dict.end());
        sort(pr_pairs.begin(), pr_pairs.end(),
             [](auto &&left, auto &&right) { return left.second > right.second; });

        auto cond_vec = vector<double>(pr_pairs.size());
        auto vol_vec = vector<size_t>(pr_pairs.size());
        auto cut_vec = vector<size_t>(pr_pairs.size());
        auto rank_map = unordered_map<size_t, size_t>();
        for (auto i = 0ul; i < pr_pairs.size(); i++) { rank_map[pr_pairs[i].first] = i; }

        //since undirectedS is not implemented for CSR, each undirected edge is counted twice
        auto vol_of_graph = num_edges(*graph_ptr_) / 2;
        auto vol_of_set = 0ul;
        auto cut_of_set = 0ul;

        for (auto i = 0ul; i < pr_pairs.size(); i++) {
            auto v = pr_pairs[i].first;
            auto deg = out_degree(v, *graph_ptr_);
            auto change = deg;

            for (auto vp = adjacent_vertices(i, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                auto neighbor_v = *vp.first;
                if (rank_map.count(neighbor_v) > 0 && rank_map[neighbor_v] < rank_map[v]) { change -= 2; }
            }

            cut_of_set += change;
            vol_of_set += deg;
            vol_vec[i] = vol_of_set;
            cut_vec[i] = cut_of_set;
            cond_vec[i] = (vol_of_set == 0 || vol_of_graph - vol_of_set == 0 ? 1 :
                           static_cast<double>(cut_of_set) / min(vol_of_set, vol_of_graph - vol_of_set));
        }

        transform(begin(pr_pairs), end(pr_pairs), back_inserter(cluster), [](auto &&ele) { return ele.first; });

        auto min_iter = min_element(begin(cond_vec), end(cond_vec));
        auto min_cond = (min_iter == end(cond_vec) ? 0.0 : *min_iter);
        auto min_cond_idx = min_iter - begin(cond_vec);
        auto status = LocalSweepCutStatus();
        status.conductance_ = min_cond;
        status.cut_ = cut_vec[min_cond_idx];
        status.volume_ = vol_vec[min_cond_idx];

        sort(begin(cluster), end(cluster));
        return make_tuple(status, cluster);
    }

    auto HKGrow::HyperCluster(const vector<size_t> &seed_set) const {
        auto seed_dict = SpareseVec();
        auto x_dict = SpareseVec();
        for (auto &seed:seed_set) { seed_dict.emplace(seed, 1.0 / seed_set.size()); }

        auto seed_iter = max_element(begin(seed_set), end(seed_set), [&](auto &&left, auto &&right) {
            return out_degree(left, *graph_ptr_) < out_degree(right, *graph_ptr_);
        });
        auto max_deg = (seed_iter == seed_set.end() ? 0ul : out_degree(*seed_iter, *graph_ptr_));

        auto step_num = DiffuseWeight(seed_dict, x_dict, static_cast<size_t >(ceil(
                pow(static_cast<double >(num_vertices(*graph_ptr_)), 1.5))));

        auto status_cluster = SweepCut(x_dict);
        auto &status = std::get<0>(status_cluster);
        status.steps_ = step_num;
        status.support_ = seed_dict.size();
        if (step_num == 0) { x_dict = seed_dict; }

        return status_cluster;
    }

    double HKGrow::GetIntersectRatio(const vector<size_t> &left_community, const vector<size_t> &right_community) {
        auto intersect_set = vector<int>(left_community.size() + right_community.size());
        auto iter_end = set_intersection(left_community.begin(), left_community.end(),
                                         right_community.begin(), right_community.end(), intersect_set.begin());
        auto intersect_set_size = iter_end - intersect_set.begin();
        auto rate = static_cast<double>(intersect_set_size) / min(left_community.size(), right_community.size());
        return rate;
    }

    vector<size_t> HKGrow::GetUnion(const vector<size_t> &left_community, const vector<size_t> &right_community) {
        auto union_set = vector<size_t>(left_community.size() + right_community.size());
        auto iter_end = set_union(left_community.begin(), left_community.end(),
                                  right_community.begin(), right_community.end(), union_set.begin());
        union_set.resize(static_cast<size_t >(iter_end - union_set.begin()));
        return union_set;
    }

    void HKGrow::MergeCommToGlobal(vector<size_t> &result_community) {
        if (overlap_community_vec_.size() == 0) {
            overlap_community_vec_.emplace_back(std::move(result_community));
        } else {
            auto is_insert = true;
            for (auto &community:overlap_community_vec_) {
                if (GetIntersectRatio(community, result_community) > 1 - DOUBLE_ACCURACY) {
                    community = GetUnion(community, result_community);
                    is_insert = false;
                    break;
                }
            }
            if (is_insert) { overlap_community_vec_.emplace_back(std::move(result_community)); }
        }
    }

    HKGrow::CommunityVec HKGrow::ExecuteHRGRow() {
        auto seeds_vec = vector<vector<size_t>>();
        seeds_vec.reserve(num_vertices(*graph_ptr_));
        auto vp = vertices(*graph_ptr_);
        transform(vp.first, vp.second, back_inserter(seeds_vec), [](size_t val) { return vector<size_t>(1, val); });

        for (auto &seeds:seeds_vec) {
            auto stats_cluster = HyperCluster(seeds);
            MergeCommToGlobal(std::get<1>(stats_cluster));
        }
        return overlap_community_vec_;
    }
}
