//
// Created by cheyulin on 1/5/17.
//

#include "hr_grow_sequential_algorithm.h"

namespace yche {
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

    vector<double> HKGrow::ComputePushCOefficientVec(size_t taylor_deg, double eps, double t, vector<double> &psi_vec) {
        auto push_coefficient_vec = vector<double>(taylor_deg + 1, 0);
        push_coefficient_vec[0] = ((exp(t) * eps) / (double) taylor_deg) / psi_vec[0];
        for (int k = 1; k <= taylor_deg; k++) {
            push_coefficient_vec[k] = push_coefficient_vec[k - 1] * (psi_vec[k - 1] / psi_vec[k]);
        }
        return push_coefficient_vec;
    }

    size_t
    HKGrow::ExpandSeed(sparse_row &graph, dict_wrapper &residual_dict, dict_wrapper &x_dict, double t, double eps,
                       size_t max_push_count) {
        auto task_queue = queue<pair<size_t, size_t>>();
        auto taylor_deg = GetTaylorDegree(t, eps);
        auto psi_vec = ComputePsiVec(taylor_deg, t);
        auto push_coefficient_vec = ComputePushCOefficientVec(taylor_deg, eps, t, psi_vec);

        auto ri = 0ul;
        auto rij = 0.0;
        auto push_num = 0ul;
        dict_wrapper residual_vec;

#define rentry(i, j, n) ((i)+(j)*(n))
        for (auto &ele: residual_dict.weight_map_) {
            tie(ri, rij) = ele;
            residual_vec.weight_map_[rentry(ri, 0, graph.n_)] += rij;
            task_queue.emplace(ri, 0);
        }

        while (push_num < max_push_count) {
            size_t i, j;
            tie(i, j) = task_queue.front();
            task_queue.pop();

            auto deg_of_i = graph.sr_degree(i);
            rij = residual_vec.weight_map_[ri];
            x_dict.weight_map_[i] += rij;
            residual_vec.weight_map_[ri] = 0;

            auto rijs = t * rij / (j + 1);
            auto ajv = 1.0 / deg_of_i;
            auto update = rijs * ajv;

            if (j == taylor_deg - 1) {
                for (auto nzi = graph.vertices_[i]; nzi < graph.vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph.edges_[nzi];
                    x_dict.weight_map_[dst_v] += update;
                }
            } else {
                for (auto nzi = graph.vertices_[i]; nzi < graph.vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph.edges_[nzi];
                    auto re = rentry(dst_v, j + 1, graph.n_);
                    auto re_old = residual_vec.get(re);
                    auto re_new = re_old + update;
                    double dv = graph.sr_degree(dst_v);
                    residual_vec.weight_map_[re] = re_new;
                    if (re_new >= dv * push_coefficient_vec[j + 1] && re_old < dv * push_coefficient_vec[j + 1]) {
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

    void HKGrow::SweepCut(sparse_row &G, dict_wrapper &x_dict, vector<size_t> &cluster,
                          double *out_cond, double *out_volume, double *out_cut) {
        auto pr_pairs = vector<pair<size_t, double>>(x_dict.weight_map_.begin(), x_dict.weight_map_.end());
        sort(pr_pairs.begin(), pr_pairs.end(), [](auto &&left, auto &&right) { return left.second > right.second; });

        auto conductance_vec = vector<double>(pr_pairs.size());
        auto volume_vec = vector<size_t>(pr_pairs.size());
        auto cut_size_vec = vector<size_t>(pr_pairs.size());
        auto rank_map = unordered_map<size_t, size_t>();
        for (auto i = 0ul; i < pr_pairs.size(); i++) { rank_map[pr_pairs[i].first] = i; }

        auto total_degree = G.vertices_[G.m_];
        auto cur_volume = 0ul;
        auto cur_cut_size = 0ul;

        for (auto i = 0ul; i < pr_pairs.size(); i++) {
            auto v = pr_pairs[i].first;
            auto deg = G.vertices_[v + 1] - G.vertices_[v];
            auto change = deg;
            for (auto nzi = G.vertices_[v]; nzi < G.vertices_[v + 1]; ++nzi) {
                auto neighbor_v = G.edges_[nzi];
                if (rank_map.count(neighbor_v) > 0 && rank_map[neighbor_v] < rank_map[v]) {
                    change -= 2;
                }
            }

            cur_cut_size += change;
            cur_volume += deg;
            volume_vec[i] = cur_volume;
            cut_size_vec[i] = cur_cut_size;
            conductance_vec[i] = (cur_volume == 0 || total_degree - cur_volume == 0 ? 1 :
                                  static_cast<double>(cur_cut_size) / min(cur_volume, total_degree - cur_volume));
        }

        transform(begin(pr_pairs), end(pr_pairs), back_inserter(cluster), [](auto &&ele) { return ele.first; });

        auto min_iter = min_element(begin(conductance_vec), end(conductance_vec));
        auto min_cond = (min_iter == end(conductance_vec) ? 0.0 : *min_iter);
        auto min_cond_idx = min_iter - begin(conductance_vec);

        if (out_cond) { *out_cond = min_cond; }
        if (out_volume) { *out_volume = volume_vec[min_cond_idx]; }
        if (out_cut) { *out_cut = cut_size_vec[min_cond_idx]; }

    }


    int HKGrow::HyperCluster(sparse_row &G, const vector<size_t> &seed_set, double t, double eps,
                             dict_wrapper &x_dict, dict_wrapper &residual_dict, vector<size_t> &cluster,
                             local_hkpr_stats *stats) {
        auto max_deg = 0ul;
        for (auto i = 0; i < seed_set.size(); ++i) {
            auto seed_vertex = seed_set[i];
            auto v_degree = G.sr_degree(seed_vertex);
            residual_dict.weight_map_[seed_vertex] = 1.0 / seed_set.size();
            max_deg = max(max_deg, v_degree);
        }

        auto step_num = ExpandSeed(G, residual_dict, x_dict, t, eps, static_cast<size_t >(ceil(pow(G.n_, 1.5))));

        if (stats) { stats->steps = step_num; }
        if (stats) { stats->support = residual_dict.weight_map_.size(); }
        double *out_cond = nullptr;
        double *out_volume = nullptr;
        double *out_cut = nullptr;
        if (stats) { out_cond = &stats->conductance; }
        if (stats) { out_volume = &stats->volume; }
        if (stats) { out_cut = &stats->cut; }

        if (step_num == 0) { x_dict = residual_dict; }
        for (auto &ele:x_dict.weight_map_) {
            ele.second *= (1.0 / max(G.sr_degree(ele.first), static_cast<size_t >(1)));
        }
        SweepCut(G, x_dict, cluster, out_cond, out_volume, out_cut);
        return 0;
    }

    void HKGrow::ExecuteHRGRow(sparse_row &G, vector<size_t> &seeds, double t, double eps, double &f_cond,
                               double &f_cut, double &f_vol, dict_wrapper &p, double &num_push) {
        dict_wrapper r;
        vector<size_t> best_cluster_vec;
        local_hkpr_stats stats;

        HyperCluster(G, seeds, t, eps, p, r, best_cluster_vec, &stats);
        seeds = best_cluster_vec;

        num_push = stats.steps;
        f_cond = stats.conductance;
        f_cut = stats.cut;
        f_vol = stats.volume;
    }
}
