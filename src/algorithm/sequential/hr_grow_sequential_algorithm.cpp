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
            r_dict[rentry(ri, 0, graph_ptr_->n_)] += rij;
            task_queue.emplace(ri, 0);
        }

        while (push_num < max_push_count) {
            size_t i, j;
            tie(i, j) = task_queue.front();
            task_queue.pop();

            auto deg_of_i = graph_ptr_->sr_degree(i);
            rij = r_dict[ri];
            x_dict[i] += rij;
            r_dict[ri] = 0;

            auto rijs = t_ * rij / (j + 1);
            auto ajv = 1.0 / deg_of_i;
            auto update = rijs * ajv;

            if (j == taylor_deg_ - 1) {
                for (auto nzi = graph_ptr_->vertices_[i]; nzi < graph_ptr_->vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph_ptr_->edges_[nzi];
                    x_dict[dst_v] += update;
                }
            } else {
                for (auto nzi = graph_ptr_->vertices_[i]; nzi < graph_ptr_->vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph_ptr_->edges_[nzi];
                    auto re = rentry(dst_v, j + 1, graph_ptr_->n_);
                    auto re_old = r_dict[re];
                    auto re_new = re_old + update;
                    r_dict[re] = re_new;
                    if (re_new >= graph_ptr_->sr_degree(dst_v) * push_coefficient_vec_[j + 1] &&
                        re_old < graph_ptr_->sr_degree(dst_v) * push_coefficient_vec_[j + 1]) {
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
        for (auto &ele:x_dict) { ele.second *= (1.0 / max(graph_ptr_->sr_degree(ele.first), 1ul)); }
        auto pr_pairs = vector<pair<size_t, double>>(x_dict.begin(), x_dict.end());
        sort(pr_pairs.begin(), pr_pairs.end(), [](auto &&left, auto &&right) { return left.second > right.second; });

        auto cond_vec = vector<double>(pr_pairs.size());
        auto vol_vec = vector<size_t>(pr_pairs.size());
        auto cut_vec = vector<size_t>(pr_pairs.size());
        auto rank_map = unordered_map<size_t, size_t>();
        for (auto i = 0ul; i < pr_pairs.size(); i++) { rank_map[pr_pairs[i].first] = i; }

        auto vol_of_graph = graph_ptr_->vertices_[graph_ptr_->m_];
        auto vol_of_set = 0ul;
        auto cut_of_set = 0ul;

        for (auto i = 0ul; i < pr_pairs.size(); i++) {
            auto v = pr_pairs[i].first;
            auto deg = graph_ptr_->vertices_[v + 1] - graph_ptr_->vertices_[v];
            auto change = deg;
            for (auto nzi = graph_ptr_->vertices_[v]; nzi < graph_ptr_->vertices_[v + 1]; ++nzi) {
                auto neighbor_v = graph_ptr_->edges_[nzi];
                if (rank_map.count(neighbor_v) > 0 && rank_map[neighbor_v] < rank_map[v]) {
                    change -= 2;
                }
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
        return make_tuple(status, cluster);
    }

    auto HKGrow::HyperCluster(const vector<size_t> &seed_set, SpareseVec &x_dict) const {
        auto seed_dict = SpareseVec();
        for (auto &seed:seed_set) { seed_dict.emplace(seed, 1.0 / seed_set.size()); }

        auto seed_iter = max_element(begin(seed_set), end(seed_set), [&](auto &&left, auto &&right) {
            return graph_ptr_->sr_degree(left) < graph_ptr_->sr_degree(right);
        });
        auto max_deg = (seed_iter == seed_set.end() ? 0ul : graph_ptr_->sr_degree(*seed_iter));

        auto step_num = DiffuseWeight(seed_dict, x_dict, static_cast<size_t >(ceil(pow(graph_ptr_->n_, 1.5))));

        auto status_cluster = SweepCut(x_dict);
        auto &status = std::get<0>(status_cluster);
        status.steps_ = step_num;
        status.support_ = seed_dict.size();
        if (step_num == 0) { x_dict = seed_dict; }

        return status_cluster;
    }

    void HKGrow::ExecuteHRGRow(vector<size_t> &seeds, SpareseVec &x_dict) {
        auto stats = HyperCluster(seeds, x_dict);
    }

    HKGrow::HKGrow(unique_ptr<SparseRow> graph_ptr, double t, double eps) : t_(t), graph_ptr_(move(graph_ptr)) {
        taylor_deg_ = GetTaylorDegree(t, eps);
        psi_vec_ = ComputePsiVec(taylor_deg_, t);
        push_coefficient_vec_ = ComputePushCoefficientVec(taylor_deg_, eps, t, psi_vec_);
    }
}
