//
// Created by cheyulin on 1/5/17.
//

#ifndef CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
#define CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H

/**
 * Implement a seeded heat-kernel clustering scheme.
 * [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
 */

#include <cmath>
#include <cassert>

#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#define mwIndex size_t
#define mwSize size_t
#define rentry(i, j) ((i)+(j)*n)

using namespace std;

struct sparse_vec {
    unordered_map<mwIndex, double> weight_map_;

    double get(mwIndex index, double default_value = 0.0) {
        auto it = weight_map_.find(index);
        return it == weight_map_.end() ? default_value : it->second;
    }

    double sum() {
        return accumulate(weight_map_.begin(), weight_map_.end(), 0,
                          [](auto &&left, auto &&right) { return left + right.second; });
    }

    mwIndex max_index() {
        return weight_map_.size() == 0 ? 0 :
               max_element(weight_map_.begin(), weight_map_.end(),
                           [](auto &&left, auto &&right) { return left.second < right.second; })->first;
    }
};

struct sparse_row {
    mwSize n_, m_;
    mwIndex *vertices_;
    mwIndex *edges_;
    double *weight_;

    mwIndex sr_degree(mwIndex u) {
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

unsigned int get_taylor_degree(double t, double eps);


mwIndex gs_qexpm_seed(sparse_row *graph, sparse_vec &set, sparse_vec &y, const double t, const double eps,
                      const mwIndex max_push_count, queue<mwIndex> &Q);


void cluster_from_sweep(sparse_row *G, sparse_vec &p, vector<mwIndex> &cluster, double *outcond,
                        double *outvolume, double *outcut);

int hyper_cluster_heat_kernel_multiple(sparse_row *G, const vector<mwIndex> &seed_set, double t, double eps,
                                       sparse_vec &p, sparse_vec &r, queue<mwIndex> &q, vector<mwIndex> &cluster,
                                       local_hkpr_stats *stats);

void hk_grow(sparse_row *G, vector<mwIndex> &seeds, double t, double eps, double *fcond, double *fcut,
             double *fvol, sparse_vec &p, double *npushes);

#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
