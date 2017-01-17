//
// Created by cheyulin on 1/5/17.
//

#include "hr_grow_sequential_algorithm.h"

unsigned int get_taylor_degree(double t, double eps) {
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


mwIndex gs_qexpm_seed(sparse_row *graph, sparse_vec &set, sparse_vec &y, const double t, const double eps,
                      const mwIndex max_push_count, queue<mwIndex> &Q) {
    mwIndex n = graph->n_;
    mwIndex taylor_deg = (mwIndex) get_taylor_degree(t, eps);

    // initialize the weights for the different residual partitions
    // r(i,j) > d(i)*exp(t)*eps/(taylor_deg*psi_j(t))
    //  since each coefficient but d(i) stays the same,
    //  we combine all coefficients except d(i)
    //  into the vector "push_coeff"
    // psi_vec[k] = psi_k(t)

    vector<double> psi_vec(taylor_deg + 1, 0);
    psi_vec[taylor_deg] = 1;
    for (int k = 1; k <= taylor_deg; k++) {
        psi_vec[taylor_deg - k] = psi_vec[taylor_deg - k + 1] * t / (double) (taylor_deg - k + 1) + 1;
    }

    vector<double> push_coeff(taylor_deg + 1, 0.);
    push_coeff[0] = ((exp(t) * eps) / (double) taylor_deg) / psi_vec[0]; // This is the correct version
    for (int k = 1; k <= taylor_deg; k++) {
        push_coeff[k] = push_coeff[k - 1] * (psi_vec[k - 1] / psi_vec[k]);
    } // push_coeff[j] = exp(t)*eps/(taylor_deg*psi_vec[j])

    mwIndex ri = 0;
    mwIndex npush = 0;
    double rij = 0;

    sparse_vec rvec;
    // i is the node index, j is the "step"
    // set the initial residual, add to the queue
    for (auto &ele: set.weight_map_) {
        tie(ri, rij) = ele;
        rvec.weight_map_[rentry(ri, 0)] += rij;
        Q.push(rentry(ri, 0));
    }

    while (npush < max_push_count) {
        // STEP 1: pop top element off of heap
        ri = Q.front();
        Q.pop();

        mwIndex i = ri % n;
        mwIndex j = ri / n;

        double degofi = (double) graph->sr_degree(i);
        rij = rvec.weight_map_[ri];

        y.weight_map_[i] += rij;

        // update r, no need to update heap here
        rvec.weight_map_[ri] = 0;

        double rijs = t * rij / (double) (j + 1);
        double ajv = 1. / degofi;
        double update = rijs * ajv;

        if (j == taylor_deg - 1) {
            // this is the terminal case, and so we add the column of A, directly to the solution vector y
            for (mwIndex nzi = graph->vertices_[i]; nzi < graph->vertices_[i + 1]; ++nzi) {
                mwIndex v = graph->edges_[nzi];
                y.weight_map_[v] += update;
            }
            npush += degofi;
        } else {
            // this is the interior case, and so we add the column of A
            // to the residual at the next time step.
            for (mwIndex nzi = graph->vertices_[i]; nzi < graph->vertices_[i + 1]; ++nzi) {
                mwIndex v = graph->edges_[nzi];
                mwIndex re = rentry(v, j + 1);
                double reold = rvec.get(re);
                double renew = reold + update;
                double dv = graph->sr_degree(v);
                rvec.weight_map_[re] = renew;
                if (renew >= dv * push_coeff[j + 1] && reold < dv * push_coeff[j + 1]) {
                    Q.push(re);
                }
            }
            npush += degofi;
        }
        // terminate when Q is empty, i.e. we've pushed all r(i,j) > eps*exp(t)*d(i)/(taylor_deg*psi_j(t))
        if (Q.size() == 0) { return npush; }
    }//end 'while'
    return npush;
}


void cluster_from_sweep(sparse_row *G, sparse_vec &p, vector<mwIndex> &cluster, double *outcond, double *outvolume,
                        double *outcut) {
    // now we have to do the sweep over p in sorted order by value
    vector<pair<size_t, double>> pr_pairs(p.weight_map_.begin(), p.weight_map_.end());
    sort(pr_pairs.begin(), pr_pairs.end(), [](auto &&left, auto &&right) { return left.second > right.second; });

    // compute cutsize, volume, and conductance
    vector<double> conductance(pr_pairs.size());
    vector<mwIndex> volume(pr_pairs.size());
    vector<mwIndex> cutsize(pr_pairs.size());

    unordered_map<size_t, size_t> rank_map;

    for (auto i = 0ul; i < pr_pairs.size(); i++) {
        rank_map[pr_pairs[i].first] = i;
    }
    mwIndex total_degree = G->vertices_[G->m_];
    mwIndex curcutsize = 0;
    mwIndex curvolume = 0;

    for (auto i = 0ul; i < pr_pairs.size(); i++) {
        mwIndex v = pr_pairs[i].first;
        mwIndex deg = G->vertices_[v + 1] - G->vertices_[v];
        mwIndex change = deg;
        for (mwIndex nzi = G->vertices_[v]; nzi < G->vertices_[v + 1]; ++nzi) {
            mwIndex nbr = G->edges_[nzi];
            if (rank_map.count(nbr) > 0 && rank_map[nbr] < rank_map[v]) {
                change -= 2;
            }
        }

        curcutsize += change;
        curvolume += deg;
        volume[i] = curvolume;
        cutsize[i] = curcutsize;
        if (curvolume == 0 || total_degree - curvolume == 0) {
            conductance[i] = 1;
        } else {
            conductance[i] = static_cast<double>(curcutsize) / min(curvolume, total_degree - curvolume);
        }
    }

    // we stopped the iteration when it finished, or when it hit target_vol
    double min_cond = numeric_limits<double>::max();
    size_t min_cond_idx = 0; // set to zero so that we only add one vertex

    for (auto i = 0ul; i < pr_pairs.size(); i++) {
        if (conductance[i] < min_cond) {
            min_cond = conductance[i];
            min_cond_idx = i;
        }
    }
    if (pr_pairs.size() == 0) {
        min_cond = 0.0;
    }
    for (auto &ele:pr_pairs) { cluster.emplace_back(ele.first); }
    if (outcond) { *outcond = min_cond; }
    if (outvolume) { *outvolume = volume[min_cond_idx]; }
    if (outcut) { *outcut = cutsize[min_cond_idx]; }
}


int hyper_cluster_heat_kernel_multiple(sparse_row *G, const vector<mwIndex> &seed_set, double t, double eps,
                                       sparse_vec &p, sparse_vec &r, queue<mwIndex> &q, vector<mwIndex> &cluster,
                                       local_hkpr_stats *stats) {
    p.weight_map_.clear();
    r.weight_map_.clear();
    q.empty();

    size_t max_deg = 0;
    for (size_t i = 0; i < seed_set.size(); ++i) {
        size_t v_degree = G->sr_degree(seed_set[i]);
        r.weight_map_[seed_set[i]] = 1.0 / seed_set.size();
        max_deg = max(max_deg, v_degree);
    }

    auto nsteps = gs_qexpm_seed(G, r, p, t, eps, static_cast<size_t >(ceil(pow(G->n_, 1.5))), q);

    if (nsteps == 0) {
        p = r; // just copy over the residual
    }

    if (stats) { stats->steps = nsteps; }
    if (stats) { stats->support = static_cast<int>(r.weight_map_.size()); }

    // scale the probabilities by their degree
    for (auto &ele:p.weight_map_) { ele.second *= (1.0 / max(G->sr_degree(ele.first), (mwIndex) 1)); }

    double *out_cond = nullptr;
    double *out_volume = nullptr;
    double *out_cut = nullptr;
    if (stats) { out_cond = &stats->conductance; }
    if (stats) { out_volume = &stats->volume; }
    if (stats) { out_cut = &stats->cut; }

    cluster_from_sweep(G, p, cluster, out_cond, out_volume, out_cut);
    return 0;
}

void hk_grow(sparse_row *G, vector<mwIndex> &seeds, double t, double eps, double *fcond, double *fcut,
             double *fvol, sparse_vec &p, double *npushes) {
    sparse_vec r;
    queue<mwIndex> q;
    local_hkpr_stats stats;
    vector<mwIndex> bestclus;
    hyper_cluster_heat_kernel_multiple(G, seeds, t, eps, p, r, q, bestclus, &stats);
    seeds = bestclus;
    *npushes = stats.steps;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
}
