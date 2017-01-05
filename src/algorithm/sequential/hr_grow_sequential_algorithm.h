//
// Created by cheyulin on 1/5/17.
//

#ifndef CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
#define CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H

/**
 * Implement a seeded heat-kernel clustering scheme.
 * [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
 * mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkgrow_mex.cpp
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

struct sparsevec {
    unordered_map<mwIndex, double> map;

    double get(mwIndex index, double default_value = 0.0) {
        auto it = map.find(index);
        return it == map.end() ? default_value : it->second;
    }

    double sum() {
        return accumulate(map.begin(), map.end(), 0,
                          [](auto &&left, auto &&right) { return left.second + right.second; });
    }

    mwIndex max_index() {
        return map.size() == 0 ? 0 :
               max_element(map.begin(), map.end(),
                           [](auto &&left, auto &&right) { return left.second < right.second; })->first;
    }
};

struct sparserow {
    mwSize n, m;
    mwIndex *ai;
    mwIndex *aj;
    double *a;

    mwIndex sr_degree(mwIndex u) {
        return ai[u + 1] - ai[u];
    }
};

//mwIndex sr_degree(sparserow *s, mwIndex u) {
//    return s->ai[u + 1] - s->ai[u];
//}

unsigned int get_taylor_degree(const double t, const double eps) {
    double eps_exp_t = eps * exp(t);
    double error = exp(t) - 1;
    double last = 1;
    double k = 0.;
    while (error > eps_exp_t) {
        k = k + 1;
        last = (last * t) / k;
        error = error - last;
    }
    return max(static_cast<unsigned int>(k), 1u);
}


/**
 *  gsqexpmseed inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms; normalized so
 *                  set[i] = 1/set.size(); )
 *  output:
 *      y = exp(tP) * set
 *              with infinity-norm accuracy of eps * e^t
 *              in the degree weighted norm
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 *      max_push_count - the total number of steps to run
 *      Q - the queue data structure
 */
template<class Queue>
int gsqexpmseed(sparserow *G, sparsevec &set, sparsevec &y, const double t, const double eps,
                const mwIndex max_push_count, Queue &Q) {
    mwIndex n = G->n;
    mwIndex N = (mwIndex) get_taylor_degree(t, eps);

    // initialize the weights for the different residual partitions
    // r(i,j) > d(i)*exp(t)*eps/(N*psi_j(t))
    //  since each coefficient but d(i) stays the same,
    //  we combine all coefficients except d(i)
    //  into the vector "pushcoeff"
    vector<double> psivec(N + 1, 0.);
    psivec[N] = 1;
    for (int k = 1; k <= N; k++) {
        psivec[N - k] = psivec[N - k + 1] * t / (double) (N - k + 1) + 1;
    } // psivec[k] = psi_k(t)
    vector<double> pushcoeff(N + 1, 0.);
    pushcoeff[0] = ((exp(t) * eps) / (double) N) / psivec[0]; // This is the correct version
    for (int k = 1; k <= N; k++) {
        pushcoeff[k] = pushcoeff[k - 1] * (psivec[k - 1] / psivec[k]);
    } // pushcoeff[j] = exp(t)*eps/(N*psivec[j])

    mwIndex ri = 0;
    mwIndex npush = 0;
    double rij = 0;
    // allocate data
    sparsevec rvec;

    // i is the node index, j is the "step"

    // set the initial residual, add to the queue
    for (auto it = set.map.begin(), itend = set.map.end(); it != itend; ++it) {
        ri = it->first;
        rij = it->second;
        rvec.map[rentry(ri, 0)] += rij;
        Q.push(rentry(ri, 0));
    }

    while (npush < max_push_count) {
        // STEP 1: pop top element off of heap
        ri = Q.front();
        Q.pop();
        mwIndex i = ri % n;
        mwIndex j = ri / n;

        double degofi = (double) G->sr_degree(i);
        rij = rvec.map[ri];

        y.map[i] += rij;

        // update r, no need to update heap here
        rvec.map[ri] = 0;

        double rijs = t * rij / (double) (j + 1);
        double ajv = 1. / degofi;
        double update = rijs * ajv;

        if (j == N - 1) {
            // this is the terminal case, and so we add the column of A, directly to the solution vector y
            for (mwIndex nzi = G->ai[i]; nzi < G->ai[i + 1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                y.map[v] += update;
            }
            npush += degofi;
        } else {
            // this is the interior case, and so we add the column of A
            // to the residual at the next time step.
            for (mwIndex nzi = G->ai[i]; nzi < G->ai[i + 1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                mwIndex re = rentry(v, j + 1);
                double reold = rvec.get(re);
                double renew = reold + update;
                double dv = G->sr_degree(v);
                rvec.map[re] = renew;
                if (renew >= dv * pushcoeff[j + 1] && reold < dv * pushcoeff[j + 1]) {
                    Q.push(re);
                }
            }
            npush += degofi;
        }
        // terminate when Q is empty, i.e. we've pushed all r(i,j) > eps*exp(t)*d(i)/(N*psi_j(t))
        if (Q.size() == 0) { return npush; }
    }//end 'while'
    return npush;
}


void cluster_from_sweep(sparserow *G, sparsevec &p, vector<mwIndex> &cluster, double *outcond, double *outvolume,
                        double *outcut) {
    // now we have to do the sweep over p in sorted order by value
    typedef vector<pair<int, double> > vertex_prob_type;
    vertex_prob_type prpairs(p.map.begin(), p.map.end());

    sort(prpairs.begin(), prpairs.end(), [](auto &&left, auto &&right) { return left.second > right.second; });

    // compute cutsize, volume, and conductance
    vector<double> conductance(prpairs.size());
    vector<mwIndex> volume(prpairs.size());
    vector<mwIndex> cutsize(prpairs.size());

    size_t i = 0;
    unordered_map<int, size_t> rank;
    for (auto it = prpairs.begin(), itend = prpairs.end();
         it != itend; ++it, ++i) {
        rank[it->first] = i;
    }
    mwIndex total_degree = G->ai[G->m];
    mwIndex curcutsize = 0;
    mwIndex curvolume = 0;
    i = 0;
    for (auto it = prpairs.begin(), itend = prpairs.end(); it != itend; ++it, ++i) {
        mwIndex v = it->first;
        mwIndex deg = G->ai[v + 1] - G->ai[v];
        mwIndex change = deg;
        for (mwIndex nzi = G->ai[v]; nzi < G->ai[v + 1]; ++nzi) {
            mwIndex nbr = G->aj[nzi];
            if (rank.count(nbr) > 0) {
                if (rank[nbr] < rank[v]) {
                    change -= 2;
                }
            }
        }
        curcutsize += change;
        //if (curvolume + deg > target_vol) {
        //break;
        //}
        curvolume += deg;
        volume[i] = curvolume;
        cutsize[i] = curcutsize;
        if (curvolume == 0 || total_degree - curvolume == 0) {
            conductance[i] = 1;
        } else {
            conductance[i] = (double) curcutsize /
                             (double) min(curvolume, total_degree - curvolume);
        }
    }

    // we stopped the iteration when it finished, or when it hit target_vol
    size_t lastind = i;
    double mincond = numeric_limits<double>::max();
    size_t mincondind = 0; // set to zero so that we only add one vertex
    for (i = 0; i < lastind; i++) {
        if (conductance[i] < mincond) {
            mincond = conductance[i];
            mincondind = i;
        }
    }
    if (lastind == 0) {
        mincond = 0.0;
    }
    for (auto &ele:prpairs) { cluster.emplace_back(ele.first); }
    if (outcond) { *outcond = mincond; }
    if (outvolume) { *outvolume = volume[mincondind]; }
    if (outcut) { *outcut = cutsize[mincondind]; }
}

struct local_hkpr_stats {
    double conductance;
    double volume;
    double support;
    double steps;
    double eps;
    double cut;
};

/** Cluster will contain a list of all the vertices in the cluster
 * @param set the set of starting vertices to use
 * @param t the value of t in the heatkernelPageRank computation
 * @param eps the solution tolerance eps
 * @param p the heatkernelpagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template<class Queue>
int hypercluster_heatkernel_multiple(sparserow *G, const vector<mwIndex> &set, double t, double eps,
                                     sparsevec &p, sparsevec &r, Queue &q,
                                     vector<mwIndex> &cluster, local_hkpr_stats *stats) {
    // reset data
    p.map.clear();
    r.map.clear();
    q.empty();

    size_t maxdeg = 0;
    for (size_t i = 0; i < set.size(); ++i) { //populate r with indices of "set"
        assert(set[i] >= 0);
        assert(set[i] < G->n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = G->sr_degree(set[i]);
        r.map[set[i]] = 1. / (double) (set.size()); // r is normalized to be stochastic
        maxdeg = max(maxdeg, setideg);
    }

    int nsteps = gsqexpmseed(G, r, p, t, eps, static_cast<int>(ceil(pow(G->n, 1.5))), q);

    if (nsteps == 0) {
        p = r; // just copy over the residual
    }
    int support = static_cast<int>(r.map.size());
    if (stats) { stats->steps = nsteps; }
    if (stats) { stats->support = support; }

    // scale the probablities by their degree
    for (auto it = p.map.begin(), itend = p.map.end();
         it != itend; ++it) {
        it->second *= (1.0 / (double) max(G->sr_degree(it->first), (mwIndex) 1));
    }

    double *outcond = NULL;
    double *outvolume = NULL;
    double *outcut = NULL;
    if (stats) { outcond = &stats->conductance; }
    if (stats) { outvolume = &stats->volume; }
    if (stats) { outcut = &stats->cut; }
    cluster_from_sweep(G, p, cluster, outcond, outvolume, outcut);
    return 0;
}

/** Grow a set of seeds via the heat-kernel.
 * @param G sparserow version of input matrix A
 * @param seeds a vector of input seeds seeds (index 0, N-1), and then
 *          updated to have the final solution nodes as well.
 * @param t the value of t in the heat-kernel
 * @param eps the solution tolerance epsilon
 * @param fcond the final conductance score of the set.
 * @param fcut the final cut score of the set
 * @param fvol the final volume score of the set
 */
void hkgrow(sparserow *G, vector<mwIndex> &seeds, double t,
            double eps, double *fcond, double *fcut,
            double *fvol, sparsevec &p, double *npushes) {
    sparsevec r;
    queue<mwIndex> q;
    local_hkpr_stats stats;
    vector<mwIndex> bestclus;
    hypercluster_heatkernel_multiple(G, seeds, t, eps, p, r, q, bestclus, &stats);
    seeds = bestclus;
    *npushes = stats.steps;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
}

/*
void copy_array_to_index_vector(const mxArray *v, vector<mwIndex> &vec) {
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);

    vec.resize(n);

    for (size_t i = 0; i < n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex) elem - 1;
    }
}


// [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 5) {
        debugflag = (int) mxGetScalar(prhs[4]);
    }

    const mxArray *mat = prhs[0];
    const mxArray *set = prhs[1];

    mxArray *cond = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *cut = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *vol = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *npushes = mxCreateDoubleMatrix(1, 1, mxREAL);

    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut; }
    if (nlhs > 3) { plhs[3] = vol; }
    if (nlhs > 5) { plhs[5] = npushes; }

    double eps = pow(10, -3);
    double t = 15.;

    if (nrhs >= 4) {
        t = mxGetScalar(prhs[2]);
        eps = mxGetScalar(prhs[3]);
    }

    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);

    vector<mwIndex> seeds;
    copy_array_to_index_vector(set, seeds);
    sparsevec hkpr;

    hkgrow(&r, seeds, t, eps, mxGetPr(cond), mxGetPr(cut), mxGetPr(vol), hkpr, mxGetPr(npushes));

    if (nlhs > 0) { // sets output "bestset" to the set of best conductance
        mxArray *cassign = mxCreateDoubleMatrix(seeds.size(), 1, mxREAL);
        plhs[0] = cassign;

        double *ci = mxGetPr(cassign);
        for (size_t i = 0; i < seeds.size(); ++i) {
            ci[i] = (double) (seeds[i] + 1);
        }
    }
    if (nlhs > 4) { // sets output "y" to the heat kernel vector computed
        mxArray *hkvec = mxCreateDoubleMatrix(r.n, 1, mxREAL);
        plhs[4] = hkvec;
        double *ci = mxGetPr(hkvec);
        for (sparsevec:: auto it = hkpr.map.begin(), itend = hkpr.map.end();
             it != itend; ++it) {
            ci[it->first] = it->second;
        }
    }
}*/

#endif //CODES_YCHE_HR_GROW_SEQUENTIAL_ALGORITHM_H
