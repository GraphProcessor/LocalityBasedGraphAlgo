##HrGrow Algo

- sparse vec

```cpp
 struct SpareseVec {
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

```

- info

```cpp
/**
 * Implement a seeded heat-kernel clustering scheme.
 * [bestset,cond,cut_,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
 */
```

- `ComputePsiVec`

```cpp
// initialize the weights for the different residual partitions
        // r(i,j) > d(i)*exp(t)*eps/(taylor_deg*psi_j(t))
        //  since each coefficient but d(i) stays the same,
        //  we combine all coefficients except d(i)
        //  into the vector "push_coeff_vec"
        // psi_vec[k] = psi_k(t)
```

- `compute threshold`

```cpp
        // push_coefficient_vec[j] = exp(t)*eps/(taylor_deg*psi_vec[j])

```

- comment info

```cpp
/**
 *  gsqexpmseed inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms; normalized so set[i] = 1/set.size(); )
 *  output:
 *      y = exp(tP) * set
 *              with infinity-norm accuracy of eps * e^t in the degree weighted norm
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 *      max_push_count - the total number of steps to run
 *      Q - the queue data structure
 */
 
// i is the node index, j is the "step"

// STEP 1: pop top element off of heap
// update r, no need to update heap here

// terminate when task_queue is empty, i.e. we've pushed all r(i,j) > eps*exp(t)*d(i)/(taylor_deg*psi_j(t))

 
 
 /** Cluster will contain a list of all the vertices in the cluster
  * @param seed_set the set of starting vertices to use
  * @param t the value of t in the heat-kernel-PageRank computation
  * @param eps the solution tolerance eps
  * @param p the heat-kernel-pagerank vector
  * @param r the residual vector
  * @param a vector which supports .push_back to add vertices for the cluster
  * @param stats a structure for statistics of the computation
  */
  
// we stopped the iteration when it finished, or when it hit target_vol
  

/** Grow a set of seeds via the heat-kernel.
 * @param G sparserow version of input matrix A
 * @param seeds a vector of input seeds seeds (index 0, N-1), and then
 *          updated to have the final solution nodes as well.
 * @param t the value of t in the heat-kernel
 * @param eps the solution tolerance epsilon
 * @param fcond the final conductance_ score of the set.
 * @param fcut the final cut_ score of the set
 * @param fvol the final volume score of the set
 */

```

- seed expand

```cpp
            if (j == taylor_deg - 1) {
                // this is the terminal case, and so we add the column of A, directly to the solution vector y
                for (size_t nzi = graph.vertices_[i]; nzi < graph.vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph.edges_[nzi];
                    y.weight_map_[dst_v] += update;
                }
                push_num += deg_of_i;
            } else {
                // this is the interior case, and so we add the column of A to the residual at the next time step.
                for (auto nzi = graph.vertices_[i]; nzi < graph.vertices_[i + 1]; ++nzi) {
                    auto dst_v = graph.edges_[nzi];
                    auto re = rentry(dst_v, j + 1, graph.n_);
                    auto re_old = residual_vec.get(re);
                    auto re_new = re_old + update;
                    double dv = graph.out_degree(dst_v);
                    residual_vec.weight_map_[re] = re_new;
                    if (re_new >= dv * push_coefficient_vec[j + 1] && re_old < dv * push_coefficient_vec[j + 1]) {
                        task_queue.emplace(dst_v, j + 1);
                    }
                }
                push_num += deg_of_i;
            }
```

- parameters

```zsh

/*
// [bestset,cond,cut_,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 5) {
        debugflag = (int) mxGetScalar(prhs[4]);
    }

    const mxArray *mat = prhs[0];
    const mxArray *set = prhs[1];

    mxArray *cond = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *cut_ = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *vol = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *npushes = mxCreateDoubleMatrix(1, 1, mxREAL);

    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut_; }
    if (nlhs > 3) { plhs[3] = vol; }
    if (nlhs > 5) { plhs[5] = npushes; }

    double eps = pow(10, -3);
    double t = 15.;

    if (nrhs >= 4) {
        t = mxGetScalar(prhs[2]);
        eps = mxGetScalar(prhs[3]);
    }

    SparseRow r;
    r.m_ = mxGetM(mat);
    r.n_ = mxGetN(mat);
    r.vertices_ = mxGetJc(mat);
    r.edges_ = mxGetIr(mat);
    r.weight_ = mxGetPr(mat);

    vector<mwIndex> seeds;
    copy_array_to_index_vector(set, seeds);
    SpareseVec hkpr;

    ExecuteHRGRow(&r, seeds, t, eps, mxGetPr(cond), mxGetPr(cut_), mxGetPr(vol), hkpr, mxGetPr(npushes));

    if (nlhs > 0) { // sets output "bestset" to the set of best conductance_
        mxArray *cassign = mxCreateDoubleMatrix(seeds.size(), 1, mxREAL);
        plhs[0] = cassign;

        double *ci = mxGetPr(cassign);
        for (size_t i = 0; i < seeds.size(); ++i) {
            ci[i] = (double) (seeds[i] + 1);
        }
    }
    if (nlhs > 4) { // sets output "y" to the heat kernel vector computed
        mxArray *hkvec = mxCreateDoubleMatrix(r.n_, 1, mxREAL);
        plhs[4] = hkvec;
        double *ci = mxGetPr(hkvec);
        for (SpareseVec:: auto it = hkpr.weight_map_.begin(), itend = hkpr.weight_map_.end();
             it != itend; ++it) {
            ci[it->first] = it->second;
        }
    }
}*/

```