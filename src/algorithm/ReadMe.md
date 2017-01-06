#Algorithm Descriptions
##OVerlapping Community Detection Algorithms
- Cis : Connected Iterative Scan
- Demon: Democratically Estimate of Modular Organization of Network

##Sequential Implementation
- CIs: [cis_sequential_algorithm.h](sequential/cis_sequential_algorithm.h), [cis_sequential_algorithm.cpp](sequential/cis_sequential_algorithm.cpp)
- Demon: [demon_sequential_algorithm.h](sequential/demon_sequential_algorithm.h), [demon_seqential_algorithm.cpp](sequential/demon_seqential_algorithm.cpp)

##HrGrow Algo

- parameters

```zsh

/*
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

    sparse_row r;
    r.m_ = mxGetM(mat);
    r.n_ = mxGetN(mat);
    r.vertices_ = mxGetJc(mat);
    r.edges_ = mxGetIr(mat);
    r.weight_ = mxGetPr(mat);

    vector<mwIndex> seeds;
    copy_array_to_index_vector(set, seeds);
    sparse_vec hkpr;

    hk_grow(&r, seeds, t, eps, mxGetPr(cond), mxGetPr(cut), mxGetPr(vol), hkpr, mxGetPr(npushes));

    if (nlhs > 0) { // sets output "bestset" to the set of best conductance
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
        for (sparse_vec:: auto it = hkpr.weight_map_.begin(), itend = hkpr.weight_map_.end();
             it != itend; ++it) {
            ci[it->first] = it->second;
        }
    }
}*/

```

