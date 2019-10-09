#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sdsl/wavelet_trees.hpp>
#include <string>
#include <sdsl/wt_int.hpp>
#include <tuple>
#include <vector>
using namespace std;
using namespace sdsl;

int MIN_LENGTH = 3;
int MAX_LENGTH = 10;

// returns the next smallest value in the wt that is bigger than x, within i and j
int range_next_value(wt_int<> wt, wt_int<>::const_iterator v, int i, int j, int x) {
    if (i > j) {
        return -1;
    }
    else if (!wt.is_leaf(*v)){
        return *v;
    }
    else {
        //auto il <- wt.rank();
    }
    return *v;
    /*
     * rnv(v, i, j, p, x)
if i > j then
return (⊥, 0, 0)
else if v is a leaf then
return (x, j − i + 1, p)
else
il ← rank0(Bv, i − 1) + 1
jl ← rank0(Bv, j)
ir ← i − il
, jr ← j − jr
nl ← jl − il + 1
if x ∈ labels(vr) then
return rnv(vr
, ir
, jr
, p + nl
, x)
else
(y, f, p
0
) ← rnv(vl
, il
, jl
, p, x)
if y ,⊥ then
return (y, f, p
0
)
else
return rnv(vr
, ir
, jr
, p + nl
,
min labels(vr))
end if
end if
end if
     */
}

int main(int argc, char *argv[]) {
    cst_sada<> cst;
    int min_reps;
    if (argc < 3) {
        min_reps = 2;
    } else {
        char *pEnd;
        min_reps = strtol(argv[2], NULL, 10);
    }
    //string file ;
    //cin >> file;
    //construct(cst, file, 1);
    construct(cst, "/home/anouk/Documents/memoria/data/testSeq2.fasta", 1);
    // a candidate is a tuple: text length, lb, rb
    vector<tuple<int, int, int>> candidate_list;
    // iterate over all nodes
    for (auto it = cst.begin(); it != cst.end(); ++it) {
        if (cst.depth(*it) > 1 && it.visit() == 1) {  // node visited for the first time
            auto v = *it;       // get the node by dereferencing the iterator
            // if depth of node is more than 1 and label has more than min_rep occs
            if (cst.size(v) >= min_reps) {
                // process node
                auto candidate = extract(cst, v);
                // if candidate length is within desired limits
                if (candidate.length() >= MIN_LENGTH && candidate.length() <= MAX_LENGTH) {
                    // string depth, leftmost leaf in SA and rightmost leaf in SA
                    // we can obtain the number of occs using rb-lb+1
                    candidate_list.emplace_back(cst.depth(v), cst.lb(v), cst.rb(v));
                }
            } else { // skip the subtree otherwise
                it.skip_subtree();
            }
        }
    }
    // the SA associated with the suffix tree
    csa_sada<> csa = cst.csa;

    // printing for debug purposes
    for (int i = 0; i < candidate_list.size(); i++) {
        cout << get<0>(candidate_list[i]) << " lb: " << get<1>(candidate_list[i]) << " rb: " << get<2>(candidate_list[i]) << endl;
    }
    wt_int<> wt;
    construct_im(wt, sdsl::util::to_string(csa), 'd');

    //cout << wt << endl;
    int a = range_next_value(wt, wt.begin()+=3, 1, 10, 3);
    cout << a << endl;
}

