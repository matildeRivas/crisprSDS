#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/construct_sa.hpp>
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

int MIN_LENGTH = 2;
int MAX_LENGTH = 3;
int MIN_SPACER_LENGTH = 2;
int MAX_SPACER_LENGTH = 3;

// Compares two tuples according to their second element
bool comp_by_second(const tuple<int, int> &a,
                    const tuple<int, int> &b) {
    return (get<1>(a) < get<1>(b));
}


int_vector<> create_sa(string infile) {
    int_vector<> seq;
    int32_t n;
    {
        load_vector_from_file(seq, infile, 1);
        n = seq.size();
        seq.resize(n + 1);
        n = seq.size();
        seq[n - 1] = 0; // Represents the symbol $
    }

    int_vector<> sa(1, 0, bits::hi(n) + 1);
    sa.resize(n);
    algorithm::calculate_sa((const unsigned char *) seq.data(), n, sa);
    return sa;
}

// @candidate: text length, lower bound in SA, right bound in SA
// returns candidates that meet the requirements for being CRISPR as a tuple containing
// initial position in text,the DR's length and the number of occs
int validate(tuple<int, int, int> candidate, int_vector<> sa, wt_int<> wt) {
    // obtain bounds
    int len = get<0>(candidate);
    auto lb = get<1>(candidate);
    auto rb = get<2>(candidate);
    cout << len << " lb: " << lb << " rb: " << rb << endl;
    vector<tuple<int, int>> cont_pairs;
    // for each candidate in group
    for (int c = lb; c <= rb; c++) {
        int sa_value = sa[c];
        cout << "searching value after " << sa_value << " in range: " << sa_value + len + MIN_SPACER_LENGTH << " "
             << sa_value + len + MAX_SPACER_LENGTH << endl;
        // search for next repeat in group that's within spacer range
        auto rs = wt.range_search_2d(lb, rb, sa_value + len + MIN_SPACER_LENGTH,
                                     sa_value + len + MAX_SPACER_LENGTH);
        cout << "values found " << get<0>(rs) << endl;
        if (get<0>(rs) > 0) {
            for (int k = 0; k < get<0>(rs); k++) {
                // tuple (position, value)
                auto c = get<1>(rs).at(k);
                cont_pairs.emplace_back(sa_value, (int) get<1>(rs).at(k).second);
                cout << sa_value << " to " << get<1>(rs).at(k).second << endl;
            }
        }

    }
    sort(cont_pairs.begin(), cont_pairs.end(), comp_by_second);
    sort(cont_pairs.begin(), cont_pairs.end());
    // TODO: dictionary for last element in group
    // TODO: obtain CRISPR

    return 0;
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
    construct(cst, "/home/anouk/Documents/memoria/data/testSeq.fasta", 1);
    //construct_im(cst, "ATCGTACGTTCGAACT", 1);

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
    //int_vector<> sa = create_sa(file);
    int_vector<> sa = create_sa("/home/anouk/Documents/memoria/data/testSeq.fasta");
    // construct the wavelet tree using the SA
    wt_int<> wt;
    construct_im(wt, sa);
    cout << sa << endl;
    cout << wt << endl;
    vector<tuple<int, int, int>> filtered_candidates;

    for (int i = 0; i < candidate_list.size(); i++) {
        //filtered_candidates.emplace_back(validate(candidate_list[i], sa, wt));
        validate(candidate_list[i], sa, wt);
        //cout << get<0>(candidate_list[i]) << " " <<get<1>(candidate_list[i]) <<" " <<get<2>(candidate_list[i]) <<endl ;
    }


}

//vector<tuple<int, int, int>>