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
#include <map>
#include <sstream>

#include <crispr.h>

using namespace std;
using namespace sdsl;


//23 47 28 37
int MIN_LENGTH = 23;
int MAX_LENGTH = 47;
int MIN_SPACER_LENGTH = 28;
int MAX_SPACER_LENGTH = 37;


void write_output(vector<tuple<int, int, int, int>> detected_crispr, string infile, string outfile_name) {
    ofstream outfile;
    outfile.open(outfile_name);
    for (auto &crispr_chain : detected_crispr) {
        int length = get<1>(crispr_chain);
        int start = get<2>(crispr_chain);
        outfile << "Number of repetitions: " << get<0>(crispr_chain) << " DR length: " << length
                << " Position of first DR: "
                << start <<
                " Position of last DR: " << get<3>(crispr_chain) << '\n' << endl;
    }
    outfile.close();
}


//selects repeats that meet basic criteria
vector<tuple<int, int, int>> select_candidates(cst_sada<> cst, int min_reps) {
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
    return candidate_list;
}


// Compares two tuples according to their second element
bool comp_by_second(const tuple<int, int> &a,
                    const tuple<int, int> &b) {
    return (get<1>(a) < get<1>(b));
}

// Compares two tuples according to their third element
bool comp_by_third(const tuple<int, int, int, int> &a,
                   const tuple<int, int, int, int> &b) {
    return (get<2>(a) < get<2>(b));
}

// Constructs suffix array from string in given file
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
// the number of repeats, their length, start and end position in text
vector<tuple<int, int, int, int>> validate(tuple<int, int, int> candidate, int_vector<> sa, wt_int<> wt) {
    // obtain bounds
    int len = get<0>(candidate);
    auto lb = get<1>(candidate);
    auto rb = get<2>(candidate);
    //cout << len << " lb: " << lb << " rb: " << rb << endl;
    vector<tuple<int, int>> cont_pairs;
    // for each candidate in group
    for (int c = lb; c <= rb; c++) {
        int sa_value = sa[c];
        // cout << "searching value after " << sa_value << " in range: " << sa_value + len + MIN_SPACER_LENGTH << " "
        //     << sa_value + len + MAX_SPACER_LENGTH << endl;
        // search for next repeat in group that's within spacer range
        auto rs = wt.range_search_2d(lb, rb, sa_value + len + MIN_SPACER_LENGTH,
                                     sa_value + len + MAX_SPACER_LENGTH);
        //cout << "values found " << get<0>(rs) << endl;
        if (get<0>(rs) > 0) {
            for (int k = 0; k < get<0>(rs); k++) {
                // tuple (position, value)
                auto c = get<1>(rs).at(k);
                cont_pairs.emplace_back(sa_value, (int) get<1>(rs).at(k).second);
                // cout << sa_value << " to " << get<1>(rs).at(k).second << endl;
            }
        }

    }
    sort(cont_pairs.begin(), cont_pairs.end(), comp_by_second);
    sort(cont_pairs.begin(), cont_pairs.end());

    // key is the position of last occurence in chain so far
    // value stores number of occs, pos of first and last found so far
    std::map<int, tuple<int, int, int>>
            chain_map;
    for (auto &p : cont_pairs) {
        int from = get<0>(p);
        int to = get<1>(p);
        // check if there's a chain that ends in "from" value
        if (chain_map.count(from)) {
            //grab the value from that chain, update number of occs and last value, and store in "to"
            auto prev = chain_map[from];
            int occs = get<0>(prev);
            int first_occ = get<1>(prev);
            chain_map[to] = tuple<int, int, int>(occs + 1, first_occ, to);
            chain_map.erase(from);
        } else {
            // add new start of chain in "to" with "from" as first occurence
            chain_map[to] = tuple<int, int, int>(2, from, to);
        }

    }

    vector<tuple<int, int, int, int>> crispr_list;
    for (auto &x : chain_map) {
        // for each CRISPR we store the number of repeats, their length, start and end position
        tuple<int, int, int, int> crispr = tuple<int, int, int, int>(get<0>(x.second), len, get<1>(x.second) + 1,
                                                                     get<2>(x.second) + 1);
        crispr_list.emplace_back(crispr);

    }

    return crispr_list;
}


vector<tuple<int, int, int, int>> find_crispr(string filename, int min_reps) {
    cst_sada<> cst;
    construct(cst, filename, 1);
    //start time
    // a candidate is a tuple: text length, lb, rb
    vector<tuple<int, int, int>> candidate_list;
    // select candidates that meet basic criteria
    candidate_list = select_candidates(cst, min_reps);
    //end time
    // get the SA associated with the suffix tree
    int_vector<> sa = create_sa(filename);
    //int_vector<> sa = create_sa("/home/anouk/Documents/memoria/data/GI326314823.fasta");
    // construct the wavelet tree using the SA
    wt_int<> wt;
    construct_im(wt, sa);
    //cout << sa << endl;
    vector<tuple<int, int, int, int>> pre_filtered_candidates;
    vector<tuple<int, int, int, int>> filtered_candidates;
    //start time
    // validate candidates and insert them in pre_filtered list
    for (int i = 0; i < candidate_list.size(); i++) {
        for (auto &cr : validate(candidate_list[i], sa, wt))
            pre_filtered_candidates.emplace_back(cr);
    }
    write_output(pre_filtered_candidates, filename, "/home/anouk/Documents/memoria/output/NC_013210_preout.txt");
    // deletes chains that are contained in another one, CHECK length
    if (pre_filtered_candidates.size()) {
        sort(pre_filtered_candidates.begin(), pre_filtered_candidates.end(), comp_by_third);
        tuple<int, int, int, int> checker = pre_filtered_candidates[0];
        for (auto i = pre_filtered_candidates.begin() + 1; i != pre_filtered_candidates.end(); i++) {
            if (get<2>(*i) < get<3>(checker)) {
                if (get<0>(*i) > get<0>(checker)) {
                    checker = *i;
                }
            } else {
                filtered_candidates.emplace_back(checker);
                checker = *i;
            }
        }
        filtered_candidates.emplace_back(checker);

    }
    //end time
    return filtered_candidates;

}

// verifies a detected crispr
void verify_crispr(Crispr &crispr, vector<tuple<int, int, int, int>> detected_crispr) {
    for (auto &crispr_chain : detected_crispr) {
        if (crispr.percentage_detected() < 100) {
            crispr.check_candidate(get<0>(crispr_chain), get<1>(crispr_chain), get<2>(crispr_chain),
                                   get<3>(crispr_chain));
        } else { break; }
    }
}

//por cada crispr voy a tener quee hacer una linea con su ID, la cantidad de reps que tiene,
// las parcialmente detectadas y la cantidad detecada. Luego puedo procesar esta informacion de otra manera

void run_test(vector<Crispr> ground_truth, string filename, int min_reps, string outfilename) {
    // run algorithm to find crispr
    vector<tuple<int, int, int, int>> detected_crispr;
    detected_crispr = find_crispr(filename, min_reps);
    std::stringstream ss;
    ss << "number of crispr \t reported crispr \n" << ground_truth.size() << "\t" << detected_crispr.size() << "\n \n";
    ss << "id \t repeats \t detected occs \t partially detected occs \n";
    // for each real crispr, check whether it was correctly detected
    for (auto &crispr : ground_truth) {
        verify_crispr(crispr, detected_crispr);
        ss << crispr.id << "\t" << crispr.positions.size() << "\t" << crispr.get_detected() << "\t"
           << crispr.get_partially_detected() << "\n";
    }
    std::string str = ss.str();
    ofstream outfile;
    outfile.open(outfilename);
    outfile << str;
    outfile.close();
}


int main(int argc, char *argv[]) {

    int min_reps;
    if (argc < 3) {
        min_reps = 2;
    } else {
        char *pEnd;
        min_reps = strtol(argv[2], NULL, 10);
    }
    string file;
    //cin >> file;

    //Clostridioides difficile  NZ https://crispr.i2bc.paris-saclay.fr/crispr/crispr_db.php?checked%5B%5D=NZ_LN614756
    //
    //GI326314823
    file = "/home/anouk/Documents/memoria/data/NC_013210.fasta";
    string outfile;
    outfile = "/home/anouk/Documents/memoria/output/NC_013210_out.txt";
    /*vector<Crispr> genome23;
    vector<int> NC_015138_1_pos{300343, 300409, 300475, 300541, 300607, 300673, 300739, 300806, 300872, 300938, 301004,
                                301070,
                                301136, 301202, 301268, 301334, 301400, 301466, 301532, 301598, 301664, 301730, 301796,
                                301862,
                                301928, 301995, 302061, 302127, 302193, 302259, 302325, 302391, 302457, 302523, 302589,
                                302655,
                                302721, 302786, 302852, 302918, 302984, 303050, 303115, 303181, 303247, 303313, 303379,
                                303445};
    Crispr NC_015138_1 = Crispr("NC_015138_1", 36, NC_015138_1_pos, "AGTCTAGATCACTGGGATATGCGCACTGGCCGGAAC");
    vector<int> NC_015138_2_pos{304387, 304453};
    Crispr NC_015138_2 = Crispr("NC_015138_2",36, NC_015138_2_pos, "AGTCTAGATCACTGGGATATGCGCACTGGCCGGAAC");
    vector<int> NC_015138_3_pos{305395, 305462, 305527, 305593, 305659, 305725, 305791, 305857, 305923, 305989, 306055,
                                306121, 306187};
    Crispr NC_015138_3 = Crispr("NC_015138_3", 36, NC_015138_3_pos, "AGTCTAGATCACTGGGATATGCGCACTGGCCGGAAC");


    genome23.emplace_back(NC_015138_1);
    */
    //tuple<int, int, int> test23 = run_test(genome23, file, min_reps, 85);
    //cout << "TP: " << get<0>(test23) << " FP: " << get<1>(test23) << " FN: " << get<2>(test23) << endl;

    // printing for debug reasons
    //file = "/home/anouk/Documents/memoria/data/Acidithiobacillus_ferrivorans_ACHa_45_pseudochromosome.fna";

    vector<Crispr> genome013210;
    vector<int> NC_013210_pos{123524, 123585, 123646, 123707, 123768, 123829, 123890, 123951, 124012, 124073, 124134,
                              124195, 124256, 124317, 124378, 124439, 124500, 124561, 124622, 124683, 124744, 124805,
                              124866, 124927
    };
    Crispr NC_013210 = Crispr("NC_013210", 29, NC_013210_pos, "GTGTTCCCCGCACACGCGGGGATGAACCG");
    genome013210.emplace_back(NC_013210);
    string test_output = "/home/anouk/Documents/memoria/output/NC_0123210_test.txt";
    run_test(genome013210, file, min_reps, test_output);


}