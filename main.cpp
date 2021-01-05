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
#include <chrono>

#include <crispr.h>

using namespace std;
using namespace sdsl;


//23 47 28 37
int MIN_LENGTH = 23;
int MAX_LENGTH = 47;
int MIN_SPACER_LENGTH = 28;
int MAX_SPACER_LENGTH = 37;


void write_output(vector<tuple<int, int, int, int>> detected_crispr, string outfile_name) {
    std::stringstream ss;
    ss << "Number of repetitions\tDR length\tPosition of first DR\tPosition of last DR\n";

    for (auto &crispr_chain : detected_crispr) {
        int length = get<1>(crispr_chain);
        int start = get<2>(crispr_chain);
        ss << get<0>(crispr_chain) << "\t" << length
           << "\t"
           << start <<
           "\t" << get<3>(crispr_chain) << "\n" << endl;
    }
    std::string str = ss.str();
    ofstream outfile;
    outfile.open(outfile_name);
    outfile << str;
    outfile.close();
}


//selects repeats that meet basic criteria
vector<tuple<int, int, int>> select_candidates(cst_sada<> cst) {
    vector<tuple<int, int, int>> candidate_list;
    // iterate over all nodes
    for (auto it = cst.begin(); it != cst.end(); ++it) {
        if (cst.depth(*it) >= MIN_LENGTH && it.visit() == 1) {  // node visited for the first time
            auto v = *it;       // get the node by dereferencing the iterator
            // if depth of node is more than 1 and label has more than min_rep occs
            if (cst.depth(v) <= MAX_LENGTH) {
                // process node
                // if candidate length is within desired limits
                //if (candidate.length() >= MIN_LENGTH && candidate.length() <= MAX_LENGTH) {
                // string depth, leftmost leaf in SA and rightmost leaf in SA
                // we can obtain the number of occs using rb-lb+1
                candidate_list.emplace_back(cst.depth(v), cst.lb(v), cst.rb(v));
                //}
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
    vector<tuple<int, int, int, int>> crispr_list;
    auto rs = wt.range_search_2d(lb, rb, 0, sa.size());
    int number_of_points = get<0>(rs);
    auto points = get<1>(rs);
    //obtener primer valor y meterlo en cadena
    tuple<int, int, int, int> current_crispr = {1, len, (int) points.at(0).second, (int) points.at(0).second};
    //ir extrayendo valores, si el valor en cadena empieza entre y-s1-l e y-s2-l entonces añadir a cadena (reps+=1 y actualizar última occ)
    for (int k = 1; k < number_of_points; k++) {
        int c = points.at(k).second;
        int latest_pos = get<3>(current_crispr);
        if (c >= latest_pos + len + MIN_SPACER_LENGTH and
            c <= latest_pos + len + MAX_SPACER_LENGTH + 2) {
            get<0>(current_crispr) = get<0>(current_crispr) + 1;
            get<3>(current_crispr) = c;
        }
        //si no, entonces si la cadena tiene largo > 1, añadir a lista final (si reps > 1, o sea primer item)
        //empezar nueva cadena
        else {
            if (get<0>(current_crispr) > 1) {
                get<2>(current_crispr) = get<2>(current_crispr) + 1;
                get<3>(current_crispr) += 1;
                crispr_list.emplace_back(current_crispr);
            }
            current_crispr = {1, len, c, c};
        }
    }
    if (get<0>(current_crispr) > 1) {
        get<2>(current_crispr) = get<2>(current_crispr) + 1;
        get<3>(current_crispr) += 1;
        crispr_list.emplace_back(current_crispr);
    }

    //remember to add 1 to positions
    return crispr_list;
}


vector<tuple<int, int, int, int>>
find_crispr(string filename, tuple<double, double, double> &time, tuple<int, int, int, int, int> &size) {
    // construct needed structures
    cst_sada<> cst;
    construct(cst, filename, 1);
    // get the SA associated with the suffix tree
    int_vector<> sa = create_sa(filename);
    // construct the wavelet tree using the SA
    wt_int<> wt;
    construct_im(wt, sa);
    std::chrono::steady_clock::time_point begin_selection = std::chrono::steady_clock::now();
    // a candidate is a tuple: text length, lb, rb
    vector<tuple<int, int, int>> candidate_list;
    // select candidates that meet basic criteria
    candidate_list = select_candidates(cst);
    std::chrono::steady_clock::time_point end_selection = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point begin_verification = std::chrono::steady_clock::now();
    //cout << sa << endl;
    vector<tuple<int, int, int, int>> pre_filtered_candidates;
    vector<tuple<int, int, int, int>> filtered_candidates;
    //start time
    // validate candidates and insert them in pre_filtered list
    for (int i = 0; i < candidate_list.size(); i++) {
        for (auto &cr : validate(candidate_list[i], sa, wt))
            pre_filtered_candidates.emplace_back(cr);
    }
    std::chrono::steady_clock::time_point end_verification = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point begin_filter = std::chrono::steady_clock::now();
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
    std::chrono::steady_clock::time_point end_filter = std::chrono::steady_clock::now();
    //end time
    get<0>(time) = std::chrono::duration_cast<std::chrono::milliseconds>(end_selection - begin_selection).count();
    get<1>(time) = std::chrono::duration_cast<std::chrono::milliseconds>(end_verification - begin_verification).count();
    get<2>(time) = std::chrono::duration_cast<std::chrono::milliseconds>(end_filter - begin_filter).count();
    get<0>(size) = size_in_bytes(cst); //suffix tree
    get<1>(size) = size_in_bytes(wt); //wavelet
    get<2>(size) = size_in_bytes(sa); //suffix array
    get<3>(size) = sizeof(candidate_list[0]) * candidate_list.size(); //candidate list after suffix
    get<4>(size) = sizeof(filtered_candidates[0]) * filtered_candidates.size(); //list of crispr

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

void run_test(vector<Crispr> ground_truth, string filename, string test_outfilename, string output_file) {
    // run algorithm to find crispr
    vector<tuple<int, int, int, int>> detected_crispr;
    tuple<double, double, double> time;
    tuple<int, int, int, int, int> size;
    detected_crispr = find_crispr(filename, time, size);
    std::stringstream ss;
    ss << "number of crispr\treported crispr\n" << ground_truth.size() << "\t" << detected_crispr.size() << "\n";
    ss << "selection time\tvalidation time\tfiltering time\n";
    ss << get<0>(time) << "\t" << get<1>(time) << "\t" << get<2>(time) << "\n";
    ss << "suffix tree size\twavelet tree size\tsuffix array size\tcandidate list size\tcrispr list size\n";
    ss << get<0>(size) << "\t" << get<1>(size) << "\t" << get<2>(size) << "\t" << get<3>(size) << "\t" << get<4>(size)
       << "\n";
    ss << "id \t repeats \t detected occs \t partially detected occs \n";
    // for each real crispr, check whether it was correctly detected
    for (auto &crispr : ground_truth) {
        verify_crispr(crispr, detected_crispr);
        ss << crispr.id << "\t" << crispr.positions.size() << "\t" << crispr.get_detected() << "\t"
           << crispr.get_partially_detected() << "\n";
    }
    std::string str = ss.str();
    ofstream outfile;
    outfile.open(test_outfilename);
    outfile << str;
    outfile.close();
    write_output(detected_crispr, output_file);
}


void memory_test(string filename, string test_outfilename) {
    vector<tuple<int, int, int, int>> detected_crispr;
    tuple<double, double, double> time;
    tuple<int, int, int, int, int> size;
    detected_crispr = find_crispr(filename, time, size);
    std::stringstream ss;
    ss << "selection time\tvalidation time\tfiltering time\n";
    ss << get<0>(time) << "\t" << get<1>(time) << "\t" << get<2>(time) << "\n";
    ss << "suffix tree size\twavelet tree size\tsuffix array size\tcandidate list size\tcrispr list size\n";
    ss << get<0>(size) << "\t" << get<1>(size) << "\t" << get<2>(size) << "\t" << get<3>(size) << "\t" << get<4>(size)
       << "\n";
    std::string str = ss.str();
    ofstream outfile;
    outfile.open(test_outfilename);
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
    //Clostridioides difficile  NZ https://crispr.i2bc.paris-saclay.fr/crispr/crispr_db.php?checked%5B%5D=NZ_LN614756
    //
    /*string tamanos[8] = {"100615", "191799", "503000", "1057280", "1967000", "4296000", "5482170", "11121000"};
    for (string s : tamanos) {
        file = "/home/anouk/Documents/memoria/data/resource_test/" + s + ".fasta";
        string test_output = "/home/anouk/Documents/memoria/data/resource_test/output/" + s + "_out.txt";
        memory_test(file, 2, test_output);
        cout << "done with " << s << endl;
    }*/
    tuple<double, double, double> time;
    tuple<int, int, int, int, int> size;
    string file = "/home/anouk/Documents/memoria/data/resource_test/503000.fasta";
    string outfile = "/home/anouk/Documents/memoria/data/resource_test/output/503000_test.txt";
    string test_output = "/home/anouk/Documents/memoria/data/output/NC_015138_test_nf.txt";
    memory_test(file, outfile);


}