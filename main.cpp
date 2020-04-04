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
vector<tuple<int, int, int>> select_candidates(cst_sada<> cst, int min_reps) {
    vector<tuple<int, int, int>> candidate_list;
    // iterate over all nodes
    for (auto it = cst.begin(); it != cst.end(); ++it) {
        if (cst.depth(*it) > 1) {  // node visited for the first time
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
    cout << len << " lb: " << lb << " rb: " << rb << endl;
    vector<tuple<int, int>> cont_pairs;
    // for each candidate in group
    for (int c = lb; c <= rb; c++) {
        int sa_value = sa[c];
        cout << "searching value after " << sa_value << " in range: " << sa_value + len + MIN_SPACER_LENGTH << " "
             << sa_value + len + MAX_SPACER_LENGTH << endl;
        // search for next repeat in group that's within spacer range
        auto rs = wt.range_search_2d(lb, rb, sa_value + len + MIN_SPACER_LENGTH,
                                     sa_value + len + MAX_SPACER_LENGTH + 2);
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


vector<tuple<int, int, int, int>>
find_crispr(string filename, int min_reps, tuple<double, double, double> &time, tuple<int, int, int, int, int> &size) {
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
    candidate_list = select_candidates(cst, min_reps);
    std::chrono::steady_clock::time_point end_selection = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point begin_prefilter = std::chrono::steady_clock::now();
    //cout << sa << endl;
    vector<tuple<int, int, int, int>> pre_filtered_candidates;
    vector<tuple<int, int, int, int>> filtered_candidates;
    //start time
    // validate candidates and insert them in pre_filtered list
    for (int i = 0; i < candidate_list.size(); i++) {
        for (auto &cr : validate(candidate_list[i], sa, wt))
            pre_filtered_candidates.emplace_back(cr);
    }
    std::chrono::steady_clock::time_point end_prefilter = std::chrono::steady_clock::now();
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
    get<1>(time) = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_prefilter - begin_prefilter).count();
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

void run_test(vector<Crispr> ground_truth, string filename, int min_reps, string test_outfilename, string output_file) {
    // run algorithm to find crispr
    vector<tuple<int, int, int, int>> detected_crispr;
    tuple<double, double, double> time;
    tuple<int, int, int, int, int> size;
    detected_crispr = find_crispr(filename, min_reps, time, size);
    std::stringstream ss;
    ss << "number of crispr\treported crispr\n" << ground_truth.size() << "\t" << detected_crispr.size() << "\n";
    ss << "selection time\tfiltering time\n";
    ss << get<0>(time) << "\t" << get<1>(time) + get<2>(time) << "\n";
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


void memory_test(string filename, int min_reps, string test_outfilename) {
    vector<tuple<int, int, int, int>> detected_crispr;
    tuple<double, double, double> time;
    tuple<int, int, int, int, int> size;
    detected_crispr = find_crispr(filename, min_reps, time, size);
    std::stringstream ss;
    ss << "selection time\tfiltering time\n";
    ss << get<0>(time) << "\t" << get<1>(time) + get<2>(time) << "\n";
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
    string file = "/home/anouk/Documents/memoria/data/NC_014392.fasta";
    string outfile = "/home/anouk/Documents/memoria/data/output/NC_014392_out2.txt";
    string test_output = "/home/anouk/Documents/memoria/data/output/NC_014392_test2.txt";

    vector<Crispr> genomeNC_014392;
    vector<int> NC_014392_1_pos{145343, 145410, 145476, 145541, 145608, 145675, 145741, 145807, 145873, 145939, 146003,
                                146068, 146134, 146200, 146266, 146332, 146398, 146464, 146531, 146599, 146665, 146732,
                                146797, 146861, 146928, 146993, 147058, 147123, 147189, 147255, 147321, 147387, 147453,
                                147518, 147584, 147650, 147716, 147780, 147845, 147911, 147977, 148041, 148107, 148174,
                                148240, 148306, 148372, 148438, 148502, 148568, 148633, 148699, 148766, 148833, 148900,
                                148966, 149033, 149099, 149165, 149232, 149298, 149363, 149428, 149494, 149559, 149624,
                                149689, 149754, 149820, 149885, 149951, 150018, 150084, 150150, 150216, 150278, 150345,
                                150411, 150477, 150543, 150610, 150677, 150743, 150809, 150876, 150939, 151006, 151071,
                                151137, 151203, 151269, 151335, 151401, 151467, 151535, 151601, 151668, 151733, 151799,
                                151867, 151933, 152001, 152067, 152133, 152200, 152266, 152332, 152399, 152465, 152531,
                                152597, 152663, 152730, 152795, 152862, 152927, 152991, 153058, 153122, 153189, 153254,
                                153319, 153385, 153450, 153516, 153581, 153646, 153711, 153777, 153842, 153907, 153972,
                                154038, 154103, 154168, 154234, 154299, 154365, 154431, 154497, 154563, 154629, 154694,
                                154760, 154825, 154891, 154956, 155022, 155087, 155153, 155218, 155283, 155349, 155415,
                                155481, 155547, 155614, 155680, 155746, 155812, 155878, 155944, 156009, 156075, 156140,
                                156206, 156271, 156336, 156402, 156468, 156535, 156602, 156668, 156734, 156801, 156866,
                                156932, 156999
    };
    Crispr NC_014392_1 = Crispr("NC_014392_1", 29, NC_014392_1_pos);
    genomeNC_014392.emplace_back(NC_014392_1);
    vector<int> NC_014392_2_pos{157108, 157176, 157241, 157307, 157373, 157438, 157503, 157569, 157635, 157701, 157768,
                                157834, 157901, 157967, 158033

    };
    Crispr NC_014392_2 = Crispr("NC_014392_2", 29, NC_014392_2_pos);
    genomeNC_014392.emplace_back(NC_014392_2);

    vector<int> NC_014392_3_pos{160079, 160145, 160211, 160277};
    Crispr NC_014392_3 = Crispr("NC_014392_3", 29, NC_014392_3_pos);
    genomeNC_014392.emplace_back(NC_014392_3);

    vector<int> NC_014392_4_pos{2433976, 2434043, 2434114, 2434180, 2434245, 2434311, 2434377, 2434444, 2434509,
                                2434578, 2434644, 2434712, 2434778, 2434844, 2434912, 2434980, 2435048
    };
    Crispr NC_014392_4 = Crispr("NC_014392_4", 30, NC_014392_4_pos);
    genomeNC_014392.emplace_back(NC_014392_4);

    vector<int> NC_014392_5_pos{2443408, 2443475};
    Crispr NC_014392_5 = Crispr("NC_014392_5", 30, NC_014392_5_pos);
    genomeNC_014392.emplace_back(NC_014392_5);

    vector<int> NC_014392_6_pos{2449154, 2449219, 2449285, 2449351, 2449418, 2449486, 2449554, 2449621, 2449688,
                                2449756, 2449822, 2449888, 2449952, 2450018, 2450084, 2450150, 2450216, 2450282,
                                2450347, 2450415, 2450482, 2450550, 2450615, 2450679
    };
    Crispr NC_014392_6 = Crispr("NC_014392_6", 30, NC_014392_6_pos);
    genomeNC_014392.emplace_back(NC_014392_6);

    vector<int> NC_014392_7_pos{2460383, 2460448, 2460514, 2460579, 2460646, 2460711, 2460777, 2460845, 2460910,
                                2460977, 2461043, 2461109, 2461175, 2461243, 2461309, 2461375, 2461442, 2461510,
                                2461576, 2461642, 2461709, 2461777, 2461849
    };
    Crispr NC_014392_7 = Crispr("NC_014392_7", 30, NC_014392_7_pos);
    genomeNC_014392.emplace_back(NC_014392_7);
    vector<int> NC_014392_8_pos{2464720, 2464789, 2464856, 2464924, 2464990, 2465058, 2465126, 2465193, 2465260,
                                2465326, 2465392, 2465459, 2465525, 2465591, 2465657, 2465723, 2465789, 2465855,
                                2465921, 2465987, 2466055, 2466120, 2466186, 2466251, 2466317, 2466383, 2466452,
                                2466520, 2466586, 2466650
    };
    Crispr NC_014392_8 = Crispr("NC_014392_8", 30, NC_014392_8_pos);
    genomeNC_014392.emplace_back(NC_014392_8);
    vector<int> NC_014392_9_pos{2477566, 2477631, 2477699, 2477767, 2477833, 2477899, 2477965, 2478033, 2478099,
                                2478165, 2478230, 2478296, 2478364, 2478430, 2478495, 2478559, 2478626, 2478694,
                                2478762, 2478828
    };
    Crispr NC_014392_9 = Crispr("NC_014392_9", 30, NC_014392_9_pos);
    genomeNC_014392.emplace_back(NC_014392_9);
    vector<int> NC_014392_10_pos{2486476, 2486542, 2486608, 2486676, 2486744, 2486810, 2486878
    };
    Crispr NC_014392_10 = Crispr("NC_014392_10", 30, NC_014392_10_pos);
    genomeNC_014392.emplace_back(NC_014392_10);

    run_test(genomeNC_014392, file, min_reps, test_output, outfile);
}