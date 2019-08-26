#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sdsl/wavelet_trees.hpp>
#include <string>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    cst_sct3<> cst;
    int min_reps;
    if (argc < 3){
         min_reps = 2;
    }
    else {
        char * pEnd;
         min_reps = strtol(argv[2], NULL, 10);
    }
    //string file ;
    //cin >> file;
    //construct(cst, file, 1);
    construct(cst, "/home/anouk/Documents/memoria/data/testSeq.fasta", 1);
    //todo: get text with more than 1 non leaf child
    for (auto it=cst.begin(); it!=cst.end(); ++it) {
        if (cst.depth(*it) > 1 && it.visit() == 1) {  // node visited the first time
            auto v = *it;       // get the node by dereferencing the iterator
            if (cst.node_depth(v) >= 2 && cst.size(v) >= min_reps) {   // if depth of node is more than 1
                // process node, e.g. output it in format d-[lb, rb]
                cout << "node: " << extract(cst, v) << endl;
                cout<<cst.depth(v)<<"-["<<cst.lb(v)<< ","<<cst.rb(v)<<"]"<<endl;
            } else { // skip the subtree otherwise
                it.skip_subtree();
            }
        }
    }
}