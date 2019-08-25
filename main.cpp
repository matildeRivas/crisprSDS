#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sdsl/wavelet_trees.hpp>
#include <string>

using namespace std;
using namespace sdsl;

int main()
{
    cst_sct3<> cst;
    //string file ;
    //cin >> file;
    //construct(cst, file, 1);
    construct(cst, "/home/anouk/Documents/memoria/data/testSeq.fasta", 1);
    cout << "inner nodes : " << cst.nodes() - cst.csa.size() << endl;
    auto v = cst.select_child(cst.child(cst.root(), 'T'),1);
    auto d = cst.depth(v);
    cout << "v : " << d << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << endl;
    cout << "extract(cst, v) : " << extract(cst, v) << endl;
}