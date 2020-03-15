//
// Created by anouk on 15-03-20.
//

#ifndef CRISPRSDS_CRISPR_H
#define CRISPRSDS_CRISPR_H


#include <string>
#include <vector>


class Crispr {
    int length;
    std::vector<int> positions;
    int detected;
    string dr;

    Crispr(int dr_length, std::vector<int> dr_positions, string dr_text);

    //tally the number of DRs detected by a candidate
    void check_candidate(int reps, int candidate_length, int start, int end);

    double percentage_detected();

};


#endif //CRISPRSDS_CRISPR_H
