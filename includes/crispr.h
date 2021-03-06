//
// Created by anouk on 15-03-20.
//

#ifndef CRISPRSDS_CRISPR_H
#define CRISPRSDS_CRISPR_H


#include <string>
#include <vector>


class Crispr {
public:
    std::string id;
    int length;
    std::vector<int> positions;
    int detected;
    int partially_detected;

public:
    Crispr(std::string id, int dr_length, std::vector<int> dr_positions);

    //tally the number of DRs detected by a candidate
    void check_candidate(int reps, int candidate_length, int start, int end);

    double percentage_detected();

    std::string get_text();

    int get_detected();

    int get_partially_detected();

};


#endif //CRISPRSDS_CRISPR_H
