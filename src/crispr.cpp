//
// Created by anouk on 15-03-20.
//

#include "../includes/crispr.h"
#include <algorithm>

Crispr::Crispr(int dr_length, std::vector<int> dr_positions, std::string dr_text) {
    detected = 0;
    length = dr_length;
    positions = dr_positions;
    dr = dr_text;
}

//TODO: test

void Crispr::check_candidate(int reps, int cand_length, int start, int end) {
    //Discard candidate if its DR length is a different length and if it is out of range of the crispr
    if (cand_length == length and end <= positions.back()) {
        //If the first and last candidate DR positions are in the array of occurences, the DR repetitions are counted as detected
        if (std::find(positions.begin(), positions.end(), start) != positions.end() and
            std::find(positions.begin(), positions.end(), end) != positions.end()) {
            detected += reps;
        }

    }
}

double Crispr::percentage_detected() {
    return 100 * detected / positions.size();
}

