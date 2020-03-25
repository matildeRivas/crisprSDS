//
// Created by anouk on 15-03-20.
//
using namespace std;

#include "../includes/crispr.h"
#include <algorithm>

Crispr::Crispr(std::string name_id, int dr_length, std::vector<int> dr_positions) {
    id = name_id;
    detected = 0;
    partially_detected = 0;
    length = dr_length;
    positions = dr_positions;
}


void Crispr::check_candidate(int reps, int cand_length, int start, int end) {
    //Discard candidate if its DR length is a different length and if it is out of range of the crispr
    if (end <= positions.back()) {
        int delta = length - cand_length;
        //If the first and last candidate DR positions are in the array of occurences, the DR repetitions are counted as detected
        for (int i = 0; i <= delta; i++) {
            if (std::find(positions.begin(), positions.end(), start - i) != positions.end() &&
                std::find(positions.begin(), positions.end(), end - i) != positions.end()) {
                if (length == cand_length and i == 0) {
                    detected += reps;
                } else {
                    partially_detected += reps;
                }
                break;
            }
        }


    }
}

double Crispr::percentage_detected() {
    return 100 * (double) detected / positions.size();
}

int Crispr::get_detected() {
    return detected;
}

int Crispr::get_partially_detected() {
    return partially_detected;
}