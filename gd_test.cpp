/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include <ctime>
#include "log_wrapper.h"
#include "default_types.h"
#include "LBBuilder.h"
#include "GeneticDualizer.h"
#include "CollFamily.h"

#include "SampleSet.h"
#include "SampleHandler.h"

void get_rand_matrix(Mat<bool>* mat);

int main(void) {
    srand(time(0));
    init_logging();

    GeneticDualizer<int, bool> test_object(5, 0.1, false, kWeightedScore, kDelayConverged, 0.001);
    LOG_(debug) << "Test of GeneticDualizer has been initiated.";

    Mat<bool> rand_mat(100, 100);
    get_rand_matrix(&rand_mat);
    LOG_(debug) << "Random L matrix:";
    LOG_(debug) << rand_mat;

    LOG_(trace) << "Processing Genetic Algorithm...";
    test_object.set_matrix(rand_mat);
    Vec< Vec<bool> > encoded_colls = test_object.execute_ga();

    return 0;
}

void get_rand_matrix(Mat<bool>* mat) {
    int sx = mat->get_sx();
    int sy = mat->get_sy();

    for (int i = 0; i < sx; i++)
        for (int j = 0; j < sy; j++)
            (*mat)[i][j] = rand() % 2;
}
