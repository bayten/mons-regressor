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

void get_rand_data(Mat<int>* mat, Vec<int>* vec);
void get_rand_matrix(Mat<bool>* mat);

int main(void) {
    srand(time(0));
    init_logging();

    // RandomLBBuilder<int> lb_builder;
    GeneticDualizer<int, bool> test_object(5, 0.1, false, kWeightedScore, kDelayConverged, 0.001);
    LOG_(debug) << "Test of GeneticDualizer has been initiated.";

    // Mat<int> rand_data(10, 5);
    // Vec<int> rand_target(10);
    Mat<bool> rand_mat(100, 100);
    get_rand_matrix(&rand_mat);
    LOG_(debug) << "Random L matrix:";
    LOG_(debug) << rand_mat;

    // SampleHandler<int> sh_obj;
    // SampleSet<int> train = sh_obj.make_samples(rand_data, rand_target);
    // LOG_(debug) << "Train_data:" << train;
    //
    // for (int i = 0; i < 2; i++) {
    //     LOG_(debug) << "Processing first iteration for " << i << " class...";
    //     ElColl<int> local_basis = lb_builder.build_lb(train, i);
    //     LOG_(debug) << "Local basis:" << local_basis;

    LOG_(trace) << "Processing Genetic Algorithm...";
    test_object.set_matrix(rand_mat);
    // test_object.set_init_data(train, local_basis, i);
    Vec< Vec<bool> > encoded_colls = test_object.execute_ga();
    // Vec< ElColl<int> > new_colls = test_object.decode_collections(encoded_colls);

    return 0;
}

void get_rand_data(Mat<int>* mat, Vec<int>* vec) {
    int sx = mat->get_sx();
    int sy = mat->get_sy();

    for (int i = 0; i < sx; i++) {
        (*vec)[i] = (2*i < sx) ? 0 : 1;
        bool is_unique = 1;
        do {
            is_unique = 1;
            for (int j = 0; j < sy; j++)
                (*mat)[i][j] = rand() % 10;

            for (int k = 0; k < i; k++) {
                if ((*mat)[i] == (*mat)[k]) {
                    is_unique = 0;
                    break;
                }
            }
        } while (!is_unique);
    }
}

void get_rand_matrix(Mat<bool>* mat) {
    int sx = mat->get_sx();
    int sy = mat->get_sy();

    for (int i = 0; i < sx; i++)
        for (int j = 0; j < sy; j++)
            (*mat)[i][j] = rand() % 2;
}
